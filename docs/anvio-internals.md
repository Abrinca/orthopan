# anvi'o pangenomics internals: research notes for OrthoPan

This document covers how anvi'o builds, scores, orders, stores and displays a pangenome. It is meant to guide the OrthoPan design (see `VISION.md`).

All paths are relative to `resources/anvio/anvio/` unless stated otherwise. anvi'o is GPL-3.0, so algorithms are described in prose and pseudocode only. Anything marked **verified** was checked by running anvi'o's own `homogeneityindex.py` against an independent reimplementation on random alignments (300 cases, matching to 1e-12).

---

## 0. Pipeline overview (`anvi-pan-genome`)

`Pangenome.process()` (`panops.py:1690-1784`) runs these steps in order:

1. Load genomes from the genomes-storage DB, then run sanity checks. More than 100 genomes is a hard error unless the user passes `--I-know-this-is-not-a-good-idea` (`panops.py:1416-1429`).
2. Get gene clusters, either de novo (DIAMOND/BLAST → MCL, `panops.py:1621-1649`) or from a user-supplied `gene-clusters-txt` (`panops.py:1652-1687`). The second route is how external clusterings (such as OrthoFinder output) get into anvi'o today.
3. Align every non-singleton cluster in parallel (`panops.py:1490-1618`).
4. Build the frequency and presence/absence matrices and per-cluster statistics, apply the `--min-occurrence` filter, and store both matrices as "views" (`panops.py:1130-1232`).
5. Store cluster membership and alignment summaries in the `gene_clusters` table (`panops.py:1444-1463`).
6. Hierarchically cluster the gene clusters, which gives the item order (`panops.py:1275-1319`).
7. Add one "forced synteny" item order per genome (`panops.py:1235-1272`).
8. Order genomes (layers) and compute layer metadata (`panops.py:1352-1396`).
9. Compute homogeneity indices and AAI (`panops.py:1322-1349`, `dbops.py:2288-2391`).

---

## 1. Geometric homogeneity

Source: `homogeneityindex.py:71-148`.

### Input
The input is the aligned amino-acid strings of **every gene in the cluster**, pooled across genomes (`homogeneityindex.py:162-169`). Only the literal character `-` counts as a gap (`homogeneityindex.py:82, 90`). `.` and lowercase characters are treated as residues.

### Binary encoding
The alignment becomes an n × L gap matrix G, with 1 for a gap and 0 for a residue (n = genes, L = columns). anvi'o packs it two ways as Python big integers (`homogeneityindex.py:71-95`):
- **By residue:** one integer per column, holding n bits (one per gene).
- **By gene:** one integer per gene, holding L bits.

### Site-level ("vertical", by residue) score: `homogeneityindex.py:109-126`
anvi'o does **not** compare each column with a consensus. It compares **every column with every other column**:

```
for each column i:
    u_i = mean over columns j≠i of  (#genes whose gap state agrees in i and j) / n
by_residue = mean_i u_i
```
Agreement is counted as the zeros in `popcount(col_i XOR col_j)`, done through `bin(...).count('0')` on a zero-filled string (`homogeneityindex.py:118-119`).

### Gene-level ("horizontal", by gene) score: `homogeneityindex.py:131-144`
This is the same idea with rows and columns swapped:
```
for each gene a:
    v_a = mean over genes b≠a of (#columns where a and b agree on gap/residue) / L
by_gene = mean_a v_a
```

### Combination
- Full mode: `geometric = (by_residue + by_gene) / 2`, an unweighted arithmetic mean (`homogeneityindex.py:146`).
- Quick mode (`--quick-homogeneity`): the score is just `by_residue`, with no gene-level term (`homogeneityindex.py:128-129`; CLI flag at `cli/compute_gene_cluster_homogeneity.py:151-153`).

No other normalisation is applied. Both terms are already in [0, 1].

### Closed form (verified): O(n·L) instead of O(L²·n + n²·L)
Summing agreements over ordered pairs separates by gene (or by column). Let `k_g` be the number of gaps in gene g, and `m_c` the number of gaps in column c:

```
by_residue = Σ_g [k_g(k_g−1) + (L−k_g)(L−k_g−1)] / (n · L · (L−1))
by_gene    = Σ_c [m_c(m_c−1) + (n−m_c)(n−m_c−1)] / (L · n · (n−1))
```

So the index depends **only on the row and column gap counts** of the alignment, which you get from one pass over the matrix. For OrthoPan this means:
- no pairwise XOR sweep is needed;
- with polars or NumPy it is two `sum` reductions over a boolean matrix per alignment, and it can be batched across all MSAs;
- a Rust extension is not needed for this metric.

Some consequences of the formula:
- `by_residue = 1` only when every gene is either all residues or all gaps.
- `by_gene = 1` only when every column is either all gaps or gap-free.
- In practice the index is 1.0 only for gap-free alignments, which matches the blog's statement (`resources/anvio-blog.md:301`).
- A truncated gene (one gene gapped over a fraction f of the columns) lowers both terms. Example: `["AAAA","AA--"]` gives `by_residue = 0.667`, `by_gene = 0.5`, `geometric = 0.583`.

### Edge cases
| Case | Behaviour | Where |
|---|---|---|
| 1 sequence | returns 1.0 | `homogeneityindex.py:106-107` |
| 0 sequences | IndexError, caught by the worker, all indices set to −1 | `homogeneityindex.py:108`, `dbops.py:2373-2389` |
| L = 1 (one column) with n ≥ 2 | division by zero (`len(differences)==0`), so −1 | `homogeneityindex.py:124` |
| Sequences of unequal length (alignment failed) | IndexError, so −1 | `dbops.py:2376-2389` |
| All-gap columns | kept. They **lower** `by_residue`, because they add a gap to every gene's k_g. They do not hurt `by_gene` (m_c = n gives full agreement). | formula above |
| Failed alignment | the raw unaligned sequences are stored as if aligned (`panops.py:1594-1600`), and homogeneity then fails, so −1 | blog line 297 |

### Complexity in anvi'o
- Site-level work is O(L² · n) bit operations. Because they run on Python big-int strings, the real cost is O(L² · n) character work.
- Gene-level work is O(n² · L).
- For a 1,000-column alignment with 2,000 genes, that is about 10⁹ string-character operations per cluster for site-level and about 4×10⁹ for gene-level. This is why `--quick-homogeneity` exists (blog line 303).

---

## 2. Functional homogeneity

Source: `homogeneityindex.py:29-68`, `utils.py:3395-3412`, `constants.py:438-481`.

### Algorithm
```
score = max = 0
for each column:
    for each unordered pair of sequences (a, b):
        if neither residue is '-':
            max += 3
            if r_a == r_b and r_a not in {X, J, B, Z}: score += 3
            elif conserved(r_a, r_b):                 score += 2
            else:                                     score += 0
functional = score / max
```

- **Gaps:** any pair that involves a gap contributes nothing to either the score or the maximum (`homogeneityindex.py:53`). Gap–gap and gap–residue pairs are therefore neutral, as the unit tests confirm (`tests/unit/test_homogeneityindex.py`, "two gaps don't reduce the score").
- **Final score:** a single ratio over the whole alignment. It is **not** an average of per-column ratios, so columns with more non-gap residues weigh more.

### Amino-acid groups ("conserved" relation)
There is no substitution matrix. `conserved(a, b)` checks whether b is in the conserved group of a (`constants.py:438-456, 460-481`):

| Residue a | Counts as conserved with |
|---|---|
| A V M C | L V I M C H A |
| L I | L V I M C H A, J |
| F W | F W Y |
| K R | K R H |
| H | H Y, plus anything in Nonpolar or Bases |
| Y | H Y, plus F W |
| S T | S T |
| N | Q N B |
| Q | Q N Z |
| D | D E B |
| E | D E Z |
| B / Z / J | {B N D} / {Z Q E} / {J L I} |
| G P X | nothing |

I checked the relation by hand and it is symmetric. Ambiguity codes behave as follows: identical X/X scores 0, while B/B, Z/Z and J/J score 2 through their group.

### Pitfall
`amino_acid_property_group[r_a]` raises KeyError for any residue not in the table: `*`, `U`, `O`, lowercase letters, `.`. The worker catches it and sets **all** indices of the cluster to −1 (`dbops.py:2373-2389`). Because only the first residue of a pair is looked up, failure depends on sequence order. OrthoPan should normalise or upper-case residues and treat unknown characters explicitly.

### Edge cases
- 1 sequence: 1.0. 0 sequences: 0.0 (`homogeneityindex.py:36-39`).
- No column has two non-gap residues (for example two non-overlapping fragments): `max = 0`, so ZeroDivisionError and then −1.

### Closed form
Per column, take the residue counts c_x over non-gap residues, with R = Σ c_x:
- `max_c = 3·C(R,2)`
- `score_c = 3·Σ_{x∉XJBZ} C(c_x,2) + 2·[Σ_{x<y} c_x c_y·cons(x,y) + Σ_{x∈JBZ} C(c_x,2)]`

This is O(L · (n + 25²)) per alignment, and it vectorises as a counts matrix (L × 25) multiplied by a 25 × 25 score matrix.

---

## 3. Combined index and other per-cluster metrics

- **Combined homogeneity** is the harmonic mean `2·f·g / (f+g)`, set to 0 when both are 0 (`homogeneityindex.py:175-178`).
- **AAI_min / AAI_max / AAI_avg** (`dbops.py:2249-2285`): over all pairs of aligned sequences, identity is matches divided by alignment length. Gap–gap positions **count as matches**, and the source carries a FIXME saying so (`dbops.py:2272`). A single-sequence cluster returns (0, 0, 0), not 1.
- Computed in `process_gene_clusters` (`panops.py:1143-1161`) and stored as items additional data (`panops.py:1222-1227`):
  - `num_genes_in_gene_cluster`: the total number of genes, paralogs included.
  - `num_genomes_gene_cluster_has_hits`: the number of genomes with at least one gene.
  - `max_num_paralogs`: the largest copy number in any single genome.
  - `SCG`: 1 when **every** genome has exactly one copy (`set(counts) == {1}`), meaning single-copy core.
- Per-genome layer data (`panops.py:1364-1396`): total length, GC, completion/redundancy if available, `num_genes`, `avg_gene_length`, `num_genes_per_kb`, `singleton_gene_clusters`, `num_gene_clusters`, and `num_gene_clusters_raw` when filtering was applied.
- When homogeneity fails, every one of these indices (AAI included) is −1. Downstream filters must exclude −1.

---

## 4. Paralogs

- **Homogeneity and AAI:** paralogs are not treated specially. All genes in a cluster are pooled into one sequence list (`homogeneityindex.py:166-169`, `dbops.py:2259-2263`).
  - Every pair counts equally, so a genome with k paralogs contributes O(k·n) of the O(n²) pairs to the functional and gene-level scores, and k rows to the site-level score.
  - Expanded families dominate, and a divergent paralog pulls the whole cluster down. There is no per-genome weighting or representative selection.
- **Presence/absence view:** binary, so paralogs are invisible (`panops.py:1149`).
- **Frequency view:** raw copy counts per genome (`panops.py:1148`).
- **Clustering:** the "frequency" item order applies `log10(x+1)` (`data/clusterconfigs/pan/frequency`, `clustering.py:353-354`).
- **Paralogs show up only in:** `max_num_paralogs`, the SCG flag, the frequency view and order, and the frequency-based genome tree.

---

## 5. Ordering

### Genomes (layers): `panops.py:1352-1360`, `clustering.py:87-161`
- Two trees are always computed and stored as newick in `layer_orders`:
  - `gene_cluster presence absence`;
  - `gene_cluster frequencies`.
- Data: the genome × cluster matrix (transposed view). Each genome vector is **L1-normalised** first (`norm='l1'`, `clustering.py:87, 155-156`), so it becomes "the share of the genome's clusters", which removes genome-size effects.
- Distance and linkage default to Euclidean and Ward (`constants.py:134-135`) and can be changed with `--distance` / `--linkage`. The clustering is scipy `hierarchy.linkage` on the raw vectors, which computes a full pdist internally (`clustering.py:253-262`).
- The blog argues that frequency-based ordering matches ANI better for isolates, while MAGs lack paralogs (`resources/anvio-blog.md:600-620`).

### Gene clusters (items): `panops.py:1275-1319`
There are two clustering configs, `presence-absence` (the default, `constants.py:112`) and `frequency` (`data/clusterconfigs/pan/*`). anvi'o rewrites them on the fly to append a weighted "occupancy" feature:
- **presence-absence:** the binary row (genomes, not normalised) plus `num_genomes_gene_cluster_has_hits` **repeated round(N/2) times** as extra columns (`panops.py:1299-1301`). Duplicated columns survive config parsing (`clusteringconfuguration.py:133`, `utils.py:1657`). In effect, occupancy gets weight ≈ N/2 in Euclidean space, so clusters sort into a core → accessory → singleton gradient first and by co-occurrence pattern second.
- **frequency:** `log10(count+1)` per genome plus `num_genes_in_gene_cluster` repeated round(√N) times, also log-scaled (`panops.py:1302-1304`).
- Matrices are concatenated (`clustering.py:332-371`), then clustered with Euclidean/Ward and stored as newick in `item_orders` (`dbops.py:6621-6652`).
- **Hard soft-limit:** when there are more than 20,000 clusters, item clustering is skipped unless `--enforce-hierarchical-clustering` is given (`constants.py:123`, `panops.py:1187-1200`).
- **Forced synteny <> genome:** one order per genome. Clusters are listed in the order of gene-caller IDs in that genome, followed by every cluster missing from it. These are stored as plain lists, not newick (`panops.py:1254-1262`). The synteny edge list is computed and thrown away (`panops.py:1264-1272`, FIXME).

### Precomputed vs on demand
- **Precomputed and stored:** both genome trees, both cluster trees, and the synteny lists.
- **On demand in the UI:** ordering items by any additional-data column, for example sorting by homogeneity (blog line 317), plus alphabetical and user orders (`interactive.py:551-619`).
- **Re-clustering:** after import this needs a CLI round trip (`anvi-import-misc-data` or a re-run).

---

## 6. Alignment

- **Aligner:** MUSCLE by default, with FAMSA available (`drivers/__init__.py:30-39`).
  - MUSCLE v5 runs as `muscle -align in.fa -output out.fa -threads N`, or with `-super5` through `MUSCLE_PARAMS` (`drivers/muscle.py:84-96`).
  - Every cluster gets its own temp directory, FASTA write and subprocess.
- **When:** during `anvi-pan-genome`, for every cluster with more than one gene, using `num_threads` worker processes and a `multiprocessing.Manager` queue (`panops.py:1499-1555`). Singletons are skipped. `--skip-alignments` disables alignment and therefore homogeneity (`panops.py:1323-1325`).
- **Storage:** the aligned sequence is **not** stored. Each gene keeps an `alignment_summary` run-length string of alternating residue and gap block lengths, starting with `.` (residue first) or `-` (gap first), for example `-|4|2|3|…` (`utils.py:1404-1441`). The full alignment is rebuilt on demand by pulling the raw AA sequence from the genomes-storage DB and replaying the runs (`utils.py:1444-1474`, `dbops.py:2104-2130`).
- **Re-computation:** `anvi-compute-gene-cluster-homogeneity` re-reads the summaries, rebuilds the alignments and recomputes the scores (`cli/compute_gene_cluster_homogeneity.py:110-111`).

---

## 7. Storage (pan DB, SQLite)

Created in `PanDatabase.create` (`dbops.py:~5519-5546`), with table definitions in `tables/__init__.py`.

| Table | Columns | Content |
|---|---|---|
| `gene_clusters` | gene_caller_id, gene_cluster_id, genome_name, alignment_summary (`tables/__init__.py:24-26`) | long-form membership with one row per gene; the only "primary" data |
| `views` | view_id, target_table (`:222`) | registry of matrix views |
| `gene_cluster_frequencies`, `gene_cluster_presence_absence` | item, layer, value (`:228`) | **dense** long-form matrices: one row for every (cluster, genome) pair, **zeros included** (`tables/views.py:64-70, 121-122`) |
| `item_additional_data` | item_name, data_key, data_value (TEXT), data_type, data_group (`:193-195`) | per-cluster EAV: stats, SCG, homogeneity, AAI |
| `item_orders` | name, type, data, additional (`:189-191`) | newick trees and comma-joined lists |
| `layer_additional_data`, `layer_orders` | EAV / newick (`:197-203`) | per-genome metadata and trees |
| `states`, `collections_*` | | UI state JSON, bins of clusters |
| `self` (meta) | key/value | num_gene_clusters, num_genes_in_gene_clusters, params, etc. |
| `pan_reaction_network_*` | | optional metabolic-network extension |

- **Stored outside the pan DB:** sequences and gene functions live in the separate genomes-storage DB (`genomestorage.py:244-247`). Cluster-level function summaries are rebuilt in memory on every load by iterating every cluster × genome × gene (`dbops.py:2846-2859`).
- **Stored:** membership, alignment run-lengths, both matrices, all scores, trees.
- **Recomputed on every load:** alignments (from run-lengths), the gene_clusters nested dict, and function summaries.

---

## 8. Scale bottlenecks (with evidence)

1. **Hard and soft limits in the code itself:**
   - more than 100 genomes is an error unless the user overrides it (`panops.py:1416-1429`);
   - more than 20k clusters means no item tree (`constants.py:123`);
   - more than 20k clusters or more than 150 genomes triggers a warning that `anvi-display-pan` may be "impossible" (`panops.py:1754-1773`).
2. **Dense long-form matrices:** both views store N × M rows including zeros. For the bundled test set (`resources/large-dataset`: 2,193 genomes × 50,515 OGs) that is about 110M rows per view, or about 220M SQLite rows with TEXT item and layer names. Rows are inserted through `executemany` (`tables/views.py:122`).
3. **EAV text tables:** every numeric per-cluster value is a TEXT row in `item_additional_data` and is re-parsed on load.
4. **Nested-dict memory:** `init_gene_clusters` builds `{cluster: {genome: [] for every genome}}` (`dbops.py:3412-3415`), which is O(M·N) Python objects even for empty cells. The function dict repeats this (`dbops.py:2846-2849`), and `load_pan_views` loads both full matrices into dicts (`dbops.py:3442-3452`).
5. **Hierarchical clustering:** scipy `linkage` on observation vectors means an O(M²) condensed distance matrix (20k items is about 1.6 GB of float64) and O(M²) Ward time. The occupancy columns duplicated round(N/2) times make the feature matrix wider still (N + N/2 columns).
6. **Homogeneity:** pure-Python O(L²·n + n²·L) string operations per cluster (§1).
   - For big clusters such as transposases this is the dominant cost, and alignment of such clusters often fails (`panops.py:1719-1733`).
   - The **whole** sequences dict is passed to every worker process (`dbops.py:2318-2319`).
   - Results come back through a `Manager().Queue()` one cluster at a time.
7. **Alignment:** one MUSCLE subprocess per cluster, with temp files. MUSCLE v5 `-align` scales roughly quadratically in n.
8. **Rendering:** the browser UI draws SVG, one element per (item, layer) cell in intensity mode (`data/interactive/js/drawer.js:1442-1470`, `svg-helpers.js:25ff`). M × N DOM nodes is unworkable at 10⁶ or more.
9. **Radial geometry:** `circlephylogram` is the default and uses `drawPie` with radii fixed per layer. A `phylogram` (linear) mode exists (`drawer.js:1246-1257`) but is secondary and still SVG.

---

## 9. Other things worth reusing or avoiding

**Reuse**
- The occupancy-weighted presence/absence clustering trick. It gives the familiar core → accessory gradient while keeping co-occurrence structure.
- L1 normalisation of genome vectors before genome clustering.
- Two parallel orderings (binary vs copy number) and per-genome synteny orders.
- The harmonic-mean combined score.
- The −1 sentinel idea, but as a null plus a reason code rather than a magic number.
- Run-length alignment summaries, if OrthoPan ever needs to store alignments compactly next to raw sequences. With OrthoFinder MSAs on disk, storing per-OG gap statistics is enough.

**Avoid or fix**
- Pairwise XOR sweeps. Use the closed form (§1, §2).
- Dense storage of zeros. Store membership long-form (OG, genome, gene) and derive counts sparsely.
- EAV TEXT columns.
- AAI counting gap–gap positions as matches.
- AAI = 0 for singletons.
- Order-dependent KeyError on non-standard residues.
- Treating raw unaligned sequences as alignments after a failure.
- Silent pooling of paralogs (see below).

**Behaviour to document if OrthoPan claims "anvi'o-compatible" scores**
- Only `-` is a gap.
- All-gap columns are kept and lower `by_residue`. This matters if OrthoPan subsets an MSA, for example after dropping paralogs: it must then decide whether to strip newly all-gap columns.
- OrthoFinder trims MSAs by default (columns that are 90% or more gaps; the `-z` flag disables trimming; check this against the installed OrthoFinder version). That systematically **raises** geometric homogeneity compared with anvi'o's untrimmed MUSCLE alignments.

---

## Implications for OrthoPan

1. **Homogeneity is cheap if computed correctly.**
   - Geometric: only per-row and per-column gap counts are needed.
   - Functional: a per-column residue-count matrix times a fixed 25 × 25 score table.
   - Both are O(n·L) per OG. 30k MSAs can be processed in seconds with NumPy (vectorised per MSA, multiprocessed across files). polars or Rust are optional.
   - Keep an "anvi'o-exact" mode, validated against `homogeneityindex.py` on fixtures, plus an optional improved mode.
2. **Paralogs: make the policy explicit.**
   - Default to anvi'o-compatible pooling.
   - Offer per-genome weighting (each genome's genes weighted 1/k so genomes count equally) or a one-representative-per-genome option (longest, or closest to the OG consensus).
   - Report `max_copies`, a single-copy flag and a single-copy-core flag as anvi'o does.
3. **MSA source.** OrthoFinder `-M msa` MSAs are trimmed by default, so either require `-z` for comparability or label scores as "on trimmed MSAs". OGs without MSAs (singletons, very large OGs) get null scores with a reason code instead of −1.
4. **Ordering.**
   - Precompute genome order from `SpeciesTree_rooted.txt`, with optional L1-normalised Ward on presence/absence and frequency matrices.
   - For OGs, use cheap orders by default: occupancy, then copy number, then a pattern hash; core → accessory gradient; synteny per reference genome; sort by any metric.
   - Run Ward with the occupancy-weight trick only when M ≤ ~20k, or use a scalable alternative (clustering unique presence patterns only, since many OGs share a pattern; or `fastcluster` / an approximate method).
5. **Storage.**
   - Write Parquet or DuckDB in long, sparse form: `membership(og, genome, gene)`, `og_stats(og, n_genes, n_genomes, max_copies, scg, geo, fun, comb, aai_*)`, `genome_stats`, `orders(name, kind, payload)`.
   - Derive the presence/frequency matrix as a bit-packed or uint8 array for the viewer: 2,193 × 50,515 ≈ 110 MB as uint8 and ≈ 14 MB bit-packed. The SQLite approach would need about 220M rows.
6. **Rendering.** The 110M-cell test set rules out SVG and DOM. Use WebGL or canvas textures with level-of-detail aggregation (for example, the density of present cells per screen pixel), which a uniform equal-area grid supports naturally.
