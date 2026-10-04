# OrthoFinder 3.1.5 outputs: inventory and what OrthoPan should use

Investigation date: 2026-10-04. Data: `resources/orthofinder-demo/` (read-only). OrthoFinder source consulted: `scratch/.venv/lib/python3.12/site-packages/orthofinder/` (3.1.5.post1.dev1). See `VISION.md`, `docs/design-decisions.md` and `docs/anvio-internals.md` for context.

## 0. The three runs

| Label | Path | What it is | Genomes | Tools (from `Log.txt`) |
|---|---|---|---|---|
| `core17` | `results/Results_Oct04` | regular run, `orthofinder -f faa-core` | 17 | DIAMOND, FAMSA, FastTree |
| `assign24` | `results/Results_Oct04_1` | `--assign faa-add --core results/Results_Oct04` | 17 + 7 | DIAMOND (profiles), FAMSA, **FastTree -fastest**, ASTRAL-Pro species tree |
| `denovo24` | `results-full-24/Results_Oct04` | regular run on **all 24** genomes | 24 | DIAMOND, FAMSA, FastTree |

**`results-full-24` is not a duplicate of the 17-genome run.** Its `Log.txt` lists 24 species (0–23, including all 7 `faa-add` genomes), it has 24 per-species files in `Putative_Xenologs/` and `Orthologues/`, and 576 (= 24²) `Blast*_*.txt.gz` files. It is a de novo run over the same 24 genomes that `assign24` covers, so it is the natural ground truth for judging `--assign`. Keep it.

Headline counts:

| | core17 | assign24 | denovo24 |
|---|---|---|---|
| genes | 40,085 | 56,760 | 56,760 |
| orthogroups (≥2 genes, `Orthogroups.tsv`) | 4,784 | 5,657 | 5,442 |
| unassigned genes (singleton OGs) | 3,519 | 4,873 | 4,663 |
| OGs in all genomes / single-copy core | 368 / 233 | 332 / 201 | 336 / 209 |
| MSAs (`MultipleSequenceAlignments/`) | 4,670 | 5,544 | 5,316 |
| resolved gene trees | 2,410 | 2,832 | 2,817 |
| duplication events (`Duplications.tsv`) | 3,975 | 5,563 | 5,724 |
| results dir size (bytes, excl. block overhead) | ~205 MB | ~185 MB | ~325 MB |

---

## 1. Inventory

Gene ID formats matter for joins:
- **Original IDs**, meaning the first word of the FASTA header, e.g. `gnl|extdb|FAM1079-i1-1.1_000590`. These appear in `Orthogroups*.tsv`, HOG `N*.tsv`, `Orthologues/`, `Putative_Xenologs/`, `Phylogenetically_Misplaced_Genes/` and `Orthogroup_Sequences/`.
- **Species-prefixed IDs**, e.g. `FAM1079-i1-1_1_gnl|extdb|FAM1079-i1-1.1_000590`. The genome name has `.` replaced by `_`, and the default option `-X` is not set. These appear in `MultipleSequenceAlignments/`, `Resolved_Gene_Trees.txt` and `Duplications.tsv`. To get the original ID back, strip `genome.replace('.', '_') + '_'`.
- **Internal IDs** `<speciesIdx>_<seqIdx>`, e.g. `0_590`. These appear only in `WorkingDirectory/`. The mappings are in `SequenceIDs.txt` and `SpeciesIDs.txt`.

### 1.1 Global files (one file for the whole run)

| File | Content | core17 size | OrthoPan use | Scaling to 10k genomes × 50k OGs (~25M genes) |
|---|---|---|---|---|
| `Orthogroups/Orthogroups.tsv` | OG × genome matrix of gene lists, OGs with ≥2 genes | 1.3 MB (4,784 rows) | **Matrix source (MVP)**. Also gives gene membership for the detail view | ~1–1.5 GB of wide text. Stream-parse into sparse Parquet |
| `Orthogroups/Orthogroups.GeneCount.tsv` | same shape, copy numbers plus a `Total` column | 227 KB | fast path for the count matrix when gene IDs are not needed | dense: 50k × 10k cells, ~1 GB |
| `Orthogroups/Orthogroups_UnassignedGenes.tsv` | singleton "OGs" (one gene each), same wide layout | 211 KB (3,519 rows) | optional: genome-specific singletons | **bad**: each row has 10k mostly-empty columns. With 100k+ singletons that is >1 GB of tabs |
| `Orthogroups/Orthogroups.txt` | `OGxxxxxxx: gene gene …`, **including** singletons | 1.4 MB (8,303 lines) | compact alternative to the two TSVs | ~ size of gene IDs, ~1 GB |
| `Orthogroups/Orthogroups_SingleCopyOrthologues.txt` | list of OG IDs, single copy in every genome | 2.3 KB (233) | SCO flag (also derivable) | tiny, tends to 0 at 10k genomes |
| `Species_Tree/SpeciesTree_rooted.txt` | rooted Newick with branch lengths | 539 B | **genome tree (MVP)** | small (~1 MB) |
| `Species_Tree/SpeciesTree_rooted_node_labels.txt` | same tree with internal labels `N0…N(n-2)` | 577 B | **preferred tree file**: labels link to HOG files and to `Duplications.tsv` | small |
| `Species_Tree/Orthogroups_for_concatenated_alignment.txt` | OGs used for the concatenated species-tree alignment | 2.7 KB (272) | not needed | absent in `--assign` runs |
| `Resolved_Gene_Trees/Resolved_Gene_Trees.txt` | one line per OG, `OGxxxxxxx: <newick>`, rooted and resolved, for OGs with ≥4 genes only | 2.1 MB (2,410 trees) | detail view gene tree; per-OG tree metrics | ~60 B/gene, so ~1.5 GB. One file, which is good, but it needs a byte-offset index |
| `Gene_Duplication_Events/Duplications.tsv` | one row per duplication: OG, species-tree node, gene-tree node, support, type, both gene lists | 1.8 MB (3,975) | detail view (where the paralogs arose) | multi-GB because it lists genes. Use the per-OG summary instead |
| `Gene_Duplication_Events/SpeciesTree_Gene_Duplications_0.5_Support.txt` | species tree with `node_<count>` labels | 683 B | optional per-clade annotation on the genome tree | small |
| `Comparative_Genomics_Statistics/Duplications_per_Orthogroup.tsv` | OG, dups (all), dups (≥50% support). Only OGs that have a tree | 36 KB (2,410) | **cheap per-OG track** | 1 row/OG |
| `Comparative_Genomics_Statistics/Duplications_per_Species_Tree_Node.tsv` | per species-tree node | 0.6 KB | genome-tree annotation | small |
| `Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv` | per-genome rows: genes, % in OGs, unassigned, species-specific OGs, histograms | 8 KB | **genome metadata track** (bottom panel) | ~100 rows × 10k columns, fine |
| `Comparative_Genomics_Statistics/Statistics_Overall.tsv` | run summary (G50, O50, counts, histograms) | 1.6 KB | run info panel | small |
| `Comparative_Genomics_Statistics/OrthologuesStats_{Totals,one-to-one,one-to-many,many-to-one,many-to-many}.tsv` | genome × genome counts of orthologue relationships | ~4 KB each | optional genome × genome heatmap | N² = 10⁸ cells, ~0.5 GB each |
| `Orthologues/<genome>.tsv` | for each OG, gene(s) of this genome mapped to orthologs in each other genome | 52 MB total (core17) | per-gene detail only | **quadratic**: ~105 MB at 24 genomes, so roughly 10–20 TB at 10k. Do not ingest |
| `Orthologues/Orthologues_<genome>/<A>__v__<B>.tsv` | the same data split per pair (17 × 16 files) | (in the 52 MB) | no | 10⁸ files. Not written by `--assign` |
| `Putative_Xenologs/<genome>.tsv` | OG, gene, "Other" gene (putative HGT partner) | 20 rows in total (core17), 138 (assign24) | optional flag | sparse |
| `Phylogenetically_Misplaced_Genes/<genome>.txt` | gene IDs removed as phylogenetically misplaced | 2 genes (core17), 10 (assign24) | optional QC | sparse |
| `Log.txt`, `Citation.txt` | command line, tools, species index | 1 KB | provenance (detect `--assign`, tools) | small |

`Statistics_Overall.tsv` mentions `Orthogroups_SpeciesOverlaps.tsv`, but 3.1.5 does not write that file. The `Gene_Trees/` directory (unresolved trees) from OrthoFinder 2 is also absent. Only `Resolved_Gene_Trees.txt` is written.

### 1.2 Per-OG files (one file per orthogroup)

| Directory | Content | Coverage | core17 | Scaling |
|---|---|---|---|---|
| `MultipleSequenceAlignments/OGxxxxxxx.fa` | aligned FASTA (FAMSA), species-prefixed IDs, 80-column lines | every OG with ≥2 genes except 114 OGs of size 2–3. Always contains **all** OG members (n_MSA = n_genes for 4,670/4,670) | 4,670 files, 16.9 MB | ~50k files, ~470 B/gene today. L grows with n, so expect 15–50 GB |
| `Orthogroup_Sequences/OGxxxxxxx.fa` | unaligned member sequences (original IDs) | every OG **including singletons** | 8,303 files, 13.5 MB | 100k+ files ≈ total proteome size. Prefer the input FASTAs |
| `Single_Copy_Orthologue_Sequences/OGxxxxxxx.fa` | unaligned sequences of SCO OGs | SCOs only | 233 files | tends to 0 |
| `Phylogenetic_Hierarchical_Orthogroups/N<k>.tsv` | one file per **internal species-tree node** (N1…N15 for 17 genomes) with HOGs rooted at that clade | **`N0.tsv` is missing**, see §3 | 15 files, 7.2 MB | ~10k files, each up to 10k columns wide. Do not ingest |

### 1.3 WorkingDirectory (internal, not a stable interface)

| Item | Content | core17 |
|---|---|---|
| `SpeciesIDs.txt` | `0: FAM1079-i1-1.1.faa` | 17 lines |
| `SequenceIDs.txt` | `0_0: gnl|extdb|FAM1079-i1-1.1_000001 IS30 family transposase, partial [...]`. **Keeps the full FASTA description**, which is a free annotation source | 40,085 lines |
| `Species<i>.fa`, `diamondDBSpecies<i>.dmnd` | renamed proteomes, DIAMOND DBs | 17 + 17 |
| `Blast<i>_<j>.txt.gz` | DIAMOND tabular, outfmt 6 (`qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore`), internal IDs | 289 files, 27 MB gz, 1.34M hits |
| `Alignments_ids/`, `Trees_ids/`, `Sequences_ids/` | pre-HOG MSAs, trees and sequences with internal IDs (the MCL-level OGs) | 4,268 `.fa` + 226 `.txt`; 2,806; 7,681 |
| `N0.ids.tsv` | **the N0 HOG table**, internal IDs | 4,798 rows |
| `OrthoFinder_graph.txt`, `clusters_OrthoFinder_I1.2.txt*` | MCL graph and clusters | 5.7 MB |
| `profile_sequences.10_km.fa(.dmnd)` | 10 k-means representative sequences per OG, used by `--assign`. **Written into the core directory by the `--assign` run** (mtime 16:03, after the core run finished) | 7.9 + 8.2 MB |

The assign run's WorkingDirectory has `previous_wd.txt`, an absolute path back to the core WorkingDirectory. So a `--assign` result directory is not self-contained. OrthoPan should read only the top-level result folders and never `WorkingDirectory/`.

### 1.4 Format samples (header + 2 rows, long cells truncated with `…`)

`Orthogroups/Orthogroups.tsv`
```
Orthogroup	FAM1079-i1-1.1	FAM13496-i1-1.1	FAM14217-p1-1.1	…	FAM3257-i1-1.1
OG0000000	gnl|extdb|FAM1079-i1-1.1_002191, gnl|extdb|FAM1079-i1-1.1_002192	gnl|extdb|FAM13496-i1-1.1_001855, …	…
OG0000001	gnl|extdb|FAM1079-i1-1.1_000680, gnl|extdb|FAM1079-i1-1.1_001255, …	…
```
`Orthogroups/Orthogroups.GeneCount.tsv`
```
Orthogroup	FAM1079-i1-1.1	FAM13496-i1-1.1	FAM14217-p1-1.1	FAM17891-i1-1.1	…	Total
OG0000000	2	2	3	0	…	162
OG0000001	7	13	3	20	…	125
```
`Orthogroups/Orthogroups_UnassignedGenes.tsv` (IDs continue after the last OG, so OG0004784 to OG0008302 are singletons)
```
Orthogroup	FAM1079-i1-1.1	…	FAM3257-i1-1.1
OG0004784	…	gnl|extdb|FAM20833-i1-1.1_002273	…
OG0004785	…	gnl|extdb|FAM19471-i1-1.1_000493	…
```
`Phylogenetic_Hierarchical_Orthogroups/N1.tsv`
```
HOG	OG	Gene Tree Parent Clade	FAM1079-i1-1.1	FAM13496-i1-1.1	FAM14217-p1-1.1	…
N1.HOG0000000	OG0000000	n0			gnl|extdb|FAM14217-p1-1.1_000081, …_000084, …_000471	…
```
`Species_Tree/SpeciesTree_rooted_node_labels.txt`
```
((FAM14217-p1-1.1:0.25276,(FAM19036-p1-1.1:0.08744,FAM19038-p1-1.1:0.09552)N3:0.19457)N1:0.46664,(…)N2:0.46664)N0;
```
`Resolved_Gene_Trees/Resolved_Gene_Trees.txt`
```
OG0000000: (((((FAM13496-i1-1_1_gnl|extdb|FAM13496-i1-1.1_001855:0.00054,FAM13496-i1-1_1_gnl|extdb|FAM13496-i1-1.1_001876:0.00054)n4:0.6411,…
```
`MultipleSequenceAlignments/OG0000100.fa`
```
>FAM1079-i1-1_1_gnl|extdb|FAM1079-i1-1.1_000590
----------------------------------------------------MTELE-KL
DA-------GEPYYF----------LDDQV-----AARKARAAKLCQEFN----SIPATE
```
`Gene_Duplication_Events/Duplications.tsv`. `Support` is the fraction of expected descendant species that retain both copies. `Type` is one of `Terminal` (2,162), `Non-Terminal` (1,635) or `Non-Terminal: STRIDE` (178).
```
Orthogroup	Species Tree Node	Gene Tree Node	Support	Type	Genes 1	Genes 2
OG0000000	N4	n1	0.3076923076923077	Non-Terminal	FAM13496-i1-1_1_gnl|extdb|…_001855, …	…
OG0000000	FAM13496-i1-1.1	n4	1.0	Terminal	FAM13496-i1-1_1_gnl|extdb|FAM13496-i1-1.1_001855	FAM13496-i1-1_1_gnl|extdb|FAM13496-i1-1.1_001876
```
`Comparative_Genomics_Statistics/Duplications_per_Orthogroup.tsv`
```
Orthogroup	Duplications (all)	Duplications (50% support)
OG0000000	96	91
OG0000001	62	55
```
`Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv` (transposed: rows are statistics, columns are genomes)
```
	FAM1079-i1-1.1	FAM13496-i1-1.1	…
Number of genes	2735	1702	…
Number of genes in orthogroups	2434	1443	…
```
`Comparative_Genomics_Statistics/OrthologuesStats_one-to-one.tsv`
```
	FAM1079-i1-1.1	FAM13496-i1-1.1	FAM14217-p1-1.1	…
FAM1079-i1-1.1	0.0	829.0	568.0	…
FAM13496-i1-1.1	829.0	0.0	517.0	…
```
`Orthologues/FAM1079-i1-1.1.tsv`
```
Orthogroup	Species	FAM1079-i1-1.1	Orthologs
OG0000000	FAM19036-p1-1.1	gnl|extdb|FAM1079-i1-1.1_002191, gnl|extdb|FAM1079-i1-1.1_002192	gnl|extdb|FAM19036-p1-1.1_002761
OG0000000	FAM14217-p1-1.1	gnl|extdb|FAM1079-i1-1.1_002191, …_002192	gnl|extdb|FAM14217-p1-1.1_000471, …
```
`Putative_Xenologs/FAM21789-i1-1.1.tsv` (assign24)
```
Orthogroup	FAM21789-i1-1.1	Other
OG0000177	gnl|extdb|FAM21789-i1-1.1_000028	gnl|extdb|FAM19471-i1-1.1_001738
OG0000177	gnl|extdb|FAM21789-i1-1.1_000028	gnl|extdb|FAM22472-i1-1.1_000686
```
`Phylogenetically_Misplaced_Genes/FAM21789-i1-1.1.txt` (assign24): one gene ID per line, e.g. `gnl|extdb|FAM21789-i1-1.1_000028`.

`WorkingDirectory/Blast0_1.txt.gz`
```
0_1	1_829	33.5	388	232	9	37	414	34	405	2.54e-52	176
0_2	1_1268	31.2	250	163	7	129	372	188	434	2.03e-34	128
```

### 1.5 Parsing pitfalls observed

- **Trailing space in MSAs.** 2,449 of 4,670 MSA files end with `" "` and no newline. This comes from `tools/trim.py`, which appends a space to test writability. A parser that does not strip whitespace sees one sequence that is one column longer than the others. Strip every line.
- **All-gap columns.** 593/4,670 MSAs (core17) contain all-gap columns, up to 908 of them. The pre-HOG alignments in `WorkingDirectory/Alignments_ids/` have none (0/4,268). The final MSAs are pruned from those alignments to the N0 HOG members (`file_updates/trees.py`: `pruned_alignments = {g: gene_dict[g] …}`), and no columns are removed afterwards. anvi'o's by-residue term counts all-gap columns against every gene, so **drop all-gap columns before computing homogeneity**.
- OG IDs are assigned in order of decreasing OG size in all three runs (0 increases between consecutive IDs). Singleton IDs follow the last multi-gene OG.

---

## 2. Core run versus `--assign` run

### 2.1 What `--assign` produces

From `run/main.py` (`qFastAdd` branch) and the outputs:
1. Builds a profile DB of 10 k-means representative sequences per core OG. It is written into the **core** `WorkingDirectory/profile_sequences.10_km.fa`.
2. Searches new genomes against the profiles (`Blast<new>_-1.txt.gz`) and assigns genes to the core **N0 HOGs**.
3. Clusters the remaining genes within new-species clades (`clusters_OrthoFinder_unassigned_clade_*`).
4. **Re-infers MSAs (FAMSA) and gene trees (FastTree `-fastest`) for every OG**, infers a new species tree with ASTRAL-Pro, then runs the full orthology step: HOGs, duplications, orthologues, xenologs, statistics.

So `assign24` is **not** only orthogroup assignments. Everything is regenerated for all 24 genomes:

| Output | assign24 covers the 7 new genomes? | Evidence |
|---|---|---|
| Orthogroups.tsv / GeneCount / Unassigned | yes | 24 genome columns. 14,717 new-genome genes in OGs, 1,958 unassigned |
| MSAs | yes, realigned | `OG0000000.fa`: 236 sequences (162 in core17), 75 from new genomes |
| Resolved gene trees | yes, re-inferred (FastTree-fastest) | 2,492 of 2,832 trees contain new-genome leaves |
| Species tree | yes, re-inferred (ASTRAL-Pro from gene trees) | no `Orthogroups_for_concatenated_alignment.txt` |
| HOGs `N1…N22.tsv` | yes | 24 genome columns |
| Duplications | yes | 1,169 terminal duplications on new genomes, 45 species-tree nodes |
| Statistics_PerSpecies / OrthologuesStats | yes | 24 columns |
| Orthologues | per-genome TSVs only (24) | the per-pair `Orthologues_<genome>/` directories are skipped for `--assign` (`main.py:828`) |
| Xenologs / misplaced genes | yes | 24 xenolog files |
| DIAMOND all-vs-all | **no** | only new-vs-profile and within-clade searches: 49 files, 3.0 MB, against 576 files and 55 MB in denovo24 |

### 2.2 OG IDs are not stable

OG IDs are reassigned by size in every run. Comparing core17 with assign24 on the 17 core genomes (gene-set comparison restricted to core genes):

- Only **19 of 4,784** core OGs keep both their ID and their core-gene set. 35,884 of 36,566 core genes end up under a different OG ID.
- Ignoring IDs, **4,280 / 4,784 (89%)** core OGs survive with exactly the same core-gene set. The remaining 504 were split, merged or lost members.
- 687 core singletons were absorbed into OGs, and 219 assign24 OGs contain only new-genome genes. These were appended within the same size-sorted ID space.
- Species tree, core17 against assign24 restricted to the 17 core taxa: RF = 2 (of 14 + 14 bipartitions). core17 against denovo24: RF = 0.

Comparing assign24 with denovo24 (same 24 genomes):

- Gene-pair co-membership: 92% of the gene pairs that assign24 groups together are also grouped in denovo24. denovo24 groups 12% more pairs (736,781 against 655,811), so assign is slightly more split. There are 5,657 assign OGs against 5,442 de novo OGs.
- 4,377 OGs are identical as gene sets. Only 41 of them also share the same ID.
- Species tree RF = 4 (of 21 + 21 bipartitions).
- Metric distributions shift slightly. For example, median gene-tree length is 5.82 in assign24 against 5.98 in denovo24, a consequence of FastTree-fastest and different alignments.

### 2.3 Implications for OrthoPan's input contract

1. **One results directory is one immutable snapshot.** Take OG IDs as opaque labels valid only within that snapshot. Store a run fingerprint (path, `Log.txt` command line, species list, file mtimes). Never join OG IDs across runs. If cross-run linking is needed later, match on gene-set Jaccard.
2. **Read the newest results directory, not the core one.** After `--assign`, every output in the core directory is stale: it covers 17 of 24 genomes and uses old IDs. Everything in the assign directory is complete and self-consistent at the top level.
3. Do not depend on `WorkingDirectory/`. In assign runs it is incomplete (`previous_wd.txt` points elsewhere, and there is no all-vs-all DIAMOND), and it is internal anyway.
4. At scale (`--core`/`--assign` is OrthoFinder's recommended path), all per-OG artefacts OrthoPan cares about are still produced: MSAs, trees, duplications and HOGs. The weaker items are DIAMOND hits (gone), per-pair orthologue files (gone) and tree quality (FastTree-fastest). The input contract should therefore require `Orthogroups.tsv` + the species tree, and treat everything else as optional and as possibly produced by different tools. Record the tool names from `Log.txt` (`MSA program`, `Tree program`) in the dataset metadata.
5. Detect `--assign` from the `Log.txt` command line (`--assign … --core …`) and show it in the UI.

---

## 3. `Orthogroups.tsv` versus `Phylogenetic_Hierarchical_Orthogroups/N0.tsv`

**In OrthoFinder 3.1.5, `N0.tsv` is not written, and `Orthogroups.tsv` *is* the N0 HOG set.** This holds for all three runs: no `N0.tsv`, while `N1.tsv…` are present.

- The source shows this directly. With the defaults `fix_files=True` and `rm_legacy=True` (`run/process_args.py:137,147`), `comparative_genomics/orthologues.py` calls `file_updates.update_output_files()`. That function converts N0 HOGs into the new orthogroups (`ogs.post_hogs_processing`), prunes the MSAs and gene trees to HOG members, renumbers the result into `Orthogroups.tsv`, and then **deletes** `Phylogenetic_Hierarchical_Orthogroups/N0.tsv` (`os.remove(files.FileHandler.HierarchicalOrthogroupsFNN0())`). `--no-fix-files` restores OrthoFinder-2 behaviour (MCL orthogroups plus a separate N0.tsv), and `-rmlg` keeps the legacy files.
- Empirical check (core17): `WorkingDirectory/N0.ids.tsv` has 4,798 N0 HOGs with ≥2 genes. **4,774 of them have an identical gene set in `Orthogroups.tsv`** (4,784 OGs). In the first five rows the HOG index even equals the OG index (`N0.HOG0000000` corresponds to `OG0000000`). The pre-HOG MCL clusters (`OG` column of N0.ids, `WorkingDirectory/Alignments_ids`) are an internal intermediate.
- Counts: Orthogroups.tsv has 4,784 / 5,657 / 5,442 OGs. N0 HOGs number 4,798 / 5,671 / 5,456, about 14 more per run (this small difference was not investigated further). N1…N(k) are nested sub-HOGs for internal clades. They are only useful for a "split this OG at clade X" view.

**Default matrix unit: `Orthogroups.tsv`.**
- For OrthoFinder ≥3.1 it is already the phylogenetically refined N0 set.
- It is the only table whose IDs match `MultipleSequenceAlignments/`, `Resolved_Gene_Trees.txt`, `Duplications*.tsv` and `Orthogroup_Sequences/`.
- For OrthoFinder 2.x inputs, prefer `N0.tsv` when it exists, because OF2's `Orthogroups.tsv` holds the coarser MCL groups. N0 IDs (`N0.HOGxxxxxxx`) will then not match the MSA/tree file names, which are keyed by the MCL OG. The `OG` column of N0.tsv provides that link, but the MSA contains extra sequences. Flag this case.
- The importer should therefore detect the version from `Log.txt` (`Started OrthoFinder version …`) and accept both layouts, which matches the existing decision in `docs/design-decisions.md`.

---

## 4. Per-OG metrics: OrthoFinder-native metrics against anvi'o homogeneity

### 4.1 What OrthoFinder already computes

| Metric | Source | Coverage | Notes |
|---|---|---|---|
| genome count, gene count, max copies, copies/genome, single-copy core flag | `Orthogroups.tsv` / GeneCount | all OGs | free |
| duplications (all, ≥50% support) | `Duplications_per_Orthogroup.tsv` | OGs with ≥4 genes (gene trees). Missing means 0 or no tree | Spearman ρ = 0.90 with copies/genome, so mostly redundant |
| duplication location on the species tree, support, terminal vs non-terminal | `Duplications.tsv` | OGs with a tree | detail view, and "expanded in clade X" filters |
| gene tree (branch lengths) | `Resolved_Gene_Trees.txt` | OGs with ≥4 genes, **50% of OGs** (2,410/4,784) | gives divergence metrics (below) |
| orthologue relations (1:1, 1:m, m:m) | `Orthologues/`, `OrthologuesStats_*` | all | per-gene and per-genome-pair. Too large for per-OG tracks at scale |
| putative xenologs (HGT) | `Putative_Xenologs/` | sparse: 2 OGs (core17), 6 OGs (assign24) | flag only |
| phylogenetically misplaced genes | `Phylogenetically_Misplaced_Genes/` | 2–12 genes per run | removed genes, QC only |

OrthoFinder computes **no** alignment quality, gap or identity metric. The anvi'o-style numbers therefore have to be derived from the MSAs or sequences.

### 4.2 Prototype metrics

The scripts were throwaway, in `/tmp/opencode/of_parse.py`, `of_metrics.py`, `of_counts.py`, `of_misc.py` and `of_aai.py`. They ran with `uv run --project . python …`. A full pass over core17 (4,784 OGs: parse every MSA, sequence file and tree, compute everything) takes **2.4 s** single-threaded in Python + NumPy. Adding the DIAMOND within-OG join (1.34M hits, polars) takes another 0.5 s.

Definitions:
- `geo`: anvi'o geometric homogeneity, closed form from per-gene and per-column gap counts (`docs/anvio-internals.md` §1), with all genes pooled as in anvi'o. All-gap columns were **not** removed here, so this matches what anvi'o would report on these files.
- `func`: anvi'o functional homogeneity, closed form from per-column residue counts with the anvi'o conserved-group table (§2).
- `msa_ident`: Σ_col Σ_x C(c_x,2) / Σ_col C(R,2), the fraction of identical residue pairs over non-gap pairs. This is `func` without the conserved-substitution credit.
- `gap_frac`: fraction of gap cells in the MSA.
- `len_cv`: coefficient of variation of unaligned member lengths. It needs no MSA.
- `tree_mpd`: mean pairwise patristic distance in the resolved gene tree. Closed form Σ_e b_e·k_e(n−k_e)/C(n,2), O(n).
- `tree_len`: total gene-tree length. `dups_50` comes from `Duplications_per_Orthogroup.tsv`.
- `diamond_pident`: mean DIAMOND %identity over within-OG gene pairs (core17 all-vs-all only).

Distributions (core17):

| metric | n | min | p25 | median | p75 | max |
|---|---|---|---|---|---|---|
| geo | 4,670 | 0.441 | 0.797 | 0.906 | 0.979 | 1.000 |
| func | 4,670 | 0.239 | 0.507 | 0.613 | 0.757 | 1.000 |
| msa_ident | 4,670 | 0.116 | 0.390 | 0.511 | 0.685 | 1.000 |
| gap_frac | 4,670 | 0.000 | 0.012 | 0.065 | 0.175 | 0.895 |
| len_cv | 4,784 | 0.000 | 0.006 | 0.024 | 0.072 | 1.472 |
| tree_mpd | 2,410 | 0.000 | 1.044 | 1.459 | 1.887 | 4.860 |
| dups_50 | 2,410 | 0 | 0 | 0 | 1 | 91 |
| diamond_pident | 4,780 | 20.7 | 42.9 | 53.2 | 69.9 | 100 |

assign24 and denovo24 have very similar distributions (median geo 0.898 / 0.894, func 0.601 / 0.598). 797 core17 MSAs (17%) have geo = 1, i.e. gap-free. 2,260 MSAs (48%) have only 2–3 sequences.

Spearman correlations, core17, OGs with a gene tree (n = 2,409):

|  | geo | func | msa_ident | gap_frac | len_cv | tree_mpd | tree_len/branch | dups_50 | copies/genome | n_genomes | diamond_pident |
|---|---|---|---|---|---|---|---|---|---|---|---|
| geo | 1 | .54 | .53 | **−.96** | **−.80** | −.53 | −.36 | −.38 | −.39 | −.01 | .39 |
| func | | 1 | **.99** | −.49 | −.45 | **−.93** | −.74 | −.34 | −.35 | −.01 | **.86** |
| msa_ident | | | 1 | −.49 | −.44 | **−.93** | −.74 | −.35 | −.36 | −.02 | **.89** |
| gap_frac | | | | 1 | .72 | .49 | .27 | .37 | .37 | .08 | −.36 |
| len_cv | | | | | 1 | .46 | .30 | .36 | .36 | .02 | −.31 |
| tree_mpd | | | | | | 1 | .71 | .29 | .29 | .12 | **−.83** |
| dups_50 | | | | | | | | 1 | **.90** | −.07 | −.29 |

All OGs with an MSA (n = 4,669, including 2–3-gene OGs): geo is −0.98 with gap_frac and −0.91 with len_cv. func is 0.99 with msa_ident and 0.94 with diamond_pident.

Genome-level: one-way best-hit DIAMOND "AAI" between genome pairs (length-weighted %id, alignments ≥50 aa) has a median of 45.6 (range 35.2–100). Its Spearman correlation with species-tree patristic distance is **−0.955**, so the species tree's branch lengths already carry that information.

### 4.3 Reading

The metrics fall on **three nearly independent axes**:
1. **Structure / completeness**: geo ≈ 1 − gap_frac ≈ −len_cv (|ρ| 0.80–0.98). Geometric homogeneity is essentially a monotone transform of the gap fraction. For flagging fragments and truncated gene calls, `len_cv` gives most of the signal **without any MSA**.
2. **Sequence divergence**: func ≈ msa_ident ≈ −tree_mpd ≈ diamond_pident (|ρ| 0.83–0.99). Functional homogeneity is almost exactly pairwise identity (ρ = 0.99). The conserved-group credit adds next to nothing at the ranking level. OrthoFinder's gene-tree MPD gives the same ordering for free, but only for OGs with ≥4 genes.
3. **Paralogy**: copies/genome ≈ dups_50 (ρ = 0.90). This is available from the matrix alone. Duplication events add *where* on the species tree a duplication happened and with what support, which matters for the detail view and clade filters. They add little as a scalar track.

Cross-axis correlations are only 0.3–0.55, so one track per axis is the right granularity.

---

## 5. Reusable outputs for computation

**Are MSAs trimmed?** Yes, partially. `tools/trim.py` runs by default (`-z` disables it). It trims only pre-HOG alignments with **>500 columns and >10 sequences** (`msa_min_seq`). It drops columns with >90% gaps while keeping at least 500 columns and at least 75% of the residues. Trimming happens **before** the HOG pruning, so small final MSAs can be trimmed too. Measured by comparing residues in the MSA with unaligned sequence length:

| | core17 | assign24 | denovo24 |
|---|---|---|---|
| MSAs with ≥1 trimmed sequence | 509 / 4,670 (10.9%) | 639 / 5,544 | 650 / 5,316 |
| …of which final n > 10 | 374 | 489 | 513 |
| sequences shortened | 4,603 / 36,280 (12.7%) | 7,394 / 51,610 | 7,792 / 51,779 |

Gap-profile heuristics cannot detect trimming reliably: trimmed alignments still contain gappy columns, and pruned ones contain all-gap columns. **Exact detection** compares each aligned sequence's residue count with its unaligned length (from `Orthogroup_Sequences/` or the input FASTA). This is cheap and should be the importer's per-MSA `trimmed` flag. The warning in `design-decisions.md` should refer to this flag rather than a dataset-wide guess.

**Which OGs have MSAs?** All OGs with ≥2 genes except a few 2–3-gene OGs (114 / 113 / 126). No upper size limit was observed. The largest OG (162 genes in core17, 306 in denovo24) has an MSA and a tree. Gene trees exist only for OGs with ≥4 genes (`min_seq=4`).

**DIAMOND hits** (`WorkingDirectory/Blast*_*.txt.gz`): yes, they are usable for AAI-like similarity, both per OG (ρ = 0.95 with MSA identity) and per genome pair (ρ = −0.955 with the species tree). But:
- they only exist for de novo runs (all-vs-all, N² files);
- `--assign` runs only have hits against OG profiles;
- they are internal-ID based and sit in `WorkingDirectory/`.

At 10k genomes they would not exist. **Not worth building on.**

**ID mappings**: `SpeciesIDs.txt` (index to file name) and `SequenceIDs.txt` (internal ID to original ID plus the full description line) are needed only for WorkingDirectory files. The description text in `SequenceIDs.txt` (e.g. `IS30 family transposase, partial [Lactobacillus parafarraginis]`) is a free per-gene annotation. The same text is in the input FASTA headers, which is the better source.

**Gene lengths** for `len_cv`: `Orthogroup_Sequences/` works, but it is one file per OG including singletons. At scale, read lengths from the input proteome FASTAs once and join them through `Orthogroups.tsv`.

---

## 6. Recommendation

### 6.1 Minimal MVP ingest

1. `Orthogroups/Orthogroups.tsv` builds the sparse OG × genome copy-number matrix and gene membership. For OrthoFinder 2.x, prefer `Phylogenetic_Hierarchical_Orthogroups/N0.tsv` if present.
2. `Species_Tree/SpeciesTree_rooted_node_labels.txt` supplies the genome tree, with node labels that match `Duplications.tsv` and the HOG files. Fall back to `SpeciesTree_rooted.txt`, then to the existing UPGMA fallback.
3. `Log.txt` supplies provenance: OrthoFinder version, `--assign` detection, MSA and tree programs.

Everything needed for MVP tracks and ordering derives from (1) and (2): genome count, copy number, single-copy core flag, frequency, MRCA and clade purity.

### 6.2 Optional extras

| Purpose | File(s) | Note |
|---|---|---|
| singleton genes | `Orthogroups_UnassignedGenes.tsv` (or `Orthogroups.txt`) | off by default. The wide format is wasteful at scale, so stream it |
| genome metadata panel | `Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv` | transpose on import |
| structure and divergence tracks | `MultipleSequenceAlignments/` (+ input FASTA or `Orthogroup_Sequences/` for lengths and trim detection) | one streaming pass |
| paralogy-origin track and detail | `Duplications_per_Orthogroup.tsv`, `Duplications.tsv` | per-OG scalar, then per-event detail on demand |
| detail viewer: gene tree | `Resolved_Gene_Trees.txt` | build a byte-offset index at import |
| detail viewer: MSA | `MultipleSequenceAlignments/OGxxxxxxx.fa` | read lazily, strip whitespace, hide all-gap columns |
| HGT/QC flags | `Putative_Xenologs/`, `Phylogenetically_Misplaced_Genes/` | tiny |
| sub-clade split of an OG | `Phylogenetic_Hierarchical_Orthogroups/N*.tsv` | later. Quadratic in size at 10k genomes |
| **skip** | `Orthologues/`, `OrthologuesStats_*`, `Single_Copy_Orthologue_Sequences/`, `WorkingDirectory/` | quadratic, redundant or internal |

### 6.3 Which per-OG metrics

Value for cost, best first:
1. **copies/genome, n_genomes, SCO flag** come from the matrix and cost nothing. They cover the paralogy axis, together with `dups_50` when present.
2. **`len_cv`** needs only gene lengths, not MSAs. It covers about 80–90% of the structural signal that geometric homogeneity carries, is immune to OrthoFinder's trimming, and works for every OG, including those without an MSA.
3. **One MSA pass** gives `geo` (computed after dropping all-gap columns, with a `trimmed` flag per OG) and `msa_ident`/`func` (both are column-count reductions, so cost is the same). This covers 2–3-gene OGs too. On this dataset it took about 2 s for 4.7k MSAs. At 10k genomes the cost is dominated by reading tens of GB of FASTA, so parallelise per file.
4. `tree_mpd` from the gene trees is a fallback for the divergence axis when MSAs are not supplied. It only covers OGs with ≥4 genes.

### 6.4 Keep, replace or drop anvi'o homogeneity?

**Keep both indices, but present them as two of three axes, and stop treating them as special.**
- **Geometric homogeneity: keep.** It is the familiar name for the structure axis and is O(n·L) in closed form. However, on OrthoFinder MSAs it is mostly a gap-fraction transform, and it is biased in two specific ways: trimming inflates it (11% of MSAs), and pruning adds all-gap columns that deflate it (13% of MSAs). Compute it on cleaned MSAs, flag trimmed OGs, and offer `len_cv` as the MSA-free alternative that always works.
- **Functional homogeneity: keep, as the divergence axis.** It is practically identical to pairwise identity (ρ = 0.99), and OrthoFinder offers no per-OG equivalent for small OGs. Gene-tree MPD agrees (ρ = −0.93) but only covers half of the OGs, and DIAMOND identity is absent from `--assign` runs. Call the track "identity / functional homogeneity" and do not spend effort on anvi'o-exact substitution groups.
- **Add from OrthoFinder:** duplications per OG (with species-tree location in the detail view) and xenolog/misplaced flags. These carry information that anvi'o does not have: where on the tree a paralog arose and with what support. They are cheap to read.
- **Drop:** DIAMOND-based AAI (it does not exist at scale, and the species tree already captures it) and per-pair orthologue tables (quadratic).
