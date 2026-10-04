# OrthoPan: Vision & Architecture Brief

## 1. Executive Summary & Purpose
**OrthoPan** is a scalable, modern, and interactive pangenome exploration tool designed as an alternative to anvi'o's pangenomic visualization workflow. 

While anvi'o introduced breakthrough metrics for evaluating cluster quality—specifically **Geometric Homogeneity** (alignment gap/residue consistency) and **Functional Homogeneity** (biochemical conservation)—its default presentation relies on a circular/radial visualization. In radial plots, genomes placed on the outermost rings receive disproportionately large visual area compared to genomes on the inner rings. Furthermore, anvi'o's backend architecture bottlenecks when handling tens of thousands of genomes.

**OrthoPan solves this by:**
1. Decoupling sequence clustering from visual exploration by natively leveraging **OrthoFinder** output (explicit orthology inference + precomputed phylogenetic MSAs).
2. Implementing anvi'o-style **Geometric & Functional Homogeneity** algorithms directly over OrthoFinder alignments using high-performance vector processing (`polars`?).
3. Providing a **linear, equal-area interactive viewer** where every genome and every orthogroup occupies an unbiased, strictly proportional footprint on screen.

---

## 2. Core Concepts & Terminology

* **Orthogroup (OG):** The set of genes descended from a single gene in the last common ancestor of the species under consideration (computed via OrthoFinder). This serves as the fundamental unit (column) of the pangenome matrix.
* **Presence / Absence Matrix:** A 2D matrix ($N \times M$) of $N$ genomes across $M$ orthogroups, indicating presence, absence, or copy number.
* **Geometric Homogeneity Index ($[0.0, 1.0]$):**
  * Evaluates structural and alignment continuity independently of amino acid identity.
  * **Site-level (Vertical):** Evaluates whether alignment columns across all genomes consistently exhibit gaps or residues, identifying spurious insertions or regional indels.
  * **Gene-level (Horizontal):** Evaluates whether member sequences share comparable lengths and domain spans, immediately flagging truncated gene calls or contig-edge fragments.
* **Functional Homogeneity Index ($[0.0, 1.0]$):**
  * Evaluates biochemical conservation across non-gap residues across aligned sites using classified physicochemical properties or substitution matrices.

---

## 3. Visual Paradigm: The Linear Canvas

Unlike circular layouts, OrthoPan treats all genomes with identical visual weight.

```
                    [ Top Bar: OG Dendrogram / Cluster Frequency ]
                    [ Track 1: Geometric Homogeneity (Bar / Profile) ]
                    [ Track 2: Functional Homogeneity (Bar / Profile) ]
                    [ Track 3: Functional Annotation / COG Category ]
                  ┌─────────────────────────────────────────────────────┐
Phylogenomic Tree │                                                     │
(Genomes / Rows)  │              Presence / Absence Heatmap             │
                  │        (Rows = Genomes, Columns = Orthogroups)      │
                  │                                                     │
                  └─────────────────────────────────────────────────────┘
                    [ Bottom: Genome Metadata (Clades, Phenotypes, QC)  ]
```

### Key UI/UX Advantages:
* **Equal Area Guarantee:** Every cell $(g_i, og_j)$ has a uniform width and height ($w \times h$).
* **Independent Scalability:** Orthogroups (X-axis, easily $20{,}000+$ items) and genomes (Y-axis, $10$ to $1{,}000+$ items) can be scaled and scrolled independently.
* **Standard Dendrogram Coupling:** Preserves canonical phylogenetic tree structures without angular distortion or radial branch-bending.

---

## 4. Technical Foundations & Tech Stack

* **Package & Environment Management:** `uv`
* **Test Suite:** `pytest`
* **Data Processing & Analytics:** `polars` (fast tabular operations, Parquet caching, multi-threaded aggregation)
* **Alignment Parsing & Matrix Operations:** Python / NumPy / Rust extensions for binary bitwise XOR sweeps over MSA matrices
* **Inputs:**
  * `Orthogroups.tsv` / `orthologs.tsv` (OrthoFinder cluster definitions)
  * OrthoFinder MSA outputs (`MultipleSequenceAlignments/*.fa`)
  * Species tree / gene trees (`SpeciesTree_rooted.txt`)
  * Metadata annotations (`orthologs-descriptions.tsv`)

---

## 5. Agent "Grill" Areas: Technical & Architecture Challenges

To finalize the design, the AI agent should interrogate the following open architectural and computational decisions:

1. **Alignment Source & Execution Strategy:**
   * Does OrthoPan mandate that OrthoFinder was run with `-M msa` (using precomputed MAFFT MSAs), or will OrthoPan run on-demand/cached MSAs on orthogroups if only the TSV is supplied?
   * How should OrthoPan handle hierarchical orthogroups (HOGs / `N0.tsv`) vs. broad orthogroups (`Orthogroups.tsv`)?
2. **Homogeneity Computation at Scale:**
   * Anvi’o computes pairwise XOR operations across aligned sequences. How do we vectorize this in NumPy or Polars to ensure sub-second calculation across $30{,}000$ MSAs?
   * How do we handle paralogs (multi-copy genes within a single genome inside an orthogroup) when computing geometric homogeneity?
3. **Frontend Rendering Engine:**
   * A linear matrix of $1{,}000 \text{ genomes} \times 30{,}000 \text{ orthogroups} = 30{,}000{,}000 \text{ cells}$. Standard DOM/SVG implementations will crash.
   * What renderer is best suited: WebGL (PIXI.js, Deck.gl), HTML5 Canvas virtual slicing, or server-side tiled image generation?
4. **Data Interchange & Intermediate Storage:**
   * Rather than an monolithic SQLite database like anvi'o, should OrthoPan serialize processed pangenomes into Apache Arrow / Parquet partitions or an embedded DuckDB database?
5. **Interactive Filtering & Sorting:**
   * What heuristic ordering schemes (Euclidean clustering, Ward's minimum variance, frequency sorting, core-to-accessory gradient) must be computed upfront vs. on-the-fly?