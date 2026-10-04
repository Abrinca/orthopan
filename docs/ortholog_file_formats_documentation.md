# Ortholog File Formats Specification

This document details the file formats used for ortholog group imports, highlighting differences between versions and layouts, edge cases (headers, compression), and required sanitization steps.

---

## 1. Overview of Formats

Depending on the pipeline version and settings (e.g., OrthoFinder v2 vs. v3, or `Orthogroups.tsv` vs. `Orthogroups.txt` / HOG outputs), input files generally arrive in one of two main structural layouts:

1. **Per-Sample Matrix (Tab-Delimited Multi-Column)**
2. **Flat Group-to-Gene Mapping (Two-Column TSV)**

Files can be supplied either uncompressed (`.tsv`, `.txt`) or gzip-compressed (`.tsv.gz`, `.txt.gz`).

---

## 2. Format Specifications

### Format A: Flat Two-Column Format (OrthoFinder v2 / HOGs)

In this format, each row maps a hierarchical orthogroup (HOG) or standard orthogroup identifier directly to a comma-separated list of all member genes across all samples.

#### Example (Raw Preview)
```text
N0.HOG0000000	gnl|extdb|FAM23217-i1-1.1_001153, gnl|extdb|FAM23217-i1-1.1_001100, gnl|extdb|FAM24227-i1-1.1_001852
N0.HOG0000001	gnl|extdb|FAM23217-i1-1.1_001457, gnl|extdb|FAM24235-i1-1.1_001685, gnl|extdb|FAM3257-i1-1.1_002790
N0.HOG0000016	gnl|extdb|FAM18356-i1-1.1_001324
```

#### Structural Details
* **Columns (2):**
  1. `HOG` / `Orthogroup`: Identifier (e.g., `N0.HOG0000000`, `OG0000000`).
  2. `Genes`: Comma-separated list of member gene identifiers (delimited by `, `).
* **Header:** Typically absent.
* **Sample Inference:** Sample/strain IDs are embedded directly in the gene identifier prefix (e.g., `FAM23217` from `FAM23217-i1-1.1_001153`).

---

### Format B: Per-Sample Matrix Format (Tab-Delimited Table)

In this format, the table layout partitions gene assignments into distinct columns corresponding to each analyzed sample/genome.

#### Example
```tsv
Orthogroup	sample-1	sample-2
OG0000000	sample-1_000342	sample-2_000093, sample-2_000127
OG0000001		sample-2_000310
```

#### Structural Details
* **Column 1:** Orthogroup identifier (e.g., `OG0000000`).
* **Remaining Columns:** One column per sample/isolate name.
* **Cells:** Empty if the sample is missing from the group, or populated with one or more comma-separated gene identifiers (`<gene1>, <gene2>`).
* **Header:** May or may not be present depending on whether the file was generated as test data or directly exported by OrthoFinder.

---

## 3. Parsing & Ingestion Logic

### A. Compression Handling
* Inspect file extension or magic bytes (`1f 8b`):
  * If `.gz`, stream through a gzip decompressor (e.g., `gzip.open()` in Python or `zcat` / `pbgzip`).
  * Otherwise, open as plain text with UTF-8 encoding.

### B. Header Detection
* Read the first line:
  * If the first line starts with `Orthogroup\t` or `HOG\t`, treat row 0 as column headers.
  * If the first line matches pattern `^(N\d+\.HOG\d+|OG\d+)\t`, treat row 0 as data (no header).

### C. Layout Detection
Inspect the number of tab-separated fields in line 1 (or line 2 if line 1 is a header):
* **2 Columns:** Flat format (`[Group_ID, Gene_List]`).
* **> 2 Columns:** Matrix format (`[Group_ID, Sample_1, Sample_2, ...]`).

### D. Gene ID Sanitization (Prefix Trimming)
Gene names from external databases or annotation pipelines often include database prefix tags that must be stripped during ingestion:

* **Prefix Rule:** Strip common prefixes matching `gnl|<source>|` or general regex `^gnl\|[^|]+\|`.
* **Example Transformations:**
  * `gnl|extdb|FAM23217-i1-1.1_001153` $\rightarrow$ `FAM23217-i1-1.1_001153`
  * `gnl|extdb|FAM1079-i1-1.1_000370` $\rightarrow$ `FAM1079-i1-1.1_000370`

---

## 4. Normalized Output Schema

Regardless of input version or format, the parser should convert rows into a standard internal representation:

| Field | Type | Description | Example |
| :--- | :--- | :--- | :--- |
| `orthogroup_id` | String | Group or HOG identifier | `N0.HOG0000000` / `OG0000000` |
| `sample_id` | String | Extracted sample/genome identifier | `FAM23217` / `sample-1` |
| `gene_id` | String | Cleaned gene/transcript identifier | `FAM23217-i1-1.1_001153` |


---


# Orthogroup Functional Annotations (`orthogroups.annotations.tsv`)

This file contains consensus functional annotations for hierarchical orthogroups (HOGs) / orthogroups, derived from majority voting across member gene descriptions.

## File Format Specification

- **Delimiter:** Tab-separated (`\t`)
- **Header:** None (raw key-value mappings)
- **Compression:** Plain text (`.tsv`) or gzipped (`.tsv.gz`)

## Structure

A 2-column table mapping each orthogroup identifier to its consensus functional description:

| Column | Field Name | Description | Example |
| :--- | :--- | :--- | :--- |
| **1** | `Orthogroup_ID` | Orthogroup or HOG identifier (key) | `N0.HOG0000000` |
| **2** | `Description` | Putative biochemical function / consensus description | `amino acid ABC transporter ATP-binding protein` |

## Example Content

```tsv
N0.HOG0000000 ABC transporter
N0.HOG0000001 methionine dehydrogenase
```
