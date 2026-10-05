#!/usr/bin/env python3
"""Generate deterministic synthetic test metadata and annotations.

Given an orthologs TSV (:mod:`orthopan.io` Format A or B, plain or gzip) and a
Newick genome tree, write seeded synthetic inputs for the test datasets:

- ``genome-metadata.tsv``  ``genome<TAB>clade<TAB>gc``; ``clade`` is a letter
  derived from the tree's leaf order, ``gc`` a float in [0.30, 0.70].
- with ``--annotations``:
  - ``descriptions.tsv``   headerless ``orthogroup<TAB>free text``;
  - ``cog.tsv``            header ``Orthogroup<TAB>cog_category`` (categorical).

A single seeded ``random.Random`` is used throughout, so repeated runs are
byte-identical. It replaces an earlier ad-hoc scratch generator.
"""

import argparse
import random
from pathlib import Path

# No public Newick parser exists; orthopan.importer._parse_newick is the one
# used by the importer and the server, so a synthetic output matches them.
from orthopan.importer import _parse_newick
from orthopan.io import open_text, read_matrix

WORDS = [
    "putative", "ABC", "transporter", "ATP-binding", "protein", "dehydrogenase",
    "kinase", "subunit", "membrane", "permease", "oxidase", "reductase", "synthase",
    "hydrolase", "domain", "family", "conserved", "hypothetical", "ferredoxin",
    "aminotransferase", "phosphatase", "ligase", "methyltransferase", "chaperone",
    "protease", "ribosomal", "flagellar", "secretion", "two-component", "response",
    "regulator", "sensor", "histidine", "transcriptional", "NADH", "iron",
]
COG = list("JKLDOUNMCGEFHIPTQRSVW")


def leaf_order(root):
    """Leaf labels in pre-order (the same order the importer numbers rows)."""
    out, stack = [], [root]
    while stack:
        node = stack.pop()
        if node.children:
            stack.extend(reversed(node.children))
        else:
            out.append(node.label)
    return out


def read_genomes_and_ogs(path):
    """Return ``(genome labels, orthogroup IDs)`` from an orthologs TSV."""
    with open_text(path) as fh:
        lines = (ln.rstrip("\r\n") for ln in fh if ln.strip())
        head = next(lines).split("\t")
        first_col = lambda ln: ln.split("\t", 1)[0]  # noqa: E731
        if head[0] in ("Orthogroup", "HOG"):  # Format B with header
            return head[1:], [first_col(ln) for ln in lines]
        if len(head) == 2:  # Format A: genome from the gene ID
            labels, rows = read_matrix(path)
            rows.close()
            with open_text(path) as fh2:
                ogs = [first_col(ln.rstrip("\r\n")) for ln in fh2 if ln.strip()]
            return list(labels), ogs
        # Headerless Format B: the first line is data, so the importer names the
        # genomes genome_1..n by the data columns (read_matrix handles this).
        labels, rows = read_matrix(path)
        rows.close()
        return list(labels), [head[0], *(first_col(ln) for ln in lines)]


def clade_letters(leaves, n):
    """Map each leaf to one of ``n`` evenly sized clade letters by leaf order."""
    size = max(1, (len(leaves) + n - 1) // n)
    return {leaf: chr(ord("A") + i // size) for i, leaf in enumerate(leaves)}


def write(path, header, rows):
    with open(path, "w") as fh:
        if header:
            fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--orthologs", required=True, type=Path,
                   help="orthologs TSV/TSV.GZ (Format A or B)")
    p.add_argument("--tree", required=True, type=Path, help="Newick genome tree")
    p.add_argument("--out", required=True, type=Path, help="output directory")
    p.add_argument("--annotations", action="store_true",
                   help="also write synthetic descriptions.tsv and cog.tsv")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--clades", type=int, default=8,
                   help="number of clade letters in the genome metadata")
    a = p.parse_args(argv)

    a.out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(a.seed)
    labels, ogs = read_genomes_and_ogs(a.orthologs)
    clade_of = clade_letters(leaf_order(_parse_newick(a.tree.read_text())), a.clades)

    write(a.out / "genome-metadata.tsv", ["genome", "clade", "gc"],
          [[g, clade_of.get(g, "?"), f"{rng.uniform(0.30, 0.70):.4f}"] for g in labels])

    if a.annotations:
        with open(a.out / "descriptions.tsv", "w") as fh:
            for og in ogs:
                phrase = " ".join(rng.choice(WORDS) for _ in range(rng.randint(3, 6)))
                fh.write(f"{og}\t{phrase}\n")
        write(a.out / "cog.tsv", ["Orthogroup", "cog_category"],
              [[og, rng.choice(COG)] for og in ogs])

    assert labels, "no genomes found in the orthologs TSV"
    print(f"{a.out}: {len(labels)} genomes, {len(ogs)} orthogroups"
          + (" (+annotations)" if a.annotations else ""))


if __name__ == "__main__":
    main()
