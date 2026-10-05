// Selected-genome state for the Members panel: a matrix click selects an
// orthogroup and a genome row. When the row is binned (several genomes per
// pixel row) the selection is a contiguous range of genome names.

import type {MemberGene} from './api';
import type {MemberGroup} from './members';

/** Genes above this count are not listed for a binned selection. */
export const SELECTION_GENE_CAP = 50;

export interface GenomeSelection {
  /** selected genome names (one for a single row, many for a binned row) */
  names: string[];
  /** true when several genomes share the clicked pixel row */
  binned: boolean;
  /** human label, e.g. "gA – gB" for a range */
  label: string;
  /** genome row range [r0, r1) in tree order, for a binned selection */
  rows?: [number, number];
}

export function singleSelection(name: string): GenomeSelection {
  return {names: [name], binned: false, label: name};
}

export function rangeSelection(names: string[], rows?: [number, number]): GenomeSelection {
  const label = names.length <= 1 ? (names[0] ?? '') : `${names[0]} – ${names[names.length - 1]}`;
  return {names, binned: names.length > 1, label, rows};
}

/** `&genome=<name>` for a single genome, `&rows=r0-r1` for a binned range. */
export function encodeSelection(sel: GenomeSelection | null): string {
  if (!sel) return '';
  if (sel.binned && sel.rows) return `&rows=${sel.rows[0]}-${sel.rows[1]}`;
  if (sel.names.length) return `&genome=${encodeURIComponent(sel.names[0])}`;
  return '';
}

/**
 * Rebuild a selection from the hash. A binned `rows` range is preferred over a
 * single `genome`; `names` are looked up in the current genome (row) order.
 */
export function decodeSelection(
  genome: string | null,
  rows: string | null,
  genomeNames: string[],
): GenomeSelection | null {
  if (rows) {
    const m = /^(\d+)-(\d+)$/.exec(rows);
    if (m) {
      const r0 = Number(m[1]);
      const r1 = Number(m[2]);
      if (r0 >= 0 && r1 > r0 && r1 <= genomeNames.length) {
        return rangeSelection(genomeNames.slice(r0, r1), [r0, r1]);
      }
    }
  }
  if (genome) return singleSelection(genome);
  return null;
}

export interface SelectionGroup {
  genome: string;
  genes: MemberGene[];
}

export interface GenomeSummary {
  label: string;
  binned: boolean;
  /** number of selected genomes */
  selected: number;
  /** selected genomes that have the orthogroup */
  withGenes: number;
  /** total genes across the selected genomes */
  genes: number;
  /** selected groups that have genes, in member order */
  groups: SelectionGroup[];
  /** too many genes to list for a binned selection */
  tooMany: boolean;
  /** single selection whose genome has no genes in this OG */
  absent: boolean;
}

/** Summarise a selection against the OG's genome-grouped members. */
export function summarizeSelection(sel: GenomeSelection, members: MemberGroup[]): GenomeSummary {
  const wanted = new Set(sel.names);
  const groups: SelectionGroup[] = [];
  let genes = 0;
  for (const m of members) {
    if (!wanted.has(m.genome)) continue;
    groups.push({genome: m.genome, genes: m.genes});
    genes += m.genes.length;
  }
  return {
    label: sel.label,
    binned: sel.binned,
    selected: sel.names.length,
    withGenes: groups.length,
    genes,
    groups,
    tooMany: genes > SELECTION_GENE_CAP,
    absent: !sel.binned && groups.length === 0,
  };
}

/** Genome names to highlight (member list and alignment rows). */
export function selectionNames(sel: GenomeSelection | null): Set<string> {
  return new Set(sel?.names ?? []);
}
