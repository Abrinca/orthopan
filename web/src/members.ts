// Member-list windowing: keep the drawer responsive for orthogroups with tens
// of thousands of members by rendering only the first `limit` genes.

import type {MemberGene} from './api';

export interface MemberGroup {
  genome: string;
  genes: MemberGene[];
}

export interface SlicedGroup {
  genome: string;
  genes: MemberGene[];
  /** genes in this genome not rendered in the slice */
  hidden: number;
  /** true when this group was force-included because its genome is selected */
  pinned: boolean;
}

export interface MemberSlice {
  groups: SlicedGroup[];
  /** total genes rendered */
  shown: number;
  /** total genes in all groups */
  total: number;
}

/**
 * Take the first `limit` genes across the genome groups, preserving order.
 * Genomes in `pinned` are always included in full (so a selected genome stays
 * visible even when it falls beyond the page).
 */
export function memberSlice(groups: MemberGroup[], limit: number, pinned: Set<string> = new Set()): MemberSlice {
  let total = 0;
  for (const g of groups) total += g.genes.length;
  const out: SlicedGroup[] = [];
  const cap = Math.max(0, limit);
  let shown = 0;
  let paged = 0;
  for (const g of groups) {
    if (pinned.has(g.genome)) {
      out.push({genome: g.genome, genes: g.genes, hidden: 0, pinned: true});
      shown += g.genes.length;
      continue;
    }
    if (paged >= cap) continue;
    const take = Math.min(g.genes.length, cap - paged);
    out.push({genome: g.genome, genes: g.genes.slice(0, take), hidden: g.genes.length - take, pinned: false});
    paged += take;
    shown += take;
  }
  return {groups: out, shown, total};
}

