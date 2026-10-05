import {describe, expect, test} from 'bun:test';
import type {MemberGene} from '../src/api';
import type {MemberGroup} from '../src/members';
import {decodeSelection, encodeSelection, rangeSelection, selectionNames, SELECTION_GENE_CAP, singleSelection, summarizeSelection} from '../src/selection';

const gene = (id: string): MemberGene => ({id, raw: `raw_${id}`});
const members: MemberGroup[] = [
  {genome: 'g0', genes: [gene('a'), gene('b')]},
  {genome: 'g1', genes: [gene('c')]},
  {genome: 'g2', genes: [gene('d'), gene('e'), gene('f')]},
];

describe('genome selection', () => {
  test('single genome present', () => {
    const s = summarizeSelection(singleSelection('g2'), members);
    expect(s.absent).toBe(false);
    expect(s.selected).toBe(1);
    expect(s.withGenes).toBe(1);
    expect(s.genes).toBe(3);
    expect(s.groups.map((g) => g.genome)).toEqual(['g2']);
  });

  test('single genome absent', () => {
    const s = summarizeSelection(singleSelection('nope'), members);
    expect(s.absent).toBe(true);
    expect(s.withGenes).toBe(0);
    expect(s.genes).toBe(0);
    expect(s.groups).toEqual([]);
  });

  test('binned range counts genomes with the OG, not all columns', () => {
    const sel = rangeSelection(['g0', 'g1', 'gX', 'g2']);
    expect(sel.binned).toBe(true);
    expect(sel.label).toBe('g0 – g2');
    const s = summarizeSelection(sel, members);
    expect(s.selected).toBe(4);
    expect(s.withGenes).toBe(3); // gX has no genes
    expect(s.genes).toBe(6);
    expect(s.tooMany).toBe(false);
    expect(s.groups.map((g) => g.genome)).toEqual(['g0', 'g1', 'g2']);
  });

  test('binned range over the cap reports tooMany', () => {
    const big: MemberGroup[] = [{genome: 'g0', genes: Array.from({length: SELECTION_GENE_CAP + 1}, (_, i) => gene(`x${i}`))}];
    const s = summarizeSelection(rangeSelection(['g0', 'g1']), big);
    expect(s.genes).toBe(SELECTION_GENE_CAP + 1);
    expect(s.tooMany).toBe(true);
  });

  test('selectionNames returns the highlighted set', () => {
    expect(selectionNames(rangeSelection(['a', 'b']))).toEqual(new Set(['a', 'b']));
    expect(selectionNames(null)).toEqual(new Set());
  });
});

describe('selection hash encoding', () => {
  const names = ['g0', 'g1', 'g2', 'g3', 'g4', 'g5'];

  test('single genome round-trips through genome=', () => {
    const sel = singleSelection('g2');
    expect(encodeSelection(sel)).toBe('&genome=g2');
    expect(decodeSelection('g2', null, names)).toEqual(singleSelection('g2'));
  });

  test('binned range round-trips through rows=', () => {
    const sel = rangeSelection(names.slice(2, 5), [2, 5]);
    expect(encodeSelection(sel)).toBe('&rows=2-5');
    const back = decodeSelection(null, '2-5', names);
    expect(back?.binned).toBe(true);
    expect(back?.rows).toEqual([2, 5]);
    expect(back?.names).toEqual(['g2', 'g3', 'g4']);
    expect(back?.label).toBe('g2 – g4');
  });

  test('rows wins over a stale genome param; invalid rows decode to null', () => {
    expect(decodeSelection('g9', '1-3', names)?.names).toEqual(['g1', 'g2']);
    expect(decodeSelection(null, '4-2', names)).toBeNull();
    expect(decodeSelection(null, '3-99', names)).toBeNull();
    expect(decodeSelection(null, null, names)).toBeNull();
  });
});
