import {describe, expect, test} from 'bun:test';
import type {MemberGene} from '../src/api';
import {memberSlice} from '../src/members';

const gene = (id: string): MemberGene => ({id, raw: `raw_${id}`});
const groups = [
  {genome: 'g0', genes: [gene('a'), gene('b'), gene('c')]},
  {genome: 'g1', genes: [gene('d'), gene('e')]},
  {genome: 'g2', genes: [gene('f')]},
];

describe('member list slicing', () => {
  test('total counts every gene, groups preserve order', () => {
    const s = memberSlice(groups, 100);
    expect(s.total).toBe(6);
    expect(s.shown).toBe(6);
    expect(s.groups.map((g) => g.genome)).toEqual(['g0', 'g1', 'g2']);
    expect(s.groups.every((g) => g.hidden === 0)).toBe(true);
  });

  test('limit inside the first group', () => {
    const s = memberSlice(groups, 2);
    expect(s.shown).toBe(2);
    expect(s.groups).toHaveLength(1);
    expect(s.groups[0].genes.map((g) => g.id)).toEqual(['a', 'b']);
    expect(s.groups[0].hidden).toBe(1);
  });

  test('limit spanning groups splits the boundary group', () => {
    const s = memberSlice(groups, 4);
    expect(s.shown).toBe(4);
    expect(s.groups.map((g) => g.genome)).toEqual(['g0', 'g1']);
    expect(s.groups[0].hidden).toBe(0);
    expect(s.groups[1].genes.map((g) => g.id)).toEqual(['d']);
    expect(s.groups[1].hidden).toBe(1);
  });

  test('zero limit renders nothing but still reports the total', () => {
    const s = memberSlice(groups, 0);
    expect(s.groups).toHaveLength(0);
    expect(s.shown).toBe(0);
    expect(s.total).toBe(6);
  });

  test('pinned genome is always included even past the page', () => {
    const s = memberSlice(groups, 1, new Set(['g2']));
    // g0 takes the one paged gene; g2 is force-included in full
    expect(s.groups.map((g) => g.genome)).toEqual(['g0', 'g2']);
    expect(s.groups[0].hidden).toBe(2);
    expect(s.groups[1].hidden).toBe(0);
    expect(s.groups[1].pinned).toBe(true);
    expect(s.shown).toBe(2); // 1 paged + 1 pinned
    expect(s.total).toBe(6);
  });

  test('pinned genome in the middle does not consume the page cap', () => {
    const s = memberSlice(groups, 2, new Set(['g1']));
    expect(s.groups.map((g) => g.genome)).toEqual(['g0', 'g1']);
    expect(s.groups[0]).toEqual({genome: 'g0', genes: [gene('a'), gene('b')], hidden: 1, pinned: false});
    expect(s.groups[1]).toEqual({genome: 'g1', genes: [gene('d'), gene('e')], hidden: 0, pinned: true});
    expect(s.shown).toBe(4);
  });

  test('pinned genome absent from the OG adds nothing', () => {
    const s = memberSlice(groups, 2, new Set(['nope']));
    expect(s.groups.map((g) => g.genome)).toEqual(['g0']);
    expect(s.shown).toBe(2);
  });
});
