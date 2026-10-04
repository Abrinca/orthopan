import {describe, expect, test} from 'bun:test';
import type {GeneTreeArrays, Genome} from '../src/api';
import {categoryColor} from '../src/annotations';
import {numericColor} from '../src/metadata';
import {layoutGeneTree, leafColor, fitTreeScale, hideLeafLabels, treeNodeX, MIN_TREE_ROW, HIDE_LABEL_ROW, MISSING_COLOR} from '../src/genetree';

const arrays: GeneTreeArrays = {
  orthogroup: 'OG1',
  // valid pre-order: 0, 1, 2(leaf), 3(leaf), 4(leaf)
  parent: [-1, 0, 1, 1, 0],
  length: [0, 0.1, 0.05, 0.05, 0.2],
  label: ['n0', 'n1', 'l1', 'l2', 'n2'],
  gene: [null, null, 'g1', 'g2', null],
  genome_row: [-1, -1, 0, 1, -1],
  duplication: [false, true, false, false, false],
  support: [null, 0.9, null, null, null],
};

describe('gene-tree array layout', () => {
  test('children, depth and leaf ranges', () => {
    const l = layoutGeneTree(arrays);
    expect(l.n).toBe(5);
    expect(l.children[0]).toEqual([1, 4]);
    expect(l.children[1]).toEqual([2, 3]);
    expect(l.nLeaves).toBe(3);
    expect(l.leafIndex[2]).toBe(0);
    expect(l.leafIndex[3]).toBe(1);
    expect(l.leafIndex[4]).toBe(2);
    expect(l.leafIndex[0]).toBe(-1);
    expect(l.lo[1]).toBe(0);
    expect(l.hi[1]).toBe(1);
    expect(l.lo[0]).toBe(0);
    expect(l.hi[0]).toBe(2);
    expect(l.depth[1]).toBeCloseTo(0.1);
    expect(l.depth[2]).toBeCloseTo(0.15);
    expect(l.maxDepth).toBeCloseTo(0.2);
  });

  test('edge depth from the root and max edges', () => {
    const l = layoutGeneTree(arrays);
    expect(l.edges[0]).toBe(0);
    expect(l.edges[1]).toBe(1);
    expect(l.edges[2]).toBe(2);
    expect(l.edges[4]).toBe(1);
    expect(l.maxEdges).toBe(2);
  });
});

describe('gene-tree view helpers', () => {
  test('fit-to-height keeps every leaf visible, with a minimum row height', () => {
    // 13 leaves in 700 px: fits at full 15 px rows
    expect(fitTreeScale(13, 700)).toBe(1);
    // huge tree: clamp so rows are at least MIN_TREE_ROW px, then scroll
    const s = fitTreeScale(1000, 700);
    expect(15 * s).toBeGreaterThanOrEqual(MIN_TREE_ROW - 1e-9);
  });

  test('leaf labels hide below the threshold', () => {
    expect(hideLeafLabels(HIDE_LABEL_ROW - 0.01)).toBe(true);
    expect(hideLeafLabels(HIDE_LABEL_ROW)).toBe(false);
  });

  test('cladogram aligns all tips at the right edge', () => {
    const width = 300;
    const pad = 10;
    const tip = treeNodeX('cladogram', 0.1, 0.5, 1, 4, true, width, pad);
    const deepTip = treeNodeX('cladogram', 0.5, 0.5, 4, 4, true, width, pad);
    expect(tip).toBe(width - pad);
    expect(deepTip).toBe(width - pad);
  });

  test('phylogram positions by branch length', () => {
    expect(treeNodeX('phylogram', 0.25, 0.5, 0, 0, false, 300, 10)).toBeCloseTo(10 + 0.5 * 280);
    expect(treeNodeX('phylogram', 0.5, 0.5, 0, 0, false, 300, 10)).toBeCloseTo(290);
  });
});

describe('leaf colour', () => {
  const cat = {name: 'clade', kind: 'categorical' as const, categories: ['A', 'B']};
  const num = {name: 'gc', kind: 'numeric' as const, min: 0, max: 1};
  const g = (metadata: Record<string, unknown>): Genome => ({row: 0, genome_idx: 0, name: 'g', n_orthogroups: 0, n_genes: 0, metadata}) as Genome;

  test('categorical uses the category index', () => {
    expect(leafColor(g({clade: 'A'}), cat)).toBe(categoryColor(0));
    expect(leafColor(g({clade: 'B'}), cat)).toBe(categoryColor(1));
  });

  test('numeric uses the gradient', () => {
    expect(leafColor(g({gc: 0.5}), num)).toBe(numericColor(0.5));
  });

  test('missing value and no field', () => {
    expect(leafColor(g({clade: ''}), cat)).toBe(MISSING_COLOR);
    const a = leafColor(g({}), undefined);
    const b = leafColor(g({}), undefined);
    expect(a).toBe(b); // stable identity colour
    expect(a).toMatch(/^#[0-9a-f]{6}$/);
  });
});
