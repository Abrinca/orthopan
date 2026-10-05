import {describe, expect, test} from 'bun:test';
import type {GeneTreeArrays, MsaResult} from '../src/api';
import {ALN_ROW_H, locateOg, orderAlignmentRows, stepPosition, ungappedPosition, visibleBand} from '../src/alignment';

const msa: MsaResult = {
  orthogroup: 'OG1',
  source: 'orthofinder',
  length: 3,
  trimmed: null,
  rows: [
    {genome: 'gc', row: 2, gene: 'g3', raw: 'r3', seq: 'ABE'},
    {genome: 'ga', row: 0, gene: 'g1', raw: 'r1', seq: 'ABD'},
    {genome: 'gb', row: 1, gene: 'g2', raw: 'r2', seq: 'ABC'},
  ],
};

const tree: GeneTreeArrays = {
  orthogroup: 'OG1',
  parent: [-1, 0, 0, 0],
  length: [0, 0.1, 0.1, 0.1],
  label: ['n0', 'g1', 'g2', 'g3'],
  gene: [null, 'g1', 'g2', 'g3'],
  genome_row: [-1, 0, 1, 2],
  duplication: [false, false, false, false],
  support: [null, null, null, null],
};

describe('alignment row order', () => {
  test('tree leaf order wins and MSA rows are matched by id', () => {
    const {rows, hasTree} = orderAlignmentRows(msa, tree);
    expect(hasTree).toBe(true);
    expect(rows.map((r) => r.gene)).toEqual(['g1', 'g2', 'g3']);
    expect(rows.map((r) => r.seq)).toEqual(['ABD', 'ABC', 'ABE']);
    expect(rows.map((r) => r.genomeRow)).toEqual([0, 1, 2]);
    expect(rows.map((r) => r.msaIndex)).toEqual([1, 2, 0]);
  });

  test('falls back to MSA order when there is no tree', () => {
    const {rows, hasTree} = orderAlignmentRows(msa, null);
    expect(hasTree).toBe(false);
    expect(rows.map((r) => r.gene)).toEqual(['g3', 'g1', 'g2']);
    expect(rows.map((r) => r.seq)).toEqual(['ABE', 'ABD', 'ABC']);
  });

  test('tree without a matching MSA still yields leaves (seq null)', () => {
    const {rows, hasTree} = orderAlignmentRows(null, tree);
    expect(hasTree).toBe(true);
    expect(rows.map((r) => r.gene)).toEqual(['g1', 'g2', 'g3']);
    expect(rows.every((r) => r.seq === null && r.msaIndex === -1)).toBe(true);
  });
});

describe('visible band', () => {
  test('top of the list includes overscan', () => {
    const b = visibleBand(0, 320, 200, ALN_ROW_H, 2);
    expect(b.startRow).toBe(0);
    expect(b.endRow).toBe(22); // ceil(320/16) + 2
    expect(b.offset).toBe(0);
  });

  test('mid-list band and offset', () => {
    const b = visibleBand(1600, 320, 200, ALN_ROW_H, 2);
    expect(b.startRow).toBe(98); // floor(1600/16) - 2
    expect(b.endRow).toBe(122); // ceil(1920/16) + 2
    expect(b.offset).toBe(98 * 16 - 1600);
  });

  test('band is clamped to the row count', () => {
    const b = visibleBand(3000, 320, 200, ALN_ROW_H, 2);
    expect(b.endRow).toBe(200);
    expect(b.startRow).toBeLessThanOrEqual(b.endRow);
  });
});

describe('prev/next position', () => {
  test('locate an OG in a known window', () => {
    expect(locateOg(['OG1', 'OG2', 'OG3'], 10, 'OG2')).toBe(11);
    expect(locateOg(['OG1', 'OG2', 'OG3'], 10, 'OG9')).toBeNull();
  });

  test('step clamps at both ends', () => {
    expect(stepPosition(11, 1, 20)).toBe(12);
    expect(stepPosition(11, -1, 20)).toBe(10);
    expect(stepPosition(0, -1, 20)).toBeNull();
    expect(stepPosition(19, 1, 20)).toBeNull();
  });
});

describe('ungapped position', () => {
  test('is 1-based for non-gap residues', () => {
    expect(ungappedPosition('A-C-D', 0)).toBe(1); // A
    expect(ungappedPosition('A-C-D', 1)).toBe(-1); // gap
    expect(ungappedPosition('A-C-D', 2)).toBe(2); // C
    expect(ungappedPosition('A-C-D', 4)).toBe(3); // D
  });

  test('out of range is a gap', () => {
    expect(ungappedPosition('ABCD', 9)).toBe(-1);
  });
});
