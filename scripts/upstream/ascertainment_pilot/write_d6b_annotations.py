#!/usr/bin/env python2
"""Write chr1-22 binary annotation files for the D6B PI x brain-cell analysis.

The script deliberately separates gene-set reconstruction from annotation writing.
Gene sets are read from the audited D6B reconstruction directory.  For each set,
gene coordinates are expanded by 100 kb, overlapping intervals are merged, and
the intervals are projected onto the SNP order in the corresponding 1000 Genomes
EUR BIM file.  Regression SNPs are taken from the official HM3/no-MHC LDSC
weight files so that the resulting annotation rows can be aligned with the
subsetted PLINK files used by LDSC.
"""

from __future__ import print_function

import argparse
import csv
import gzip
import os
import sys


CELL_ORDER = [
    'ASC1', 'ASC2', 'END', 'exCA1', 'exCA3', 'exDG', 'exPFC1',
    'exPFC2', 'GABA1', 'GABA2', 'MG', 'NSC', 'ODC', 'OPC'
]
ANNOTATION_ORDER = ['PI'] + CELL_ORDER + ['PI_x_' + x for x in CELL_ORDER]


def open_text(path, mode='r'):
    if path.endswith('.gz'):
        return gzip.open(path, mode)
    return open(path, mode)


def mkdir(path):
    if not os.path.isdir(path):
        os.makedirs(path)


def read_gene_sets(gene_set_dir):
    result = {}
    for name in ANNOTATION_ORDER:
        path = os.path.join(gene_set_dir, name + '.genes.txt')
        if not os.path.isfile(path):
            raise IOError('missing gene-set file: %s' % path)
        genes = set()
        with open(path, 'r') as handle:
            for line in handle:
                symbol = line.strip()
                if symbol:
                    genes.add(symbol)
        if not genes:
            raise ValueError('empty gene-set file: %s' % path)
        result[name] = genes
    return result


def read_gene_locations(loc_path):
    locations = {}
    with open(loc_path, 'r') as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            fields = line.split()
            if len(fields) < 6:
                continue
            chrom = fields[1].replace('chr', '')
            if chrom not in [str(i) for i in range(1, 23)]:
                continue
            try:
                start = int(fields[2])
                end = int(fields[3])
            except ValueError:
                continue
            symbol = fields[5]
            locations.setdefault(symbol, []).append((chrom, start, end))
    return locations


def merge_intervals(intervals):
    if not intervals:
        return []
    ordered = sorted(intervals)
    merged = [list(ordered[0])]
    for chrom, start, end in ordered[1:]:
        last = merged[-1]
        if start <= last[2] + 1:
            if end > last[2]:
                last[2] = end
        else:
            merged.append([chrom, start, end])
    return [tuple(x) for x in merged]


def build_intervals(gene_sets, locations, window_bp, interval_path):
    intervals_by_annotation = {}
    with open(interval_path, 'w') as out:
        out.write('annotation\tchrom\tstart\tend\n')
        for name in ANNOTATION_ORDER:
            raw_intervals = []
            mapped_genes = set()
            for symbol in gene_sets[name]:
                for chrom, start, end in locations.get(symbol, []):
                    mapped_genes.add(symbol)
                    raw_intervals.append((chrom, max(1, start - window_bp), end + window_bp))
            merged = {}
            for chrom in [str(i) for i in range(1, 23)]:
                merged[chrom] = merge_intervals([x for x in raw_intervals if x[0] == chrom])
                for chrom2, start, end in merged[chrom]:
                    out.write('%s\t%s\t%d\t%d\n' % (name, chrom2, start, end))
            intervals_by_annotation[name] = merged
    return intervals_by_annotation


def read_regression_snps(weights_path):
    snps = set()
    with open_text(weights_path, 'r') as handle:
        first = True
        for raw in handle:
            line = raw.decode('utf-8') if not isinstance(raw, str) else raw
            line = line.strip()
            if not line:
                continue
            fields = line.split()
            if first:
                first = False
                continue
            if len(fields) >= 2:
                snps.add(fields[1])
    return snps


def read_bim(bim_path, regression_snps):
    records = []
    with open(bim_path, 'r') as handle:
        for raw in handle:
            fields = raw.strip().split()
            if len(fields) < 4:
                continue
            snp = fields[1]
            if snp not in regression_snps:
                continue
            chrom = fields[0].replace('chr', '')
            try:
                bp = int(float(fields[3]))
            except ValueError:
                continue
            cm = fields[2]
            records.append((chrom, bp, snp, cm))
    return records


def annotate_records(records, intervals_by_annotation):
    # A pointer per annotation keeps the projection linear in SNP and interval
    # count.  Intervals are sorted and merged within each chromosome.
    pointers = dict((name, 0) for name in ANNOTATION_ORDER)
    rows = []
    for chrom, bp, snp, cm in records:
        values = []
        for name in ANNOTATION_ORDER:
            intervals = intervals_by_annotation[name].get(chrom, [])
            idx = pointers[name]
            while idx < len(intervals) and intervals[idx][2] < bp:
                idx += 1
            pointers[name] = idx
            hit = 0
            if idx < len(intervals):
                if intervals[idx][1] <= bp <= intervals[idx][2]:
                    hit = 1
            values.append(str(hit))
        rows.append((chrom, bp, snp, cm, values))
    return rows


def write_annotation(path, rows):
    with gzip.open(path, 'wb') as out:
        header = ['CHR', 'BP', 'SNP', 'CM'] + ANNOTATION_ORDER
        out.write(('\t'.join(header) + '\n').encode('utf-8'))
        for chrom, bp, snp, cm, values in rows:
            fields = [chrom, str(bp), snp, cm] + values
            out.write(('\t'.join(fields) + '\n').encode('utf-8'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gene-set-dir', required=True)
    parser.add_argument('--loc', required=True)
    parser.add_argument('--bim-dir', required=True)
    parser.add_argument('--weights-dir', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--window-bp', type=int, default=100000)
    args = parser.parse_args()

    annotation_dir = os.path.join(args.out, 'annotations')
    regression_dir = os.path.join(args.out, 'regression_snps')
    mkdir(args.out)
    mkdir(annotation_dir)
    mkdir(regression_dir)

    gene_sets = read_gene_sets(args.gene_set_dir)
    locations = read_gene_locations(args.loc)
    interval_path = os.path.join(args.out, 'annotation_intervals.tsv')
    intervals = build_intervals(gene_sets, locations, args.window_bp, interval_path)

    inventory_path = os.path.join(args.out, 'annotation_inventory.tsv')
    with open(inventory_path, 'w') as inventory:
        inventory.write('chrom\tn_snp\tannotated_sum\tannotation\tn_intervals\n')
        for chrom_num in range(1, 23):
            chrom = str(chrom_num)
            weights = os.path.join(args.weights_dir, 'weights.hm3_noMHC.%s.l2.ldscore.gz' % chrom)
            bim = os.path.join(args.bim_dir, '1000G.EUR.QC.%s.bim' % chrom)
            if not os.path.isfile(weights):
                raise IOError('missing weights file: %s' % weights)
            if not os.path.isfile(bim):
                raise IOError('missing BIM file: %s' % bim)
            regression_snps = read_regression_snps(weights)
            regression_file = os.path.join(regression_dir, 'regression_snps_chr%s.txt' % chrom)
            records = read_bim(bim, regression_snps)
            with open(regression_file, 'w') as out:
                for record in records:
                    out.write(record[2] + '\n')
            rows = annotate_records(records, intervals)
            annot_path = os.path.join(annotation_dir, 'd6b.%s.annot.gz' % chrom)
            write_annotation(annot_path, rows)
            for index, name in enumerate(ANNOTATION_ORDER):
                total = sum(int(row[4][index]) for row in rows)
                n_intervals = len(intervals[name].get(chrom, []))
                inventory.write('%s\t%d\t%d\t%s\t%d\n' %
                                (chrom, len(rows), total, name, n_intervals))
            print('chr%s: %d regression SNPs' % (chrom, len(records)))

    audit = os.path.join(args.out, 'D6B_ANNOTATION_WRITER_AUDIT.txt')
    with open(audit, 'w') as out:
        out.write('analysis=D6B published PI x human brain-cell annotation writer\n')
        out.write('annotation_count=%d\n' % len(ANNOTATION_ORDER))
        out.write('annotation_definition=PI + 14 cell sets + 14 PI intersections\n')
        out.write('window_bp=%d\n' % args.window_bp)
        out.write('gene_set_dir=%s\n' % args.gene_set_dir)
        out.write('coordinate_reference=%s\n' % args.loc)
        out.write('regression_reference=%s\n' % args.weights_dir)
        out.write('snp_reference=%s\n' % args.bim_dir)
        out.write('note=Gene sets were reconstructed and audited separately; annotation projection uses merged gene windows.\n')


if __name__ == '__main__':
    main()
