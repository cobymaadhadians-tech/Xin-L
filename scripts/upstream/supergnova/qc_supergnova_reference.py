#!/usr/bin/env python3
import argparse
import bisect
import csv
import os
import sys


def read_bim(path):
    records = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            chrom, snp, cm, bp, a1, a2 = line.split()
            records.append((int(chrom), snp, float(cm), int(bp), a1, a2))
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True)
    args = parser.parse_args()

    root = os.path.abspath(args.root)
    ref = os.path.join(root, "reference", "1000G_phase3_EUR503_GRCh37")
    partition_path = os.path.join(root, "reference", "ldetect_eur", "eur_all.bed")
    plink_dir = os.path.join(ref, "plink")
    qc_dir = os.path.join(ref, "qc")

    variants = {}
    for chrom in range(1, 23):
        bim = read_bim(os.path.join(plink_dir, f"eur_chr{chrom}.bim"))
        variants[chrom] = [record[3] for record in bim]

    rows = []
    with open(partition_path, encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for block_index, row in enumerate(reader, start=1):
            chrom = int(row["CHR"])
            start = int(row["START"])
            end = int(row["END"])
            positions = variants[chrom]
            count = bisect.bisect_right(positions, end) - bisect.bisect_left(positions, start)
            rows.append((block_index, chrom, start, end, count))

    out_path = os.path.join(qc_dir, "partition_reference_coverage.tsv")
    with open(out_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["BLOCK", "CHR", "START", "END", "N_REFERENCE_SNPS"])
        writer.writerows(rows)

    counts = [row[4] for row in rows]
    counts_sorted = sorted(counts)
    median = counts_sorted[len(counts_sorted) // 2]
    summary = {
        "N_BLOCKS": len(rows),
        "N_BLOCKS_LT120_REFERENCE_SNPS": sum(x < 120 for x in counts),
        "MIN_REFERENCE_SNPS": min(counts),
        "MEDIAN_REFERENCE_SNPS": median,
        "MAX_REFERENCE_SNPS": max(counts),
    }
    with open(os.path.join(qc_dir, "partition_reference_summary.tsv"), "w", encoding="utf-8") as handle:
        handle.write("METRIC\tVALUE\n")
        for key, value in summary.items():
            handle.write(f"{key}\t{value}\n")

    if len(rows) != 1703:
        raise SystemExit(f"Expected 1703 LDetect blocks, found {len(rows)}")

    code = os.path.join(root, "software", "SUPERGNOVA")
    sys.path.insert(0, code)
    import ld.ldscore as ld
    import ld.parse as ps

    prefix = os.path.join(plink_dir, "eur_chr22")
    snps = ps.PlinkBIMFile(prefix + ".bim")
    indivs = ps.PlinkFAMFile(prefix + ".fam")
    geno = ld.PlinkBEDFile(prefix + ".bed", len(indivs.IDList), snps)
    if geno.n != 503 or geno.m != len(snps.IDList):
        raise SystemExit("SUPERGNOVA BED parser dimensions do not match BIM/FAM")

    with open(os.path.join(qc_dir, "supergnova_parser_smoke_test.tsv"), "w", encoding="utf-8") as handle:
        handle.write("CHECK\tVALUE\n")
        handle.write("CHROMOSOME\t22\n")
        handle.write(f"N_INDIVIDUALS\t{geno.n}\n")
        handle.write(f"N_VARIANTS\t{geno.m}\n")
        handle.write("BED_PARSER\tPASS\n")


if __name__ == "__main__":
    main()
