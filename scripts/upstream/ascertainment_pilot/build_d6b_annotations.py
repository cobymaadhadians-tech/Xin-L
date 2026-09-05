#!/usr/bin/env python3
"""Reconstruct the published PI × human brain-cell gene sets for D6B.

The source data are the public gnomAD pLI table, the GTEx hippocampus/PFC
DroNc-seq UMI matrix, and the published cell-cluster mapping workbook.
This script deliberately writes an audit trail for every deterministic
choice. It does not create LDSC annotation files; that is handled by the
separate annotation/LD-score step after the gene-set audit passes.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd


CELL_MAP = {
    1: "exPFC1",
    2: "exPFC2",
    3: "exCA1",
    4: "exCA3",
    5: "GABA1",
    6: "GABA2",
    7: "exDG",
    8: "ASC1",
    9: "ASC2",
    10: "ODC",  # ODC1 and ODC2 are the published Oligodendrocytes set.
    11: "ODC",
    12: "OPC",
    13: "MG",
    14: "NSC",
    15: "END",
}

CELL_ORDER = [
    "ASC1", "ASC2", "END", "exCA1", "exCA3", "exDG", "exPFC1",
    "exPFC2", "GABA1", "GABA2", "MG", "NSC", "ODC", "OPC",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--raw-dir", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--loc", required=True)
    p.add_argument("--min-detection-fraction", type=float, default=0.01)
    p.add_argument("--top-k", type=int, default=1600)
    p.add_argument("--window-bp", type=int, default=100_000)
    return p.parse_args()


def read_cluster_map(path: Path) -> dict[str, int]:
    out: dict[str, int] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) < 2:
                continue
            try:
                cluster = int(fields[1])
            except ValueError:
                continue
            out[fields[0]] = cluster
    return out


def read_gnomad_pi(path: Path) -> tuple[set[str], int, int]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = 0
        names: list[str] = []
        for row in reader:
            try:
                pli = float(row.get("pLI", "nan"))
            except ValueError:
                pli = math.nan
            if math.isfinite(pli) and pli > 0.9:
                rows += 1
                gene = (row.get("gene") or "").strip()
                if gene:
                    names.append(gene)
    return set(names), rows, len(names)


def read_gene_coordinates(path: Path) -> dict[str, tuple[str, int, int]]:
    result: dict[str, tuple[str, int, int]] = {}
    with path.open() as handle:
        for line in handle:
            f = line.split()
            if len(f) < 6:
                continue
            chrom = f[1]
            if chrom not in {str(i) for i in range(1, 23)}:
                continue
            try:
                start, end = int(f[2]), int(f[3])
            except ValueError:
                continue
            if start > 0 and end >= start:
                result.setdefault(f[5], (chrom, start, end))
    return result


def matrix_path(raw_dir: Path) -> Path:
    candidates = [
        raw_dir / "GTEx_droncseq_hip_pcf.umi_counts.txt.gz",
        raw_dir / "GTEx_droncseq_hip_pcf" / "GTEx_droncseq_hip_pcf.umi_counts.txt.gz",
        raw_dir / "umi_counts.txt.gz",
    ]
    for p in candidates:
        if p.exists():
            return p
    raise FileNotFoundError("GTEx UMI matrix not found")


def cluster_path(raw_dir: Path) -> Path:
    candidates = [
        raw_dir / "GTEx_droncseq_hip_pcf.clusters.txt.gz",
        raw_dir / "GTEx_droncseq_hip_pcf" / "GTEx_droncseq_hip_pcf.clusters.txt.gz",
        raw_dir / "clusters.txt.gz",
    ]
    for p in candidates:
        if p.exists():
            return p
    raise FileNotFoundError("GTEx cluster file not found")


def main() -> None:
    args = parse_args()
    raw_dir = Path(args.raw_dir)
    out_dir = Path(args.out_dir)
    genes_dir = out_dir / "gene_sets"
    out_dir.mkdir(parents=True, exist_ok=True)
    genes_dir.mkdir(parents=True, exist_ok=True)

    matrix = matrix_path(raw_dir)
    clusters = read_cluster_map(cluster_path(raw_dir))
    pi_genes, pi_rows, pi_names = read_gnomad_pi(
        raw_dir / "gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"
    )
    coords = read_gene_coordinates(Path(args.loc))

    header = pd.read_csv(matrix, sep="\t", compression="gzip", nrows=0)
    cell_names = [str(x) for x in header.columns[1:]]
    cluster_ids = np.array([clusters.get(x, -1) for x in cell_names], dtype=np.int16)
    keep = np.array([x in CELL_MAP for x in cluster_ids], dtype=bool)
    if not keep.all():
        cell_names = [x for x, k in zip(cell_names, keep) if k]
        cluster_ids = cluster_ids[keep]
    if len(cell_names) < 10_000:
        raise RuntimeError(f"Unexpected retained cell count: {len(cell_names)}")

    categories = {name: np.flatnonzero(np.array([CELL_MAP[c] == name for c in cluster_ids]))
                  for name in CELL_ORDER}
    if any(len(v) == 0 for v in categories.values()):
        raise RuntimeError(f"At least one published category has no cells: {categories}")

    # Pass 1: library sizes, detection counts, and stable gene order.
    libsize = np.zeros(len(cell_names), dtype=np.float64)
    detection: list[np.ndarray] = []
    gene_names: list[str] = []
    for chunk in pd.read_csv(matrix, sep="\t", compression="gzip", chunksize=512,
                             index_col=0):
        arr = chunk.to_numpy(dtype=np.float64, copy=False)
        if arr.shape[1] != keep.size:
            raise RuntimeError("Matrix column count does not match cluster map")
        arr = arr[:, keep]
        if arr.shape[1] != len(cell_names):
            raise RuntimeError("Cell columns changed between matrix passes")
        gene_names.extend([str(x) for x in chunk.index])
        libsize += arr.sum(axis=0)
        detection.append((arr > 0).sum(axis=1).astype(np.int32))
    detection_all = np.concatenate(detection)
    if len(gene_names) != len(detection_all):
        raise RuntimeError("Gene and detection vector lengths differ")
    min_detected = int(math.ceil(args.min_detection_fraction * len(cell_names)))
    expressed = detection_all >= min_detected
    expressed_index = np.flatnonzero(expressed)
    if len(expressed_index) < args.top_k:
        raise RuntimeError("Low-expression filter left fewer genes than top-k")

    # Pass 2: mean library-size-normalized expression in each published set.
    means = np.zeros((len(expressed_index), len(CELL_ORDER)), dtype=np.float64)
    write_pos = 0
    gene_pos = 0
    scale = 1_000_000.0 / np.maximum(libsize, 1.0)
    for chunk in pd.read_csv(matrix, sep="\t", compression="gzip", chunksize=512,
                             index_col=0):
        arr = chunk.to_numpy(dtype=np.float64, copy=False)
        if arr.shape[1] != keep.size:
            raise RuntimeError("Matrix column count does not match cluster map")
        arr = arr[:, keep]
        idx = np.flatnonzero(expressed[gene_pos:gene_pos + arr.shape[0]])
        if len(idx):
            norm = arr[idx, :] * scale[None, :]
            for j, name in enumerate(CELL_ORDER):
                means[write_pos:write_pos + len(idx), j] = norm[:, categories[name]].mean(axis=1)
            write_pos += len(idx)
        gene_pos += arr.shape[0]
    if write_pos != len(expressed_index):
        raise RuntimeError("Expression pass did not match filtered gene count")

    expressed_genes = np.array(gene_names, dtype=object)[expressed_index]
    # Specificity is each category's mean normalized expression divided by
    # the equal-weight mean across the 14 published categories.
    specificity = means / np.maximum(means.mean(axis=1, keepdims=True), 1e-12)

    manifest_rows = []
    for j, name in enumerate(CELL_ORDER):
        order = np.argsort(-specificity[:, j], kind="mergesort")
        selected = [str(x) for x in expressed_genes[order[:args.top_k]]]
        selected = list(dict.fromkeys(selected))
        if len(selected) != args.top_k:
            raise RuntimeError(f"Non-unique top genes in {name}: {len(selected)}")
        path = genes_dir / f"{name}.genes.txt"
        path.write_text("\n".join(selected) + "\n")
        mapped = sum(g in coords for g in selected)
        manifest_rows.append({
            "annotation": name,
            "source_cluster_ids": ",".join(str(k) for k, v in CELL_MAP.items() if v == name),
            "n_cells": len(categories[name]),
            "n_genes_after_low_expression_filter": len(expressed_genes),
            "min_detected_cells": min_detected,
            "top_k": args.top_k,
            "selected_genes": len(selected),
            "selected_genes_with_hg19_coordinates": mapped,
            "gene_set_file": str(path),
        })

    pi_path = genes_dir / "PI.genes.txt"
    pi_path.write_text("\n".join(sorted(pi_genes)) + "\n")
    manifest_rows.insert(0, {
        "annotation": "PI",
        "source_cluster_ids": "",
        "n_cells": "",
        "n_genes_after_low_expression_filter": "",
        "min_detected_cells": "",
        "top_k": "",
        "selected_genes": len(pi_genes),
        "selected_genes_with_hg19_coordinates": sum(g in coords for g in pi_genes),
        "gene_set_file": str(pi_path),
    })
    for row in list(manifest_rows[1:]):
        cell_path = Path(row["gene_set_file"])
        cell_genes = {x.strip() for x in cell_path.read_text().splitlines() if x.strip()}
        inter = sorted(cell_genes & pi_genes)
        inter_path = genes_dir / f"PI_x_{row['annotation']}.genes.txt"
        inter_path.write_text("\n".join(inter) + ("\n" if inter else ""))
        manifest_rows.append({
            "annotation": f"PI_x_{row['annotation']}",
            "source_cluster_ids": row["source_cluster_ids"],
            "n_cells": row["n_cells"],
            "n_genes_after_low_expression_filter": row["n_genes_after_low_expression_filter"],
            "min_detected_cells": row["min_detected_cells"],
            "top_k": row["top_k"],
            "selected_genes": len(inter),
            "selected_genes_with_hg19_coordinates": sum(g in coords for g in inter),
            "gene_set_file": str(inter_path),
        })

    fieldnames = list(manifest_rows[0].keys())
    with (out_dir / "gene_set_manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(manifest_rows)

    audit = [
        "analysis=D6B published PI x human brain-cell annotation reconstruction",
        "protocol=Grotzinger_Nature_Genetics_2022",
        f"matrix={matrix}",
        f"n_cells_total_in_matrix={len(clusters)}",
        f"n_cells_retained_published_clusters={len(cell_names)}",
        "retained_cluster_ids=1-15; unclassified cluster IDs excluded",
        "published_cell_categories=14; ODC1 and ODC2 combined as ODC",
        f"n_genes_in_matrix={len(gene_names)}",
        f"min_detection_fraction={args.min_detection_fraction}",
        f"min_detected_cells={min_detected}",
        f"n_genes_after_low_expression_filter={len(expressed_genes)}",
        f"top_k={args.top_k}",
        f"pLI_threshold=0.9",
        f"pLI_rows_gt_0.9={pi_rows}",
        f"pLI_gene_name_rows={pi_names}",
        f"pLI_unique_gene_symbols={len(pi_genes)}",
        f"n_annotations={len(manifest_rows)}",
        "annotation_count_definition=PI + 14 cell sets + 14 PI intersections",
        f"window_bp={args.window_bp}",
        f"gene_coordinate_reference={args.loc}",
    ]
    (out_dir / "D6B_GENESET_AUDIT.txt").write_text("\n".join(audit) + "\n")
    print("D6B_GENESETS_COMPLETE")
    print(f"n_annotations={len(manifest_rows)}")
    print(f"n_filtered_genes={len(expressed_genes)}")
    print(f"n_pi_rows={pi_rows}; n_pi_unique_symbols={len(pi_genes)}")


if __name__ == "__main__":
    main()
