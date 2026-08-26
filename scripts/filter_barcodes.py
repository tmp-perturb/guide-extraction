#!/usr/bin/env python3
"""Extract cell barcode whitelist from alevin-fry H5 count matrix via hard thresholds.

Usage:
    filter_barcodes.py --h5 <af_counts_mat.h5> \\
                       --out-wl <barcode_whitelist.csv> \\
                       [--min-umi 1000] [--min-genes 500] \\
                       [--out-noheader <barcode_whitelist_noheader.txt>]

Output:
    barcode_whitelist.csv         — CSV with header "barcode"
    barcode_whitelist_noheader.txt — bare barcode list (for simpleaf --explicit-pl)
"""

import argparse
import sys

import h5py
import numpy as np
from scipy.sparse import csc_matrix


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Extract cell barcode whitelist from GEX H5 matrix.")
    p.add_argument("--h5", required=True, help="Path to af_counts_mat.h5 from scprocess/alevin-fry")
    p.add_argument("--out-wl", required=True, help="Output CSV whitelist file path")
    p.add_argument("--out-noheader", default=None, help="Optional: output whitelist without header (for simpleaf)")
    p.add_argument("--min-umi", type=int, default=1000, help="Minimum UMI count per cell (default: 1000)")
    p.add_argument("--min-genes", type=int, default=500, help="Minimum genes detected per cell (default: 500)")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    # ---- 1. Load mRNA matrix (CSC: genes × cells) ----
    print(f"Loading: {args.h5}", file=sys.stderr)
    with h5py.File(args.h5, "r") as f:
        g = f["matrix"]
        shape = tuple(g["shape"][:])  # type: ignore[arg-type]
        mat = csc_matrix(
            (g["data"][:], g["indices"][:], g["indptr"][:]),
            shape=shape,
        )
        all_barcodes = g["barcodes"][:].astype(str)

    n_genes, n_cells = shape
    print(f"  Matrix: {n_genes} genes × {n_cells:,} cells, {mat.nnz:,} non-zero", file=sys.stderr)

    # ---- 2. Per-cell QC metrics ----
    umi = np.array(mat.sum(axis=0)).flatten()
    genes_detected = np.array((mat > 0).sum(axis=0)).flatten()
    print(
        f"  UMI:   median={np.median(umi):.0f}  min={umi.min():.0f}  max={umi.max():.0f}",
        file=sys.stderr,
    )
    print(
        f"  Genes: median={np.median(genes_detected):.0f}  min={genes_detected.min():.0f}  max={genes_detected.max():.0f}",
        file=sys.stderr,
    )

    # ---- 3. Hard threshold filter ----
    mask = (umi > args.min_umi) & (genes_detected > args.min_genes)
    mrna_bc = all_barcodes[mask]
    print(
        f"  Pass:  {mrna_bc.shape[0]:,} / {n_cells:,}  ({mrna_bc.shape[0] / n_cells * 100:.1f}%)",
        file=sys.stderr,
    )

    # ---- 4. Export ----
    np.savetxt(args.out_wl, mrna_bc, fmt="%s", header="barcode", comments="")
    print(f"  Saved: {args.out_wl}", file=sys.stderr)

    if args.out_noheader:
        np.savetxt(args.out_noheader, mrna_bc, fmt="%s", comments="")
        print(f"  Saved (no header): {args.out_noheader}", file=sys.stderr)

    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
