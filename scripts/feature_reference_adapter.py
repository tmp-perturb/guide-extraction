#!/usr/bin/env python3
"""Feature Reference Adapter: wide-to-long decomposition of dual-sgRNA pair library
into single-guide indices compatible with simpleaf (guides.fasta + t2g.tsv).

Usage:
    feature_reference_adapter.py --csv <input.csv> \\
                                 --out-fasta <guides.fasta> \\
                                 --out-t2g <t2g.tsv> \\
                                 [--out-feature-ref <feature_ref.csv>]

The script:
  1. Loads raw_guides_k562_essential.csv (wide: sgID_A, sgID_B per row)
  2. Filters duplicated guide pairs
  3. Decomposes into single-guide long format (4,536 guides × 20bp)
  4. Collapses non-targeting/control guides to "Non-Targeting"
  5. Outputs FASTA + 2-column t2g (guide_id \\t guide_id, for alevin-fry)
"""

import argparse
import re
import sys

import pandas as pd

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
CONST_READ = "R2"
CONST_PATTERN = "(BC)"
CONST_FEATURE_TYPE = "CRISPR Guide Capture"
NONTARGETING_LABEL = "Non-Targeting"

SANITIZE_RE = re.compile(r"[^a-zA-Z0-9_-]")
NONTARGETING_RE = re.compile(r"non-targeting|ntc|control|intergenic", re.IGNORECASE)

FEATURE_REF_COLS = [
    "id", "name", "read", "pattern", "sequence",
    "feature_type", "target_gene_id", "target_gene_name",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Decompose dual-sgRNA library CSV into single-guide FASTA + t2g for simpleaf."
    )
    p.add_argument("--csv", required=True, help="Path to raw_guides_k562_essential.csv")
    p.add_argument("--out-fasta", required=True, help="Output FASTA file path")
    p.add_argument("--out-t2g", required=True, help="Output 2-column t2g file path")
    p.add_argument("--out-feature-ref", default=None, help="Optional: output Cell Ranger feature_ref.csv")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    # ---- 1. Load ----
    raw = pd.read_csv(args.csv, dtype=str, keep_default_na=False)
    print(f"Loaded {len(raw)} rows from {args.csv}", file=sys.stderr)

    # ---- 2. Pre-filter duplicates ----
    dup_col = "duplicated guide pair?"
    mask_dup = raw[dup_col].str.lower().isin(["true", "yes"])
    raw = raw[~mask_dup].copy()
    print(f"  After dedup filter: {len(raw)} rows", file=sys.stderr)

    # ---- 3. Wide-to-long decomposition ----
    guide_a = raw[["sgID_A", "targeting sequence A", "ensembl gene id", "gene"]].rename(
        columns={
            "sgID_A": "id",
            "targeting sequence A": "sequence",
            "ensembl gene id": "target_gene_id",
            "gene": "target_gene_name",
        }
    )
    guide_b = raw[["sgID_B", "targeting sequence B", "ensembl gene id", "gene"]].rename(
        columns={
            "sgID_B": "id",
            "targeting sequence B": "sequence",
            "ensembl gene id": "target_gene_id",
            "gene": "target_gene_name",
        }
    )
    long_df = pd.concat([guide_a, guide_b], ignore_index=True).fillna("")
    long_df = long_df[
        long_df["id"].str.strip().ne("") & long_df["sequence"].str.strip().ne("")
    ]
    print(f"  Long-format: {len(long_df)} single guides", file=sys.stderr)

    # ---- 4. Inject constants ----
    long_df["name"] = long_df["id"]
    long_df["read"] = CONST_READ
    long_df["pattern"] = CONST_PATTERN
    long_df["feature_type"] = CONST_FEATURE_TYPE

    # ---- 5. Sanitise ----
    long_df["id"] = long_df["id"].str.replace(SANITIZE_RE, "_", regex=True)
    long_df["name"] = long_df["name"].str.replace(SANITIZE_RE, "_", regex=True)

    # ---- 6. Collapse non-targeting ----
    is_ntc = long_df["target_gene_name"].str.contains(NONTARGETING_RE, na=False, regex=True)
    long_df.loc[is_ntc, "target_gene_name"] = NONTARGETING_LABEL
    long_df.loc[is_ntc, "target_gene_id"] = NONTARGETING_LABEL

    # ---- 7. Global dedup by id ----
    long_df = long_df.drop_duplicates(subset="id", keep="first")

    # ---- 8. Output FASTA ----
    with open(args.out_fasta, "w") as fh:
        for _, row in long_df.iterrows():
            fh.write(f">{row['id']}\n{row['sequence']}\n")
    print(f"  FASTA: {args.out_fasta}", file=sys.stderr)

    # ---- 9. Output t2g (2-column: guide_id \\t guide_id, no header) ----
    with open(args.out_t2g, "w") as fh:
        for _, row in long_df.iterrows():
            fh.write(f"{row['id']}\t{row['id']}\n")
    print(f"  T2G:   {args.out_t2g}", file=sys.stderr)

    # ---- 10. Optional feature_ref.csv ----
    if args.out_feature_ref:
        long_df[FEATURE_REF_COLS].to_csv(args.out_feature_ref, index=False)
        print(f"  Feature ref: {args.out_feature_ref}", file=sys.stderr)

    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
