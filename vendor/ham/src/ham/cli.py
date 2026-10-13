"""HAM CLI — Hash Accelerated Matcher command-line interface."""

import sys
import argparse

from .hash_builder import build_guide_hash
from .matcher import match_reads, load_whitelist, load_guide_hash
from .dedup import build_count_matrix
from .merge import merge_all


def cmd_build_hash(args):
    """Build guide sequence hash table."""
    build_guide_hash(args.fasta, args.output)


def cmd_match(args):
    """Match sgRNA reads to guide reference."""
    wl = load_whitelist(args.whitelist)
    gh = load_guide_hash(args.guide_hash)

    # Build custom chem_cfg kwargs when --chemistry custom
    match_kwargs = {}
    if args.chemistry == "custom":
        match_kwargs["chem_cfg"] = {
            "cb_start": args.cb_start,
            "cb_end": args.cb_end,
            "umi_start": args.umi_start,
            "umi_end": args.umi_end,
            "window_start": args.window_start,
            "window_end": args.window_end,
            "guide_len": args.guide_len,
        }

    match_reads(
        r1_path=args.read1,
        r2_path=args.read2,
        whitelist=wl,
        guide_hash=gh,
        output_path=args.output,
        max_reads=args.max_reads,
        threads=args.threads,
        low_memory=args.low_memory,
        cb_max_hamming=args.cb_max_hamming,
        chemistry=args.chemistry,
        **match_kwargs,
    )


def cmd_dedup(args):
    """UMI deduplication + MEX generation."""
    build_count_matrix(
        hits_path=args.input,
        output_dir=args.output_dir,
        umi_threshold=args.umi_threshold,
        umi_len=args.umi_len,
    )


def cmd_merge(args):
    """Merge per-lane count matrices."""
    lanes = []
    with open(args.lanes) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                lanes.append((parts[0], parts[1], parts[2]))
    merge_all(lanes, args.out, args.prefix)


def main():
    parser = argparse.ArgumentParser(
        description='HAM: Hash Accelerated Matcher for Perturb-seq sgRNA quantification')
    sub = parser.add_subparsers(dest='command', required=True)

    # ── build-hash ──
    p_build = sub.add_parser('build-hash', help='Build guide sequence hash table')
    p_build.add_argument('--fasta', required=True, help='Guide FASTA file')
    p_build.add_argument('--output', required=True, help='Output hash table (.pkl)')
    p_build.set_defaults(func=cmd_build_hash)

    # ── match ──
    p_match = sub.add_parser('match', help='Match sgRNA reads to guide reference')
    p_match.add_argument('-1', '--read1', required=True,
                         help='Read1 FASTQ (CB+UMI), comma-separated for multi-file')
    p_match.add_argument('-2', '--read2', required=True,
                         help='Read2 FASTQ (guide), comma-separated for multi-file')
    p_match.add_argument('-w', '--whitelist', required=True,
                         help='Cell barcode whitelist')
    p_match.add_argument('-g', '--guide-hash', required=True,
                         help='Guide hash table (.pkl)')
    p_match.add_argument('-o', '--output', required=True,
                         help='Output hits (.npz)')
    p_match.add_argument('-n', '--max-reads', type=int, default=None,
                         help='Max reads to process (for testing)')
    p_match.add_argument('-t', '--threads', type=int, default=1,
                         help='Parallel workers (one per FASTQ pair)')
    p_match.add_argument('--low-memory', action='store_true',
                         help='Use binary-search CB hash (~110 MB vs ~700 MB, '
                              'slightly slower)')
    p_match.add_argument('--cb-max-hamming', type=int, default=1,
                         help='Max Hamming distance for CB correction '
                              '(default: 1; use 2 if whitelist chemistry '
                              'differs from sequencing chemistry)')
    p_match.add_argument('--chemistry', type=str, default='10xv3',
                         choices=['10xv3', '10xv2-5p', '10xv2-5p-12umi', 'custom'],
                         help='10x chemistry: 10xv3 (3-prime, 12bp UMI), '
                              '10xv2-5p (5-prime v1/v2, 10bp UMI), '
                              '10xv2-5p-12umi (5-prime v3 GEM-X, 12bp UMI), '
                              'or custom (use --cb-start/--umi-start/etc.) '
                              '[default: 10xv3]')
    # Custom chemistry flags (used when --chemistry custom)
    p_match.add_argument('--cb-start', type=int, default=0,
                         help='CB start position in R1 [default: 0]')
    p_match.add_argument('--cb-end', type=int, default=16,
                         help='CB end position in R1 [default: 16]')
    p_match.add_argument('--umi-start', type=int, default=16,
                         help='UMI start position in R1 [default: 16]')
    p_match.add_argument('--umi-end', type=int, default=28,
                         help='UMI end position in R1 [default: 28]')
    p_match.add_argument('--window-start', type=int, default=28,
                         help='Guide window start in R2 [default: 28]')
    p_match.add_argument('--window-end', type=int, default=54,
                         help='Guide window end in R2 [default: 54]')
    p_match.add_argument('--guide-len', type=int, default=20,
                         help='Guide protospacer length in bp [default: 20]')
    p_match.set_defaults(func=cmd_match)

    # ── dedup ──
    p_dedup = sub.add_parser('dedup', help='UMI dedup + MEX matrix generation')
    p_dedup.add_argument('-i', '--input', required=True,
                         help='Assignments .npz from ham match')
    p_dedup.add_argument('-o', '--output-dir', required=True,
                         help='Output directory for MEX files')
    p_dedup.add_argument('-t', '--umi-threshold', type=int, default=1,
                         help='UMI Hamming distance threshold (default: 1)')
    p_dedup.add_argument('--umi-len', type=int, default=12,
                         help='UMI length in bp (12 for v3, 10 for 5prime v1)')
    p_dedup.set_defaults(func=cmd_dedup)

    # ── merge ──
    p_merge = sub.add_parser('merge', help='Merge per-lane count matrices')
    p_merge.add_argument('--lanes', required=True,
                         help='TSV file: lane_id, matrix_dir, suffix')
    p_merge.add_argument('--out', required=True, help='Output directory')
    p_merge.add_argument('--prefix', default='merged', help='Output file prefix')
    p_merge.set_defaults(func=cmd_merge)

    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
