#!/usr/bin/env python3
"""Translate cell barcodes between 10x sequencer format and 3v3 whitelist standard.

10x 3' v3 chemistry barcodes exist in two equivalent representations that differ
by a complementary base swap at positions 7–8:
  - FROM (sequencer format): the barcode as read by the sequencer.
  - TO   (3v3 whitelist):    the barcode as defined in Cell Ranger's 3v3 whitelist.

Cell Ranger / alevin-fry internally applies a FROM→TO translation using the
official translation table (cellranger_whitelist_translation_3v3.txt) during
barcode correction.  If you use a whitelist derived from Cell Ranger output but
quantify with a tool that does NOT perform this translation (e.g., HAM), you
must apply it yourself so that the whitelist and FASTQ barcodes are in the same
representation.

This script applies FROM→TO by default; --direction to_from applies TO→FROM.
A backup copy with "_from_backup" suffix is created before overwriting the input.

Usage:
    python3 translate_barcodes.py merged_barcodes.tsv.gz \\
        --trans-table cellranger_whitelist_translation_3v3.txt \\
        --direction from_to
"""

import argparse
import gzip
import os
import shutil
import sys


def load_translation(path: str, direction: str = "from_to") -> dict:
    """Load the FROM↔TO translation table.

    The file format is two columns (tab-separated):
        FROM_barcode  TO_barcode
    """
    trans = {}
    with open(path) as f:
        for line in f:
            from_bc, to_bc = line.strip().split()
            if direction == "from_to":
                trans[from_bc] = to_bc   # FROM → TO
            else:
                trans[to_bc] = from_bc   # TO → FROM
    return trans


def translate_barcodes(input_path: str, trans: dict) -> tuple:
    """Translate barcodes in a gzipped TSV.

    Barcodes may have a lane suffix (e.g. '-L01'); only the 16 bp prefix
    is translated.  Returns (translated_lines, untranslated_count).
    """
    lines = []
    untranslated = 0

    opener = gzip.open if input_path.endswith('.gz') else open
    with opener(input_path, 'rt') as f:
        for line in f:
            bc = line.strip()
            parts = bc.rsplit('-', 1)
            base = parts[0]
            suffix = f'-{parts[1]}' if len(parts) > 1 else ''

            translated = trans.get(base)
            if translated is not None:
                lines.append(f'{translated}{suffix}\n')
            else:
                lines.append(line)
                untranslated += 1

    return lines, untranslated


def main():
    p = argparse.ArgumentParser(
        description="Translate 10x cell barcodes between sequencer format "
                    "and 3v3 whitelist standard")
    p.add_argument('input', help='Gzipped barcode TSV file to translate')
    p.add_argument('--trans-table', required=True,
                   help='Path to cellranger_whitelist_translation_3v3.txt')
    p.add_argument('--direction', default='from_to',
                   choices=['from_to', 'to_from'],
                   help='Direction: from_to (sequencer→whitelist, default) '
                        'or to_from (whitelist→sequencer)')
    args = p.parse_args()

    trans = load_translation(args.trans_table, args.direction)
    print(f"Loaded {len(trans):,} translation entries ({args.direction})",
          file=sys.stderr)

    lines, untranslated = translate_barcodes(args.input, trans)
    print(f"Translated {len(lines):,} barcodes, {untranslated} untranslated",
          file=sys.stderr)

    # Backup original. Derive a distinct backup path that works whether or not
    # the input is gzipped (the whitelist is a plain .txt; merged barcodes are
    # .tsv.gz). A blind .replace('.gz', ...) is a no-op for non-.gz inputs and
    # would collide with the source file.
    if args.input.endswith('.gz'):
        backup = args.input[:-3] + '_from_backup.gz'
    else:
        backup = args.input + '.from_backup'
    shutil.copy2(args.input, backup)
    print(f"Backup: {backup}", file=sys.stderr)

    # Write translated
    opener = gzip.open if args.input.endswith('.gz') else open
    with opener(args.input, 'wt') as f:
        f.writelines(lines)
    print(f"Written: {args.input}", file=sys.stderr)


if __name__ == '__main__':
    main()
