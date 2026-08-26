#!/usr/bin/env python3
"""
build_guide_hash.py — Build guide sequence hash table

Design rationale:
  Based on fba diagnostics: (1) >99% of Illumina errors are substitutions,
  with no insertions/deletions; (2) all guides are exactly 20bp.
  
  Therefore: pre-compute all Hamming=1 single-base substitution variants
  for O(1) hash lookup, replacing fba's BK-tree O(log N) edit-distance query.

Output: Python pickle format hash table (v3 — dual string + integer keys)
  - seq_to_idx:     {20bp_seq_str: guide_index}
  - seq_to_idx_int: {uint64_encoded: guide_index}
  - idx_to_id:      [guide_id_string]  (index -> ID mapping)
  
Size: 4,536 x 61 ~ 276,696 entries ~ 10 MB
"""

import sys
import os
import pickle
import time
from pathlib import Path

BASES = ['A', 'C', 'G', 'T']

# DNA → 2-bit encoding
_BASE2BITS = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
# Byte-level lookup (256 entries) for speed
_BYTE2BITS = [-1] * 256
_BYTE2BITS[ord('A')] = 0
_BYTE2BITS[ord('C')] = 1
_BYTE2BITS[ord('G')] = 2
_BYTE2BITS[ord('T')] = 3


def _encode_seq(seq: str) -> int:
    """Encode DNA sequence as uint64 (2 bits/base, max 32bp)."""
    val = 0
    for ch in seq:
        val = (val << 2) | _BASE2BITS[ch]
    return val


def _encode_seq_bytes(b: bytes, start: int, length: int) -> int:
    """Encode DNA subsequence from bytes using pre-computed lookup table."""
    val = 0
    for i in range(start, start + length):
        val = (val << 2) | _BYTE2BITS[b[i]]
    return val


def build_guide_hash(fasta_path: str, hash_path: str) -> dict:
    """
    Build two-level guide hash table from FASTA:
      Level 1: exact sequence -> guide index
      Level 2: all Hamming=1 neighbors -> parent guide index
    """
    t0 = time.time()

    # ── 1. Load guides ──
    guides = {}
    with open(fasta_path) as f:
        current_id = None
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                current_id = line[1:]
            elif current_id and line:
                guides[current_id] = line
                current_id = None

    n_guides = len(guides)
    guide_lengths = set(len(s) for s in guides.values())
    print(f"[1/3] Loaded {n_guides} guides, lengths: {guide_lengths}")

    # ── 2. Build exact sequence hash (dual string + integer keys) ──
    seq_to_idx: dict[str, int] = {}
    seq_to_idx_int: dict[int, int] = {}
    idx_to_id: list[str] = []

    for gid, seq in guides.items():
        idx = len(idx_to_id)
        idx_to_id.append(gid)
        seq_to_idx[seq] = idx
        seq_to_idx_int[_encode_seq(seq)] = idx

    exact_count = len(seq_to_idx)
    duplicate_count = n_guides - exact_count
    print(f"[2/3] Exact sequences: {exact_count} "
          f"(duplicates removed: {duplicate_count})")

    # ── 3. Generate Hamming=1 neighbors (dual string + integer keys) ──
    hamming1_added = 0
    collision_count = 0

    for gid, seq in guides.items():
        parent_idx = seq_to_idx[seq]
        for pos in range(len(seq)):
            orig = seq[pos]
            for alt in BASES:
                if alt == orig:
                    continue
                variant = seq[:pos] + alt + seq[pos + 1:]
                variant_int = _encode_seq(variant)
                if variant not in seq_to_idx:
                    seq_to_idx[variant] = parent_idx
                    seq_to_idx_int[variant_int] = parent_idx
                    hamming1_added += 1
                else:
                    collision_count += 1

    total = len(seq_to_idx)
    print(f"[3/3] Hamming=1 variants added: {hamming1_added:,}")
    print(f"      Collisions (variant == existing guide): {collision_count}")
    print(f"      Total hash entries: {total:,}")

    # ── 4. Save ──
    data = {
        'seq_to_idx': seq_to_idx,
        'seq_to_idx_int': seq_to_idx_int,
        'idx_to_id': idx_to_id,
        'n_guides': n_guides,
        'guide_length': list(guide_lengths)[0],
    }

    os.makedirs(os.path.dirname(hash_path) if os.path.dirname(hash_path) else '.',
                exist_ok=True)

    with open(hash_path, 'wb') as f:
        pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)

    size_mb = os.path.getsize(hash_path) / (1024 * 1024)
    elapsed = time.time() - t0
    print(f"\nSaved: {hash_path} ({size_mb:.1f} MB)")
    print(f"Build time: {elapsed:.1f}s")

    return data


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit("Usage: build_guide_hash.py <guides.fasta> <guide_hash.pkl>")
    build_guide_hash(sys.argv[1], sys.argv[2])
