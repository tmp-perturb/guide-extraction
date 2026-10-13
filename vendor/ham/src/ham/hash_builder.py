"""HAM: Guide hash table construction."""

import os
import pickle
import time

from .encoding import BASES, encode_seq


def build_guide_hash(fasta_path: str, hash_path: str) -> dict:
    """Build two-level guide hash table from FASTA.

    Level 1: exact sequence -> guide index
    Level 2: all Hamming=1 neighbours -> parent guide index

    Returns dict with keys: seq_to_idx, seq_to_idx_int, idx_to_id, n_guides.
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
        seq_to_idx_int[encode_seq(seq)] = idx

    exact_count = len(seq_to_idx)
    print(f"[2/3] Exact sequences: {exact_count} "
          f"(duplicates removed: {n_guides - exact_count})")

    # ── 3. Generate Hamming=1 neighbours (dual string + integer keys) ──
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
                variant_int = encode_seq(variant)
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
