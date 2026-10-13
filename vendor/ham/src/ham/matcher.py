"""HAM: Core matching engine — read pair processing and guide assignment."""

import os
import gzip
import pickle
import time
import numpy as np
from typing import Optional

from .encoding import (
    _BYTE2BITS, _BASE2BITS, BASES,
    encode_barcode, encode_seq, encode_umi_bytes,
    encode_window_bigint, extract_guides_from_bigint,
    generate_cb_variants,
)

# ── Read layout constants per chemistry ──
CHEMISTRY_CONFIGS = {
    "10xv3": {
        "cb_start": 0, "cb_end": 16,           # 16 bp cell barcode
        "umi_start": 16, "umi_end": 28,         # 12 bp UMI
        "window_start": 28, "window_end": 54,   # 26 bp guide window (R2[28:54])
        "guide_len": 20,                         # 20 bp protospacer
        "default_whitelist": "3M-february-2018.txt",
    },
    "10xv2-5p": {
        "cb_start": 0, "cb_end": 16,            # 16 bp cell barcode
        "umi_start": 16, "umi_end": 26,         # 10 bp UMI (5' v1/v2)
        "window_start": 16, "window_end": 35,   # 19 bp guide window at R2[16:35]
        "guide_len": 19,                         # 19 bp protospacer
        "default_whitelist": "737K-august-2016.txt",
    },
    "10xv2-5p-12umi": {
        "cb_start": 0, "cb_end": 16,            # 16 bp cell barcode
        "umi_start": 16, "umi_end": 28,         # 12 bp UMI (5' v3, GEM-X)
        "window_start": 16, "window_end": 35,   # 19 bp guide window at R2[16:35]
        "guide_len": 19,                         # 19 bp protospacer
        "default_whitelist": "3M-5pgex-jan-2023.txt",
    },
}
DEFAULT_CHEMISTRY = "10xv3"


def load_whitelist(path: str) -> set:
    """Load cell barcode whitelist as a Python set."""
    with open(path) as f:
        return set(line.strip() for line in f if line.strip())


def load_guide_hash(path: str) -> dict:
    """Load pre-built guide hash table from pickle."""
    with open(path, 'rb') as f:
        return pickle.load(f)


def build_cb_hash(whitelist: set, low_memory: bool = False,
                  cb_max_hamming: int = 2) -> dict:
    """Build compact CB correction hash table.

    Plan H (default, low_memory=False): Python dict[int→int], O(1) lookup, ~700 MB.
    Plan G (low_memory=True): sorted numpy int32 arrays + binary search,
                               O(log N) lookup, ~110 MB.

    cb_max_hamming: max Hamming distance for CB correction (1 or 2).
                    Hamming≤1 is stricter, reducing ambient false matches.
    """
    if low_memory:
        return _build_cb_hash_numpy(whitelist, cb_max_hamming)
    return _build_cb_hash_dict(whitelist, cb_max_hamming)


def _build_cb_hash_dict(whitelist: set, cb_max_hamming: int = 2) -> dict:
    """Plan H: Python int dict for O(1) CB lookup."""
    t0 = time.time()
    barcode_list = list(whitelist)
    bc_to_idx = {bc: i for i, bc in enumerate(barcode_list)}
    cb_map: dict[int, int] = {}
    total_theoretical = 0

    for bc in barcode_list:
        parent_idx = bc_to_idx[bc]
        for variant_str, _ in generate_cb_variants(bc, max_hamming=cb_max_hamming):
            total_theoretical += 1
            encoded = encode_barcode(variant_str)
            if encoded not in cb_map:
                cb_map[encoded] = parent_idx

    elapsed = time.time() - t0
    collision_rate = (total_theoretical - len(cb_map)) / total_theoretical * 100
    est_mb = len(cb_map) * 52 / (1024 * 1024)
    print(f"  CB hash (Plan H, dict, Hamming≤{cb_max_hamming}): {len(cb_map):,} entries "
          f"(collision {collision_rate:.1f}%) | "
          f"~{est_mb:.0f}MB [{elapsed:.1f}s]")

    return {'cb_mode': 'dict', 'cb_map': cb_map, 'barcode_list': barcode_list}


def _build_cb_hash_numpy(whitelist: set, cb_max_hamming: int = 2) -> dict:
    """Plan G: sorted numpy int32 arrays + binary search for low-memory CB lookup.

    Stores encoded CB variants and parent indices as two parallel sorted
    int32 arrays.  Lookup via np.searchsorted (C-level binary search).
    Memory: ~110 MB (vs ~700 MB for dict).
    """
    t0 = time.time()
    barcode_list = list(whitelist)
    bc_to_idx = {bc: i for i, bc in enumerate(barcode_list)}

    # Build dict first (fast dedup), then convert to numpy arrays
    cb_map: dict[int, int] = {}
    total_theoretical = 0

    for bc in barcode_list:
        parent_idx = bc_to_idx[bc]
        for variant_str, _ in generate_cb_variants(bc, max_hamming=cb_max_hamming):
            total_theoretical += 1
            encoded = encode_barcode(variant_str)
            if encoded not in cb_map:
                cb_map[encoded] = parent_idx

    # Convert dict to sorted numpy arrays (C-level sort, fast)
    n = len(cb_map)
    cb_keys = np.fromiter(cb_map.keys(), dtype=np.uint32, count=n)
    cb_vals = np.fromiter(cb_map.values(), dtype=np.int32, count=n)
    order = np.argsort(cb_keys, kind='quicksort')
    cb_keys = cb_keys[order]
    cb_vals = cb_vals[order]

    elapsed = time.time() - t0
    collision_rate = (total_theoretical - n) / total_theoretical * 100
    est_mb = (cb_keys.nbytes + cb_vals.nbytes) / (1024 * 1024)
    print(f"  CB hash (Plan G, numpy): {n:,} entries "
          f"(collision {collision_rate:.1f}%) | "
          f"~{est_mb:.0f}MB [{elapsed:.1f}s]")

    return {'cb_mode': 'numpy', 'cb_keys': cb_keys, 'cb_vals': cb_vals,
            'barcode_list': barcode_list}


# ── Worker function for multi-process I/O ──

def _process_one_fastq_pair(args: tuple) -> tuple:
    """Worker: process one R1+R2 FASTQ pair, return hits array.

    Args: (r1_path, r2_path, cb_hash, guide_hash_int, max_reads, chem_cfg)
    chem_cfg is a dict from CHEMISTRY_CONFIGS with constants for this chemistry.
    Returns: (np.ndarray shape (N,3) int32, stats_dict)
    """
    r1_path, r2_path, cb_hash, guide_hash_int, max_reads, chem_cfg = args

    cb_start, cb_end = chem_cfg["cb_start"], chem_cfg["cb_end"]
    umi_start, umi_end = chem_cfg["umi_start"], chem_cfg["umi_end"]
    win_start, win_end = chem_cfg["window_start"], chem_cfg["window_end"]
    guide_len = chem_cfg["guide_len"]

    cb_mode = cb_hash['cb_mode']
    barcode_list = cb_hash['barcode_list']
    bc_to_idx = {bc: i for i, bc in enumerate(barcode_list)}

    if cb_mode == 'numpy':
        cb_keys = cb_hash['cb_keys']
        cb_vals = cb_hash['cb_vals']

    MAX_ESTIMATE = 30_000_000
    hits = np.zeros((MAX_ESTIMATE, 3), dtype=np.int32)
    write_ptr = 0

    stats = {'total': 0, 'valid_cb': 0, 'cb_exact': 0, 'cb_corrected': 0,
             'matched': 0, 'exact_hit': 0, 'hamming1_hit': 0, 'no_guide': 0,
             'short_read2': 0}

    opener = gzip.open if r1_path.endswith('.gz') else open
    with opener(r1_path, 'rb') as f1, opener(r2_path, 'rb') as f2:
        while True:
            # Read1 record
            h1 = f1.readline()
            if not h1: break
            seq1 = f1.readline().strip()
            f1.readline(); f1.readline()

            # Read2 record
            h2 = f2.readline()
            if not h2: break
            seq2 = f2.readline().strip()
            f2.readline(); f2.readline()

            stats['total'] += 1

            # ── Step 1: CB correction ──
            if len(seq1) < umi_end:
                continue
            raw_cb_bytes = seq1[cb_start:cb_end]
            if any(_BYTE2BITS[b] < 0 for b in raw_cb_bytes):
                continue
            cb_int = encode_barcode(raw_cb_bytes.decode())

            if cb_mode == 'dict':
                cb_parent = cb_hash['cb_map'].get(cb_int)
                if cb_parent is None:
                    continue
            else:  # numpy
                pos = np.searchsorted(cb_keys, cb_int)
                if pos >= len(cb_keys) or cb_keys[pos] != cb_int:
                    continue
                cb_parent = int(cb_vals[pos])

            stats['valid_cb'] += 1
            cb_idx = bc_to_idx.get(barcode_list[cb_parent], -1)
            if cb_idx < 0:
                continue

            parent_cb = barcode_list[cb_parent]
            if raw_cb_bytes.decode() == parent_cb:
                stats['cb_exact'] += 1
            else:
                stats['cb_corrected'] += 1

            # ── Step 2: UMI encoding ──
            umi_len = umi_end - umi_start
            umi_int = encode_umi_bytes(seq1, umi_start, umi_len)

            # ── Step 3: Guide matching ──
            if len(seq2) < win_end:
                stats['short_read2'] += 1
                continue

            window_bytes = seq2[win_start:win_end]
            big_int = encode_window_bigint(window_bytes)
            guides = extract_guides_from_bigint(big_int, guide_len)

            found_idx = -1
            # Fast path: exact match
            for g in guides:
                idx = guide_hash_int.get(g)
                if idx is not None:
                    found_idx = idx
                    stats['exact_hit'] += 1
                    break

            # Slow path: Hamming=1
            if found_idx < 0:
                for offset, g in enumerate(guides):
                    for pos in range(guide_len):
                        orig_bits = (g >> ((guide_len - 1 - pos) * 2)) & 0x3
                        for alt_bits in range(4):
                            if alt_bits == orig_bits:
                                continue
                            shift = (guide_len - 1 - pos) * 2
                            variant = (g & ~(0x3 << shift)) | (alt_bits << shift)
                            idx = guide_hash_int.get(variant)
                            if idx is not None:
                                found_idx = idx
                                stats['hamming1_hit'] += 1
                                break
                        if found_idx >= 0:
                            break
                    if found_idx >= 0:
                        break

            if found_idx >= 0:
                stats['matched'] += 1
                hits[write_ptr, 0] = cb_idx
                hits[write_ptr, 1] = umi_int
                hits[write_ptr, 2] = found_idx
                write_ptr += 1
            else:
                stats['no_guide'] += 1

            if max_reads and stats['total'] >= max_reads:
                break

    return hits[:write_ptr], stats


# ── Main entry point ──

def match_reads(
    r1_path: str,
    r2_path: str,
    whitelist: set,
    guide_hash: dict,
    output_path: str,
    max_reads: Optional[int] = None,
    threads: int = 1,
    low_memory: bool = False,
    cb_max_hamming: int = 2,
    chemistry: str = "10xv3",
    report_interval: int = 1_000_000,
    chem_cfg: Optional[dict] = None,
) -> dict:
    """Core matching loop (HAM: integer encoding + numpy + multi-process I/O).

    When chemistry="custom", chem_cfg must be provided with keys:
    cb_start, cb_end, umi_start, umi_end, window_start, window_end, guide_len.
    """
    if chemistry == "custom":
        if not chem_cfg:
            raise ValueError(
                "chemistry='custom' requires a chem_cfg dict with keys: "
                "cb_start, cb_end, umi_start, umi_end, "
                "window_start, window_end, guide_len")
    else:
        chem_cfg = CHEMISTRY_CONFIGS.get(
            chemistry, CHEMISTRY_CONFIGS[DEFAULT_CHEMISTRY])
    mode_label = "Plan G (numpy, low-memory)" if low_memory else "Plan H (dict)"
    seq_to_idx = guide_hash['seq_to_idx']
    idx_to_id = guide_hash['idx_to_id']

    # Use integer-encoded guide hash if available
    guide_hash_int: dict = guide_hash.get('seq_to_idx_int')
    if guide_hash_int is None:
        print("Note: integer guide hash not found, converting...")
        guide_hash_int = {}
        for seq_str, idx in seq_to_idx.items():
            guide_hash_int[encode_seq(seq_str)] = idx

    # Build CB hash (Plan H dict or Plan G numpy)
    mode_label = "Plan G (numpy, low-memory)" if low_memory else "Plan H (dict)"
    print(f"Building CB hash ({mode_label}, Hamming≤{cb_max_hamming})...")
    cb_hash = build_cb_hash(whitelist, low_memory=low_memory,
                            cb_max_hamming=cb_max_hamming)
    print()

    # Split FASTQ paths
    r1_list = [p.strip() for p in r1_path.split(',') if p.strip()]
    r2_list = [p.strip() for p in r2_path.split(',') if p.strip()]
    assert len(r1_list) == len(r2_list), \
        f"Mismatched R1/R2 file counts: {len(r1_list)} vs {len(r2_list)}"

    n_files = len(r1_list)
    use_mp = threads > 1 and n_files > 1

    t0 = time.time()

    if use_mp:
        import multiprocessing as mp
        mp.set_start_method('fork', force=True)

        worker_args = [
            (r1, r2, cb_hash, guide_hash_int, max_reads, chem_cfg)
            for r1, r2 in zip(r1_list, r2_list)
        ]

        print(f"Processing {n_files} FASTQ pairs with "
              f"{min(threads, n_files)} workers...")
        with mp.Pool(processes=min(threads, n_files)) as pool:
            results = pool.map(_process_one_fastq_pair, worker_args)

        all_arrays = []
        merged_stats: dict = {}
        for arr, stats in results:
            if arr.size > 0:
                all_arrays.append(arr)
            for k, v in stats.items():
                merged_stats[k] = merged_stats.get(k, 0) + v
        all_hits = np.vstack(all_arrays) if all_arrays else \
            np.zeros((0, 3), dtype=np.int32)
    else:
        all_arrays = []
        merged_stats = {'total': 0, 'valid_cb': 0, 'cb_exact': 0, 'cb_corrected': 0,
                        'matched': 0, 'exact_hit': 0, 'hamming1_hit': 0, 'no_guide': 0,
                        'short_read2': 0}

        for r1, r2 in zip(r1_list, r2_list):
            print(f"Processing: {os.path.basename(r1)} + {os.path.basename(r2)}")
            arr, stats = _process_one_fastq_pair(
                (r1, r2, cb_hash, guide_hash_int, max_reads, chem_cfg))
            if arr.size > 0:
                all_arrays.append(arr)
            for k, v in stats.items():
                merged_stats[k] = merged_stats.get(k, 0) + v

        all_hits = np.vstack(all_arrays) if all_arrays else \
            np.zeros((0, 3), dtype=np.int32)

    elapsed = time.time() - t0

    # ── Save as .npz ──
    os.makedirs(os.path.dirname(output_path) if os.path.dirname(output_path) else '.',
                exist_ok=True)
    print(f"\nWriting {all_hits.shape[0]:,} hits...")
    barcode_list = cb_hash['barcode_list']
    np.savez_compressed(output_path,
                        hits=all_hits,
                        barcode_list=np.array(barcode_list),
                        idx_to_id=np.array(idx_to_id))

    # ── Summary ──
    stats = merged_stats
    print(f"\n{'='*60}")
    print("MATCHING COMPLETE (HAM v1.0)")
    print(f"{'='*60}")
    print(f"  Total reads:       {stats['total']:>12,}")
    valid_pct = stats['valid_cb'] / max(stats['total'], 1) * 100
    print(f"  Valid CB:          {stats['valid_cb']:>12,}  ({valid_pct:.1f}%)")
    print(f"    Exact CB:        {stats['cb_exact']:>12,}")
    print(f"    Corrected CB:    {stats['cb_corrected']:>12,}")
    print(f"  Short Read2:       {stats['short_read2']:>12,}")
    matched_pct = stats['matched'] / max(stats['valid_cb'], 1) * 100
    print(f"  Guide matched:     {stats['matched']:>12,}  ({matched_pct:.1f}%)")
    print(f"    Exact:           {stats['exact_hit']:>12,}")
    print(f"    Hamming=1:       {stats['hamming1_hit']:>12,}")
    print(f"  No guide detected: {stats['no_guide']:>12,}")
    print(f"  Time: {elapsed:.1f}s  ({stats['total']/elapsed:.0f} reads/s)")
    print(f"  Output: {output_path}")
    print(f"{'='*60}")

    return stats
