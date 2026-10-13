"""HAM: Integer encoding utilities for DNA sequences.

All DNA sequences are encoded as 2-bit-per-base integers:
  A=00, C=01, G=10, T=11

This module provides the shared encoding/decoding functions used across
matcher, dedup, and hash_builder modules.
"""

# ── Base-to-bits mapping ──
_BASE2BITS = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
_BITS2BASE = 'ACGT'
BASES = ('A', 'C', 'G', 'T')

# ── Byte-level lookup table (256 entries) ──
_BYTE2BITS = [-1] * 256
_BYTE2BITS[ord('A')] = 0
_BYTE2BITS[ord('C')] = 1
_BYTE2BITS[ord('G')] = 2
_BYTE2BITS[ord('T')] = 3

# ── 5bp chunk LUT (P0 optimisation) ──
# Index: b0<<8 | b1<<6 | b2<<4 | b3<<2 | b4.  Identity function since
# the index IS the correctly-ordered 10-bit encoding.
_5MER_ENC = list(range(1024))


def encode_barcode(bc: str) -> int:
    """Encode a 16bp DNA barcode into a 32-bit integer (2 bits/base)."""
    val = 0
    for ch in bc:
        val = (val << 2) | _BASE2BITS[ch]
    return val


def encode_seq(seq: str) -> int:
    """Encode a DNA sequence string as uint64 (2 bits/base)."""
    val = 0
    for ch in seq:
        val = (val << 2) | _BASE2BITS[ch]
    return val


def encode_seq_bytes(b: bytes, start: int, length: int) -> int:
    """Encode a DNA subsequence from bytes using pre-computed LUT."""
    val = 0
    for i in range(start, start + length):
        val = (val << 2) | _BYTE2BITS[b[i]]
    return val


def encode_umi_bytes(b: bytes, start: int, umi_len: int = 12) -> int:
    """Encode UMI from bytes as uint32 (2 bits/base), loop-unrolled.

    For 12bp UMI (10xv3): umi_len=12 → 24 bits.
    For 10bp UMI (10xv2-5p): umi_len=10 → 20 bits.
    """
    v = _BYTE2BITS
    if umi_len == 12:
        return (v[b[start]]   << 22 | v[b[start+1]]  << 20 |
                v[b[start+2]] << 18 | v[b[start+3]]  << 16 |
                v[b[start+4]] << 14 | v[b[start+5]]  << 12 |
                v[b[start+6]] << 10 | v[b[start+7]]  << 8  |
                v[b[start+8]] << 6  | v[b[start+9]]  << 4  |
                v[b[start+10]] << 2 | v[b[start+11]])
    elif umi_len == 10:
        return (v[b[start]]   << 18 | v[b[start+1]]  << 16 |
                v[b[start+2]] << 14 | v[b[start+3]]  << 12 |
                v[b[start+4]] << 10 | v[b[start+5]]  << 8  |
                v[b[start+6]] << 6  | v[b[start+7]]  << 4  |
                v[b[start+8]] << 2  | v[b[start+9]])
    else:
        raise ValueError(f"Unsupported UMI length: {umi_len}")


def decode_umi(val: int, umi_len: int = 12) -> str:
    """Decode uint32 UMI back to DNA string."""
    chars = []
    for _ in range(umi_len):
        chars.append(_BITS2BASE[val & 0x3])
        val >>= 2
    return ''.join(reversed(chars))


def encode_window_bigint(window_bytes: bytes) -> int:
    """Encode a DNA window as a bigint (up to 52-bit, optimised).

    For 26bp windows (10xv3): 5×5bp chunk LUT + 1 trailing base.
    For 19bp windows (10xv2-5p): 3×5bp chunk LUT + 4 trailing bases.
    """
    v = _BYTE2BITS
    b = window_bytes
    n = len(b)

    if n == 26:
        i0 = (v[b[0]]  << 8 | v[b[1]]  << 6 | v[b[2]]  << 4 |
              v[b[3]]  << 2 | v[b[4]])
        i1 = (v[b[5]]  << 8 | v[b[6]]  << 6 | v[b[7]]  << 4 |
              v[b[8]]  << 2 | v[b[9]])
        i2 = (v[b[10]] << 8 | v[b[11]] << 6 | v[b[12]] << 4 |
              v[b[13]] << 2 | v[b[14]])
        i3 = (v[b[15]] << 8 | v[b[16]] << 6 | v[b[17]] << 4 |
              v[b[18]] << 2 | v[b[19]])
        i4 = (v[b[20]] << 8 | v[b[21]] << 6 | v[b[22]] << 4 |
              v[b[23]] << 2 | v[b[24]])
        last = v[b[25]]
        return (_5MER_ENC[i0] << 42 |
                _5MER_ENC[i1] << 32 |
                _5MER_ENC[i2] << 22 |
                _5MER_ENC[i3] << 12 |
                _5MER_ENC[i4] << 2  |
                last)
    elif n == 19:
        i0 = (v[b[0]]  << 8 | v[b[1]]  << 6 | v[b[2]]  << 4 |
              v[b[3]]  << 2 | v[b[4]])
        i1 = (v[b[5]]  << 8 | v[b[6]]  << 6 | v[b[7]]  << 4 |
              v[b[8]]  << 2 | v[b[9]])
        i2 = (v[b[10]] << 8 | v[b[11]] << 6 | v[b[12]] << 4 |
              v[b[13]] << 2 | v[b[14]])
        last4 = (v[b[15]] << 6 | v[b[16]] << 4 | v[b[17]] << 2 | v[b[18]])
        return (_5MER_ENC[i0] << 28 |
                _5MER_ENC[i1] << 18 |
                _5MER_ENC[i2] << 8  |
                last4)
    else:
        raise ValueError(f"Unsupported window size: {n}")


def extract_guides_from_bigint(big_int: int, guide_len: int = 20) -> list:
    """Extract sliding guide codes from a window int via shift+mask.

    For 20bp guide in 26bp window (10xv3): offset=0..6, 7 guides.
    For 19bp guide in 19bp window (10xv2-5p): offset=0 only, 1 guide.
    """
    n_offsets = 7 if guide_len == 20 else 1
    mask = (1 << (guide_len * 2)) - 1
    guides = []
    for offset in range(n_offsets):
        shift = ((n_offsets + guide_len - 1 - offset) * 2) if guide_len == 20 else 0
        guides.append((big_int >> shift) & mask)
    return guides


def generate_cb_variants(bc: str, max_hamming: int = 2):
    """Generate all Hamming <= max_hamming variants of a 16bp barcode.

    Yields (variant_str, hamming_distance).
    """
    yield bc, 0
    n = len(bc)
    # Hamming = 1
    for p1 in range(n):
        o1 = bc[p1]
        for a1 in BASES:
            if a1 == o1:
                continue
            yield bc[:p1] + a1 + bc[p1+1:], 1
    # Hamming = 2
    if max_hamming >= 2:
        for p1 in range(n):
            o1 = bc[p1]
            for a1 in BASES:
                if a1 == o1:
                    continue
                v1 = bc[:p1] + a1 + bc[p1+1:]
                for p2 in range(p1+1, n):
                    o2 = bc[p2]
                    for a2 in BASES:
                        if a2 == o2:
                            continue
                        yield v1[:p2] + a2 + v1[p2+1:], 2
