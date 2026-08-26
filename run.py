#!/usr/bin/env python3
"""Omnibenchmark module entrypoint: guide_extraction.

Single-lane guide extraction: from an sgRNA FASTQ pair (R1=CB+UMI, R2=protospacer),
a GEX matrix (for the cell-barcode whitelist), and a guide-library CSV, produce one
merged (cells x guides) UMI count matrix in MEX format.

This is a re-orchestration of the standalone Snakemake guide_extraction workflow.
The scientific logic is unchanged: it calls the SAME helper scripts
(feature_reference_adapter.py, filter_barcodes.py, translate_barcodes.py,
build_guide_hash.py) and the SAME tools (simpleaf / piscem / alevin-fry, ham).
Only the outer orchestration + machine-specific plumbing (Snakemake driver, conda
activation, hard-coded paths) is replaced by this entrypoint.

Pipeline (per lane):
    guide_csv -> feature_reference_adapter -> guides.fasta + t2g.tsv
      simpleaf method:  fasta -> simpleaf index (piscem)
      hash_matcher:     fasta -> build_guide_hash -> guide_hash.pkl
    gex_h5 -> whitelist (UMI/genes QC) -> barcode_whitelist_noheader.txt
      simpleaf method:  simpleaf quant (--explicit-pl wl | --knee) -> alevin MEX
      hash_matcher:     ham match -> ham dedup -> MEX
    per-lane MEX -> ham merge -> merged_{matrix,barcodes,features}

Omnibenchmark CLI contract (flags injected at runtime):
    --output_dir <dir>       directory to write outputs into (REQUIRED)
    --name <node_id>         current node id, used for the lane suffix (REQUIRED)
    --data.guide_csv <path>  guide-library CSV
    --data.gex_h5 <path>     GEX matrix (.h5 / .h5ad / .h5mu)
    --data.sgRNA_r1 <path>   sgRNA FASTQ R1 (comma-separated for multi-file lane)
    --data.sgRNA_r2 <path>   sgRNA FASTQ R2 (comma-separated for multi-file lane)
    parameters: --method --tenx_chemistry --kmer_length --minimizer_length
                --resolution --use_knee --min_umi --min_genes
                --umi_threshold --cb_max_hamming [--translation_table]

Outputs written into <output_dir>:
    merged_matrix.mtx.gz     cells x guides, deduplicated UMI counts
    merged_barcodes.tsv.gz   cell barcodes (with -L<NN> lane suffix)
    merged_features.tsv.gz   guide feature ids
"""
import argparse
import os
import re
import shutil
import subprocess
import sys

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(HERE, "scripts")
CHEMISTRY_SPEC_PATH = os.path.join(HERE, "config", "chemistry_spec.yaml")


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _run(cmd, env=None):
    """Run a subprocess, echoing the command; raise on failure."""
    print("+ " + " ".join(str(c) for c in cmd), flush=True)
    subprocess.run([str(c) for c in cmd], check=True, env=env)


def _bool(v):
    if isinstance(v, bool):
        return v
    return str(v).strip().lower() in ("1", "true", "yes", "on")


def load_chemistry(chem_name):
    """Resolve tenx_chemistry -> downstream parameters via chemistry_spec.yaml.

    Ported from the standalone Snakefile's _resolve_chemistry (custom chemistry
    branch omitted for this minimal module)."""
    with open(CHEMISTRY_SPEC_PATH) as fh:
        spec_all = yaml.safe_load(fh)
    if chem_name not in spec_all:
        sys.exit(f"ERROR: unknown tenx_chemistry '{chem_name}'. "
                 f"Known: {', '.join(sorted(spec_all))}")
    return spec_all[chem_name]


def lane_suffix(name):
    """Trailing digits -> '-L<NN>', else '-<name>' (matches merge.smk logic)."""
    m = re.search(r"(\d+)$", name)
    return f"-L{m.group(1)}" if m else f"-{name}"


# ---------------------------------------------------------------------------
# pipeline steps
# ---------------------------------------------------------------------------
def build_reference(guide_csv, ref_dir):
    fasta = os.path.join(ref_dir, "guides.fasta")
    t2g = os.path.join(ref_dir, "t2g.tsv")
    _run([sys.executable, os.path.join(SCRIPTS, "feature_reference_adapter.py"),
          "--csv", guide_csv, "--out-fasta", fasta, "--out-t2g", t2g])
    return fasta, t2g


def build_index(fasta, ref_dir, kmer, minimizer, threads, af_home):
    idx_dir = os.path.join(ref_dir, "piscem_index")
    os.makedirs(idx_dir, exist_ok=True)
    env = dict(os.environ, ALEVIN_FRY_HOME=af_home)
    # simpleaf needs its info JSON in ALEVIN_FRY_HOME before index/quant.
    _run(["simpleaf", "set-paths"], env=env)
    _run(["simpleaf", "index", "--output", idx_dir, "--threads", threads,
          "--ref-seq", fasta, "--kmer-length", kmer,
          "--minimizer-length", minimizer], env=env)
    # simpleaf writes to idx_dir/index/; quant references the piscem prefix there.
    return os.path.join(idx_dir, "index", "piscem_idx")


def build_hash(fasta, ref_dir):
    hash_file = os.path.join(ref_dir, "guide_hash.pkl")
    _run([sys.executable, os.path.join(SCRIPTS, "build_guide_hash.py"),
          fasta, hash_file])
    return hash_file


def extract_whitelist(gex_h5, lane_dir, min_umi, min_genes, chem,
                      translation_table):
    """Derive cell-barcode whitelist from GEX matrix (.h5/.h5ad/.h5mu).

    Ported from whitelist.smk. Applies TruSeq->Nextera translation (to_from)
    for dual-oligo (3') chemistries when translation is enabled."""
    wl_csv = os.path.join(lane_dir, "barcode_whitelist.csv")
    wl_noheader = os.path.join(lane_dir, "barcode_whitelist_noheader.txt")
    ext = gex_h5.rsplit(".", 1)[-1].lower()
    # inline (h5ad/h5mu) branches compare numerically; filter_barcodes.py takes
    # CLI strings itself, so cast here only for the inline paths.
    min_umi_i, min_genes_i = int(min_umi), int(min_genes)

    if ext == "h5ad":
        import anndata
        from scipy.sparse import issparse
        ad = anndata.read_h5ad(gex_h5)
        X = ad.X
        total = X.sum(axis=1).A1 if issparse(X) else X.sum(axis=1)
        ngenes = (X > 0).sum(axis=1).A1 if issparse(X) else (X > 0).sum(axis=1)
        qc = (total >= min_umi_i) & (ngenes >= min_genes_i)
        _write_wl(_clean_barcodes(ad.obs_names[qc]), wl_csv, wl_noheader)
    elif ext == "h5mu":
        import mudata as md
        from scipy.sparse import issparse
        ad = md.read_h5mu(gex_h5).mod["rna"]
        X = ad.X
        total = X.sum(axis=1).A1 if issparse(X) else X.sum(axis=1)
        ngenes = (X > 0).sum(axis=1).A1 if issparse(X) else (X > 0).sum(axis=1)
        qc = (total >= min_umi_i) & (ngenes >= min_genes_i)
        _write_wl(_clean_barcodes(ad.obs_names[qc]), wl_csv, wl_noheader)
    else:  # scprocess raw .h5 (CSC) via filter_barcodes.py
        _run([sys.executable, os.path.join(SCRIPTS, "filter_barcodes.py"),
              "--h5", gex_h5, "--out-wl", wl_csv,
              "--out-noheader", wl_noheader,
              "--min-umi", min_umi, "--min-genes", min_genes])

    if chem.get("translation"):
        if not translation_table or not os.path.exists(translation_table):
            sys.exit("ERROR: chemistry requires barcode translation but no valid "
                     "--translation_table was provided. Provide the 10x RNA<->Feature "
                     "translation table, or use a 5' chemistry (translation:false).")
        _run([sys.executable, os.path.join(SCRIPTS, "translate_barcodes.py"),
              wl_noheader, "--trans-table", translation_table,
              "--direction", "to_from"])
    return wl_noheader


def _clean_barcodes(obs_names):
    out = []
    for b in obs_names:
        b = b.decode() if isinstance(b, bytes) else b
        out.append(b.split("_")[1] if "_" in b else b)
    return out


def _write_wl(barcodes, wl_csv, wl_noheader):
    with open(wl_csv, "w") as f:
        f.write("barcode\n")
        f.writelines(b + "\n" for b in barcodes)
    with open(wl_noheader, "w") as f:
        f.writelines(b + "\n" for b in barcodes)
    print(f"Whitelist: {len(barcodes)} cells", flush=True)


def simpleaf_quant(reads1, reads2, wl, idx_prefix, t2g, lane_dir, chem,
                   resolution, use_knee, threads, af_home):
    out_dir = os.path.join(lane_dir, "simpleaf_quant")
    if os.path.isdir(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    chemistry = chem.get("geometry_override") or chem["af_chemistry"]
    env = dict(os.environ, ALEVIN_FRY_HOME=af_home)
    _run(["simpleaf", "set-paths"], env=env)
    cmd = ["simpleaf", "quant", "--chemistry", chemistry,
           "--output", out_dir, "--threads", threads,
           "--index", idx_prefix, "--reads1", reads1, "--reads2", reads2,
           "--t2g-map", t2g, "--resolution", resolution]
    cmd += ["--knee"] if use_knee else ["--explicit-pl", wl]
    _run(cmd, env=env)
    # merge (ham merge) reads alevin's quants_mat.{mtx,rows,cols}
    return os.path.join(out_dir, "af_quant", "alevin")


def hash_quant(reads1, reads2, wl, hash_file, lane_dir, chem,
               umi_threshold, cb_max_hamming, threads):
    out_dir = os.path.join(lane_dir, "guide_quant")
    os.makedirs(out_dir, exist_ok=True)
    hits = os.path.join(out_dir, "hits.npz")
    _run(["ham", "match", "-1", reads1, "-2", reads2, "-w", wl,
          "-g", hash_file, "-o", hits, "-t", threads,
          "--cb-max-hamming", cb_max_hamming,
          "--chemistry", chem.get("ham_chemistry", "10xv3")])
    umi_len = chem.get("umi_len", 12)
    _run(["ham", "dedup", "-i", hits, "-o", os.path.join(out_dir, "matrix"),
          "-t", umi_threshold, "--umi-len", umi_len])
    return os.path.join(out_dir, "matrix")


def merge(quant_dir, name, output_dir, chem, translation_table):
    """ham merge over a single lane -> merged MEX trio (ported from merge.smk)."""
    lane_list = os.path.join(output_dir, ".lane_list.tsv")
    with open(lane_list, "w") as f:
        f.write(f"{name}\t{quant_dir}\t{lane_suffix(name)}\n")
    _run(["ham", "merge", "--lanes", lane_list, "--out", output_dir,
          "--prefix", "merged"])
    if chem.get("translation") and translation_table:
        _run([sys.executable, os.path.join(SCRIPTS, "translate_barcodes.py"),
              os.path.join(output_dir, "merged_barcodes.tsv.gz"),
              "--trans-table", translation_table, "--direction", "from_to"])


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description="Omnibenchmark module: guide_extraction")
    # contract: framework-injected
    p.add_argument("--output_dir", required=True)
    p.add_argument("--name", default="lane")
    # contract: stage inputs (dotted ids -> read via getattr)
    p.add_argument("--data.guide_csv", required=True)
    p.add_argument("--data.gex_h5", required=True)
    p.add_argument("--data.sgRNA_r1", required=True)
    p.add_argument("--data.sgRNA_r2", required=True)
    # contract: parameters
    p.add_argument("--method", default="simpleaf",
                   choices=["simpleaf", "hash_matcher"])
    p.add_argument("--tenx_chemistry", default="3v3")
    p.add_argument("--kmer_length", default="15")
    p.add_argument("--minimizer_length", default="11")
    p.add_argument("--resolution", default="parsimony-gene")
    p.add_argument("--use_knee", default="false")
    p.add_argument("--min_umi", default="1000")
    p.add_argument("--min_genes", default="500")
    p.add_argument("--umi_threshold", default="1")
    p.add_argument("--cb_max_hamming", default="1")
    # optional: required only for translating (3') chemistries
    p.add_argument("--translation_table", default="")
    # optional: threads (execution, not a benchmark parameter)
    p.add_argument("--threads", default="4")
    args = p.parse_args()

    guide_csv = getattr(args, "data.guide_csv")
    gex_h5 = getattr(args, "data.gex_h5")
    reads1 = getattr(args, "data.sgRNA_r1")
    reads2 = getattr(args, "data.sgRNA_r2")

    output_dir = os.path.abspath(args.output_dir)
    ref_dir = os.path.join(output_dir, "refs")
    lane_dir = os.path.join(output_dir, "lanes", args.name)
    af_home = os.path.join(ref_dir, "af_home")
    for d in (ref_dir, lane_dir, af_home):
        os.makedirs(d, exist_ok=True)

    chem = load_chemistry(args.tenx_chemistry)
    trans_table = args.translation_table or ""

    fasta, t2g = build_reference(guide_csv, ref_dir)
    wl = extract_whitelist(gex_h5, lane_dir, args.min_umi, args.min_genes,
                           chem, trans_table)

    if args.method == "simpleaf":
        idx_prefix = build_index(fasta, ref_dir, args.kmer_length,
                                 args.minimizer_length, args.threads, af_home)
        quant_dir = simpleaf_quant(
            reads1, reads2, wl, idx_prefix, t2g, lane_dir, chem,
            args.resolution, _bool(args.use_knee), args.threads, af_home)
    else:
        hash_file = build_hash(fasta, ref_dir)
        quant_dir = hash_quant(
            reads1, reads2, wl, hash_file, lane_dir, chem,
            args.umi_threshold, args.cb_max_hamming, args.threads)

    merge(quant_dir, args.name, output_dir, chem, trans_table)

    expected = [os.path.join(output_dir, f"merged_{x}")
                for x in ("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz")]
    missing = [e for e in expected if not os.path.exists(e)]
    if missing:
        sys.exit("ERROR: expected outputs not produced: " + ", ".join(missing))
    print("guide_extraction: wrote " + ", ".join(os.path.basename(e) for e in expected))


if __name__ == "__main__":
    main()
