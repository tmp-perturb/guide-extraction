# guide-extraction

Omnibenchmark module for the K562 lane01 guide-extraction benchmark. It
re-orchestrates the existing simpleaf/HAM workflow and writes the merged MEX
trio expected by the benchmark plan.

The entrypoint accepts `--output_dir`, `--name`, the four `data.*` inputs, and
the extraction parameters declared in the plan. Outputs are
`merged_matrix.mtx.gz`, `merged_barcodes.tsv.gz`, and
`merged_features.tsv.gz`.
