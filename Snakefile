from Configs import *

#configfile: "config.yaml"
cfg = ParsedConfig(config)
wildcard_constraints:
    feature_type="all_cell_peaks|by_cluster_peaks|tiles|peaks|default",
    method="[^/]+"

include: "scripts/data_cleaning/Snakefile"
include: "scripts/feature_engineering/Snakefile"
include: "scripts/clustering/Snakefile"
include: "scripts/evaluation/Snakefile"
# include: "scripts/visualization/Snakefile"

rule all:
    input:
        rules.feature_engineering.input,
        # rules.metrics.output,
        # rules.embeddings.output

# ------------------------------------------------------------------------------
# Merge benchmark files
#
# Run this after the main pipeline using:
# snakemake --configfile config.yaml --cores 1 benchmarks
# ------------------------------------------------------------------------------

# rule benchmarks:
#     input:
#         script = "scripts/merge_benchmarks.py"
#     output:
#         cfg.ROOT / "benchmarks.csv"
#     message: "Merge all benchmarks"
#     params:
#         cmd = f"conda run -n {cfg.py_env} python"
#     shell: "{params.cmd} {input.script} -o {output} --root {cfg.ROOT}"