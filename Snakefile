from Configs import *

cfg = ParsedConfig(config)
wildcard_constraints:
    feature_type="all_cell_peaks|by_cluster_peaks|tiles|peaks|default",
    method="[^/]+"

include: "scripts/feature_engineering/Snakefile"
include: "scripts/clustering/Snakefile"
include: "scripts/evaluation/Snakefile"

rule all:
    input:
        rules.feature_engineering.input,
        rules.clustering.input,
        rules.evaluation.input