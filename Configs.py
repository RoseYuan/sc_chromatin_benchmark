from pathlib import Path
import snakemake.io
from collections import defaultdict
import itertools


class ParsedConfig:
    FEATURE_TYPES = ["all_cell_peaks", "by_cluster_peaks", "tiles", "peaks", "default"]

    def __init__(self, config):

        # TODO: define and check schema of config

        self.ROOT = Path(config["ROOT"]).resolve()
        self.DATA_SCENARIOS = config["DATA_SCENARIOS"]
        self.SCALING = config["SCALING"]
        self.FEATURE_SELECTION = config["FEATURE_SELECTION"]
        self.METHODS = config["METHODS"]
        self.GRAPH_CONSTRUCTION = config["GRAPH_CONSTRUCTION"]
        self.r_env = config["r_env"]
        self.py_env = config["py_env"]
        self.MACS2_PATH = config["MACS2_PATH"]
        self.reticulate_py = config["reticulate_py"]


    def get_all_feature_selections(self):
        return list(self.FEATURE_SELECTION.keys())

    # --------------------------------------------------------------------------
    # Gets all available methods. filter for framework (R/python) if needed.
    #
    # @param framework Only methods based on the framework will be retrieved.
    # one of ["python", "R", "both"], default: both
    # --------------------------------------------------------------------------
    def get_all_methods(self, framework="both"):
        all_methods = []
        for method in self.METHODS:
            is_r = self.get_from_method(method, "R")
            if framework == "both":
                all_methods.append(method)
            elif (framework == "python") and (not is_r):
                all_methods.append(method)
            elif (framework == "R") and (is_r):
                all_methods.append(method)

        return all_methods

    def get_all_scenarios(self):
        return list(self.DATA_SCENARIOS.keys())

    def get_feature_selection(self, key):
        if key not in self.FEATURE_SELECTION:
            raise ValueError(f"{key} not a valid key for feature selection")
        value = self.FEATURE_SELECTION[key]
        if key == 'ndim':
            return value if isinstance(value, list) else [value]
        return value
    
    def get_graph_construction(self, key):
        if key not in self.GRAPH_CONSTRUCTION:
            raise ValueError(f"{key} not a valid key for feature selection")
        value = self.GRAPH_CONSTRUCTION[key]
        return value
    
    def get_from_method(self, method, key):
        if method not in self.METHODS:
            raise ValueError(f"{method} not defined as method")
        if key not in self.METHODS[method]:
            return False
            # raise ValueError(f"{key} not a valid attribute of scenario {scenario}")
        value = self.METHODS[method][key]
        if key == 'distance':
            return value if isinstance(value, list) else [value]
        if key == 'feature_type':
            return value if isinstance(value, list) else [value]
        if key == 'tile_size':
            return value if isinstance(value, list) else [value]
        return value

    def get_from_method_language(self, method):
        if self.get_from_method(method, key="R"):
            return "R"
        else:
            return "python"
    # def get_hvg(
    #         self,
    #         wildcards: snakemake.io.Wildcards,
    #         output_pattern: str = None,
    #         **kwargs
    # ) -> str:
    #     """
    #     Get hvg parameter for integration run scripts
    #     :param wildcards: wildcards passed by Snakemake containing at least 'hvg' key
    #     :param output_pattern: file pattern with placeholders 'hvg'
    #         Only needed for R integration methods to get path of separate HVG file
    #         minimal example: 'output_dir/{hvg}.h5ad'
    #     :param kwargs: additional wildcards that are not contained in output_pattern
    #     :return: empty string for full-feature (n_hvg = 0), otherwise '-v <option>'
    #         with <option> specific to python or R methods
    #     """
    #     n_hvgs = self.get_feature_selection(wildcards.hvg)
    #     if n_hvgs == 0:
    #         return ""
    #     if output_pattern is not None:
    #         p = Path(snakemake.io.expand(output_pattern, **wildcards, **kwargs)[0])
    #         hvg_path = (p.parent / f'{p.stem}_hvg').with_suffix(p.suffix)
    #         return f'-v "{hvg_path}"'
    #     return f"-v {n_hvgs}"

    def get_from_scenario(self, scenario, key):
        if scenario not in self.DATA_SCENARIOS:
            raise ValueError(f"{scenario} not defined as scenario")
        if key not in self.DATA_SCENARIOS[scenario]:
            return False
        value = self.DATA_SCENARIOS[scenario][key]
        if key == 'resolution':
            return value if isinstance(value, list) else [value]
        return value

    def get_all_python_methods(self):
        return [
            method for method in self.METHODS
            if not self.get_from_method(method, "R")
        ]

    def get_all_R_methods(self):
        return [
            method for method in self.METHODS
            if self.get_from_method(method, "R")
        ]

    def get_all_wildcards(self, methods=None, feature_types=False): #agg_resolution=False
        """
        TODO: include method subsetting
        Collect all wildcards for wildcard-dependent rule
        :param methods: subset of methods, default: None, using all methods defined in config.
        :param feature_types: feature type or list of feature types to be considered.
            If feature_types == None, feature types set to default.
            Useful if a certain metric is examined on a specific feature type.
            Feature types are ["all_cell_peaks", "by_cluster_peaks", "tiles", "peaks"]
        :return: (comb_func, wildcards)
            comb_func: function for combining wildcards in snakemake.io.expand function
            wildcards: dictionary containing wildcards
        """
        wildcards = defaultdict(list)

        if methods is None:
            methods = self.METHODS

        if feature_types is False:   # if do not specify a feature type subset to enable, enable all feature types
            feature_types = ParsedConfig.FEATURE_TYPES
            # feature_types = "default"
        elif isinstance(feature_types, list):
            for ot in feature_types:
                if ot not in ParsedConfig.FEATURE_TYPES:
                    raise ValueError(f"{feature_types} not a valid feature type")

        comb_func = zip

        for scenario in self.get_all_scenarios():
            resolution = self.get_from_scenario(scenario, key = "resolution")

            # if agg_resolution:  # if this is necessary? Will snakemake automatically manage it?
            #     resolution = ["all"]
            cp_scenario = scenario # to avoid zipping scenario in the inner loop

            if resolution is False:
                resolution = [0.2]
            cp_resolution = resolution

            for method in methods:
                resolution = cp_resolution

                distance = self.get_from_method(method, key = "distance")
                tile_size = self.get_from_method(method, key = "tile_size")
                ndim = self.get_feature_selection(key = "ndim")
                
                if distance is False:
                    distance = ["default"]
                
                if tile_size is False:
                    tile_size = [0]

                def reshape_wildcards(*lists):
                    cart_prod = itertools.product(*lists)
                    return tuple(zip(*cart_prod))

                if isinstance(feature_types, list):
                    # feature type wildcard included
                    ot = self.get_from_method(method, "feature_type")

                    if not ot:
                        continue  # skip if method feature type is not defined in feature_types     
                    ot = set(feature_types).intersection(ot)

                    ot, method, scenario, distance, ndim, tile_size, resolution = reshape_wildcards(
                        ot,
                        [method],
                        [cp_scenario],
                        distance,
                        ndim,
                        tile_size,
                        resolution
                    )

                    wildcards["feature_type"].extend(ot)
                    wildcards["method"].extend(method)
                    wildcards["distance"].extend(distance)
                    wildcards["ndim"].extend(ndim)
                    wildcards["tile_size"].extend(tile_size)
                    wildcards["resolution"].extend(resolution)
                    wildcards["scenario"].extend(scenario)

                else:
                    ot = self.get_from_method(method, "feature_type")
                    ot, method, scenario, distance, ndim, tile_size, resolution = reshape_wildcards(
                        ot,
                        [method],
                        [cp_scenario],
                        distance,
                        ndim,
                        tile_size,
                        resolution
                    )

                    wildcards["feature_type"].extend(ot)
                    wildcards["method"].extend(method)
                    wildcards["distance"].extend(distance)
                    wildcards["ndim"].extend(ndim)
                    wildcards["tile_size"].extend(tile_size)
                    wildcards["resolution"].extend(resolution)
                    wildcards["scenario"].extend(scenario)

        # print(wildcards)
        return comb_func, wildcards

    # def get_integrated_for_metrics(self, rules, method):
    #     if method == "unintegrated":
    #         return Path(rules.integration_prepare.output[0]).with_suffix(".h5ad")
    #     elif self.get_from_method(method, "R"):
    #         return rules.convert_RDS_h5ad.output
    #     else:
    #         return rules.integration_run_python.output


    def get_peak_option_for_feature_engineering(self, wildcards):
        macs2_path = self.MACS2_PATH
        if self.get_from_method(wildcards.method, "peak_filtering"):
            min_width = self.get_from_scenario(wildcards.scenario, key="min_width")
            max_width = self.get_from_scenario(wildcards.scenario, key="max_width")
            return f"-z {macs2_path} -a {min_width} -b {max_width}"
        return f"-z {macs2_path}"


    def get_archr_option(self, wildcards):
        if wildcards.method == "ArchR":
            resolutions = self.get_from_scenario(wildcards.scenario, key="iter_resolution")
            return f"-r {resolutions}"
        return ""
    
    def get_agg_option(self, wildcards):
        if wildcards.method == "aggregation" or wildcards.method == "aggregation2" :
            n_meta_features = self.get_from_method(wildcards.method, key="n_meta_features")
            n_cells = self.get_from_method(wildcards.method, key="n_cells")
            norm_method = self.get_from_method(wildcards.method, key="norm_method")

            reduce = self.get_from_method(wildcards.method, key="reduce")
            feature_method = self.get_from_method(wildcards.method, key="feature_method")
            resolution = self.get_from_method(wildcards.method, key="resolution")
            cmd = f"-q {feature_method}"
            if reduce:
                cmd = cmd + f" -e {reduce}"
            if n_meta_features:
                cmd = cmd + f" -k {n_meta_features}"
            if n_cells:
                cmd = cmd + f" -u {n_cells}"
            if norm_method:
                cmd = cmd + f" -j {norm_method}"
            # if resolution:
            #     cmd = cmd + f" -v {resolution}"
            return cmd
        return ""