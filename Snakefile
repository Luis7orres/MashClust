# ============================================================================
# MashClust Snakemake Pipeline
# ============================================================================
# Universal bacterial genome clustering using Mash distance-based approach
# 
# Usage:
#   snakemake -s Snakefile --cores 4 -c "qsub {cluster.qsub}"
#   snakemake -s Snakefile --cores 4 --use-conda
#   snakemake -s Snakefile --dry-run

configfile: "config.yaml"
import os

# Helper to get full paths
BASE_OUT = config["directories"]["output_base"]
SCRIPTS = config["directories"]["scripts_dir"]
LOGS = config["directories"].get("logs_dir", "logs")

# Ensure logs directory exists
os.makedirs(LOGS, exist_ok=True)

# Define output directories
DIR_DL = os.path.join(BASE_OUT, config["output_structure"]["download"])
DIR_SK = os.path.join(BASE_OUT, config["output_structure"]["sketch"])
DIR_DIST = os.path.join(BASE_OUT, config["output_structure"]["distances"])
DIR_CLUST = os.path.join(BASE_OUT, config["output_structure"]["cluster"])
DIR_VIS = os.path.join(BASE_OUT, config["output_structure"]["visualize"])
DIR_RES = os.path.join(BASE_OUT, config["results"]["base_dir"])

rule all:
    input:
        os.path.join(DIR_RES, config["results"]["accessions_file"]),
        os.path.join(DIR_VIS, "cluster_distribution.html")

rule download_genomes:
    input:
        accessions = config["params"]["accessions"]
    output:
        genomes_dir = directory(os.path.join(DIR_DL, "genomes")),
        genome_list = os.path.join(DIR_DL, "genome_list.txt")
    log:
        os.path.join(LOGS, "0-download.log")
    params:
        keep_zip = "--keep-zip" if config["params"]["download"]["keep_zip"] else "",
        batch_size = config["params"]["download"].get("batch_size", 50),
        delay = config["params"]["download"].get("delay", 45),
        max_retries = config["params"]["download"].get("max_retries", 5),
        api_key_arg = f"--api-key {config['params']['download']['api_key']}" if config["params"]["download"].get("use_api_key", False) else ""
    threads: 1
    shell:
        """
        python3 {SCRIPTS}/0-mashclust-download.py \
            {input.accessions} \
            -o {DIR_DL} \
            {params.keep_zip} \
            --batch-size {params.batch_size} \
            --delay {params.delay} \
            --max-retries {params.max_retries} \
            {params.api_key_arg} 2>&1 | tee {log}
        """

rule sketch_and_filter:
    input:
        genomes_dir = rules.download_genomes.output.genomes_dir
    output:
        sketch = os.path.join(DIR_SK, "genomes_sketch.msh"),
        target_list = os.path.join(DIR_SK, "targets.txt"),
        non_target_list = os.path.join(DIR_SK, "non_targets.txt")
    log:
        os.path.join(LOGS, "1-sketch.log")
    params:
        filter_arg = f"-f {config['params']['bacteria_name']}" if config["params"].get("filter_enabled", True) else "--no-filter",
        kmer = config["params"]["mash"]["kmer_size"],
        sketch_size = config["params"]["mash"]["sketch_size"]
    threads: config["params"]["threads"]
    shell:
        """
        python3 {SCRIPTS}/1-mashclust-sketch.py \
            {input.genomes_dir} \
            -o {DIR_SK} \
            {params.filter_arg} \
            -k {params.kmer} \
            -s {params.sketch_size} \
            --threads {threads} 2>&1 | tee {log}
        """

rule calculate_distances:
    input:
        sketch = rules.sketch_and_filter.output.sketch
    output:
        matrix = os.path.join(DIR_DIST, "distances.txt")
    log:
        os.path.join(LOGS, "2-distances.log")
    threads: config["params"]["threads"]
    shell:
        """
        python3 {SCRIPTS}/2-mashclust-distances.py \
            {DIR_SK} \
            -o {DIR_DIST} \
            --threads {threads} 2>&1 | tee {log}
        """

rule cluster_genomes:
    input:
        matrix = rules.calculate_distances.output.matrix
    output:
        representatives = os.path.join(DIR_CLUST, "representatives.txt"),
        clusters_json = os.path.join(DIR_CLUST, "clustering_data.json")
    log:
        os.path.join(LOGS, "3-cluster.log")
    params:
        threshold = config["params"]["id_threshold"],
        n_reps = config["params"]["n_representatives"],
        ref_arg = ("--reference " + " ".join(config["reference"]["accessions"])) if config["reference"]["enabled"] else ""
    shell:
        """
        python3 {SCRIPTS}/3-mashclust-cluster.py \
            {DIR_DIST} \
            -o {DIR_CLUST} \
            -t {params.threshold} \
            -n {params.n_reps} \
            {params.ref_arg} 2>&1 | tee {log}
        """

rule visualize_results:
    input:
        json_data = rules.cluster_genomes.output.clusters_json
    output:
        os.path.join(DIR_VIS, "cluster_distribution.html")
    log:
        os.path.join(LOGS, "4-visualize.log")
    shell:
        """
        python3 {SCRIPTS}/4-mashclust-visualize.py \
            {DIR_CLUST} \
            -o {DIR_VIS} 2>&1 | tee {log}
        """

rule finalize_results:
    input:
        non_targets = rules.sketch_and_filter.output.non_target_list,
        representatives = rules.cluster_genomes.output.representatives
    output:
        acc_file = os.path.join(DIR_RES, config["results"]["accessions_file"]),
        genomes_dir = directory(os.path.join(DIR_RES, config["results"]["genomes_dir"]))
    log:
        os.path.join(LOGS, "5-finalize.log")
    shell:
        """
        python3 {SCRIPTS}/5-mashclust-finalize.py \
            --non-targets {input.non_targets} \
            --representatives {input.representatives} \
            --out-dir {DIR_RES} \
            --genomes-subdir {config[results][genomes_dir]} \
            --acc-file {output.acc_file} \
            --dataset-manager /home/sandb8/scripts/pipelines/pipeline_1/resources/scripts/dataset-manager.py 2>&1 | tee {log}
        """