import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.1.2")

input_folder = "/IMCR_shares/Moorlab/simo/170714_NB501465_0132_AH25MKBGX3/dhd1"

##### load config and sample sheets #####

configfile: "config.yaml"
validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
validate(samples, schema="schemas/samples.schema.yaml")


units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
validate(units, schema="schemas/units.schema.yaml")

# my_add, mostly 4 fastqcs
fq_files = pd.read_table(config["fq_files"]).set_index(["fastq_name_base"], drop=False)



##### target rules #####
# I actually use the pipeline until final_multiqc, then use R interactively to define outliers and compute DE
rule all:
    input:
        expand(["results/diffexp/{contrast}.diffexp.tsv",
                "results/diffexp/{contrast}.ma-plot.svg"],
               contrast=config["diffexp"]["contrasts"]),
        "results/pca.svg",
        "qc/final_multiqc_report.html"

rule my_all:
	input:
		"qc/final_multiqc_report.html",
		"feature_counts/sample.counts.matrix.txt"
	shell:
		"echo -e '\n\n~~~~~~~ For DESeq2 plots DEGs and enrichment, see scripts/my_deseq2_enrichment.R custom script ~~~~~~~\n\n'"

##### setup singularity #####

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"


##### setup report #####

report: "report/workflow.rst"


##### load rules #####

include: "rules/common.smk"
include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/diffexp.smk"
include: "rules/qc.smk"
include: "rules/my_enrichment.smk"


#### NOTES
# cat UNtrimmed/counts/all.tsv_with_gene_symbol | filter_1col -v 2 dups.gene_symbols | cut -f2- > star.counts.unique
