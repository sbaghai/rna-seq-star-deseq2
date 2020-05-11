###### COUNTS

def get_strandness(units):
    if "strandedness" in units.columns:
        return units["strandedness"].tolist()
    else:
        strand_list=["none"]
        return strand_list*units.shape[0]

rule count_matrix:
    input:
        expand("star/{unit.sample}-{unit.unit}/ReadsPerGene.out.tab", unit=units.itertuples())
    output:
        "counts/all.tsv"
    params:
        samples=units["sample"].tolist(),
        strand=get_strandness(units)
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/count-matrix.py"

rule counts_add_symbols:
	input:
		"counts/all.tsv"
	output:
		"counts/all.tsv_with_gene_symbol"
	shell:
		"cat {input} | translate -a -v -e gene_symbol /data/references/gencode/mmusculus/ens2gene_symbol 1 > {output}"

rule subread_fc:
	input:
		"star/{sample}-{unit}/Aligned.out.bam"
	output:
		"feature_counts/{sample}-{unit}.counts",
                "feature_counts/{sample}-{unit}.counts.summary"
	log:
		"feature_counts/{sample}-{unit}.counts.log"
	params:
		threads=4,
		# TODO add a way to substitute 0, yes and reverse with 0,1,2; and also in shell bit, replace that and
		strand=2,
		# strand=get_strandness(units)
		# strand=strand.replace("none","0")
		# strand=strand.replace("yes","1")
		# strand=strand.replace("reverse","2")
		if_paired="-p",
		# TODO conditional also for PE/SE
		annotation=config["ref"]["annotation"]
	conda:
		"../envs/my_featcounts.yaml"
	shell:
		"featureCounts -T {params.threads} {params.if_paired} -C -s {params.strand} -t exon -F GTF -a {params.annotation} -g 'gene_name' -o {output} {input}  2> {log}"

# params not recognised in this format
#		shell("featureCounts "
#		"-T {params.threads} "
#		"{params.if_paired} "
#		"-C "
#		"-s {params.strand} "
#		"-t exon "
#		"-F GTF "
#		"-a {snakemake.params.annotation} "
#		"-g 'gene_name' "
#		"-o {output} "
#		"{input}  2> {log}")

rule fc_matrix:
	input:
		expand("feature_counts/{unit.sample}-{unit.unit}.counts", unit=units.itertuples())
	params:
		samples=units["sample"].tolist()
	output:
		"feature_counts/sample.counts.matrix.txt"
	run:
	# TODO eventually do this in py with pandas of numpy!
	#shell: "(for i in {input}; do cat $i | append_each_row -B $i; done) | tab2matrix -r 'gene_name' > {output} && touch {output}"
		import pandas as pd
		counts = [pd.read_table(f, index_col=0, usecols=[0, 6], header=0, comment="#")
		for f in input]
		for t, sample in zip(counts, params.samples):
    			t.columns = [sample]
		matrix = pd.concat(counts, axis=1)
		matrix.index.name = "geneID"
		matrix.to_csv(output[0], sep="\t")


############# DE ## did not use this, see below

def get_deseq2_threads(wildcards=None):
    # https://twitter.com/mikelove/status/918770188568363008
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6


rule deseq2_init:
    input:
        counts="counts/all.tsv"
    output:
        "deseq2/all.rds"
    params:
        samples=config["samples"]
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/init.log"
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2-init.R"


rule pca:
    input:
        "deseq2/all.rds"
    output:
        report("results/pca.svg", "../report/pca.rst")
    params:
        pca_labels=config["pca"]["labels"]
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/pca.log"
    script:
        "../scripts/plot-pca.R"


def get_contrast(wildcards):
    return config["diffexp"]["contrasts"][wildcards.contrast]


rule deseq2:
    input:
        "deseq2/all.rds"
    output:
        table=report("results/diffexp/{contrast}.diffexp.tsv", "../report/diffexp.rst"),
        ma_plot=report("results/diffexp/{contrast}.ma-plot.svg", "../report/ma.rst"),
    params:
        contrast=get_contrast
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/{contrast}.diffexp.log"
    threads: get_deseq2_threads
    script:
        "../scripts/deseq2.R"

