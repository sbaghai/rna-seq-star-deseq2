def get_fq(wildcards):
    if config["trimming"]["skip"]:
        # no trimming, use raw reads
        return units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    else:
        # yes trimming, use trimmed data
        if not is_single_end(**wildcards):
            # paired-end sample
            return expand("trimmed/{sample}-{unit}.{group}.fastq.gz",
                          group=[1, 2], **wildcards)
        # single end sample
        return "trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)
            

rule align:
    input:
        sample=get_fq
    output:
        # see STAR manual for additional output files
        "star/{sample}-{unit}/Aligned.out.bam",
        "star/{sample}-{unit}/ReadsPerGene.out.tab"
    log:
        "logs/star/{sample}-{unit}.log"
    params:
        # path to STAR reference genome index
        index=config["ref"]["index"],
        # optional parameters
        extra="--quantMode GeneCounts --sjdbGTFfile {} {}".format(
              config["ref"]["annotation"], config["params"]["star"])
    threads: 1
    conda:
        "../envs/my_star.yaml"
    wrapper:
        "0.19.4/bio/star/align"


rule test_align:
    input:
        "pippo"
    output:
        "pluto"
    shell:
        "STAR --version > {output} "


rule all_align:
	input:
		expand("star/{unit.sample}-{unit.unit}/Aligned.out.bam", unit=units.itertuples()),
		expand("star/{unit.sample}-{unit.unit}/ReadsPerGene.out.tab", unit=units.itertuples())
	output:
		"all_align.done"
	shell:
		"touch {output}"
