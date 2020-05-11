rule enrichr:
	input: "d2_enrichr/deseq2.{contrast}.full.txt"   # contrast = exp_vs_hem
	output:
		"d2_enrichr/{contrast}.enrichr_pos",
		"d2_enrichr/{contrast}.enrichr_neg",
		"d2_enrichr/{contrast}.enrichr_pos.lfc1",
		"d2_enrichr/{contrast}.enrichr_neg.lfc1"
	log:
		"logs/{contrast}.enrichr.log"
	conda:
		"../envs/my_enrichr.yaml"
	script:
		"../scripts/my_enrichr.R"

#META: enrichr
#1	Term
#2	Overlap
#3	P.value
#4	Adjusted.P.value
#5	Old.P.value
#6	Old.Adjusted.P.value
#7	Odds.Ratio
#8	Combined.Score
#9	Genes


rule enrichr_contrast_table:
	wildcard_constraints:
		contrast="exp_vs_hem|ko_vs_wt|exp_vs_wt"
	input:
		"d2_enrichr/{contrast}.enrichr_pos",
		"d2_enrichr/{contrast}.enrichr_neg",
		"d2_enrichr/{contrast}.enrichr_pos.lfc1",
		"d2_enrichr/{contrast}.enrichr_neg.lfc1"
	output:
		"d2_enrichr/{contrast}.enrichr_results.gz"
	shell:
		"(for i in {input} ; do " 
			" awk -F '\\t' -v i=$i '(NF>=1){{if(NF<2){{a=$0; next}}; print i FS a FS $0}}' FS='\\t' $i "  # skip empty lines and mimick fasta2tab
			#" | append_each_row -B $i "                        # or do this within the previous awk e.g.: # awk -v i=$i '/^>/{a=$0; next}{print i FS a FS $0}' FS='\t'
		"; done) "
		" | sed 's|\\t |\\t|g' "                               # one column from enrichr output starts with a space if the number is <10; I remove it
		" | sed 's|^d2_enrichr/{wildcards.contrast}.enrichr_|{wildcards.contrast}\t|' "
		" | awk -F '\\t' '{{$5=$5\"\\t\"$5; print}}' OFS='\\t' | awk -F '\\t' '{{gsub(\"/\",\"\\t\",$6); print}}' OFS='\\t' "
		" | gzip > {output} "

#META: enrichr_contrast_table  enrichr_contrast_table_top_res enrichr_contrast_table_top_res_multiDE
#1	contrast
#2	pos|neg|pos.lfc1|neg.lfc1
#3	Collection/Library
#4	Term
#5	Overlap
#6	overlap_numerator
#7	overlap_denominator
#8	P.value
#9	Adjusted.P.value
#10	Old.P.value
#11	Old.Adjusted.P.value
#12	Odds.Ratio
#13	Combined.Score
#14	Genes


# KEEP only enriched terms with at least 4 DE genes and with FDR<0.05
rule enrichr_contrast_table_top_res:
	wildcard_constraints:
		contrast="exp_vs_hem|ko_vs_wt|exp_vs_wt"
	input:
		"d2_enrichr/{contrast}.enrichr_results.gz"
	output:
		"d2_enrichr/{contrast}.enrichr_results.top.gz"
	shell:
		"zcat {input} | awk -F '\\t' '$9<0.05 && $6>3' OFS='\\t' | sort -t$'\t' -k9,9g -k6,6nr -k7,7n | gzip > {output} "


# TODO add last contrast if gender bias is acceptable
rule enrichr_contrast_table_top_res_multiDE:
	wildcard_constraints:
		contrast="exp_vs_hem|ko_vs_wt|exp_vs_wt"
	input:
		"d2_enrichr/exp_vs_hem.enrichr_results.top.gz",
		"d2_enrichr/exp_vs_wt.enrichr_results.top.gz"
	output:
		"d2_enrichr/top.enrichr_results.all_contrasts.gz"
	shell:
		"cat {input} > {output}"


# cat <(grep -A 14 "^#META.*enrichr_contrast_table_top_res_multiDE" rules/my_enrichment.smk | cut -f2 | unhead | tr "\n" "\t" |  sed 's|\t$|\n|') <(zcat d2_enrichr/top.enrichr_results.all_contrasts.gz) | tab2xls > d2_enrichr/top.enrichr_results.all_contrasts.xls 

########################
rule myex:
	input:
		"units.tsv"
	conda:
		"../envs/my_enrichr.yaml"
	shell:
		"R --version "
		#"awk -F '\t' '(NF>=1){{if(NF<2){{$1=\">\"$1}}; print}}' {input} | append_each_row -B pippo "
