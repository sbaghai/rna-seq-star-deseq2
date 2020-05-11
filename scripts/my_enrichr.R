log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(enrichR)

#dir.create("/data/projects/tbx3/d2_enrichr")
#setwd("/data/projects/tbx3/d2_enrichr")

all_res=read.table(snakemake@input[[1]], check.names=F, row.names=1, h=T)


my_FDR=0.05

pos=rownames(all_res[which(all_res$padj<my_FDR & all_res$log2FoldChange>0),])
neg=rownames(all_res[which(all_res$padj<my_FDR & all_res$log2FoldChange<0),])

# more stringent
pos.lfc1=rownames(all_res[which(all_res$padj<my_FDR & all_res$log2FoldChange > 1),])
neg.lfc1=rownames(all_res[which(all_res$padj<my_FDR & all_res$log2FoldChange < -1),])

enrichr.db=enrichR::listEnrichrDbs()

# edit enrichR::printEnrichr function so that it prints all results, not just top 10 foreach library, and all columns, not just 2,3,6:
my_printEnrich <-function (data, file, sep = "\t", columns = c(1:9)) { enrich <- file(file, "w"); for (i in 1:length(data)) { writeLines(names(data)[i], enrich); n <- nrow(data[[i]]); if (n > 0) { writeLines(paste(apply(data[[i]][1:n, columns, drop = FALSE], 1, function(x) paste(x, collapse = sep)), collapse = "\n"), enrich); writeLines("\n", enrich)}}; close(enrich)}


e.pos=enrichR::enrichr(genes=pos, databases=enrichr.db$libraryName)
my_printEnrich(e.pos, file=snakemake@output[[1]])
#snakemake@output[[paste0(snakemake@wildcards.contrast,".enrichr_pos")]]

e.neg=enrichR::enrichr(genes=neg, databases=enrichr.db$libraryName)
my_printEnrich(e.neg, file=snakemake@output[[2]])

e.pos.lfc1=enrichR::enrichr(genes=pos.lfc1, databases=enrichr.db$libraryName)
my_printEnrich(e.pos.lfc1, file=snakemake@output[[3]])

e.neg.lfc1=enrichR::enrichr(genes=neg.lfc1, databases=enrichr.db$libraryName)
my_printEnrich(e.neg.lfc1, file=snakemake@output[[4]])

