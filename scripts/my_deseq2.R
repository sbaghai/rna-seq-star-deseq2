library(DESeq2)     # in environment deseq2.1.18, see envs/my_deseq2.1.18.yml
library(ggplot2)
library(ggrepel)
#library(cowplot)
#library("vsn")
library(ggpubr)
##library(EnhancedVolcano)    # could not install it in deseq2.1.18 for compatibility issues and postpone it
#library(enrichR)   # set in other environment, enrichr2.1

# before all this
# conda activate deseq2.1.18

dir.create("/data/projects/tbx3/d2_enrichr")

setwd("/data/projects/tbx3/d2_enrichr")

all.counts=read.table("../feature_counts/sample.counts.matrix.txt", h=T, row.names=1, check.names=F)
all.meta=read.table("../local/share/data/metadata_dhd1_samples",h=T, check.names=F)
row.names(all.meta)=all.meta$sample
dds<-DESeqDataSetFromMatrix(countData=all.counts, colData=all.meta, design=~condition)
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
plotPCA(vsd, intgroup="condition")
plotPCA(rld, intgroup="condition")
#plotCounts(dds, gene="Ctnnb1", ingroup="condition")

# PCA table - would be better done on rpkms or on normalized counts extracted by dds, e.g. from vsd, see:
#edit(DESeq2:::plotPCA.DESeqTransform) ---> 'assay(vsd)'
#project.pca <- prcomp(t(all.counts)) # use instead -----> project.pca <- prcomp(t(assay(vsd)))
#summary(project.pca)
#project.pca.proportionvariances <- ((project.pca$sdev^2) / (sum(project.pca$sdev^2)))*100
#barplot(project.pca.proportionvariances, cex.names=1, xlab=paste("Principal component (PC), 1-", length(project.pca$sdev)), ylab="Proportion of variation (%)", main="Scree plot", ylim=c(0,100))
#par(cex=1.0, cex.axis=0.8, cex.main=0.8)
#pairs(project.pca$x[,1:5], col="black", main="Principal components analysis bi-plot\nPCs 1-5", pch=16)
#pairs(project.pca$x[,1:5], col=all.meta$condition, main="Principal components analysis bi-plot\nPCs 1-5")
#pairs(project.pca$x[,6:10], col="black", main="Principal components analysis bi-plot\nPCs 6-10", pch=16)

#edit(DESeq2:::plotPCA.DESeqTransform)
DESeq2::plotPCA(vsd, intgroup="condition")
#DESeq2::plotPCA(vsd, intgroup="condition") + geom_text(aes(label=all.meta$sample),hjust=0.25, vjust=-0.5, show_guide = F)   # possibly even better with geom_text_repel 

## s1: i.e. subset1 
counts.s1=all.counts[,!(colnames(all.counts) %in% "wt_dhd1_17-3d")]
meta.s1=all.meta[!(row.names(all.meta) %in% "wt_dhd1_17-3d"),]
dds.s1<-DESeqDataSetFromMatrix(countData=counts.s1, colData=meta.s1, design=~condition)
vsd.s1 <- vst(dds.s1, blind=FALSE)
plotPCA(vsd.s1, intgroup="condition")
#vsd.s1b <- vst(dds.s1, blind=TRUE)
#plotPCA(vsd.s1b, intgroup="condition")
# compare if there are neat changes using blind or not on groups of experimental design
#p.blind<-plotPCA(vsd.s1b, intgroup="condition")
#p.notblind<-plotPCA(vsd.s1, intgroup="condition")
#library(cowplot)
#plot_grid(p.blind, p.notblind)
#pdf("plotPCA_no.outlier_vst_blind_vs_not.pdf", height=4, width=10)
#plot_grid(p.blind, p.notblind)
#dev.off()
pdf("plotCounts_std_ctrl.pdf",width=8, height=10)
par(mfrow = c(3, 2)) 
plotCounts(dds.s1, gene="Actb", ingroup="condition")
plotCounts(dds.s1, gene="Gapdh", ingroup="condition")
plotCounts(dds.s1, gene="Gusb", ingroup="condition")
plotCounts(dds.s1, gene="B2m", ingroup="condition")
plotCounts(dds.s1, gene="Rpl27", ingroup="condition")
plotCounts(dds.s1, gene="Rps13", ingroup="condition")
plotCounts(dds.s1, gene="Xist", ingroup="condition")
plotCounts(dds.s1, gene="Ddx3y", ingroup="condition")
plotCounts(dds.s1, gene="Eif2s3y", ingroup="condition")
plotCounts(dds.s1, gene="Eif2s3x", ingroup="condition")
plotCounts(dds.s1, gene="Rps27a", ingroup="condition")
plotCounts(dds.s1, gene="Vim", ingroup="condition")
dev.off()

pdf("plotCounts_prj_ctrl.pdf",width=8, height=10)
#plotCounts(dds.s1, gene="Tbx3", ingroup="condition")
#plotCounts(dds.s1, gene="Ctnnb1", ingroup="condition")
#plotCounts(dds.s1, gene="Pygo1", ingroup="condition")
#plotCounts(dds.s1, gene="Bcl9", ingroup="condition")
dev.off()

### DEGs
# DE exp vs hem; "eh"= subset for exp vs hem (dhd1)
meta.eh=meta.s1[meta.s1$condition == "exp_dhd1" | meta.s1$condition == "hem_dhd1",]
counts.eh=counts.s1[,rownames(meta.eh)]
dds.eh=DESeqDataSetFromMatrix(countData=counts.eh, colData=meta.eh, design=~condition)
# VST proved to be the best normalization, so I kee using that one:
vsd.eh <- vst(dds.eh, blind=FALSE)
pdf("pca_exp_vs_hem.pdf", width=5, height=4)
plotPCA(vsd.eh, intgroup="condition") + geom_text_repel(aes(label=colnames(counts.eh)), size=2, color="darkgray")
dev.off()


deseq.eh=DESeq(dds.eh)
res.exp_vs_hem=results(deseq.eh, contrast=c("condition", "exp_dhd1", "hem_dhd1"))
#res.exp_vs_hem[order(res.exp_vs_hem$log2FoldChange, decreasing=T),]
write.table(res.exp_vs_hem[order(res.exp_vs_hem$log2FoldChange, decreasing=T),], sep="\t", quote=F,"deseq2.exp_vs_hem.full.txt")
#bawk '$3!="NA"' deseq2.exp_vs_hem.full.txt | sed 's|^baseMean|GeneID\tbaseMean|' | tab2xls > deseq2.exp_vs_hem.xls

pdf("plotMA_exp_vs_hem.pdf", width=5, height=5)
DESeq2::plotMA(res.exp_vs_hem, ylim=c(-8,8))
ggpubr::ggmaplot(res.exp_vs_hem, main = expression("dhd1: exp vs hem (FDR<0.05 and |logFC|>1)"), fdr = 0.05, fc = 2, size = 0.4, palette = c("#B31B21", "#1465AC", "darkgray"), genenames = as.vector(rownames(res.exp_vs_hem)), legend = "top", font.label = c("plain", 11), font.legend = "bold", font.main = "plain", ggtheme = ggplot2::theme_minimal(), top=20, select.top.method="padj")
dev.off()


## DE exp vs wt
meta.ew=meta.s1[meta.s1$condition == "exp_dhd1" | meta.s1$condition == "wt_dhd1",]
counts.ew=counts.s1[,rownames(meta.ew)]
dds.ew=DESeqDataSetFromMatrix(countData=counts.ew, colData=meta.ew, design=~condition)
vsd.ew <- vst(dds.ew, blind=FALSE)
pdf("pca_exp_vs_wt.pdf", width=5, height=4)
plotPCA(vsd.ew, intgroup="condition") + geom_text_repel(aes(label=colnames(counts.ew)), size=2, color="darkgray")
dev.off()

deseq.ew=DESeq(dds.ew)
res.exp_vs_wt=results(deseq.ew, contrast=c("condition", "exp_dhd1", "wt_dhd1"))
write.table(res.exp_vs_wt[order(res.exp_vs_wt$log2FoldChange, decreasing=T),], sep="\t", quote=F,"deseq2.exp_vs_wt.full.txt")
#bawk '$3!="NA"' deseq2.exp_vs_wt.full.txt | sed 's|^baseMean|GeneID\tbaseMean|' | tab2xls > deseq2.exp_vs_wt.xls


pdf("plotMA_exp_vs_wt.pdf", width=5, height=5)
DESeq2::plotMA(res.exp_vs_wt, ylim=c(-8,8))
ggpubr::ggmaplot(res.exp_vs_wt, main = expression("dhd1: exp vs wt (FDR<0.05 and |logFC|>1)"), fdr = 0.05, fc = 2, size = 0.4, palette = c("#B31B21", "#1465AC", "darkgray"), genenames = as.vector(rownames(res.exp_vs_wt)), legend = "top", font.label = c("plain", 11), font.legend = "bold", font.main = "plain", ggtheme = ggplot2::theme_minimal(), top=20, select.top.method="padj")
dev.off()

####
# there's a male wt, see:
# library(ggrepel)
dy<-plotCounts(dds.ew, gene="Ddx3y", ingroup="condition",returnData=T)
dy$label<-rownames(dy)
dx<-plotCounts(dds.ew, gene="Xist", ingroup="condition",returnData=T)
dx$label<-rownames(dx)
py <- ggplot(dy, aes(x = condition, y = count, color = condition)) + geom_point(position=position_jitter(w = 0.1,h = 0), size=4) + geom_text_repel(data=subset(dy, count > 50), aes(label=label), color="black") + theme_bw() + ggtitle("Ddx3y") + theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
px <- ggplot(dx, aes(x = condition, y = count, color = condition)) + geom_point(position=position_jitter(w = 0.1,h = 0), size=4) + geom_text_repel(data=subset(dx, count < 50), aes(label=label), color="black") + theme_bw() + ggtitle("Xist") + theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
pdf("exp_vs_wt.Ddx3y_Xist.plotCounts.pdf", height=2, width=5)
#par(mfrow=c(1,2))   # if using a 3x7 pdf
#plotCounts(dds.ew, gene="Ddx3y", ingroup="condition")
#plotCounts(dds.ew, gene="Xist", ingroup="condition")
#par(mfrow=c(1,1))
plot_grid(py,px, ncol=2)
dev.off()
######


exit()
## see scripts/my_enrichr.R
##### conda activate enrichr2.1
library(enrichR)

setwd("/data/projects/tbx3/d2_enrichr")
all_res=read.table("deseq2.exp_vs_hem.full.txt", check.names=F, row.names=1, h=T)

#contrast=EXPvsHEM
my_FDR=0.05

# FDR signif, no cutoff on logFC
pos=rownames(all_res[which(all_res$padj<my_FDR & all_res$log2FoldChange>0),])
neg=rownames(all_res[which(all_res$padj<my_FDR & all_res$log2FoldChange<0),])

# more stringent
pos.lfc1=rownames(all_res[which(all_res$padj<my_FDR & all_res$log2FoldChange > 1),])
neg.lfc1=rownames(all_res[which(all_res$padj<my_FDR & all_res$log2FoldChange < -1),])

enrichr.db=enrichR::listEnrichrDbs()

# edit enrichR::printEnrichr function so that it prints all results, not just top 10 foreach library, and all columns, not just 2,3,6:
my_printEnrich <-function (data, file, sep = "\t", columns = c(1:9)) { enrich <- file(file, "w"); for (i in 1:length(data)) { writeLines(names(data)[i], enrich); n <- nrow(data[[i]]); if (n > 0) { writeLines(paste(apply(data[[i]][1:n, columns, drop = FALSE], 1, function(x) paste(x, collapse = sep)), collapse = "\n"), enrich); writeLines("\n", enrich)}}; close(enrich)}


e.pos=enrichR::enrichr(genes=pos, databases=enrichr.db$libraryName)
my_printEnrich(e.pos, file="enrichr_pos")

e.neg=enrichR::enrichr(genes=neg, databases=enrichr.db$libraryName)
my_printEnrich(e.neg, file="enrichr_neg")

e.pos.lfc1=enrichR::enrichr(genes=pos.lfc1, databases=enrichr.db$libraryName)
my_printEnrich(e.pos.lfc1, file="enrichr_pos.lfc1")

e.neg.lfc1=enrichR::enrichr(genes=neg.lfc1, databases=enrichr.db$libraryName)
my_printEnrich(e.pos.lfc1, file="enrichr_neg.lfc1")

### repeat, with custom signatures


