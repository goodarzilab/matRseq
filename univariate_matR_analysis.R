#!/usr/bin/env Rscript
suppressMessages(library(tidyverse))
suppressMessages(library(DESeq2))
suppressMessages(library(ggpubr))
if("rstatix" %in% rownames(installed.packages()) == FALSE) {install.packages("rstatix",repos = "http://cran.us.r-project.org")}
suppressMessages(library(rstatix))

args = commandArgs(trailingOnly=TRUE)
metadata <- args[[1]]
#metadata <- '/rasis/Projects/Tools/matRseq/data/metadata_univariate.txt'

meta <- read_tsv(metadata)

indir <- dirname(metadata)
setwd(indir)

files <- paste0(meta[1] %>% pull(),".cnt")
print(files)

## Design formula
design <- args[[2]]
#design <- "~cell.line"
ref <- args[[3]]
#ref <- "SW480-Par"

outfile <- args[[4]]

## Generate the count table
datalist <- lapply(files, function(x){read.table(file=x,header=FALSE,col.names=c("tRNA", sub(".cnt", "", x)))})
m <- Reduce(function(...) merge(..., by=1, all = TRUE), datalist)
rownames(m) <- m[,1]
m <- m[,-1]
m[is.na(m)]<-0
colnames(m) <- meta[1] %>% pull()
write.table(m,file=gsub(".txt","_tRNA_raw_counts.txt", outfile), sep="\t", quote=FALSE, col.names=NA )

## Filter the count table
f <- m[rowSums(m)>2,]
print(sprintf("Count matrix was filtered down to %d rows (initially %d).", dim(f)[[1]], dim(m)[[1]]))

## DESeq2 analysis
meta <- as.data.frame(meta)
covariate <- strsplit(design,'~')[[1]][[2]]
meta[[covariate]] <- factor(meta[[covariate]])
meta[[covariate]] <- relevel(meta[[covariate]], ref=ref)

dds <- DESeqDataSetFromMatrix(countData = m, colData = meta, design = formula("~cell.line"))
dds <- DESeq(dds)
res <- results(dds)

ncu <- counts(dds, normalized=TRUE)
write.table(ncu, file=gsub(".txt","_tRNA_DE2norm.txt", outfile), sep="\t", quote=FALSE, col.names=NA )

write.table(res, outfile, quote = F, sep="\t", row.names = T)

fit.s <- as.data.frame(res)
fit.s$t <- as.factor(abs(fit.s$log2FoldChange)>0.5 & fit.s$padj<0.1)
g <- ggplot(data=as.data.frame(fit.s), aes(x=log2FoldChange, y=-log(pvalue), colour=t), alpha=0.75, size=1) +
  geom_vline(aes(xintercept=0), colour="black",linetype="dashed", size=0.5)+
  xlab(expression(log[2]~'fold-change ribosomal to total tRNA')) + ylab(expression(-log[10]~paste(italic('p'),'-value')))+
  geom_point(alpha=0.75, size=1) +
  scale_colour_manual(values=c("#11008f", "#FF9100", "FF0000"))+theme_bw(30) +
  theme(panel.background = element_rect(colour = "black",size=2), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(legend.position="none") +
  theme(text = element_text(size=12))
  #xlim(c(-5, 5)) + ylim(c(0, 50))
ggsave(gsub(".txt","_volvano.pdf", outfile), width=6, height=7)

res <- as.data.frame(res)
res$tRNA.family <-  gsub("(tRNA-\\S+-\\S+)-\\S+-\\d+$", "\\1", rownames(res))
res <- res %>% group_by(tRNA.family) %>% filter(n()>3) %>% ungroup()
stat.test <- res %>% group_by(tRNA.family) %>% wilcox_test(log2FoldChange ~ 1, mu = 0, alternative = "two.sided", p.adjust.method = "none", conf.level = 0.9) %>% adjust_pvalue() %>% mutate(y.position = 5)
stat.test$star <- ifelse(stat.test$p<0.05, ifelse(stat.test$p<0.01,"**","*"), "")
res %>% ggplot(aes(x=reorder(tRNA.family, log2FoldChange, FUN = median), y=log2FoldChange)) + geom_hline(yintercept=0, linetype="dashed") + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + stat_pvalue_manual(stat.test, label = "star", xmin = "tRNA.family", xmax = NULL) #+ylim(c(-5,5))
ggsave(gsub(".txt","_tRNA-family.pdf", outfile))
