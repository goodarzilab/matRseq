#!/usr/bin/env Rscript
suppressMessages(library(tidyverse))
suppressMessages(library(DESeq2))
suppressMessages(library(ggpubr))
suppressMessages(library(rstatix))

logit.all<-function(x, design, model, n, feature.list){
  logit.all.gene<-function(m){
    prep<-data.frame(design,m)
    fit<-glm(model, data=prep, family="binomial"(link="logit"), weights = m)
    sfit<<-summary(fit)
    return(c(t(sfit$coefficients)))
  }
  logit.all.gene_try<-function(m){
    out<-tryCatch(
      {
        result_logistic<-logit.all.gene(m)
      },
      error=function(cond){
        return(rep(NA, times=4*n))
      }
    )
    return(out)
  }
  logit.x<-t(apply(x, 1, logit.all.gene_try))
  colnames(logit.x)<-apply(expand.grid(colnames(sfit$coefficients), rownames(sfit$coefficients)), 1, paste, collapse=".")
  rownames(logit.x)<-feature.list
  return(logit.x)
}

normalize.median.of.ratios <- function (m){
  gm_mean_z = function(x){
    exp(sum(log(x)) / length(x))
  }
  edf <- m
  id.names <- rownames(edf)
  geo.mean.vec <- apply(edf, 1, function(x) gm_mean_z(x))
  #gm_means_z is used because in the original paper, genes with gm_mean 0 are excluded when calculating median.
  ratios.df <- edf/geo.mean.vec # Division by 0 gm_mean will create NAs here.
  ratios.df <- as.matrix(ratios.df)
  ratios.df[is.infinite(ratios.df)] <- NA # in case division by 0 creates Inf
  normalization.factors <- apply(ratios.df, 2, function(x) median(x, na.rm=TRUE)) # NAs will be removed from calculation of median here.
  print("Normalization factors:")
  print(normalization.factors)
  normalized.edf <- t(t(edf)/normalization.factors)
  return(normalized.edf)
}

args = commandArgs(trailingOnly=TRUE)
metadata <- args[[1]]

#metadata <- '/rasis/Projects/Tools/matRseq/data/metadata.txt'
meta <- read_tsv(metadata)

indir <- dirname(metadata)
setwd(indir)

files <- paste0(meta[1] %>% pull(),".cnt")
print(files)

## Design formula
design <- args[[2]]

outfile <- args[[3]]

## Generate the count table
datalist <- lapply(files, function(x){read.table(file=x,header=FALSE,col.names=c("tRNA", sub(".cnt", "", x)))})
m <- Reduce(function(...) merge(..., by=1, all = TRUE), datalist)
rownames(m) <- m[,1]
m <- m[,-1]
m[is.na(m)]<-0
colnames(m) <- meta[1] %>% pull()
write.table(m,file=gsub(".txt","tRNA_raw_counts.txt", outfile), sep="\t", quote=FALSE, col.names=NA )

## Filter the count table
response <- strsplit(design,'~')[[1]][[1]]
covariate <- strsplit(design,'~')[[1]][[2]]
sample.types <- unique(meta[[response]])
print(sample.types)
set1 <- (meta[[response]] == sample.types[[1]])
set2 <- (meta[[response]] == sample.types[[2]])
f <- m[rowSums(m[,set1])>1 & rowSums(m[,set2])>1,]
print(sprintf("Count matrix was filtered down to %d rows (initially %d).", dim(f)[[1]], dim(m)[[1]]))

## Logit analysis
meta <- as.data.frame(meta)
meta[[response]] <- factor(meta[[response]])
meta[[covariate]] <- factor(meta[[covariate]])

df <- cbind(normalize.median.of.ratios(f[,set1]),normalize.median.of.ratios(f[,set2]))
write.table(df,file=gsub(".txt","tRNA_normalized_counts.txt", outfile), sep="\t", quote=FALSE, col.names=NA )

fit <- logit.all(df, as.data.frame(meta), formula(design), 2,rownames(df))

fit.s <- fit[,c(5,8)]
colnames(fit.s) <- c("logit_estimate", "pvalue")
fit.s <- cbind(fit.s, p.adjust(fit.s[,'pvalue'], method="fdr"))
colnames(fit.s) <- c("logit_estimate", "pvalue", "qvalue")

write.table(fit.s, outfile, quote = F, sep="\t", row.names = T)
head(fit.s)

fit.s <- data.frame(fit.s)
fit.s$t <- as.factor(abs(fit.s$logit_estimate)>0.5 & fit.s$qvalue<0.1)
g <- ggplot(data=as.data.frame(fit.s), aes(x=logit_estimate, y=-log(pvalue), colour=t), alpha=0.75, size=1) +
  geom_vline(aes(xintercept=0), colour="black",linetype="dashed", size=0.5)+
  xlab(expression(log[2]~'fold-change ribosomal to total tRNA')) + ylab(expression(-log[10]~paste(italic('p'),'-value')))+
  geom_point(alpha=0.75, size=1) +
  scale_colour_manual(values=c("#11008f", "#FF9100", "FF0000"))+theme_bw(30) +
  theme(panel.background = element_rect(colour = "black",size=2), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(legend.position="none") +
  theme(text = element_text(size=12)) + 
  xlim(c(-5, 5)) + ylim(c(0, 50))
ggsave(gsub(".txt","_volvano.pdf", outfile), width=6, height=7)


res <- fit.s
res$tRNA.family <-  gsub("(tRNA-\\S+-\\S+)-\\S+-\\d+$", "\\1", rownames(res))
res <- res %>% group_by(tRNA.family) %>% filter(n()>3) %>% ungroup()
stat.test <- res %>% group_by(tRNA.family) %>% wilcox_test(logit_estimate ~ 1, mu = 0, alternative = "two.sided", p.adjust.method = "none", conf.level = 0.9) %>% adjust_pvalue() %>% mutate(y.position = 5)
stat.test$star <- ifelse(stat.test$p<0.05, ifelse(stat.test$p<0.01,"**","*"), "")
res %>% ggplot(aes(x=reorder(tRNA.family, logit_estimate, FUN = median), y=logit_estimate)) + geom_hline(yintercept=0, linetype="dashed") + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + stat_pvalue_manual(stat.test, label = "star", xmin = "tRNA.family", xmax = NULL) +ylim(c(-5,5))
ggsave(gsub(".txt","_tRNA-family.pdf", outfile))
