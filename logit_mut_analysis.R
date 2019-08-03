#!/usr/bin/env Rscript
if("tidyverse" %in% rownames(installed.packages()) == FALSE) {install.packages("tidyverse",repos = "http://cran.us.r-project.org")}
suppressMessages(library(tidyverse))
suppressMessages(library(DESeq2))
suppressMessages(library(ggpubr))
if("rstatix" %in% rownames(installed.packages()) == FALSE) {install.packages("rstatix",repos = "http://cran.us.r-project.org")}
suppressMessages(library(rstatix))

logit.all<-function(x, design, model, n, feature.list){
  logit.all.gene<-function(m){
    prep<-data.frame(design,m)
    fit<-glm(model, data=prep, family="binomial"(link="logit"), weights = m)
    sfit<<-summary(fit)
    if(length(c(t(sfit$coefficients)))<8){
      print(prep)
      print(c(t(sfit$coefficients)))
    }
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
  #colnames(logit.x)<-apply(expand.grid(colnames(sfit$coefficients), rownames(sfit$coefficients)), 1, paste, collapse=".")
  #rownames(logit.x)<-feature.list
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
ountfile <- args[[1]]
#countfile <- '/rasis/Projects/Tools/matRseq/data/pileup.parsed.txt'


metadata <- args[[2]]
meta <- read_tsv(metadata)
meta.rD <- meta
meta.mD <- meta
meta.rD[[1]] <- paste0(meta.rD[[1]], ".rDA")
meta.rD$status <- "rDA"
meta.mD[[1]] <- paste0(meta.mD[[1]], ".mDA")
meta.mD$status <- "mDA"
meta <- bind_rows(meta.rD, meta.mD)

indir <- dirname(metadata)
setwd(indir)

## variable
covariate <- args[[3]]


outfile <- args[[4]]


## Load the count table
m <- read.table(countfile, sep="\t", head=T, row.names=1)
m <- m[meta[[1]]]

## Filter the count table
#mutation count ~ sample.type (e.g. ribo vs. total)
design<-paste('status~',covariate)
response <- 'status'
set1 <- (meta$status == 'rDA')
set2 <- (meta$status == 'mDA')
levels <- unique(meta[[covariate]])
set3 <- (meta[[covariate]] == levels[1])
set4 <- (meta[[covariate]] == levels[2])
tmp<-m
tmp[tmp>0] <- 1
f <- m[rowSums(tmp[,set1])>2 & rowSums(tmp[,set2])>2 & rowSums(tmp[,set3])>2 & rowSums(tmp[,set4])>2,]
print(sprintf("Count matrix was filtered down to %d rows (initially %d).", dim(f)[[1]], dim(m)[[1]]))

## Logit analysis
meta <- as.data.frame(meta)
meta[[response]] <- factor(meta[[response]])
meta[[response]] <- relevel(meta[[response]], ref='rDA')
meta[[covariate]] <- factor(meta[[covariate]])

df <- cbind(normalize.median.of.ratios(f[,set1]),normalize.median.of.ratios(f[,set2]))
write.table(df,file=gsub(".txt","_mut_normalized_counts.txt", outfile), sep="\t", quote=FALSE, col.names=NA )

fit <- logit.all(df, as.data.frame(meta), formula(design), 2,rownames(df))

fit.s <- fit[,c(5,8)]
colnames(fit.s) <- c("logit_estimate", "pvalue")
fit.s <- cbind(fit.s, p.adjust(fit.s[,'pvalue'], method="fdr"))
colnames(fit.s) <- c("logit_estimate", "pvalue", "qvalue")

write.table(fit.s, outfile, quote = F, sep="\t", row.names = T)
fit.sig <- fit.s[fit.s[,'qvalue']<0.1,]
print(sprintf("%d significant mutational signatures identified.", dim(fit.sig)[[1]]))
write.table(fit.sig[order(fit.sig[,'pvalue']),], file=gsub(".txt", ".sig.txt", outfile), quote = F, sep="\t", row.names = T, col.names=NA )

dir.create("plots", showWarnings = FALSE)
setwd('/plots')

for (t in rownames(fit.sig)){
  meta$test <- df[t,]
  meta %>% ggplot(aes(x=status,y=test, color=!!as.symbol(covariate), fill=!!as.symbol(covariate)), size=1)+geom_boxplot(width=0.5, alpha=0.25)+geom_point(position = position_dodge(width=0.5))+
    theme_bw(30) +  theme(panel.background = element_rect(colour = "black",size=2), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +theme(text = element_text(size=12))
  ggsave(gsub(".txt",paste0("_",t,".pdf"), outfile), device = "pdf" )}
