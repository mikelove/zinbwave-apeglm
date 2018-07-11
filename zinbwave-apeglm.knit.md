---
title: "ZINB-WaVE + apeglm for single-cell RNA-seq"
author: "Michael Love"
output: html_document
---

This document closely follows a similar
one, [zinbwave-deseq2](https://github.com/mikelove/zinbwave-deseq2),
showing integration of *zinbwave* with *DESeq2*. The difference is
that here we additionally assess when *apeglm* effect size estimation
may provide a benefit for zero-inflated Negative Binomial count data,
as may be produced in a single cell RNA-seq experiment (scRNA-seq).
*apeglm* effect sizes appear useful in the case that the groups of
cells we want to compare have small sample size, e.g. here the two
groups being compared have ~20 cells per group. It outperforms the
Normal distribution based shrinkage estimator when the ratio of DE
genes is low.

We use the *splatter* package to simulate single-cell RNA-seq data.

* Zappia, Phipson, and Oshlack "Splatter: simulation of single-cell RNA
sequencing data" *Genome Biology* (2017)
[doi: 10.1186/s13059-017-1305-0](https://doi.org/10.1186/s13059-017-1305-0)

We then use the methods defined in the following paper to combine
*zinbwave* observation weights with *DESeq2* modeling of negative
binomial counts. Finally we use *apeglm* to estimate LFC between two
groups that have few cells.

* Van den Berge & Perraudeau *et al* "Observation weights unlock bulk
RNA-seq tools for zero inflation and single-cell applications" *Genome Biology* (2018)
[doi: 10.1186/s13059-018-1406-4](https://doi.org/10.1186/s13059-018-1406-4)

> It is important to note that while methods such as ZINB-WaVE and
> ZINGER can successfully identify excess zeros, they cannot, however,
> readily discriminate between their underlying causes, i.e., between
> technical (e.g., dropout) and biological (e.g., bursting) zeros. 

We will try two approaches to integrating with *apeglm*:

1. Using *zinbwave* to define weights to account for the excess zeros,
   then using the Negative Binomial likelihood with *apeglm*.
2. Using *zinbwave* to define weights, then using these as a
   parameter for the probability of a zero coming from the zero
   component in a Zero-Inflated Negative Binomial likelihood with
   *apeglm*.
   
These two approaches end up having nearly identical results, while the
first is already built in to *DESeq2*'s *lfcShrink* function and has
been optimized for speed.



We construct simulated single cell data using *splatter*. The custom
parameters here define 1% of genes to identify each sub-group (so
comparing two groups we would expect 2% of genes to be DE), and the
non-null LFCs are set to be centered at a log2 fold change of 2. The
low percent of DE genes and relatively high LFC values help show the
difference between the Normal- and t-distributed priors in *DESeq2* and
*apeglm*. We specify three groups, where two of the groups will have
20 cells on average, and we will see how to compare the LFC between
these two groups.


```r
suppressPackageStartupMessages(library(splatter))
params <- newSplatParams()
# LFC centered at 2
params <- setParam(params, "de.facLoc", log(2) * 2)
params <- setParam(params, "de.facScale", log(2) * 1)
params <- setParam(params, "de.prob", .01)
# include drop-out 
params <- setParam(params, "dropout.type", "experiment")
params <- setParam(params, "dropout.mid", 2)
# three groups, we will show apeglm improvements on the two small groups
sim <- splatSimulate(params, group.prob=c(.2, .2, .6),
                     method="groups", seed=1)
```

```
## Getting parameters...
```

```
## Creating simulation object...
```

```
## Simulating library sizes...
```

```
## Simulating gene means...
```

```
## Simulating group DE...
```

```
## Simulating cell means...
```

```
## Simulating BCV...
```

```
## Simulating counts..
```

```
## Simulating dropout (if needed)...
```

```
## Done!
```

```r
table(sim$Group)
```

```
## 
## Group1 Group2 Group3 
##     17     26     57
```

```r
# define the true LFC
rowData(sim)$log2FC <- with(rowData(sim), log2(DEFacGroup2/DEFacGroup1))
```

As in the 
[zinbwave-deseq2](https://github.com/mikelove/zinbwave-deseq2)
repo, we run *zinbwave* to estimate the excess zeros and to define a
weight matrix for these.


```r
library(zinbwave)
keep <- rowSums(counts(sim) >= 10) >= 5
table(keep)
```

```
## keep
## FALSE  TRUE 
##  7879  2121
```

```r
zinb <- sim[keep,]
zinb$condition <- factor(zinb$Group)
nms <- c("counts", setdiff(assayNames(zinb), "counts"))
assays(zinb) <- assays(zinb)[nms]
# epsilon as recommended from Van den Berge and Perraudeau
zinb <- zinbwave(zinb, K=0, BPPARAM=SerialParam(), epsilon=1e12)
```

We additionally run *DESeq2* to estimate dispersions and maximum
likelihood estimates of log2 fold change. This again is the
recommended code from the
[zinbwave-deseq2](https://github.com/mikelove/zinbwave-deseq2) repo.


```r
suppressPackageStartupMessages(library(DESeq2))
dds <- DESeqDataSet(zinb, design=~condition)
# arguments as recommended from Van den Berge and Perraudeau
dds <- DESeq(dds, test="LRT", reduced=~1,
             sfType="poscounts", minmu=1e-6, minRep=Inf)
```

```
## estimating size factors
```

```
## estimating dispersions
```

```
## gene-wise dispersion estimates
```

```
## mean-dispersion relationship
```

```
## final dispersion estimates
```

```
## fitting model and testing
```

We plot the dispersion estimates. If this plot fails, see the
suggested code in the
[zinbwave-deseq2](https://github.com/mikelove/zinbwave-deseq2) repo.


```r
plotDispEsts(dds)
```

<img src="zinbwave-apeglm_files/figure-html/plotdisp-1.png" width="672" />

We extract and compute a number of estimators for the LFC. 

First, we compute two simple pseudocount based LFC estimates, one without
using the excess zero weights, and the other using the weights to
account for excess zeros as identified by *zinbwave*.


```r
ncts <- counts(dds, normalized=TRUE)
idx1 <- dds$condition == "Group1"
idx2 <- dds$condition == "Group2"
pc <- .1
simple.lfc <- log2(rowMeans(ncts[,idx2])+pc) -
  log2(rowMeans(ncts[,idx1])+pc)
wtd.lfc <- log2(rowSums(ncts[,idx2]*wts[,idx2])/rowSums(wts[,idx2])+pc) -
  log2(rowSums(ncts[,idx1]*wts[,idx1])/rowSums(wts[,idx1])+pc)
```

The maximum likelihood estimate from *DESeq2*. All of the methods
using *DESeq2* will take into account the excess zero weights.


```r
res <- results(dds, name="condition_Group2_vs_Group1",
               independentFiltering=FALSE)
```

The "old" shrunken LFC from the original *DESeq2* paper, using the
Normal distribution:


```r
res.norm <- lfcShrink(dds, coef=2, type="normal")
```

The new shrunken LFC available in *DESeq2*, using the *apeglm* method
and a Negative Binomial likelihood, with excess zeros accounted for by
weights. In this case, the *apeglm* estimator is 10x faster to compute
than the normal distributed prior estimator above.


```r
ape.nb <- lfcShrink(dds, coef=2, type="apeglm")
```

```
## using 'apeglm' for LFC shrinkage. If used in published research, please cite:
##     Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for
##     sequence count data: removing the noise and preserving large differences.
##     bioRxiv. https://doi.org/10.1101/303255
```

Finally, we will create a likelihood function for the Zero-Inflated
Negative Binomial, and set up a number of parameters to pass to the
*apeglm* function itself. This approach ends up with nearly identical
results to the use of *lfcShrink* to call *apeglm* above, while
requiring additional code and slower than the above to compute.


```r
library(apeglm)
library(ZIM)
Y <- counts(dds)
design <- model.matrix(design(dds), data=colData(dds))
disps <- dispersions(dds)
wts <- assays(dds)[["weights"]]
# combine dispersion and wts into a parameter matrix,
# which will be passed to apeglm
param <- cbind(disps, 1 - wts)
offset <- matrix(log(sizeFactors(dds)), nrow=nrow(dds),
                 ncol=ncol(dds), byrow=TRUE)
# need to put to natural log scale for apeglm
mle <- log(2) * cbind(res$log2FoldChange, res$lfcSE)
logLikZINB <- function (y, x, beta, param, offset) {
    xbeta <- x %*% beta + offset
    mean.hat <- exp(xbeta)
    k <- 1/param[1]
    omega <- param[-1]
    ZIM::dzinb(y, k=k, lambda=mean.hat, omega=omega, log=TRUE)
}
# run apeglm with a ZINB likelihood and zinbwave weights
# used to define the probability of an excess zero
fit <- apeglm(Y=Y, x=design, log.lik=logLikZINB, param=param,
              coef=2, mle=mle, offset=offset)
# need to put back to log2 scale
ape.zinb <- log2(exp(1)) * fit$map[,2]
```

Finally, we construct plots to compare our methods. We compute the
correlation over all genes, the median absolute error (MAE) of the top
30 genes, and over all genes.


```r
myplot <- function(x,y,n=30,...) {
  plot(x,y,ylim=c(-6,6),xlab="",ylab="",cex=.5,...)
  idx <- rank(-abs(y)) < n
  points(x[idx], y[idx], col="blue", cex=2)
  abline(0,1,col="red")
  lgd <- c(paste("cor =", round(cor(x,y,use="complete"),3)),
           paste("MAE top =", round(mean(abs(x[idx]-y[idx]), na.rm=TRUE),2)),
           paste("MAE =", round(mean(abs(x-y), na.rm=TRUE),2)))
  legend("topleft", legend=lgd)
}
```



```r
par(mfrow=c(2,3), mar=c(2,3,2,1))
myplot(mcols(dds)$log2FC, simple.lfc, main="pseudocount")
myplot(mcols(dds)$log2FC, wtd.lfc, main="wtd pseudocount")
myplot(mcols(dds)$log2FC, res$log2FoldChange, main="MLE")
myplot(mcols(dds)$log2FC, res.norm$log2FoldChange, main="normal + wtd NB")
myplot(mcols(dds)$log2FC, ape.nb$log2FoldChange, main="apeglm + wtd NB")
myplot(mcols(dds)$log2FC, ape.zinb, main="apeglm + ZINB lik")
```

<img src="zinbwave-apeglm_files/figure-html/plotlfc-1.png" width="864" />


```r
session_info()
```

```
## Session info -------------------------------------------------------------
```

```
##  setting  value                       
##  version  R version 3.5.0 (2018-04-23)
##  system   x86_64, darwin15.6.0        
##  ui       X11                         
##  language (EN)                        
##  collate  en_US.UTF-8                 
##  tz       Europe/Rome                 
##  date     2018-07-11
```

```
## Packages -----------------------------------------------------------------
```

```
##  package              * version   date       source         
##  acepack                1.4.1     2016-10-29 CRAN (R 3.5.0) 
##  ADGofTest              0.3       2011-12-28 CRAN (R 3.5.0) 
##  annotate               1.58.0    2018-05-01 Bioconductor   
##  AnnotationDbi          1.42.1    2018-05-08 Bioconductor   
##  apeglm               * 1.2.0     2018-05-01 Bioconductor   
##  assertthat             0.2.0     2017-04-11 CRAN (R 3.5.0) 
##  backports              1.1.2     2017-12-13 cran (@1.1.2)  
##  base                 * 3.5.0     2018-04-24 local          
##  base64enc              0.1-3     2015-07-28 CRAN (R 3.5.0) 
##  bbmle                  1.0.20    2017-10-30 CRAN (R 3.5.0) 
##  bindr                  0.1.1     2018-03-13 CRAN (R 3.5.0) 
##  bindrcpp               0.2.2     2018-03-29 CRAN (R 3.5.0) 
##  Biobase              * 2.40.0    2018-05-01 Bioconductor   
##  BiocGenerics         * 0.26.0    2018-05-01 Bioconductor   
##  BiocInstaller        * 1.30.0    2018-05-04 Bioconductor   
##  BiocParallel         * 1.14.1    2018-05-06 Bioconductor   
##  bit                    1.1-14    2018-05-29 CRAN (R 3.5.0) 
##  bit64                  0.9-7     2017-05-08 CRAN (R 3.5.0) 
##  bitops                 1.0-6     2013-08-17 CRAN (R 3.5.0) 
##  blob                   1.1.1     2018-03-25 CRAN (R 3.5.0) 
##  checkmate              1.8.5     2017-10-24 CRAN (R 3.5.0) 
##  cluster                2.0.7-1   2018-04-13 CRAN (R 3.5.0) 
##  coda                   0.19-1    2016-12-08 CRAN (R 3.5.0) 
##  codetools              0.2-15    2016-10-05 CRAN (R 3.5.0) 
##  colorspace             1.3-2     2016-12-14 CRAN (R 3.5.0) 
##  compiler               3.5.0     2018-04-24 local          
##  copula                 0.999-18  2017-09-01 CRAN (R 3.5.0) 
##  data.table             1.11.4    2018-05-27 CRAN (R 3.5.0) 
##  datasets             * 3.5.0     2018-04-24 local          
##  DBI                    1.0.0     2018-05-02 CRAN (R 3.5.0) 
##  DelayedArray         * 0.6.1     2018-06-15 Bioconductor   
##  DESeq2               * 1.20.0    2018-05-01 Bioconductor   
##  devtools             * 1.13.6    2018-06-27 CRAN (R 3.5.0) 
##  digest                 0.6.15    2018-01-28 cran (@0.6.15) 
##  dplyr                  0.7.5     2018-05-19 cran (@0.7.5)  
##  edgeR                  3.22.3    2018-06-21 Bioconductor   
##  emdbook                1.3.9     2016-02-11 CRAN (R 3.5.0) 
##  evaluate               0.10.1    2017-06-24 CRAN (R 3.5.0) 
##  foreach                1.4.4     2017-12-12 CRAN (R 3.5.0) 
##  foreign                0.8-70    2017-11-28 CRAN (R 3.5.0) 
##  Formula                1.2-3     2018-05-03 CRAN (R 3.5.0) 
##  genefilter             1.62.0    2018-05-01 Bioconductor   
##  geneplotter            1.58.0    2018-05-01 Bioconductor   
##  GenomeInfoDb         * 1.16.0    2018-05-01 Bioconductor   
##  GenomeInfoDbData       1.1.0     2018-01-10 Bioconductor   
##  GenomicRanges        * 1.32.3    2018-05-16 Bioconductor   
##  ggplot2                3.0.0     2018-07-03 CRAN (R 3.5.0) 
##  glmnet                 2.0-16    2018-04-02 CRAN (R 3.5.0) 
##  glue                   1.2.0     2017-10-29 CRAN (R 3.5.0) 
##  graphics             * 3.5.0     2018-04-24 local          
##  grDevices            * 3.5.0     2018-04-24 local          
##  grid                   3.5.0     2018-04-24 local          
##  gridExtra              2.3       2017-09-09 CRAN (R 3.5.0) 
##  gsl                    1.9-10.3  2017-01-05 CRAN (R 3.5.0) 
##  gtable                 0.2.0     2016-02-26 CRAN (R 3.5.0) 
##  Hmisc                  4.1-1     2018-01-03 CRAN (R 3.5.0) 
##  htmlTable              1.12      2018-05-26 CRAN (R 3.5.0) 
##  htmltools              0.3.6     2017-04-28 CRAN (R 3.5.0) 
##  htmlwidgets            1.2       2018-04-19 CRAN (R 3.5.0) 
##  IRanges              * 2.14.10   2018-05-16 Bioconductor   
##  iterators              1.0.9     2017-12-12 CRAN (R 3.5.0) 
##  knitr                  1.20      2018-02-20 CRAN (R 3.5.0) 
##  lattice                0.20-35   2017-03-25 CRAN (R 3.5.0) 
##  latticeExtra           0.6-28    2016-02-09 CRAN (R 3.5.0) 
##  lazyeval               0.2.1     2017-10-29 CRAN (R 3.5.0) 
##  limma                  3.36.2    2018-06-21 Bioconductor   
##  locfit                 1.5-9.1   2013-04-20 CRAN (R 3.5.0) 
##  magrittr               1.5       2014-11-22 CRAN (R 3.5.0) 
##  MASS                 * 7.3-50    2018-04-30 CRAN (R 3.5.0) 
##  Matrix                 1.2-14    2018-04-13 CRAN (R 3.5.0) 
##  matrixStats          * 0.53.1    2018-02-11 CRAN (R 3.5.0) 
##  memoise                1.1.0     2017-04-21 CRAN (R 3.5.0) 
##  methods              * 3.5.0     2018-04-24 local          
##  munsell                0.5.0     2018-06-12 CRAN (R 3.5.0) 
##  mvtnorm                1.0-7     2018-01-26 CRAN (R 3.5.0) 
##  nnet                   7.3-12    2016-02-02 CRAN (R 3.5.0) 
##  numDeriv               2016.8-1  2016-08-27 CRAN (R 3.5.0) 
##  parallel             * 3.5.0     2018-04-24 local          
##  pcaPP                  1.9-73    2018-01-14 CRAN (R 3.5.0) 
##  pillar                 1.2.3     2018-05-25 CRAN (R 3.5.0) 
##  pkgconfig              2.0.1     2017-03-21 CRAN (R 3.5.0) 
##  plyr                   1.8.4     2016-06-08 CRAN (R 3.5.0) 
##  pspline                1.0-18    2017-06-12 CRAN (R 3.5.0) 
##  purrr                  0.2.5     2018-05-29 cran (@0.2.5)  
##  R6                     2.2.2     2017-06-17 CRAN (R 3.5.0) 
##  RColorBrewer           1.1-2     2014-12-07 CRAN (R 3.5.0) 
##  Rcpp                   0.12.17   2018-05-18 cran (@0.12.17)
##  RCurl                  1.95-4.10 2018-01-04 CRAN (R 3.5.0) 
##  rlang                  0.2.1     2018-05-30 cran (@0.2.1)  
##  rmarkdown            * 1.9       2018-03-01 CRAN (R 3.5.0) 
##  rpart                  4.1-13    2018-02-23 CRAN (R 3.5.0) 
##  rprojroot              1.3-2     2018-01-03 cran (@1.3-2)  
##  RSQLite                2.1.1     2018-05-06 CRAN (R 3.5.0) 
##  rstudioapi             0.7       2017-09-07 CRAN (R 3.5.0) 
##  S4Vectors            * 0.18.3    2018-06-08 Bioconductor   
##  scales                 0.5.0     2017-08-24 CRAN (R 3.5.0) 
##  SingleCellExperiment * 1.2.0     2018-05-01 Bioconductor   
##  softImpute             1.4       2015-04-08 CRAN (R 3.5.0) 
##  splatter             * 1.4.0     2018-05-01 Bioconductor   
##  splines                3.5.0     2018-04-24 local          
##  stabledist             0.7-1     2016-09-12 CRAN (R 3.5.0) 
##  stats                * 3.5.0     2018-04-24 local          
##  stats4               * 3.5.0     2018-04-24 local          
##  stringi                1.2.3     2018-06-12 CRAN (R 3.5.0) 
##  stringr                1.3.1     2018-05-10 CRAN (R 3.5.0) 
##  SummarizedExperiment * 1.10.1    2018-05-11 Bioconductor   
##  survival               2.42-3    2018-04-16 CRAN (R 3.5.0) 
##  testthat             * 2.0.0     2017-12-13 CRAN (R 3.5.0) 
##  tibble                 1.4.2     2018-01-22 CRAN (R 3.5.0) 
##  tidyselect             0.2.4     2018-02-26 CRAN (R 3.5.0) 
##  tools                  3.5.0     2018-04-24 local          
##  utils                * 3.5.0     2018-04-24 local          
##  withr                  2.1.2     2018-03-15 CRAN (R 3.5.0) 
##  XML                    3.98-1.11 2018-04-16 CRAN (R 3.5.0) 
##  xtable                 1.8-2     2016-02-05 CRAN (R 3.5.0) 
##  XVector                0.20.0    2018-05-01 Bioconductor   
##  yaml                   2.1.19    2018-05-01 CRAN (R 3.5.0) 
##  ZIM                  * 1.0.3     2017-02-06 CRAN (R 3.5.0) 
##  zinbwave             * 1.2.0     2018-05-01 Bioconductor   
##  zlibbioc               1.26.0    2018-05-01 Bioconductor
```
