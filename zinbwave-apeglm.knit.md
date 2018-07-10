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
groups being compared have ~20 cells per group.

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

suppressPackageStartupMessages(library(DESeq2))
dds <- DESeqDataSet(zinb, design=~condition)
# stabilize the weights to avoid errors
wts <- assays(dds)[["weights"]]
wts[wts < 1e-6] <- 1e-6
assays(dds)[["weights"]] <- wts

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



```r
plotDispEsts(dds)
```

<img src="zinbwave-apeglm_files/figure-html/plotdisp-1.png" width="672" />


```r
res <- results(dds, name="condition_Group2_vs_Group1",
               independentFiltering=FALSE)

# compare to "old" shrunken LFC
lfc1 <- lfcShrink(dds, coef=2, type="normal")

# the apeglm method using weights from zinbwave
lfc2 <- lfcShrink(dds, coef=2, type="apeglm")
```

```
## using 'apeglm' for LFC shrinkage. If used in published research, please cite:
##     Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for
##     sequence count data: removing the noise and preserving large differences.
##     bioRxiv. https://doi.org/10.1101/303255
```

```r
library(apeglm)
library(ZIM)

Y <- counts(dds)
design <- model.matrix(design(dds), data=colData(dds))
disps <- dispersions(dds)
wts <- assays(dds)[["weights"]]
param <- cbind(disps, 1 - wts)
offset <- matrix(log(sizeFactors(dds)), nrow=nrow(dds), ncol=ncol(dds), byrow=TRUE)
mle <- log(2) * cbind(res$log2FoldChange, res$lfcSE)

logLikZINB <- function (y, x, beta, param, offset) {
    xbeta <- x %*% beta + offset
    mean.hat <- exp(xbeta)
    k <- 1/param[1]
    omega <- param[-1]
    ZIM::dzinb(y, k=k, lambda=mean.hat, omega=omega, log=TRUE)
}

# run apeglm with a ZINB likelihood and zinbwave weights
# used to define the probability of a zero
fit <- apeglm(Y=Y, x=design, log.lik=logLikZINB, param=param,
              coef=2, mle=mle, offset=offset)
lfc3 <- log2(exp(1)) * fit$map[,2]

# simple LFC using pseudocount
ncts <- counts(dds, normalized=TRUE)
simple.lfc <- log2(rowMeans(ncts[,dds$condition == "Group2"]) + .1) -
              log2(rowMeans(ncts[,dds$condition == "Group1"]) + .1)
```



```r
myplot <- function(x,y,n=20,...) {
  plot(x,y,ylim=c(-6,6),xlab="",ylab="",cex=.5,...)
  idx <- rank(-abs(y)) < n
  points(x[idx], y[idx], col="blue", cex=2)
  abline(0,1,col="red")
  legend("topleft",
      legend=c(
      paste("cor =",
            round(cor(x,y,use="complete"),3)),
      paste("MAE top =",
            round(mean(abs(x[idx]-y[idx]), na.rm=TRUE),3)),
      paste("MAE =",
            round(mean(abs(x-y), na.rm=TRUE),3))))
}
par(mfrow=c(2,2), mar=c(2,3,2,1))
myplot(mcols(dds)$log2FC, lfc1$log2FoldChange, main="normal")
myplot(mcols(dds)$log2FC, simple.lfc, main="pseudocount")
myplot(mcols(dds)$log2FC, lfc2$log2FoldChange, main="apeglm + weight NB")
myplot(mcols(dds)$log2FC, lfc3, main="apeglm + ZINB lik")
```

<img src="zinbwave-apeglm_files/figure-html/plotlfc-1.png" width="672" />


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
##  date     2018-07-10
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
