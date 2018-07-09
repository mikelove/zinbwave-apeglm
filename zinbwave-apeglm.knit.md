---
title: "ZINB-WaVE + apeglm integration for single-cell RNA-seq"
author: "Michael Love"
output: html_document
---

Here we use the *splatter* package to simulate single-cell RNA-seq
data.

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
sim <- splatSimulate(params, group.prob=c(.15, .15, .7),
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
##     13     19     68
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
##  7845  2155
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
```

```
## Loading required package: MASS
```

```r
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
    legend=c(paste("cor =",
      round(cor(x,y,use="complete"),3)),
      paste("MAE =",
      round(mean(abs(x[idx]-y[idx]), na.rm=TRUE),3))))
}
par(mfrow=c(2,2), mar=c(2,3,2,1))
myplot(mcols(dds)$log2FC, lfc1$log2FoldChange, main="normal")
myplot(mcols(dds)$log2FC, lfc2$log2FoldChange, main="apeglm wt NB")
myplot(mcols(dds)$log2FC, lfc3, main="apeglm ZINB")
myplot(mcols(dds)$log2FC, simple.lfc, main="pseudocount")
```

<img src="zinbwave-apeglm_files/figure-html/plotlfc-1.png" width="672" />