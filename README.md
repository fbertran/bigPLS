<!-- README.md is generated from README.Rmd. Please edit that file -->




# bigPLS <img src="man/figures/logo.png" align="right" width="200"/>

# bigPLS, PLS models and their extension for big data in R
## Frédéric Bertrand and Myriam Maumy-Bertrand

<!-- badges: start -->
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/fbertran/bigPLS/workflows/R-CMD-check/badge.svg)](https://github.com/fbertran/bigPLS/actions)
[![Codecov test coverage](https://codecov.io/gh/fbertran/bigPLS/branch/master/graph/badge.svg)](https://codecov.io/gh/fbertran/bigPLS?branch=master)
[![CRAN status](https://www.r-pkg.org/badges/version/bigPLS)](https://cran.r-project.org/package=bigPLS)
[![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/bigPLS)](https://cran.r-project.org/package=bigPLS)
[![GitHub Repo stars](https://img.shields.io/github/stars/fbertran/bigPLS?style=social)](https://github.com/fbertran/bigPLS)
[![DOI](https://zenodo.org/badge/18454102.svg)](https://zenodo.org/badge/latestdoi/18454102)
<!-- badges: end -->


The goal of bigPLS is provide Cox models in a high dimensional setting in R.



Support for parallel computation and GPU is being developed.


This website and these examples were created by F. Bertrand and M. Maumy-Bertrand.

## Installation

You can install the released version of bigPLS from [CRAN](https://CRAN.R-project.org) with:


```r
install.packages("bigPLS")
```

You can install the development version of bigPLS from [github](https://github.com) with:


```r
devtools::install_github("fbertran/bigPLS")
```

# Example

## Allelotyping real dataset

### The dataset


```r
library(plsRcox)
data(micro.censure)
Y_train_micro <- micro.censure$survyear[1:80]
C_train_micro <- micro.censure$DC[1:80]
Y_test_micro <- micro.censure$survyear[81:117]
C_test_micro <- micro.censure$DC[81:117]

data(Xmicro.censure_compl_imp)
X_train_micro <- apply((as.matrix(Xmicro.censure_compl_imp)),FUN="as.numeric",MARGIN=2)[1:80,]
X_train_micro_df <- data.frame(X_train_micro)

X_train_micro_orig <- Xmicro.censure_compl_imp[1:80,]
X_train_micro_orig_df <- data.frame(X_train_micro_orig)
```

### Compute deviance residuals with some options.


```r
head(computeDR(Y_train_micro,C_train_micro,plot=TRUE))
```

## Simulated data

### Generate dataset


```r
set.seed(4669)
library(bigPLS)
x_sim <- matrix(sample(0:1, size = 20000, replace = TRUE), ncol = 2)
dCox_sim <- dataCox(10^4, lambda = 3, rho = 2, x_sim,
beta = c(1,3), cens.rate = 5)
data(dCox_sim)
```

### Compute deviance residuals with some options.


```r
with(dCox_sim,head(computeDR(time,status,plot=TRUE)))
```

### Model Matrix


```r
coxgpls(~.,Y_train_micro,C_train_micro,ncomp=6,trace=TRUE,model_matrix=TRUE,dataXplan = X_train_micro_orig_df,ind.block.x=c(3,10,20))
```

### coxgpls


```r
(cox_gpls_fit=coxgpls(X_train_micro,Y_train_micro,C_train_micro,ncomp=6,ind.block.x=c(3,10,15)))

(cox_gpls_fit2=coxgpls(~X_train_micro,Y_train_micro,C_train_micro,ncomp=6,ind.block.x=c(3,10,15)))

(cox_gpls_fit3=coxgpls(~.,Y_train_micro,C_train_micro,ncomp=6,
dataXplan=X_train_micro_df,ind.block.x=c(3,10,15)))

rm(cox_gpls_fit,cox_gpls_fit2,cox_gpls_fit3)
```

### cv.coxgpls

```r
set.seed(123456)
(cv.coxgpls.res=cv.coxgpls(list(x=X_train_micro,time=Y_train_micro, status=C_train_micro),nt=10,ind.block.x=c(3,10,15)))
```


### coxgplsDR


```r
(cox_gplsDR_fit=coxgplsDR(X_train_micro,Y_train_micro,C_train_micro,ncomp=6,ind.block.x=c(3,10,15)))

(cox_gplsDR_fit2=coxgplsDR(~X_train_micro,Y_train_micro,C_train_micro,ncomp=6,ind.block.x=c(3,10,15)))

(cox_gplsDR_fit3=coxgplsDR(~.,Y_train_micro,C_train_micro,ncomp=6,
dataXplan=X_train_micro_df,ind.block.x=c(3,10,15)))

rm(cox_gplsDR_fit,cox_gplsDR_fit2,cox_gplsDR_fit3)
```

### cv.coxgplsDR


```r
set.seed(123456)

(cv.coxsplsDR.res=cv.coxgplsDR(list(x=X_train_micro,time=Y_train_micro, status=C_train_micro),nt=10,ind.block.x=c(3,10,15)))
```


### coxDKgplsDR


```r
(cox_DKsplsDR_fit=coxDKgplsDR(X_train_micro,Y_train_micro,C_train_micro,ncomp=6, validation="CV",ind.block.x=c(3,10,15),verbose=TRUE))

(cox_DKsplsDR_fit=coxDKgplsDR(~X_train_micro,Y_train_micro,C_train_micro,ncomp=6, validation="CV",ind.block.x=c(3,10,15)))

(cox_DKsplsDR_fit=coxDKgplsDR(~.,Y_train_micro,C_train_micro,ncomp=6,
validation="CV",dataXplan=data.frame(X_train_micro),ind.block.x=c(3,10,15)))

rm(cox_DKsplsDR_fit)
```

### cv.coxDKgPLSDR


```r
set.seed(123456)

(cv.coxDKgplsDR.res=cv.coxDKgplsDR(list(x=X_train_micro,time=Y_train_micro, status=C_train_micro),nt=10,ind.block.x=c(3,10,15)))
```

### coxsgpls

```r
(cox_sgpls_fit=coxsgpls(X_train_micro,Y_train_micro,C_train_micro,ncomp=6,ind.block.x=c(3,10,15), alpha.x = rep(0.95, 6)))

(cox_sgpls_fit2=coxsgpls(~X_train_micro,Y_train_micro,C_train_micro,ncomp=6,ind.block.x=c(3,10,15), alpha.x = rep(0.95, 6)))

(cox_sgpls_fit3=coxsgpls(~.,Y_train_micro,C_train_micro,ncomp=6,
dataXplan=X_train_micro_df,ind.block.x=c(3,10,15), alpha.x = rep(0.95, 6)))

rm(cox_sgpls_fit,cox_sgpls_fit2,cox_sgpls_fit3)
```

### cv.coxsgpls

```r
set.seed(123456)
(cv.coxsgpls.res=cv.coxsgpls(list(x=X_train_micro,time=Y_train_micro, status=C_train_micro),nt=10,ind.block.x=c(3,10,15), alpha.x = rep(0.95, 10)))
```

### coxsgplsDR


```r
(cox_sgplsDR_fit=coxsgplsDR(X_train_micro,Y_train_micro,C_train_micro,ncomp=6,ind.block.x=c(3,10,15), alpha.x = rep(0.95, 10)))

(cox_sgplsDR_fit2=coxsgplsDR(~X_train_micro,Y_train_micro,C_train_micro,ncomp=6,ind.block.x=c(3,10,15), alpha.x = rep(0.95, 10)))

(cox_sgplsDR_fit3=coxsgplsDR(~.,Y_train_micro,C_train_micro,ncomp=6,
                           dataXplan=X_train_micro_df,ind.block.x=c(3,10,15), alpha.x = rep(0.95, 10)))

rm(cox_sgplsDR_fit,cox_sgplsDR_fit2,cox_sgplsDR_fit3)
```

### cv.coxsgplsDR


```r
set.seed(4669)

(cv.coxsgplsDR.res=cv.coxsgplsDR(list(x=X_train_micro,time=Y_train_micro, status=C_train_micro),nt=10,ind.block.x=c(3,10,15), alpha.x = rep(0.95, 10)))
```


### coxDKsgplsDR


```r
(cox_DKsgplsDR_fit=coxDKsgplsDR(X_train_micro,Y_train_micro,C_train_micro,ncomp=6, validation="CV",ind.block.x=c(3,10,15), alpha.x = rep(0.95, 10),verbose=TRUE))

(cox_DKsgplsDR_fit=coxDKsgplsDR(~X_train_micro,Y_train_micro,C_train_micro,ncomp=6, validation="CV",ind.block.x=c(3,10,15), alpha.x = rep(0.95, 10)))

(cox_DKsgplsDR_fit=coxDKsgplsDR(~.,Y_train_micro,C_train_micro,ncomp=6,
                              validation="CV",dataXplan=data.frame(X_train_micro),ind.block.x=c(3,10,15), alpha.x = rep(0.95, 10)))

rm(cox_DKsgplsDR_fit)
```

### cv.coxDKgplsDR


```r
set.seed(123456)

(cv.coxDKgplsDR.res=cv.coxDKgplsDR(list(x=X_train_micro,time=Y_train_micro, status=C_train_micro),nt=10,ind.block.x=c(3,10,15), alpha.x = rep(0.95, 10)))
```




### coxspls_sgpls

```r
(cox_coxspls_sgpls_fit=coxspls_sgpls(X_train_micro,Y_train_micro,C_train_micro,ncomp=6,keepX = rep(5, 6)))

(cox_coxspls_sgpls_fit2=coxspls_sgpls(~X_train_micro,Y_train_micro,C_train_micro,ncomp=6,ind.keepX = rep(5, 6)))

(cox_coxspls_sgpls_fit3=coxspls_sgpls(~.,Y_train_micro,C_train_micro,ncomp=6,
dataXplan=X_train_micro_df,keepX = rep(5, 6)))

rm(cox_coxspls_sgpls_fit,cox_coxspls_sgpls_fit2,cox_coxspls_sgpls_fit3)
```

### cv.coxspls_sgpls

```r
set.seed(123456)
(cv.coxspls_sgpls.res=cv.coxspls_sgpls(list(x=X_train_micro,time=Y_train_micro, status=C_train_micro),nt=10, keepX = rep(5, 10)))
```


### coxspls_sgplsDR


```r
(cox_spls_sgplsDR_fit=coxspls_sgplsDR(X_train_micro,Y_train_micro,C_train_micro,ncomp=6,ind.block.x=c(3,10,15)))

(cox_spls_sgplsDR_fit2=coxspls_sgplsDR(~X_train_micro,Y_train_micro,C_train_micro,ncomp=6,ind.block.x=c(3,10,15)))

(cox_spls_sgplsDR_fit3=coxspls_sgplsDR(~.,Y_train_micro,C_train_micro,ncomp=6,
                           dataXplan=X_train_micro_df,ind.block.x=c(3,10,15)))

rm(cox_spls_sgplsDR_fit,cox_spls_sgplsDR_fit2,cox_spls_sgplsDR_fit3)
```

### cv.coxspls_sgplsDR


```r
set.seed(123456)

(cv.coxspls_sgplsDR.res=cv.coxspls_sgplsDR(list(x=X_train_micro,time=Y_train_micro, status=C_train_micro),nt=10,ind.block.x=c(3,10,15)))
```


### coxDKspls_sgplsDR


```r
(cox_DKspls_sgplsDR_fit=coxDKspls_sgplsDR(X_train_micro,Y_train_micro,C_train_micro,ncomp=6, validation="CV",ind.block.x=c(3,10,15),verbose=TRUE))

(cox_DKspls_sgplsDR_fit=coxDKspls_sgplsDR(~X_train_micro,Y_train_micro,C_train_micro,ncomp=6, validation="CV",ind.block.x=c(3,10,15)))

(cox_DKspls_sgplsDR_fit=coxDKspls_sgplsDR(~.,Y_train_micro,C_train_micro,ncomp=6,
                              validation="CV",dataXplan=data.frame(X_train_micro),ind.block.x=c(3,10,15)))

rm(cox_DKspls_sgplsDR_fit)
```

### cv.coxDKspls_sgplsDR


```r
set.seed(123456)

(cv.coxDKspls_sgplsDR.res=cv.coxDKspls_sgplsDR(list(x=X_train_micro,time=Y_train_micro, status=C_train_micro),nt=10,ind.block.x=c(3,10,15)))
```




```r
data(survData, package="bigSurvSGD")
resultsBigscale <- bigscale(formula=Surv(time, status)~.,data=survData, parallel.flag=TRUE, num.cores=2)
resultsBigscale

# Simulated survival data to be read off the memory
data(survData) # a dataset with 1000 observations (rows) and 10 features (columns)
# Save dataset survSGD as bigSurvSGD to be read chunk-by-chunk off the memory 
write.csv(survData, file.path(tempdir(), "bigSurvData.csv"), row.names = FALSE) 
dataPath <- file.path(tempdir(), "bigSurvData.csv") # path to where data is
resultsBigscaleOffMemory <- bigscale(formula=Surv(time, status)~., data=dataPath, 
bigmemory.flag=TRUE, parallel.flag=TRUE, num.cores=2)
resultsBigscaleOffMemory


# Simulated sparse survival data
data(sparseSurvData,package="bigSurvSGD") # a sparse data with 100 observations (rows) and 150 features (columns)
resultsBigscaleSparse <- bigscale(formula=Surv(time, status)~.,data=sparseSurvData, parallel.flag=TRUE, num.cores=2)
resultsBigscaleSparse



# Simulated survival data - just estimation and no confidence interval
data(survData) # a dataset with 1000 observations (rows) and 10 features (columns)
resultsBig <- bigSurvSGD::bigSurvSGD(formula=Surv(time, status)~.,data=survData, inference.method="none",
parallel.flag=TRUE, num.cores=2, features.mean = resultsBigscale$features.mean, features.sd = resultsBigscale$features.sd)
(resultsBig$coef)

# Simulated survival data to be read off the memory
data(survData) # a dataset with 1000 observations (rows) and 10 features (columns)
# Save dataset survSGD as bigSurvSGD to be read chunk-by-chunk off the memory 
write.csv(survData, file.path(tempdir(), "bigSurvData.csv"), row.names = FALSE) 
dataPath <- file.path(tempdir(), "bigSurvData.csv") # path to where data is
resultsBigOffMemory <- bigSurvSGD::bigSurvSGD(formula=Surv(time, status)~., data=dataPath, 
bigmemory.flag=TRUE, parallel.flag=TRUE, num.cores=2, features.mean = resultsBigscale$features.mean, features.sd = resultsBigscale$features.sd)  #much faster without tests, inference.method="none")
(resultsBigOffMemory$coef)


# Simulated sparse survival data
data(sparseSurvData) # a sparse data with 100 observations (rows) and 150 features (columns)
resultsBigSparse <- bigSurvSGD::bigSurvSGD(formula=Surv(time, status)~.,data=sparseSurvData, 
alpha=0.9, lambda=0.1, features.mean = resultsBigscale$features.mean, features.sd = resultsBigscale$features.sd)
(resultsBigSparse$coef)
```



```r
#data(survData, package = "bigSurvSGD")
#survData[2,3] <- NA
#write.csv(survData, "~/Documents/GitHub/bigPLS/bigSurvData.csv", row.names = FALSE) 
datapath0_NA = path.expand("~/Documents/GitHub/bigPLS/add_data/bigSurvData_NA.csv")

library(bigPLS)
if(FALSE){
resultsBigscale <- bigscale(formula=Surv(time, status)~., data=datapath0_NA,bigmemory.flag=TRUE, parallel.flag=TRUE, num.cores=2)
resultsBigscale
}

# First PLS step -> compute weights

ind.col=1
name.col.all <-(colnames(read.csv(datapath0_NA,nrows=2,header=TRUE))[-c(1,2)])
name.col0 = name.col.all[ind.col]
partialbigSurvSGDv0(datapath = datapath0_NA, resBigscale=resultsBigscale, name.col=name.col0)

# Need to filter for missing values pairwise

simplify2array(lapply(name.col.all,partialbigSurvSGDv0, datapath=datapath0_NA, resBigscale=resultsBigscale, bigmemory.flag=FALSE))
simplify2array(lapply(name.col.all,partialbigSurvSGDv0, datapath=datapath0_NA, resBigscale=resultsBigscale, bigmemory.flag=TRUE))

#install.packages("bigalgebra")
#library("bigalgebra")

datapath0 = path.expand("~/Documents/GitHub/bigPLS/add_data/bigSurvData.csv")

library(bigPLS)

debug(bigPLS:::bigplsRcoxmodel.default)



bigPLS::bigplsRcox(formula=Surv(time, status)~.,data=datapath0,verbose=TRUE)

bigPLS::bigplsRcox(formula=Surv(time, status)~.,data=datapath0, backingfile="temp_plsRcox_file2.bin", backingpath=path.expand("~/Documents/GitHub/bigPLS/add_data/"),  descriptorfile="temp_plsRcox_file2.desc",verbose=TRUE)

bigPLS::bigplsRcox(formula=Surv(time, status)~.,data=NULL,descriptorfile="temp_plsRcox_file.desc", backingpath=path.expand("~/Documents/GitHub/bigPLS/add_data/"),verbose=TRUE)
```

https://stackoverflow.com/questions/49959260/removing-columns-from-big-matrix-which-have-only-one-value





```r

# Simulated survival data - just estimation and no confidence interval
data(survData) # a dataset with 1000 observations (rows) and 10 features (columns)


resultsBig <- bigSurvSGD::bigSurvSGD(formula=Surv(time, status)~.,data=survData, inference.method="none",
parallel.flag=TRUE, num.cores=2, features.mean = resultsBigscale$features.mean, features.sd = resultsBigscale$features.sd)
(resultsBig$coef)

```










```r
#install.packages("bigstatsr")
library(bigstatsr)
library(bigassertr)
library(bigparallelr)

X <- big_attachExtdata()

# No scaling
big_noscale <- big_scale(center = FALSE, scale = FALSE)
class(big_noscale) # big_scale returns a new function
str(big_noscale(X))
big_noscale2 <- big_scale(center = FALSE)
str(big_noscale2(X)) # you can't scale without centering

# Centering
big_center <- big_scale(scale = FALSE)
str(big_center(X))
str(big_scale()(X))

dim(X)
rows_along(X)
cols_along(X)

seq_range(c(3, 10))

X <- big_attachExtdata()

test <- big_parallelize(X, p.FUN = partialbigSurvSGD, p.combine = 'rbind', ncores = 2)
```



# TODO

## bigPLS
PLS regression (linear) or generalized linear regression (glm) for big data


## Cox

plsRcox

## regular PLS
  
## glm PLS

  
imputePLS


Final fit of the model.

# big data Cox (via SGD)
bigSurvSGD


# SGD cox
coxphSGD
