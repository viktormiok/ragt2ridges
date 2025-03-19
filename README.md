<img src="https://github.com/viktormiok/viktormiok.me/blob/main/software/ragt2ridges.png" align="right" height="200" width="200">

[![CRAN status](https://www.r-pkg.org/badges/version/ragt2ridges)](https://cran.r-project.org/package=ragt2ridges) ![](https://img.shields.io/badge/languages-R_,_C++_and_C-orange.svg) ![version](https://img.shields.io/badge/GiHub_version-1.1.0-519dd9) ![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/viktormiok/ragt2ridges) ![GitHub issues](https://img.shields.io/github/issues/viktormiok/ragt2ridges)

![dependencies](https://img.shields.io/badge/dependencies-up%20to%20date-orange)  	![commit](https://img.shields.io/github/last-commit/viktormiok/ragt2ridges) ![GitHub](https://img.shields.io/github/license/viktormiok/ragt2ridges)

[![Edit with Gitpod](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/viktormiok/ragt2ridges) 


# ragt2ridges
This repository contains code, benchmarking study, and illustration analysis for the R-package __`ragt2ridges`__.

- [Introduction](#introduction)
- [Code](#code)
- [Benchmarking](#benchmarking)
- [Illustration](#illustration)
- [Installation](#installation)
- [Data](#data)
- [Tutorials](#tutorials)
- [License](#license)
- [References](#references)

## Introduction

The R-package __`ragt2ridges`__ performs ridge maximum likelihood estimation of vector auto-regressive processes: the VAR(1) model (more to be added). Prior knowledge may be incorporated in the estimation through a) specification of the edges believed to be absent in the time series chain graph, and b) a shrinkage target towards which the parameter estimate is shrunken for large penalty parameter values. Estimation functionality is accompanied by methodology for penalty parameter selection.

In addition, the package offers supporting functionality for exploiting estimated models. Among others, i) a procedure to infer the support of the non-sparse ridge estimate (and thereby of the time series chain graph) is implemented, ii) a table of node-wise network summary statistics, iii) mutual information analysis, and iv) impulse response analysis.

## Code

First, the VAR(1) model, its properties, and its associated time-series chain graph are recapitulated. With this knowledge refreshed, the ridge penalized full ML estimator of the VAR(1) model is presented. The estimator is extended to allow the incorporation of prior knowledge to support both temporal and contemporaneous interactions. In both cases, memory efficient evaluation of the estimator is outlined. Cross-validation (which requires minor changes to the estimator) is described to guide the choice of the penalty parameters. Then, several strategies (e.g. selection of temporal and contemporaneous relationships, mutual information, and path analysis) for downstream exploitation of the estimated model are discussed. 

The R-package __`ragt2ridges`__ performs ridge maximum likelihood estimation of vector auto-regressive processes: the VAR(1), VAR(2), fused VAR(1), and VARX(1). The estimator is extended to allow the incorporation of prior knowledge to support both temporal and contemporaneous interactions. In both cases, memory efficient evaluation of the estimator is outlined. Cross-validation (which requires minor changes to the estimator) is described to guide the choice of the penalty parameters. 
In addition, the package offers supporting functionality for exploiting estimated models. Among others, i) a procedure to infer the support of the non-sparse ridge estimate is implemented, a table of node-wise network summary statistics, path analysis, mutual information analysis, and impulse response analysis.

## Benchmarking 

The ridge ML estimator of the VAR(1) model is compared to its SCAD counterpart (Abegaz and Wit, 2013), which has been implemented in the R-package [__`SparseTSCGM`__](https://cran.r-project.org/web/packages/SparseTSCGM/index.html), employing simulation. The two methods are compared in terms of squared Frobenius loss of the estimates and for sensitivity and specificity of their edge selection of the time-series chain graph.

In all, the proposed ridge estimator of the VAR(1) model is a worthy competitor to the SCAD estimator of Abegaz and Wit (2013). Concerning the Frobenius loss the ridge estimator seems even preferable, while for more high-dimensional settings the edge selection properties of the ridge estimator are not inferior to that of its SCAD counterpart.

## Illustration

Illustration of the time-series chain graphs underlying the various vector autoregressive models estimated using ridge penalized maximum likelihood.

<img src="https://user-images.githubusercontent.com/22052679/124133786-79437600-da82-11eb-878c-bf6405d1b4c7.jpeg" align="top" height="540" width="670">

### VAR(1)
This technique aims to unravel the dynamic interrelatedness of the variates (e.g., mRNA genes) of a single molecular level (e.g., mRNA gene expression). The model thus explains the temporal dependencies among the genes and captures the contemporaneous ones (through the inverse of the error covariance matrix).
### VAR(2)
The previous technique is extended to assess dynamic dependencies over a longer time range than that implied by the VAR(1) model. This is done through the VAR(2) model, which includes an additional explanatory time point, that is the two time points directly preceding the current one may both contribute to the observed variation in the latter. 
### Fused VAR(1)
Using the fused VAR(1) model differences among the groups' interaction networks may be identified. Hereto a group-wise VAR(1) model is assumed but fitted jointly to facilitate the borrowing of information when they share network features.
### VARX(1)
When information on additional molecular levels (e.g., DNA copy number or microRNA gene expression) is available, those levels may be incorporated into the network. The VARX(1) model integrates time-varying covariates from other molecular levels (corresponding to the “X” in VARX) into the VAR(1) model.

## Installation

The R-package __`ragt2ridges`__ depends on [__`rags2ridges`__](https://github.com/markvdwiel/ShrinkBayes) and on [R >= 3.0.0](https://cran.r-project.org/) and is also available from [__CRAN__](https://cran.r-project.org/web/packages/ragt2ridges/ragt2ridges.pdf). This requires the package [__`devtools`__](https://cran.r-project.org/web/packages/devtools/index.html):

``` r
devtools::install_github("viktormiok/ragt2ridges", build_vignettes=TRUE)
```

Please restart R before loading the package and its documentation:

``` r
library(ragt2ridges)
utils::help(ragt2ridges)
utils::vignette("ragt2ridges")
```

## Data
All the data required for performing temporal integrative genomics analysis and published in the reference articles have been deposited in the National Center for Biotechnology Information Gene Expression Omnibus (GEO) and are accessible through the GEO Series accession numbers:
| Data type     | GEO number |
| ------------- | ------------- |
| CGH Arrays  | [__`GSE138724`__](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4117045)  |
| mRNA Arrays  | [__`GSE138079`__](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138079)  |
| miRNA Arrays  | [__`GSE78279`__](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78279)  |

To access one of the data sets for instance GSE78279 you need to run the code below. Unpacking the data requires tar and gunzip, which should already be available on most systems.

```
cd ../  #To get to the main GitHub repo folder
mkdir -p data/tigaR_data_analysis/
cd data/tigaR_data_analysis/
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78279/suppl/GSE78279_RAW.tar
mkdir GSE78279_RAW
tar -C GSE78279_RAW -xvf GSE78279_RAW.tar
gunzip GSE78279_RAW/*_Regional_*
```

## Tutorials

Please see the following tutorials for detailed examples of how to use __`ragt2ridges`__: 

### ragt2ridges walkthrough:
* [Reference manual](https://cran.r-project.org/web/packages/ragt2ridges/ragt2ridges.pdf)
* [Comparison](https://github.com/viktormiok/ragt2ridges/tree/main/comparison)
* [Illustration](https://github.com/viktormiok/ragt2ridges/tree/main/illustration)

## License

__`ragt2ridges`__ is distributed under the GPL-3.0 License. Please read the license before using __`ragt2ridges`__, which it is distributed in the `LICENSE` file.


## References

Publications related to __`ragt2ridges`__ include:

 - **Miok, V.**, Wilting, S.M., van Wieringen, W.N. (2017),
   "[Ridge estimation of the VAR(1) model and its time series chain graph from multivariate time-course omics data](http://onlinelibrary.wiley.com/doi/10.1002/bimj.201500269/abstract)".       
    *Biometrical Journal*, 59(1): 172-191.
 - Babion, I., **Miok, V.**, Jaspers, A., Huseinovic, A., Steenbergen, R.D., van Wieringen, W.N., Wilting, S.M. (2018), "[Comprehensive molecular profiling of HPV-induced transformation over time](https://aacrjournals.org/cancerres/article/78/13_Supplement/5059/629528/Abstract-5059-Comprehensive-molecular-profiling-of)",  
  *Cancer Research*, 78, (13 Supplement), 5059-5059
 - **Miok, V.**, Wilting, S.M., van Wieringen, W.N. (2019),
   "[Ridge estimation of network models from time-course omics data](https://onlinelibrary.wiley.com/doi/abs/10.1002/bimj.201700195)",
    *Biometrical Journal*, 61(2):391-405.
 - Babion, I., **Miok, V.**, Jaspers, A., Huseinovic, A., Steenbergen, R.D., van Wieringen, W.N., Wilting, S.M. (2020), "[Identification of Deregulated Pathways, Key Regulators, and Novel miRNA-mRNA Interactions in HPV-Mediated Transformation](https://doi.org/10.3390/cancers12030700)",    
 *Cancers*, 12(3), 700.
 - van Wieringen, W.N. (2018), "[ragt2ridges: Ridge Estimation of Vector Auto-Regressive (VAR) Processes](https://cran.r-project.org/web/packages/ragt2ridges/index.html)". 
    *R package*, version 0.3.2

Please cite the relevant publications if you use __`ragt2ridges`__.
