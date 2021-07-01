# ragt2ridges

This repository contains code, benchmarking study and illustration analysis for the R-package ragt2ridges.

that performs ridge maximum likelihood estimation of vector auto-regressive processes: the VAR(1), VAR(2) and VARX(1) models.
# Code

First, the VAR(1) model and its properties along with its associated time-series chain graph are recapitulated. With this knowledge refreshed, the ridge penalized full ML estimator of the VAR(1) model is presented. The estimator is extended to allow the incorporation of prior knowledge on the support of both temporal and contemporaneous interactions. In both cases memory efficient evaluation of the estimator is outlined. Cross-validation (which requires minor changes to the estimator) is described to guide the choice of the penalty parameters. Then, several strategies (e.g. selection of temporal and contemporaneous relationships, mutual information, and path analysis) for down-stream exploitation of the estimated model are discussed. 

The R-package ragt2ridges performs ridge maximum likelihood estimation of vector auto-regressive processes: the VAR(1), VAR(2), fused VAR(1) and VARX(1). The estimiator is extended to allow the incorporation of prior knowledge on the support of both temporal and contemporaneous interactions. In both cases memory efficient evaluation of the estimator is outlined. Cross-validation (which requires minor changes to the estimator) is described to guide the choice of the penalty parameters. 
In addition, the package offers supporting functionality for the exploitation of estimated models. Among others, i) a procedure to infer the support of the non-sparse ridge estimate is implemented, a table of node-wise network summary statistics, path analysis, mutual information analysis, and impulse response analysis.

# Benchmarking 

The ridge ML estimator of the VAR(1) model is compared to its SCAD counterpart (Abegaz and Wit, 2013), which has been implemented in the R-package SparseTSCGM, by means of simulation. The two methods are compared in terms of squared Frobenius loss of the estimates, as well as, for sensitivity and specificity of their edge selection of the time-series chain graph.

In all, the proposed ridge estimator of the VAR(1) model is a worthy competitor to the SCAD estimator of Abegaz and Wit (2013). With respect to the Frobenius loss the ridge estimator seems even preferable, while for more high-dimensional settings the edge selection properties of the ridge estimator not inferior to that of its SCAD counterpart.

# Illustration

Illustration of the time-series chain graphs underlying the various vector autoregressive models, estimated by means of ridge penalized maximum likelihood.

Models.jpeg![Models](https://user-images.githubusercontent.com/22052679/124133786-79437600-da82-11eb-878c-bf6405d1b4c7.jpeg){:height="16px" width="16px"}



### VAR(1)
This technique aims to unravel the dynamic interrelatedness of the variates (e.g., mRNA genes) of a single molecular level (e.g., mRNA gene expression). The model thus explicates the temporal dependencies among the genes, but also captures the contemporaneous ones (through the inverse of the error covariance matrix).
### VAR(2)
The previous technique is extended to assess the presence of dynamic dependencies over a longer time range than that implied by the VAR(1) model. This is done through the VAR(2) model, which includes an additional explanatory time point, that is the two time points directly preceding the current one may both contribute to the observed variation in the latter. 
### Fused VAR(1)
Using fused VAR(1) model differences among the groups' interaction networks may be identified. Hereto a group-wise VAR(1) model is assumed but fitted jointly to facilitate the borrowing of information when they share network features.
### VARX(1)
When information on additional molecular levels (e.g., DNA copy number or microRNA gene expression) is available, those levels may be incorporated into the network. The VARX(1) model integrates time-varying covariates from other molecular levels (corresponding the “X” in VARX) into the VAR(1) model.

# References

Publications related to ```ragt2ridges``` include:

 - Miok, V., Wilting, S.M., & van Wieringen, W.N. (2017),
   "Ridge estimation of the VAR(1) model and its time series chain graph from multivariate time-course omics data".
    _Biometrical Journal_, 59(1): 172-191.
    ([doi:10.1002/bimj.201500269](http://onlinelibrary.wiley.com/doi/10.1002/bimj.201500269/abstract)). 
 - Miok, V., Wilting, S.M., & van Wieringen, W.N. (2018),
   "Ridge estimation of network models from time-course omics data",
    _Biometrical Journal_, ([doi.org/10.1002/bimj.201700195](https://onlinelibrary.wiley.com/doi/abs/10.1002/bimj.201700195)).
 - van Wieringen, W.N. (2018), 
   "ragt2ridges: Ridge Estimation of Vector Auto-Regressive (VAR) Processes". 
    _R package_, version 0.3.2
