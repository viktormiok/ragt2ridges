# ragt2ridges

This repository contains code, benchmarking study and illustration analysis for the R-package ragt2ridges that performs ridge maximum likelihood estimation of vector auto-regressive processes: the VAR(1), VAR(2) and VARX(1) models.

# Benchmarking 

The ridge ML estimator of the VAR(1) model is compared to its SCAD counterpart (Abegaz and Wit, 2013), which has been implemented in the R-package SparseTSCGM, by means of simulation. The two methods are compared in terms of squared Frobenius loss of the estimates, as well as, for sensitivity and specificity of their edge selection of the time-series chain graph.

In all, the proposed ridge estimator of the VAR(1) model is a worthy competitor to the SCAD estimator of Abegaz and Wit (2013). With respect to the Frobenius loss the ridge estimator seems even preferable, while for more high-dimensional settings the edge selection properties of the ridge estimator not inferior to that of its SCAD counterpart.

# Illustration

Illustration of the time-series chain graphs underlying the various vector autoregressive models, estimated by means of ridge penalized maximum likelihood.
https://onlinelibrary.wiley.com/doi/full/10.1002/bimj.201700195
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
