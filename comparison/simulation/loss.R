################################################################################
#                                                                              #
#                                                                              #       
#   Filename  :	  loss.R     							      				                       #
#                                                                              #       
#   Project   :   BiomJ article "Ridge estimation of VAR(1) model and its time #
#                 series chain graph from multivariate time-course omics data" #
#   Date      :   11.08.2016                                                   #
#   Purpose   :   As described in BiomJ article                                #
#																				                                       #
#   R Version :   R-3.2.2                                                      #          
#                                                                              #
#                                                                              #
#   Input data files  :    ---                                                 #                           
#   Output data files :    ---                                                 #
#                                                                              #
#   Required R packages :  longitudinal, Biobase, Matrix, SparseTSCGM,         #
#                          ragt2ridges                                         #
#                                                                              #
#                                                                              #
################################################################################


rm(list=ls())
set.seed(321)
# load libraries
library(longitudinal)
library(Biobase)
library(Matrix)
library(SparseTSCGM)
library(ragt2ridges)

# set the directory
getDir = getwd()
setwd("./simulation")
source("functions.R")

p=25  # numbers of variables (genes)
n=5  # number of samples (cell lines)
T=20  # number of time points

# generate autoregressive coefficient and precision matrix
sim <- simDataGen(p=p,
                  n=n,
                  data="sim5",
                  topologyA="clique", 
                  topologyP="banded")
trueA = sim$A; trueP = sim$P

# set the directory for output 
setwd(getDir)  
setwd("./results")

# set number of iteration itTotal
itTotal = 1

# define output values
lossPr = lossAr = lossPs = lossAs = numeric()
for(it in 1:itTotal){
    # sample data from a VAR(1) model
    Y <- dataVAR1(n, 
                  T, 
                  trueA,
                  solve(trueP))
    # covariate-wise zero centering of the data
    Yridge <- centerVAR1data(Y)
    
    # optimal lambdas for ragt2ridges
    lamR = seq(4, 0.01, length.out=20)
    LOOCVrPrev <- loglikLOOCVVAR1sim(c(lamR[1], 
                                       lamR[1]), 
                                     Yridge, 
                                     penalty="ridge"
    )
    LOOCVr = numeric()
    for (i in 1:length(lamR)) {
        LOOCVr = cbind(LOOCVr, loglikLOOCVVAR1sim(c(lamR[i],
                                                     lamR[i]), 
                                                   Yridge,
                                                   penalty="ridge")
        )
        if (LOOCVr[i] - LOOCVrPrev > 5) {
            k <- which.min(LOOCVr)
            break
        } else {
            LOOCVrPrev=LOOCVr[i]
            k <- which.min(LOOCVr)
        }
    }
    optLr <- nlminb(c(lamR[k], lamR[k]), 
                    loglikLOOCVVAR1sim, 
                    gradient=NULL,
                    lower=c(10^(-10), 10^(-10)), 
                    Y=Yridge, 
                    penalty="ridge",
                    control=list(rel.tol=10^(-10)))$par
    # Ridge ML eximtion of the VAR(1) model
    ridgeEst <- ridgeVAR1(Y=Yridge,
                          lambdaA=optLr[1],
                          lambdaP=optLr[2]
    )
    ridgeA = ridgeEst$A
    ridgeP = ridgeEst$P
    
    # convert a time-series array to a longitudinal object required for
    # SparseTSCGM
    Yscad <- array2longitudinal(Yridge)
    Yscad1 <- as.longitudinal(Yscad, 
                              repeats=n
    )
    
    # optimal lambdas for ridge method
    lamS = seq(2, 0.01, length.out=20)
    LOOCVs = numeric()
    LOOCVsPrev <- loglikLOOCVVAR1sim(c(lamS[1], lamS[1]),
                                     Yridge, 
                                     penalty="scad")
    for (j in 1:length(lamS)) {
        LOOCVs <- cbind(LOOCVs, loglikLOOCVVAR1sim(c(lamS[j], lamS[j]), 
                                                   Yridge,
                                                   penalty="scad")
        )
        if (LOOCVs[j] - LOOCVsPrev > 5) {
            l <- which.min(LOOCVs)
            break
        } else {
            LOOCVsPrev=LOOCVs[j]
            l <- which.min(LOOCVs)
        }
    }
    optLl <- nlminb(c(lamS[l], lamS[l]), 
                    loglikLOOCVVAR1sim, 
                    gradient=NULL,
                    lower=c(10^(-10),10^(-10)),
                    Y=Yridge, 
                    penalty="scad",
                    control=list(rel.tol=10^(-10)))$par
    # estimation of the autregressive coefficient matrix and precision matrix
    # using SCAD penalties
    scadEst <- sparse.tscgm(data=Yscad1, 
                            lam1=optLl[1],
                            lam2=optLl[2],
                            optimality=NULL
    )
    scadA <- t(scadEst$gamma)
    scadP = scadEst$theta
    
    # Frobenius loss for A and Sigma
    lossPr = cbind(lossPr, sqrt(sum((as.numeric(ridgeP) - as.numeric(trueP))^2)))
    lossPs = cbind(lossPs, sqrt(sum((as.numeric(scadP) - as.numeric(trueP))^2)))
    lossAr = cbind(lossAr, sqrt(sum((as.numeric(ridgeA) - as.numeric(trueA))^2)))
    lossAs = cbind(lossAs, sqrt(sum((as.numeric(scadA) - as.numeric(trueA))^2)))
    LOSSbandedClique25 = list(lossPr=lossPr, 
                               lossPs=lossPs, 
                               lossAr=lossAr, 
                               lossAs=lossAs
    )
    
    # save the output
    save(LOSSbandedClique25,
         file="LOSSbandedClique25.Rdata"
    )
}

