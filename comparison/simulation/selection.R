################################################################################
#                                                                              #
#                                                                              #       
#   Filename  :	  loss.R											                                 #
#                                                                              #       
#   Project   :   BiomJ article "Ridge estimation of VAR(1) model and its time #
#                 series chain graph from multivariate time-course omics data" #
#   Date      :   11.08.2016                                                   #
#   Purpose   :   As described in BiomJ article                                #
#																				                                       #
#   R Version :   R-3.2.2                                                      #          
#                                                                              #
#                                                                              #
#   Input data files  :                                                        #                           
#   Output data files :                                                        #
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
library(ROCR)
library(Matrix)
library(SparseTSCGM)
library(ragt2ridges)

# set the directory
getDir <- getwd()
setwd("./simulation")
source("functions.R")

p=25  # numbers of variables (genes)
n=15  # number of samples (cell lines)
T=20  # number of time points

# generate autoregressive coefficient and precision matrix
sim <- simDataGen(p=p,
                  n=n,
                  data="sim5",
                  topologyA="clique", 
                  topologyP="banded"
)
trueA1 <- sim$A; trueP <- sim$P

# true parameter matrix where non-zero element are 1
trueA <- trueA1
trueA[trueA != 0] <- 1

# set the directory for output
setwd(getDir)
setwd("./results")

# set number of iteration itTotal
itTotal=1

# define output values
scadT <- ridgeTa <- ridgeTb <- scadF <- ridgeFa <- ridgeFb <- numeric()

for(it in 1:itTotal){
    # sample data from a VAR(1) model
    Y <- dataVAR1(n, 
                  T,
                  trueA1,
                  solve(trueP)
    )
    # covariate-wise zero centering of the data
    Yridge <- centerVAR1data(Y)
    
    # convert a time-series array to a longitudinal object required for
    # SparseTSCGM
    Yscad <- array2longitudinal(Yridge)
    Yscad1 <- as.longitudinal(Yscad, 
                              repeats=n
    )
    
    # optimal penalty parameter ridge method
    lamS <- seq(2, 0.01, length.out=20)
    LOOCVs <- numeric()
    LOOCVsPrev <- loglikLOOCVVAR1sim(c(lamS[1],
                                       lamS[1]),
                                     Yridge, 
                                     penalty="scad")
    for (j in 1:length(lamS)) {
        LOOCVs <- cbind(LOOCVs, 
                        loglikLOOCVVAR1sim(c(lamS[j], 
                                             lamS[j]), 
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
    optLs <- nlminb(c(lamS[l], lamS[l]), 
                    loglikLOOCVVAR1sim,
                    gradient=NULL,
                    lower=c(10^(-10), 10^(-10)),
                    Y=Yridge, 
                    penalty="scad",
                    control=list(rel.tol=0.01))$par
    
    # estimation of the autregressive coefficient matrix and precision matrix
    # using SCAD penalties
    scadA <- t(sparse.tscgm(data=Yscad1,
                            lam1=optLs[1],
                            lam2=optLs[2],
                            optimality=NULL)$gamma)
    scadA[scadA != 0] <- 1
    
    # calculate false positive and true postive rate
    if (sum(scadA == 1) != 0) {
        PRDs <- prediction(as.numeric(scadA),
                           as.numeric(trueA)
        )
        TPRs <- PRDs@tp[[1]][2]/max(PRDs@tp[[1]])
        FPRs <- PRDs@fp[[1]][2]/max(PRDs@fp[[1]])
    } else {
        TPRs <- FPRs <- 0
    }
    
    # optimal penalty parameter for SCAD method
    lamR <- seq(4, 0.01, length.out=20)
    LOOCVrPrev <- loglikLOOCVVAR1sim(c(lamR[1], 
                                       lamR[1]),
                                     Yridge,
                                     penalty="ridge")
    LOOCVr <- numeric()
    for (i in 1:length(lamR)) {
        LOOCVr <- cbind(LOOCVr,
                        loglikLOOCVVAR1sim(c(lamR[i], lamR[i]), 
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
                    control=list(rel.tol=0.01))$par
    # Ridge ML eximtion of the VAR(1) model
    ridgeEst <- ridgeVAR1(Y=Yridge, 
                          lambdaA=optLr[1], 
                          lambdaP=optLr[2]
    )
    # calculate false positive and true postive rate
    if (sum(scadA == 1) != 0) {
        ridgeAa <- matrix(0, p, p)
        # determine null and non-null elements of A
        ridgeAa[sparsifyVAR1(A=ridgeEst$A, 
                             SigmaE=symm(solve(ridgeEst$P)),
                             threshold="top",
                             top=as.numeric(sum(scadA == 1)))$nonzeros] <- 1
        if (sum(ridgeAa == 1) != 0) {
            PRDra <- prediction(as.numeric(ridgeAa),
                                as.numeric(trueA)
            )
            TPRra <- PRDra@tp[[1]][2]/max(PRDra@tp[[1]])
            FPRra <- PRDra@fp[[1]][2]/max(PRDra@fp[[1]])
        } else {
            TPRra <- FPRra <- 0
        }
    } else {
        TPRra <- FPRra <- 0
    }
    ridgeAb <- matrix(0, p, p)
    ridgeAb[sparsifyVAR1(A=ridgeEst$A, 
                         SigmaE=symm(solve(ridgeEst$P)),
                         threshold="localFDR",
                         FDRcut=0.8)$nonzeros] <- 1
    #  Calculate false positive and true postive rate
    if(sum(ridgeAb == 1)!=0){
         PRDrb <- prediction(as.numeric(ridgeAb),as.numeric(trueA))
         TPRrb <- PRDrb@tp[[1]][2]/max(PRDrb@tp[[1]])
         FPRrb <- PRDrb@fp[[1]][2]/max(PRDrb@fp[[1]])
    } else {
        TPRrb <- FPRrb <- 0
    }

    # true postive and false positive rate of edge selection comparison
    scadT <- c(scadT, TPRs)
    ridgeTa <- c(ridgeTa, TPRra)
    ridgeTb <- c(ridgeTb, TPRrb)
    scadF <- c(scadF, FPRs)
    ridgeFa <- c(ridgeFa, FPRra)
    ridgeFb <- c(ridgeFb, FPRrb)


    ROCbandedClique25 <- list(scadT=scadT,
                              ridgeTa=ridgeTa,
                              ridgeTb=ridgeTb, 
                              scadF=scadF,
                              ridgeFa=ridgeFa,
                              ridgeFb=ridgeFb
    )
    # save the output
    save(ROCbandedClique25,
         file="ROCbandedClique25.RData"
    )
}


