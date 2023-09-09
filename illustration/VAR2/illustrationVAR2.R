################################################################################
#                                                                              #
#                                                                              #
#   Filename  :	  illustrationVAR2.R 					       #
#                                                                              #
#   Project   :   BiomJ article "Ridge estimation of network models            #
#                                from time-course omics data"                  #
#   Date      :   05.07.2018                                                   #
#   Purpose   :   As described in BiomJ article                                #
#									       #
#   R Version :   R-3.2.2                                                      #
#                                                                              #
#                                                                              #
#   Input data files  :    hpvP53 {ragt2ridges}                                #
#   Output data files :    ---                                                 #
#                                                                              #
#   Required R packages :  ragt2ridges, Biobase,                               #
#                                                                              #
#                                                                              #
################################################################################

#remove all objects from the environment
rm(list=ls())

# load libraries
library(ragt2ridges)
library(Biobase)

# load and reformat data
data(hpvP53)
Y <- longitudinal2array(t(exprs(hpvP53rna)))

################################################################################
#   VAR(2) model
################################################################################

# set the directory for results
setwd("./results")

# search for optimal penalty parameters
optLambda <- optPenaltyVAR2(Y, 
                            lambdaMin=c(0.00001, 11, 0.001), 
                            lambdaMax=c(0.08, 16, 0.01), 
                            lambdaInit=c(0.03, 13, 0.00758)
)

# specify grid for contour
lambdaA1grid = seq(0.00001, 0.08, length.out=20)
lambdaA2grid = seq(11, 16, length.out=20)
LOOCVres <- loglikLOOCVcontourVAR2(lambdaA1grid, 
                                   lambdaA2grid,
                                   Y, 
                                   lambdaP=0.00758
)

# contur plot of the full support - Figure 5a
setEPS()
postscript(file="Figure_SM1a.eps")
op = par(pty="s")
contour(lambdaA1grid, 
        lambdaA2grid,
        LOOCVres$llLOOCV,
        xlab="lambdaA1",
        ylab="lambdaA2",
        main="LOOCV log-likelihood",
        nlevels=50
)
points(optLambda[1], 
       optLambda[2],
       pch=20,
       cex=2, 
       col="red"
)
par(op)
dev.off()

# fit VAR(2) model
VAR2hat <- ridgeVAR2(Y,
                     lambdaA1=optLambda[1],
                     lambdaA2=optLambda[2],
                     lambdaP=optLambda[3]
)
Ahat1 = VAR2hat$A1
Ahat2 = VAR2hat$A2
Phat = VAR2hat$P
rownames(Ahat1) = colnames(Ahat1) = rownames(Ahat2) = 
                  colnames(Ahat2) = rownames(Phat)  = 
                  colnames(Phat) = rownames(hpvP53rna)

# determine support for A and O
zeros    <- sparsifyVAR2(Ahat1,
                         Ahat2, 
                         SigmaE=symm(solve(Phat)), 
                         threshold="localFDR",
                         FDRcut=c(0.95, 0.95), 
                         statistics=F
)
zerosA1  = zeros$zerosA1
zerosA2  = zeros$zerosA2
supportP <- sparsify(Phat, 
                     threshold="localFDR", 
                     FDRcut=0.95,
                     output="light")$zeros
zerosP   <- support4ridgeP(nNodes=64,
                           zeros=supportP
)

# determine optimal lambda's with inferred support
optLambda1 <- optPenaltyVAR2(Y,
                             lambdaMin=c(0, 0.01, 0.001), 
                             lambdaMax=c(1, 10, 0.01),
                             lambdaInit=c(0.3, 2, 0.006390177), 
                             zerosA1=zerosA1, 
                             zerosA2=zerosA2,
                             zerosP=zerosP$zeros, 
                             cliquesP=zerosP$cliques, 
                             separatorsP=zerosP$separators, 
                             zerosA1fit="sparse"
)
# specify grid for contour
lambdaA1grid = seq(0.01, 1, length.out=20)
lambdaA2grid = seq(0.1, 10, length.out=20)
LOOCVres1 <- loglikLOOCVcontourVAR2(lambdaA1grid,
                                    lambdaA2grid, 
                                    Y, 
                                    lambdaP=0.006390177, 
                                    zerosA1=zerosA1, 
                                    zerosA2=zerosA2, 
                                    zerosP=zerosP$zeros, 
                                    cliquesP=zerosP$cliques, 
                                    separatorsP=zerosP$separators
)
# contur plot of the spars support - Figure 5b
setEPS()
postscript(file="Figure_SM1b.eps")
op = par(pty="s")
contour(lambdaA1grid,
        lambdaA2grid, 
        LOOCVres1$llLOOCV,
        xlab="lambdaA1", 
        ylab="lambdaA2",
        main="cross-validated log-likelihood",
        nlevels=100
)
points(optLambda1[1],
       optLambda1[2], 
       pch=20, 
       cex=2, 
       col="red"
)
par(op)
dev.off()

# re-fit of the VAR(1) model including the prior knowledge
VAR1hat <- ridgeVAR2(Y=Y,
                     lambdaA1=optLambda1[1],
                     lambdaA2=optLambda1[2],
                     lambdaP=optLambda1[3],
                     zerosA1=zerosA1,
                     zerosA2=zerosA2, 
                     zerosP=zerosP$zeros,
                     cliquesP=zerosP$cliques, 
                     separatorsP=zerosP$separators,
                     zerosA1fit="sparse"
)
Ahat1 = VAR1hat$A1
Ahat2 = VAR1hat$A2
Phat = VAR1hat$P
rownames(Ahat1) = colnames(Ahat1) = rownames(Ahat2) = colnames(Ahat2) <-
    rownames(Phat) = colnames(Phat) = rownames(hpvP53rna)

# graph of the interaction among the genes - Figure 6
setEPS()
postscript(file="Figure_SM2.eps", fonts=c("serif", "Palatino","sans"))
op = par(pty="s")
graphVAR2(Ahat1, 
          Ahat2, 
          Phat, 
          nNames=rownames(Ahat1), 
          type="TSCG",
          vertex.label.cex=0.5,
          vertex.label.font=1,
          vertex.size=4,
          vertex.label.color.T0="darkblue",
          vertex.label.color.T1="darkblue",
          vertex.frame.color="steelblue", 
          vertex.color.T0="lightblue",
          vertex.color.T1="lightblue", 
          edge.width=1.5,
          main=""
)
par(op)
dev.off()

# node statistics table - Table 2
stats <- nodeStatsVAR2(Ahat1,
                       Ahat2,
                       Phat, 
                       as.table=TRUE
)
rownames(stats) = fData(hpvP53)[, 1]
stats[rownames(stats)%in%c("IGFBP3","IGF1","RPRM","CCND2","THBS1","SFN","CCNB1",
    "TP73","DDB2","SESN2"),  ]

# Histogram of the correlation between the fit and the observation
nCovariates = dim(Y)[1]
nTimes = dim(Y)[2]
nSamples = dim(Y)[3]
Yhat = Y[, , ]

for (i in 1:nSamples) Yhat[, -1, i] = cbind(Ahat1 %*% Y[, -nTimes, i][,1],
    Ahat1 %*% Y[, -nTimes, i][,2:7] + Ahat2 %*% Y[, -c(7,8), i])

corFit2 = numeric()
for (j in 1:nCovariates) {
    slh = numeric()
    for (i in 1:4) slh = c(slh, cor(Yhat[j, -1, i], Y[j, -1, i], m="s"))
    corFit2 = rbind(corFit2, slh)
}

################################################################################
#   VAR(1) model
################################################################################

# search for optimal penalty parameters
optLambdas1 <- optPenaltyVAR1(Y,
                              lambdaMin=c(120, 0.0024), 
                              lambdaMax=c(930,0.0026),
                              lambdaInit=c(250, 0.0025)
)
# fit VAR(1) model
VAR1hat <- ridgeVAR1(Y=Y, 
                     lambdaA=optLambdas1[1],
                     lambdaP=optLambdas1[2]
)
Ahat = VAR1hat$A
Phat = VAR1hat$P
rownames(Ahat) = colnames(Ahat) = rownames(Phat) = colnames(Phat) <-
     rownames(hpvP53rna)

# determine support for A and O
zerosA <- sparsifyVAR1(A=Ahat,
                       SigmaE=symm(solve(Phat)),
                       threshold="localFDR", 
                       FDRcut=0.95,
                       statistics=F)$zeros
supportP <- sparsify(Phat,
                     threshold="localFDR",
                     FDRcut=0.95, 
                     output="light")$zeros
zerosP <- support4ridgeP(nNodes=64,
                         zeros=supportP
)

# determine optimal lambda's with inferred support
optLambdas2 <- optPenaltyVAR1(Y, 
                              lambdaMin=c(10^(-5), 10^(-5)),
                              lambdaMax=c(10, 0.1),
                              lambdaInit=c(5, 0.01),
                              zerosA=zerosA,
                              zerosP=zerosP$zeros,
                              cliquesP=zerosP$cliques, 
                              separatorsP=zerosP$separators,
                              zerosAfit="sparse"
)
# re-fit of the VAR(1) model including the prior knowledge
VAR1hat <- ridgeVAR1(Y=Y, 
                     lambdaA=optLambdas2[1],
                     lambdaP=optLambdas2[2],
                     zerosA=zerosA,
                     zerosP=zerosP$zeros,
                     cliquesP=zerosP$cliques, 
                     separatorsP=zerosP$separators, 
                     zerosAfit="sparse"
)
Ahat = VAR1hat$A
Phat = VAR1hat$P
rownames(Ahat) = colnames(Ahat) = rownames(Phat) = colnames(Phat) <-
    rownames(hpvP53rna)

# Histogram of the correlation between the fit and the observation
nCovariates = dim(Y)[1]
nTimes = dim(Y)[2]
nSamples = dim(Y)[3]
Yhat = Y[, , ]
for (i in 1:nSamples) Yhat[, -1, i] = Ahat %*% Y[, -nTimes, i]
corFit1 = numeric()
for (j in 1:nCovariates) {
    slh = numeric()
    for (i in 1:4) slh = c(slh, cor(Yhat[j, -1, i], Y[j, -1, i], m="s"))
    corFit1 = rbind(corFit1, slh)
}
    
# make histograms and QQplot to compare VAR1 and VAR2 model - Figure 7
setEPS()
postscript("Figure_2a.eps")
op = par(pty="s")
hist(corFit1, 
     xlab="Correlation",
     ylab="Frequency", 
     main="Histogram of the correlation fit vs. observation in VAR(1) model",
     n=15,
     col="firebrick",
     border="red",
     xlim=c(-1, 1),
     ylim=c(0, 30)
) 
par(op)
dev.off()

setEPS()
postscript("Figure_2b.eps")
op = par(pty="s")
hist(corFit2, 
     xlab="Correlation", 
     ylab="Frequency", 
     main="Histogram of the correlation fit vs. observation in VAR(2) model",
     n=15, 
     col="blue", 
     border="lightblue",
     xlim=c(-1, 1),
     ylim=c(0, 30)
)
par(op)
dev.off()

setEPS()
postscript("Figure_2c.eps")
op = par(pty="s")
qqplot(corFit1, 
       corFit2, 
       main="QQ-plot: VAR(1) vs. VAR(2)", 
       xlab="VAR(1)", 
       ylab="VAR(2)"
)
abline(0, 
       1, 
       col="red", 
       lwd=2
)
par(op)
dev.off()
