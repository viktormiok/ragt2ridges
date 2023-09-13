################################################################################
#                                                                              #
#                                                                              #
#   Filename  :	  illustrationVARX.R 					       #
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


# clear all.
rm(list=ls())

# load libraries
library(Biobase)
library(CGHcall)
library(ragt2ridges)
library(lattice)
library(abind)

# set the directory for results
setwd("./results")

# load and reformat data
data(hpvP53)
Y <- longitudinal2array(t(exprs(hpvP53rna)))
X1 <- longitudinal2array(t(copynumber(hpvP53cn)))
X2 <- longitudinal2array(t(exprs(hpvP53mir)))

# load the constrained elements to zero for B
zerosB = cbind(cn2rna, mir2rna)
zerosB = which(zerosB==0, arr.ind=TRUE)

# merge X1 and X2 data into the one array
X = abind(X1, X2, along=1)

# search for optimal penalty parameters
optLambda1 <- optPenaltyVARX1(Y,
                              X,
                              lambdaMin=c(250, 100, 0.001),
                              lambdaMax=c(400,200,0.01), 
                              lambdaInit=c(255, 130, 0.00347),
                              lagX=0, 
                              zerosB=zerosB,
                              unbalanced=matrix(c(4, 5, 3, 3), 2, 2)
)

# specify grid for contour
lambdaAgrid = seq(250, 350, length.out=20)
lambdaBgrid = seq(100, 200, length.out=20)
LOOCVres1 <- loglikLOOCVcontourVARX1(lambdaAgrid,
                                     lambdaBgrid,
                                     Y,
                                     X,
                                     lagX=0,
                                     zerosB=zerosB, 
                                     lambdaP=0.00347,
                                     unbalanced=matrix(c(4, 5, 3, 3), 2, 2)
)
# make a contour plot - Figure 15a
setEPS()
postscript(file="Figure_SM6leftPanel.eps")
op = par(pty="s")
contour(lambdaAgrid, 
        lambdaBgrid,
        LOOCVres1$llLOOCV,
        xlab="lambdaA",
        ylab="lambdaB",
        main="cross-validated log-likelihood", 
        nlevels=5
)
points(optLambda1[1],
       optLambda1[2],
       pch=20,
       cex=2,
       col="red"
)
par(op)
dev.off()

# fit VAR(1) model
VAR1hat <- ridgeVARX1(Y=Y,
                      X=X,
                      lambdaA=optLambda1[1],
                      lambdaB=optLambda1[2], 
                      lambdaP=optLambda1[3],
                      zerosB=zerosB, 
                      lagX=0,
                      unbalanced=matrix(c(4, 5, 3, 3), 2, 2)
)
Ahat = VAR1hat$A
B1hat = VAR1hat$B[, 1:64]
B2hat = VAR1hat$B[, 65:170]
Phat = VAR1hat$P
rownames(Ahat) = colnames(Ahat) = rownames(Phat) = colnames(Phat) <-
rownames(B1hat) = rownames(B2hat) = colnames(B1hat) = rownames(hpvP53rna)
colnames (B2hat) = rownames(hpvP53mir)

# determine support for A and O
spar <- sparsifyVARX1(X=X, 
                      A=Ahat,
                      B=VAR1hat$B,
                      SigmaE=symm(solve(Phat)),
                      threshold="localFDR",
                      FDRcut=c(0.95, 0.1),
                      zerosB=zerosB,
                      statistics=F
)
zerosA = spar$zerosA
supportP <- sparsify(Phat, 
                     threshold="localFDR",
                     FDRcut=0.95,
                     output="light")$zeros
zerosP <- support4ridgeP(nNodes=64, 
                         zeros=supportP
)
# determine optimal lambda's with inferred support
optLambda2 <- optPenaltyVARX1(Y=Y,
                              X=X, 
                              lambdaMin=c(0.6, 0.4, 0.00001),
                              lambdaMax=c(5, 1.5, 0.0001), 
                              lambdaInit=c(1.5, 0.7, 0.0000894),
                              lagX=0,
                              zerosA=zerosA, 
                              zerosB=zerosB,
                              zerosP=zerosP$zeros,
                              cliquesP=zerosP$cliques,
                              separatorsP=zerosP$separators,
                              unbalanced=matrix(c(4, 5, 3, 3), 2, 2)
)
# determine contour
lambdaAgrid = seq(0.1, 5, length.out=20)
lambdaBgrid = seq(0.1, 1.5, length.out=20)
LOOCVres1 <- loglikLOOCVcontourVARX1(lambdaAgrid, 
                                     lambdaBgrid,
                                     Y,
                                     X,
                                     lagX=0,
                                     zerosA=zerosA,
                                     zerosB=zerosB,
                                     zerosP=zerosP$zeros,
                                     cliquesP=zerosP$cliques,
                                     separatorsP=zerosP$separators, 
                                     lambdaP=0.0000894,
                                     unbalanced=matrix(c(4, 5, 3, 3), 2, 2)
)

# plot contour
setEPS()
postscript(file="Figure_SM6rightPanel.eps")
op = par(pty="s")
contour(lambdaAgrid, 
        lambdaBgrid,
        LOOCVres1$llLOOCV,
        xlab="lambdaA",
        ylab="lambdaB",
        main="cross-validated log-likelihood", 
        nlevels=10
)
points(optLambda2[1],
       optLambda2[2],
       pch=20, 
       cex=2,
       col="red"
)
par(op)
dev.off()

# re-fit of the VARX(1) model including the prior knowledge
AhatNonsparse = Ahat
B1hatNonsparse = B1hat
B2hatNonsparse = B2hat
PhatNonsparse = Phat                       #
VAR1hat <- ridgeVARX1(Y=Y,
                      X=X,
                      lambdaA=optLambda2[1],
                      lambdaB=optLambda2[2], 
                      lambdaP=optLambda2[3],
                      zerosA=zerosA,
                      zerosB=zerosB,
                      zerosP=zerosP$zeros, 
                      cliquesP=zerosP$cliques,
                      separatorsP=zerosP$separators,
                      lagX=0,
                      unbalanced=matrix(c(4, 5, 3, 3), 2, 2)
)
Ahat = VAR1hat$A
B1hat = VAR1hat$B[,1:64]
B2hat = VAR1hat$B[,65:170]
Phat = VAR1hat$P
rownames(Ahat) = colnames(Ahat) = rownames(Phat) = colnames(Phat) <-
    rownames(B1hat) = rownames(B2hat) = colnames(B1hat) = rownames(hpvP53rna)
    colnames (B2hat) = rownames(hpvP53mir)

# heatmaps of the ridge re-estimates of A, B - Figure 16
setEPS()
postscript(file="Figure_SM7leftPanel.eps")
op = par(pty="s")
    edgeHeat(Ahat, 
             main="ridge re-estimate of A, with inferred support"
)
par(op)
dev.off()
setEPS()
postscript(file="Figure_SM7middlePanel.eps")
op = par(pty="s")
edgeHeat(B1hat, 
         main= "ridge re-estimate of B (DNA copy number),
                  with inferred support"
)
par(op)
dev.off()
setEPS()
postscript(file="Figure_SM7rightPanel.eps")
op = par(pty="s")
edgeHeat(B2hat, 
        main="ridge re-estimate of B (miRNA gene expression),
                 with inferred support"
)
par(op)
dev.off()

# graph of the temporal interaction among the genes - Figure 17
setEPS()
postscript(file="Figure_5.eps", fonts=c("serif", "Palatino","sans"))
op = par(pty="s")
graphVARX1(Ahat,
           VAR1hat$B, 
           Phat,
           nNamesY=rownames(hpvP53rna), 
           nNamesX=c(rownames(hpvP53rna),rownames(hpvP53mir)),
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

#  Histogram of the correlation between the fit and the observation
nCovariates = dim(Y)[1]
nTimes = dim(Y)[2]
nSamples = dim(Y)[3]
Yhat = Y[,,]
Xhat = X[,,]
for (i in 1:nSamples) {
    Yhat[, -1, i] = Ahat %*% Y[, -nTimes, i] + VAR1hat$B %*% Xhat[, -1, i]
}
corFit1 = numeric()
for(j in 1:nCovariates){
    slh = numeric()
    for (i in 1:4){
        slh = c(slh, cor(Yhat[j, -1, i], 
                          Y[j, -1, i],
                          m="s")
        )
    }
    corFit1 = rbind(corFit1, slh)
}

################################################################################
#   VAR(1) model
################################################################################

# search for optimal penalty parameters
optLambd1 <- optPenaltyVAR1(Y, 
                            lambdaMin=c(120, 0.0024),
                            lambdaMax=c(930, 0.0026), 
                            lambdaInit=c(250, 0.0025)
)

# fit VAR(1) model
VAR1hat1 <- ridgeVAR1(Y=Y,
                      lambdaA=optLambd1[1], 
                      lambdaP=optLambd1[2]
)
Ahat1 = VAR1hat1$A
Phat1 = VAR1hat1$P
rownames(Ahat1) = colnames(Ahat1) = rownames(Phat1) = colnames(Phat1) <-
                   rownames(hpvP53rna)

# determine support for A and O
zerosA1 <- sparsifyVAR1(A=Ahat1,
                        SigmaE=symm(solve(Phat1)),
                        threshold="localFDR", 
                        FDRcut=0.95, 
                        statistics=F)$zeros
supportP <- sparsify(Phat1, 
                     threshold="localFDR",
                     FDRcut=0.95,
                     output="light")$zeros
zerosP1 <- support4ridgeP(nNodes=64,
                          zeros=supportP)


# determine optimal lambda's with inferred support
optLambd2 <- optPenaltyVAR1(Y,
                            lambdaMin=c(10^(-5), 10^(-5)),
                            lambdaMax=c(10, 0.1), 
                            lambdaInit=c(5, 0.01),
                            zerosA=zerosA1, 
                            zerosP=zerosP1$zeros, 
                            cliquesP=zerosP1$cliques,
                            separatorsP=zerosP1$separators,
                            zerosAfit="sparse"
)
# re-fit of the VAR(1) model including the prior knowledge
VAR1hat <- ridgeVAR1(Y=Y, 
                     lambdaA=optLambd2[1],
                     lambdaP=optLambd2[2],
                     zerosA=zerosA1,
                     zerosP=zerosP1$zeros,
                     cliquesP=zerosP1$cliques,
                     separatorsP=zerosP1$separators,
                     zerosAfit="sparse"
)
Ahat1 = VAR1hat$A
Phat1 = VAR1hat$P
rownames(Ahat1) = colnames(Ahat1) = rownames(Phat1) = colnames(Phat1) <-
                   rownames(hpvP53rna)

# Histogram of the correlation between the fit and the observation
nCovariates = dim(Y)[1]
nTimes = dim(Y)[2]
nSamples = dim(Y)[3]
Yhat1 = Y[, , ]
for (i in 1:nSamples) Yhat1[, -1, i] = Ahat1 %*% Y[, -nTimes, i]
corFit2 = numeric()
for (j in 1:nCovariates) {
    slh = numeric()
    for (i in 1:4) slh = c(slh, cor(Yhat1[j, -1, i],
                                     Y[j, -1, i],
                                     m="s")
    )
    corFit2 = rbind(corFit2, slh)
}

# make the plot of the histograms to compare VAR1 and VAR2 model
setEPS()
postscript("Figure_SM8leftPanel.eps")
op = par(pty="s")
hist(corFit2, 
     xlab="Correlation",
     ylab="Frequency",
     main="Histogram of the correlation fit vs. observation in VAR(1) model",
     n=10,
     col="firebrick",
     border="red",
     xlim=c(-1, 1),
     ylim=c(0, 50)
)
par(op)
dev.off()

setEPS()
postscript("Figure_SM8rightPanel.eps")
op = par(pty="s")
hist(corFit1, 
     xlab="Correlation",
     ylab="Frequency", 
     main="Histogram of the correlation fit vs. observation in VARX(1) model",
     n=10, 
     col="firebrick", 
     border="red", 
     xlim=c(-1, 1),
     ylim=c(0, 50)
)
par(op)
dev.off()

# comparison plot of the fit with VAR1 and VARX1 - Figure 14
id = 18
setEPS()
postscript("Figure_6.eps")
op = par(pty="s")
label = c(rep("CellLine 1", 7),
          rep("CellLine 2", 7),
          rep("CellLine 3", 7),
          rep("CellLine 4", 7)
)
timefac = rep(2:8, 4)
p = data.frame(as.numeric(Y[id, -1,]),
               timefac, 
               label
)
colnames(p) = c("Expression", "TimePoints", "CellLine")
xyplot(Expression ~ TimePoints|CellLine, 
                                data=p, 
                                layout=c(4, 1),
                                col="blue", 
                                pch=20, 
                                key=list(text=list(c("VAR(1) model","VARX(1)")),
    lines=list(lty=c(2,1),
               col=c("blue","red"),
               lwd=c(3,3)), 
               columns=2,
               space="top"), aspect=1
    )
for(i in 1:4){
         trellis.focus("panel", i, 1)
         llines(Yhat[id, -1, i] ~ 2:8, 
                col="red", 
                lwd=3
         )
         llines(Yhat1[id, -1, i] ~ 2:8,
                col="blue",
                lwd=3,
                lty=2
          )
}
trellis.unfocus()
par(op)
dev.off()



