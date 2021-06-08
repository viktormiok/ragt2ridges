################################################################################
#                                                                              #
#                                                                              #       
#   Filename  :	  illustration.R   												                     #
#                                                                              #       
#   Project   :   BiomJ article "Ridge estimation of VAR(1) model and its time #
#                 series chain graph from multivariate time-course omics data" #
#   Date      :   11.08.2016                                                   #
#   Purpose   :   As described in BiomJ article                                #
#																				                                       #
#   R Version :   R-3.2.2                                                      #          
#                                                                              #
#                                                                              #
#   Input data files  :    hpvP53 {ragt2ridges}                                #                           
#   Output data files :    ---                                                 #
#                                                                              #
#   Required R packages :  ragt2ridges, Biobase, lattice                       #
#                                                                              #
#                                                                              #
################################################################################


#remove all objects from the environment
rm(list = ls())

# load libraries
library(ragt2ridges)
library(Biobase)
library(lattice)

# set the directory
setwd("./results")

# load and reformat data
data(hpvP53)
Y <- longitudinal2array(t(exprs(hpvP53)))

# search for optimal penalty parameters
optLambdas1 <- optPenaltyVAR1(Y, 
                              lambdaMin=c(120, 0.0024),
                              lambdaMax=c(930, 0.0026),
                              lambdaInit=c(250, 0.0025)
)
# specify grid for contour
lambdaAgrid <- seq(120, 930, length.out = 20)
lambdaPgrid <- sort(c(seq(0.002557, 0.0026, length.out=100),
                      seq(0.001557, 0.0036, length.out=20))
)
LOOCVres1 <- loglikLOOCVcontourVAR1(lambdaAgrid, 
                                    lambdaPgrid, 
                                    Y
)

# plot contour - Figure 30, left panel, SM VIII        
postscript(file="Figure_30a_SM-VIII.eps")
contour(lambdaAgrid, 
        lambdaPgrid,
        LOOCVres1$llLOOCV,
        xlab="lambdaA",
        ylab="lambdaP", 
        main="cross-validated log-likelihood",
        nlevels=25
)
points(optLambdas1[1], 
       optLambdas1[2], 
       pch=20,
       cex=2,
       col="red"
)
dev.off()

# fit VAR(1) model
VAR1hat <- ridgeVAR1(Y = Y,
                     lambdaA=optLambdas1[1],
                     lambdaP=optLambdas1[2]
)
Ahat <- VAR1hat$A
Phat <- VAR1hat$P
rownames(Ahat) <- colnames(Ahat) <- rownames(Phat) <- colnames(Phat) <-
                  rownames(hpvP53)
# heatmaps of the ridge estimates of A and Omega, Figure 31, SM VIII
postscript(file="Figure_31a_SM-VIII.eps")
edgeHeat(Ahat, 
         main="ridge estimate of A"
)
dev.off()
postscript(file="Figure_31b_SM-VIII.eps")
edgeHeat(Phat,
         main="ridge precision estimate"
)
dev.off()

# determine support for A and O
zerosA <- sparsifyVAR1(A=Ahat,
                       SigmaE=symm(solve(Phat)),
                       threshold="localFDR",
                       FDRcut=0.95,
                       statistics=F)$zeros

zerosP <- sparsify(Phat,
                   threshold="localFDR",
                   FDRcut=0.95,
                   output="light")$zeros

# determine optimal lambda's with inferred support
optLambdas2 <- optPenaltyVAR1(Y, 
                              lambdaMin=c(10^(-5), 10^(-5)),
                              lambdaMax=c(10, 0.1),
                              lambdaInit=c(5, 0.01),
                              zerosA=zerosA,
                              zerosP=zerosP, 
                              zerosAfit="sparse"
)
# determine contour
lambdaAgrid <- seq(-1, 1, length.out = 20) + optLambdas2[1]
lambdaPgrid <- seq(-0.001, 0.001, length.out = 20) + optLambdas2[2]
LOOCVres2 <- loglikLOOCVcontourVAR1(lambdaAgrid, 
                                    lambdaPgrid, 
                                    Y,
                                    zerosA=zerosA,
                                    zerosP=zerosP, 
                                    zerosAfit="sparse"
)
# plot contour - Figure 30, right panel, SM VIII
postscript(file="Figure_30b_SM-VIII.eps")
contour(lambdaAgrid,
        lambdaPgrid, 
        LOOCVres2$llLOOCV, 
        xlab="lambdaA",
        ylab="lambdaP",
        main="cross-validated log-likelihood", 
        nlevels=25
)
points(optLambdas2[1], 
       optLambdas2[2],
       pch=20, 
       cex=2, 
       col="red")
dev.off()

# re-fit of the VAR(1) model including the prior knowledge
AhatNonsparse <- Ahat
PhatNonsparse <- Phat
VAR1hat <- ridgeVAR1(Y = Y,
                     lambdaA=optLambdas2[1], 
                     lambdaP=optLambdas2[2],
                     zerosA=zerosA,
                     zerosP=zerosP,
                     zerosAfit="sparse"
)
Ahat <- VAR1hat$A
Phat <- VAR1hat$P
rownames(Ahat) <- colnames(Ahat) <- rownames(Phat) <- colnames(Phat) <-
                  rownames(hpvP53)
    
# heatmaps of the ridge re-estimates of A and Omega, Figure 32, SM VIII
postscript(file="Figure_32a_SM-VIII.eps")
edgeHeat(Ahat, main = "ridge re-estimate of A, with inferred support")
dev.off()
postscript(file="Figure_32b_SM-VIII.eps")
edgeHeat(Phat, main = "ridge precision re-estimate, with inferred support")
dev.off()

# Figure 1, top row, left
postscript(file="Figure_1a.eps")
edgeHeat(adjacentMat(Ahat), 
         legend=FALSE,
         main="inferred support of A"
)
dev.off()
postscript(file="Figure_33a_SM-VIII.eps")
edgeHeat(adjacentMat(Ahat), 
         legend=FALSE,
         main="inferred support of A"
)
dev.off()
postscript(file="Figure_33b_SM-VIII.eps")
edgeHeat(adjacentMat(Phat),
         legend=FALSE,
         main="inferred support of precision matrix"
)
dev.off()

# graph of the temporal interaction among the genes - Figure 1, bottom panel
postscript(file="Figure_1c.eps", fonts=c("serif", "Palatino","sans"))
graphVAR1(Ahat, 
          Phat,
          nNames=rownames(Ahat),
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
dev.off()

# node statistics table - Table 2, SM VIII
stats <- nodeStatsVAR1(Ahat,
                       Phat,
                       as.table=TRUE
)
rownames(stats) <- fData(hpvP53)[, 1]
stats[rownames(stats)%in%c("BBC3", "CCND2", "IGF1", "IGFBP3", "THBS1", "CCNG1",
                           "CDKN2A", "SERPINE1", "SESN2", "STEAP3"),  ]

# Histogram of the correlation between the fit and the observation
nCovariates <- dim(Y)[1]
nTimes <- dim(Y)[2]
nSamples <- dim(Y)[3]
Yhat <- Y[, , ]
for (i in 1:nSamples) Yhat[, -1, i] <- Ahat %*% Y[, -nTimes, i]
corFit <- numeric()
for (j in 1:nCovariates) {
    slh <- numeric()
    for (i in 1:4) slh <- c(slh, cor(Yhat[j, -1, i], Y[j, -1, i], m="s"))
    corFit <- rbind(corFit, slh)
}
# Figure 1, top row, right
postscript(file="Figure_1b.eps")
hist(corFit, 
     xlab="Correlation", 
     ylab="Frequency", 
     main="Histogram of the correlation fit vs. observation",
     n=20, 
     col="blue", 
     border="lightblue",
     xlim=c(-1, 1)
)
dev.off()

#  Fit of the expression levels of the CCNG1, EI24, SFN, STEAP3, Figure 34, SM VIII
#  CCNG1 - 18, EI24 - 30, SFN - 54, STEAP3 - 57:
label <- c("CellLine 1", "CellLine 1", "CellLine 1", "CellLine 1", 
           "CellLine 1", "CellLine 1", "CellLine 1", "CellLine 2", "CellLine 2", 
           "CellLine 2", "CellLine 2", "CellLine 2", "CellLine 2", "CellLine 2",
           "CellLine 3", "CellLine 3", "CellLine 3", "CellLine 3", "CellLine 3",
           "CellLine 3", "CellLine 3", "CellLine 4", "CellLine 4", "CellLine 4",
           "CellLine 4", "CellLine 4", "CellLine 4", "CellLine 4"
)
timefac <- rep(2:8, 4)
# CCNG1 - Figure_34a_SM-VIII
p <- data.frame(as.numeric(Y[18, -1, ]), timefac, label)
colnames(p) <- c("Expression", "TimePoints", "CellLine")
xyplot(Expression ~ TimePoints | CellLine, 
       data=p, 
       layout=c(4, 1), 
       col="blue",
       pch=20, 
       aspect=1
)
for (i in 1:4) {
    trellis.focus("panel", i, 1)
    llines(Yhat[18, -1, i] ~ 2:8, col="red", lwd=2)
}
trellis.unfocus()


# EI24 - Figure_34b_SM-VIII
p <- data.frame(as.numeric(Y[30, -1, ]), timefac, label)
colnames(p) <- c("Expression", "TimePoints", "CellLine")
xyplot(Expression ~ TimePoints | CellLine,
       data=p,
       layout=c(4, 1), 
       col="blue", 
       pch=20,
       aspect=1
)
for (i in 1:4) {
    trellis.focus("panel", i, 1)
    llines(Yhat[30, -1, i] ~ 2:8, col="red", lwd=2)
}
trellis.unfocus()


# SFN - Figure_34c_SM-VIII
p <- data.frame(as.numeric(Y[54, -1, ]), timefac, label)
colnames(p) <- c("Expression", "TimePoints", "CellLine")
xyplot(Expression ~ TimePoints | CellLine, 
       data=p,
       layout=c(4, 1), 
       col="blue", 
       pch=20, 
       aspect=1
)
for (i in 1:4) {
    trellis.focus("panel", i, 1)
    llines(Yhat[54, -1, i] ~ 2:8, col = "red", lwd = 2)
}
trellis.unfocus()


# STEAP3 - Figure_34d_SM-VIII
p <- data.frame(as.numeric(Y[57, -1, ]), timefac, label)
colnames(p) <- c("Expression", "TimePoints", "CellLine")
xyplot(Expression ~ TimePoints | CellLine,
       data=p,
       layout=c(4, 1), 
       col="blue",
       pch=20,
       aspect=1
)
for (i in 1:4) {
    trellis.focus("panel", i, 1)
    llines(Yhat[57, -1, i] ~ 2:8, col="red", lwd=2)
}
trellis.unfocus()


# plot of gene epxression data over time - Figure 2, right panel
p <- data.frame(as.numeric(Y[37, -1, ]), timefac, label)
colnames(p) <- c("Expression", "TimePoints", "CellLine")
xyplot(Expression ~ TimePoints | CellLine, 
       data = p, 
       layout = c(4, 1), 
       col="orange", 
       pch=20,
       aspect=1,
       type="l"
)
for (i in 1:4) {
    trellis.focus("panel", i, 1)
    llines(Y[46, -1, i] ~ 2:8, col="red", lty=2)
    llines(Y[54, -1, i] ~ 2:8, col="purple", lty=3)
}
trellis.unfocus()

