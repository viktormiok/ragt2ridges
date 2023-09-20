################################################################################
#                                                                              #
#                                                                              #
#   Filename  :	  illustrationVAR1fused.R				       #
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

# set the directory for results
setwd("./results")

# load and reformat data
data(hpvP53)
Y <- longitudinal2array(t(exprs(hpvP53rna)))

# group indicies
id <- c(rep(1, 2), rep(2, 2))-1

# search for optimal penalty parameters
optLambda <- optPenaltyVAR1fused(Y, 
                                 id, 
                                 rep(10^(-10), 3),
                                 rep(1000, 3),
                                 c(19, 21, 0.0024))

# specify grid for contour
lambdaA1grid = seq(0.1, 35, length.out=20)
lambdaA2grid = seq(0.1, 35, length.out=20)
LOOCVres1 <- loglikLOOCVcontourVAR1fused(lambdaAgrid, 
                                         lambdaFgrid, 
                                         Y, 
                                         lambdaP=0.0024
)
# make a contour plot
contour(lambdaAgrid, 
        lambdaFgrid,
        LOOCVres1$llLOOCV,
        xlab="lambdaA", 
        ylab="lambdaF",
        main="cross-validated log-likelihood", 
        nlevels=100
)
points(optLambda[1],
       optLambda[2], 
       pch=20, 
       cex=2, 
       col="red"
)
# fit fused VAR(1) model
VAR1hat <- ridgeVAR1fused(Y=Y, 
                          id=id, 
                          lambdaA=optLambda[1],
                          lambdaF=optLambda[2], 
                          lambdaP=optLambda[3]
)
Ahat16 = VAR1hat$As[1:64,]
Ahat18 = VAR1hat$As[65:128,]
Phat = VAR1hat$P
rownames(Ahat16) = colnames(Ahat16) = rownames(Ahat18) = colnames(Ahat18) <-
     rownames(Phat) = colnames(Phat) = rownames(hpvP53rna)

# determine support for A1, A2 and O
zerosA1 <- sparsifyVAR1(A=Ahat16,
                        SigmaE=symm(solve(Phat)),
                        threshold="localFDR", 
                        FDRcut=0.95,
                        statistics=F)$zeros
zerosA2 <- sparsifyVAR1(A=Ahat18,
                        SigmaE=symm(solve(Phat)),
                        threshold="localFDR",
                        FDRcut=0.95,
                        statistics=F)$zeros
zerosP <- sparsify(Phat, 
                   threshold="localFDR",
                   FDRcut=0.95,
                   output="light")$zeros

################################################################################
#    data from HPV16 affected cell lines                                       #
################################################################################

# determine optimal lambda's with inferred support for A1
optLambdasA1 <- optPenaltyVAR1(Y[, , 1:2],
                               lambdaMin=rep(10^(-10), 2),
                               lambdaMax=rep(1000, 2),
                               lambdaInit=c(1, 1),
                               zerosA=zerosA1, 
                               zerosP=zerosP
)

# specify grid for contour
lambdaAgrid = seq(0.5, 2.5, length.out=20)
lambdaPgrid = seq(0.01, 0.03, length.out=20)
LOOCVres1 <- loglikLOOCVcontourVAR1(lambdaAgrid,
                                    lambdaPgrid, 
                                    Y[, , 1:2],
                                    zerosA=zerosA1, 
                                    zerosP=zerosP
)

# plot contour - Figure 9a, left panel      
setEPS()  
postscript(file="Figure_SM3a.eps")
op = par(pty="s")
contour(lambdaAgrid, 
        lambdaPgrid,
        LOOCVres1$llLOOCV,
        xlab="lambdaA",
        ylab="lambdaP",
        main="cross-validated log-likelihood", 
        nlevels=50)
points(optLambdasA1[1],
       optLambdasA1[2],
       pch=20,
       cex=2,
       col="red"
)
par(op)
dev.off()

# re-fit of the VAR(1) model including the prior knowledge
VAR1hatA16 <- ridgeVAR1(Y[,,1:2], 
                        lambdaA=optLambdasA1[1], 
                        lambdaP=optLambdasA1[2],
                        zerosA=zerosA1,
                        zerosP=zerosP
)
Ahat16 = VAR1hatA16$A
rownames(Ahat16) = colnames(Ahat16) = rownames(hpvP53rna)
    
# het map of the model parameters with sparsified support-Figure 10a, left panel
setEPS()
postscript(file="Figure_SM4a.eps")
op = par(pty="s")
edgeHeat(Ahat16, main="ridge re-estimate of A from HPV16 affected cell line,
    with inferred support")
par(op)
dev.off()

# graph of interaction among the genes - Figure 11a, top panel
setEPS()
postscript(file="Figure_SM5a.eps", fonts=c("serif", "Palatino","sans"))
op = par(pty="s")
graphVAR1(Ahat16, 
          VAR1hatA16$P,
          nNames=rownames(Ahat16), 
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
                                                                                
################################################################################
#    data from HPV18 affected cell lines                                       #
################################################################################

# determine optimal lambda's with inferred support for A2
optLambdasA2 <- optPenaltyVAR1(Y[,,3:4], 
                               lambdaMin=rep(10^(-10), 2),
                               lambdaMax=rep(100, 2), 
                               lambdaInit=rep(1,2),
                               zerosA=zerosA2, 
                               zerosP=zerosP
)
# specify grid for contour
lambdaAgrid = seq(1, 2.5, length.out=20)
lambdaPgrid = seq(0.001, 0.006, length.out=20)
LOOCVres1 <- loglikLOOCVcontourVAR1(lambdaAgrid,
                                    lambdaPgrid,
                                    Y[,,3:4],
                                    zerosA=zerosA2,
                                    zerosP=zerosP
)

# plot contour - Figure 9b, right panel     
setEPS()   
postscript(file="Figure_SM3b.eps")
op = par(pty="s")
contour(lambdaAgrid, 
        lambdaPgrid,
        LOOCVres1$llLOOCV,
        xlab="lambdaA",
        ylab="lambdaP",
        main="cross-validated log-likelihood",
        nlevels=100
)
points(optLambdasA2[1],
       optLambdasA2[2],
       pch=20,
       cex=2,
       col="red"
)
par(op)
dev.off()

# re-fit of the VAR(1) model including the prior knowledge
VAR1hatA18 <- ridgeVAR1(Y[,,3:4], 
                        lambdaA=optLambdasA2[1],
                        lambdaP=optLambdasA2[2],
                        zerosA=zerosA2,
                        zerosP=zerosP
)
Ahat18 = VAR1hatA18$A
rownames(Ahat18) = colnames(Ahat18) = rownames(hpvP53rna)
    
# het map of the model parameters with sparsified support-Figure 10b,right panel
setEPS()
postscript(file="Figure_SM4b.eps")
op = par(pty="s")
edgeHeat(Ahat18, main="ridge re-estimate of A from HPV18 affected cell line,
    with inferred support")
par(op)
dev.off()

# graph of the temporal interaction among the genes - Figure 11b, bottom panel
setEPS()
postscript(file="Figure_SM5b.eps", fonts=c("serif", "Palatino","sans"))
op = par(pty="s")
graphVAR1(Ahat18,
          VAR1hatA18$P,
          nNames=rownames(Ahat18), 
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

################################################################################
# overlap and differential graph: load function
# the function below is a modication of the graphVAR1 function from the 
# ragt2-ridges-package to highlight differences among time-series chain graphs
# of two VAR(1) models
################################################################################

graphVAR1diff <- function(sparseA1,
			                    sparseA2, 
			                    type="Aonly", 
			                    prune=TRUE, 
			                    nNames=NULL, 
			                    main=NULL, 
			                    vertex.color.T0="lightcyan2", 
			                    vertex.color.T1="lightcyan2", 
			                    vertex.frame.color="steelblue", 
			                    vertex.label.cex=-1, 
			                    vertex.label.color.T0="black", 
			                    vertex.label.color.T1="black", 
			                    vertex.label.font=1.5, 
			                    vertex.size=-1, 
			                    edge.arrow.size=-1, 
			                    edge.width=-1, ...){



	# if no covariate names are specified the columns 
	# and row names of A are given names 1, 2, et cetera
	if (is.null(nNames)){ 
		nNames = as.character(1:nrow(sparseA1)) 
	}
    
	if (type=="Aonly"){


		Ajoint = (sparseA1 != 0) * (sparseA2 != 0)
		sparseA1[which(Ajoint != 0)] = 0
		sparseA2[which(Ajoint != 0)] = 0
		sparseA1[Ajoint] = 0

		# store number of nodes
		nNodes = nrow(sparseA1)

		# default node size and font size if not provided
		if (vertex.label.cex <= 0 ){ 
			vertex.label.cex = max(6*(nrow(sparseA1))^(-0.8), 0.1) 
		}
		if (vertex.size <= 0){ 
			vertex.size = max(75*(nrow(sparseA1))^(-0.7), 1) 
		}

		# default plot title if not provided
		if (is.null(main)){ 
			main = "Cross-time relations of VAR(1) model" 
		}

		# generate grid
		grid = rbind(cbind(-1, 1:nrow(sparseA1)), cbind(1, 1:nrow(sparseA1)))
    
		# adjacency matrix 
		sparseA = sparseA1 + sparseA2
		adjMat1 = rbind(cbind(0*sparseA1, t(sparseA1)),
		                 cbind(0*sparseA1, 0*t(sparseA1))
		)
		adjMat1[adjMat1 != 0] = 1
		adjMat2 = rbind(cbind(0*sparseA2, t(sparseA2)), 
		                 cbind(0*sparseA2, 0*t(sparseA2)))
		adjMat2[adjMat2 != 0] = 1
		adjMat = adjMat1 + adjMat2

		# convert adjacency matrix to graph-object
		gObj <- graph.adjacency(adjMat)

		# select pos and neg edges
		negEdges <- which(sparseA[which(sparseA != 0)] < 0)
		posEdges <- which(sparseA[which(sparseA != 0)] > 0)
		igraph::E(gObj)[negEdges]$style = 2
		igraph::E(gObj)[posEdges]$style = 1

		adjMatContrast = (sparseA1 != 0) - (sparseA2 != 0)
		A1edges <- which(adjMatContrast[adjMatContrast != 0] > 1)
		A2edges <- which(adjMatContrast[adjMatContrast != 0] < 1)
		igraph::E(gObj)[negEdges]$color = "blue"
		igraph::E(gObj)[posEdges]$color = "firebrick1"

      		# size of the edges       
		edge.widthA <- abs(t(sparseA)[which(t(sparseA) != 0, arr.ind=TRUE)]) 
		edge.widthA  <- edge.widthA  / max(edge.widthA)
		if (edge.width <= 0){ 
			edge.width = 40 * (nrow(sparseA)^(-0.8)) * edge.widthA  
		} else { 
			edge.width = edge.width * edge.widthA  
		}
		if (edge.arrow.size <= 0){ 
			edge.arrow.size = 10 * (nrow(sparseA)^(-0.8)) * edge.widthA 
		} else { 
			edge.arrow.size = edge.arrow.size * edge.widthA  
		}

		# make plot
		plot(gObj, 
			   edge.lty=0, 
			   edge.arrow.size=0, 
			   layout=grid, 
			   vertex.shape="none", 
			   vertex.label=c(nNames, nNames), 
			   vertex.label.family="sans", 
			   vertex.color="white",
			   vertex.label.color="white",
			   main=main, 
			   vertex.label.cex=vertex.label.cex, ...)
		for (e in seq_len(ecount(gObj))){
			graph2 <- delete.edges(gObj, igraph::E(gObj)[(1:ecount(gObj))[-e]])
			if (e != seq_len(ecount(gObj))[length(seq_len(ecount(gObj)))]){
			plot(graph2, 
				   edge.arrow.size=edge.arrow.size[e], 
				   layout=grid, 
				   vertex.size=rep(vertex.size, 2*nNodes), 
				   vertex.label.font=vertex.label.font, 
				   vertex.color=c(rep(vertex.color.T0, nNodes), 
					 rep(vertex.color.T1, nNodes)), 
				   vertex.frame.color=vertex.frame.color,
				   layout=grid, 
				   vertex.label.color=c(rep(vertex.label.color.T0, nNodes), 
					 rep(vertex.label.color.T1, nNodes)),  
				   edge.color=igraph::E(graph2)$color, 
				   vertex.label.cex=vertex.label.cex,
				   vertex.label.cex=0,  
				   edge.width=edge.width[e], 
				   vertex.label.family="sans", 
				   edge.lty=igraph::E(graph2)$style,     
				   # vertex.label=c(nNames, nNames), 
				   vertex.label="",
				   add=TRUE, ...)        
			}
			if (e == seq_len(ecount(gObj))[length(seq_len(ecount(gObj)))]){
			plot(graph2, 
				   edge.arrow.size=edge.arrow.size[e], 
				   layout=grid, 
				   vertex.size=rep(vertex.size, 2*nNodes), 
				   vertex.label.font=vertex.label.font, 
				   vertex.color=c(rep(vertex.color.T0, nNodes), 
					 rep(vertex.color.T1, nNodes)), 
				   vertex.frame.color=vertex.frame.color,
				   layout=grid, 
				   vertex.label.color=c(rep(vertex.label.color.T0, nNodes), 
					 rep(vertex.label.color.T1, nNodes)),  
				   edge.color=igraph::E(graph2)$color, 
				   vertex.label.cex=vertex.label.cex,
				   vertex.label.cex=0,  
				   edge.width=edge.width[e], 
				   vertex.label.family="sans", 
				   edge.lty=igraph::E(graph2)$style,     
				   vertex.label=c(nNames, nNames), 
				   add=TRUE, ...)        
			}
		}        
		text(-1,
		     -1.2,
		     expression(paste(Y[italic(t)], "", sep="")),
		     cex=1.2,
		     font=3
	  )
		text(1,
		     -1.2,
		     expression(paste(Y[italic(t)+1], "", sep="")),
		     cex=1.2
		)
	}
}


################################################################################
#  overlap and differential graph: plot graph
################################################################################

# make differential graphs
A16 = VAR1hatA16$A
P16 = VAR1hatA16$P
diag(P16) = 0
A18 = VAR1hatA18$A
P18 = VAR1hatA18$P
diag(P18) = 0
A16rm <- intersect(which(apply(A16, 
                               1,
                               function(Z){ all(Z == 0 )})),
                   which(apply(A16,
                               2,
                               function(Z){ all(Z == 0)}))
)
P16rm <- intersect(which(apply(P16, 
                               1,
                               function(Z){ all(Z == 0 )})), 
                   which(apply(P16,
                               2,
                               function(Z){ all(Z == 0)}))
)
A18rm <- intersect(which(apply(A18,
                               1,
                               function(Z){ all(Z == 0 )})),
                   which(apply(A18,
                               2,
                               function(Z){ all(Z == 0)})))
P18rm <- intersect(which(apply(P18, 
                               1, 
                               function(Z){ all(Z ==0 )})), 
                   which(apply(P18,
                               2,
                               function(Z){ all(Z == 0)}))
)
idRemove <- intersect(intersect(intersect(A16rm, A18rm), P16rm), P18rm)
Ajoint = (A16 != 0) * (A18 != 0)
Adiff  = (1-(A16 != 0)) * (A18 != 0) + (A16 != 0) * (1-(A18 != 0))
Pjoint = (P16 != 0) * (P18 != 0)
Pdiff  = (1-(P16 != 0)) * (P18 != 0) + (P16 != 0) * (1-(P18 != 0))
Aboth = (A16 + A18)
Aboth[Adiff != 0] = 0
Pboth = P16 + P18
rownames(Aboth) = colnames(Aboth) = rownames(hpvP53rna)
rownames(Pboth) = colnames(Pboth) = rownames(hpvP53rna)


setEPS()
postscript(file="Figure_4leftPanel.eps",
           fonts=c("serif", "Palatino","sans")
)
op = par(pty="s")
graphVAR1(Aboth, 
          Pboth,
          prune=FALSE,
          nNames=rownames(hpvP53rna),
          vertex.color.T0="lightgrey",
          vertex.color.T1="lightgrey",
          vertex.frame.color="dimgrey",
          vertex.label.cex=0.4,
          vertex.label.color.T0="black",
          vertex.label.color.T1="black",
          main="HPV16 and HPV18 time-series chain graph: overlap")
par(op)
dev.off()

setEPS()
postscript(file="Figure_4rightPanel.eps",
           fonts=c("serif", "Palatino","sans")
)
op = par(pty="s")
graphVAR1diff(A16,
              A18,
              type="Aonly",
              prune=FALSE,
              nNames=rownames(hpvP53rna),
              vertex.color.T0="lightgrey", 
              vertex.color.T1="lightgrey",
              vertex.frame.color="grey", 
              vertex.label.cex=0.4,
              vertex.label.color.T0="black",
              vertex.label.color.T1="black",
              main="HPV16 and HPV18 time-series chain graph: contrast"
)
legend("bottom",
       c("HPV16", "HPV18"),
       col=c("blue", "firebrick1"),
       lwd=2
)
par(op)
dev.off()


