
loglikLOOCVVAR1sim <- function(lambdas, Y, penalty = "ridge", ...) {
## ----------------------------------------------------------------------------
## Title: loglikLOOCVVAR1sim
## ----------------------------------------------------------------------------
## Authors: Wessel N. van Wieringen and Viktorian Miok
## ----------------------------------------------------------------------------
## Description: Leave-k-fold-out (minus) cross-validation log-likelihood of
## VAR(1) model
## ----------------------------------------------------------------------------
## Required Packages: ragt2ridges
## ----------------------------------------------------------------------------
## Usage: loglikLOOCVVAR1sim (lambdas, Y, penalty='ridge', ...)  
##      lambdas: ridge penalities for auto-regressive and precision parameters 
##            Y: Three-dimensional array containing the data 
##      penalty: selection of the
##
## ----------------------------------------------------------------------------
## Value: the minus k-fold cross-validated log-likelihood
## ----------------------------------------------------------------------------
    
    if (as.character(class(Y)) != "array") {
        stop("Input (Y) is of wrong class.")
    }
    if (length(dim(Y)) != 3) {
        stop("Input (Y) is of wrong dimensions: either covariate, time or sample dimension is missing.")
    }
    if (as.character(class(lambdas)) != "numeric") {
        stop("Input (lambdas) is of wrong class.")
    }
    if (length(lambdas) != 2) {
        stop("Input (lambdas) is of wrong length.")
    }
    if (any(is.na(lambdas))) {
        stop("Input (lambdas) is not a vector of non-negative numbers.")
    }
    if (any(lambdas < 0)) {
        stop("Input (lambdas) is not a vector of non-negative numbers.")
    }
    loglik = 0
    if (penalty == "ridge") {
        for (k in 1:dim(Y)[3]) {
            VAR1hat <- ridgeVAR1(Y[, , -k], 
                                 lambdas[1],
                                 lambdas[2], ...)
            for (l in 2:dim(Y)[2]) {
                res = Y[, l, k] - VAR1hat$A %*% Y[, l - 1, k]
                loglik = loglik - t(res) %*% VAR1hat$P %*% res/2 + log(det(VAR1hat$P))/2
            }
        }
        return(-loglik)
    }
    if (penalty == "scad") {
        for (k in 1:dim(Y)[3]) {
            Ylaso = array2longitudinal(Y[, , -k])
            Ylaso1 = as.longitudinal(Ylaso, repeats = dim(Y)[3] - 1)  #
            VAR1hat <- sparse.tscgm(data=Ylaso1,
                                    lam1=lambdas[1], 
                                    lam2=lambdas[2], 
                                    model="ar1",
                                    optimality=NULL,
                                    control=list(maxit.out=5, 
                                                 maxit.in = 10)
            )
            Omega = VAR1hat$theta
            
            for (l in 2:dim(Y)[2]) {
                res = Y[, l, k] - t(VAR1hat$gamma) %*% Y[, l - 1, k]  #as.numeric(determinant(Omega,logarithm=TRUE)$modulus)/2
                loglik = loglik - t(res) %*% Omega %*% res/2 + log(det(Omega))/2
            }
        }
        return(-loglik)
    }
}


simDataGen <- function(p=25,
                       n=15,
                       data="sim5",
                       topologyA="clique", 
                       topologyP="banded") {
## ----------------------------------------------------------------------------
## Title: simDataGen
## ----------------------------------------------------------------------------
## Authors:  Viktorian Miok
## ----------------------------------------------------------------------------
## Description: Generation of autoregressive coeficient and precision matrix 
## from VAR(1) model
## ----------------------------------------------------------------------------
## Required Packages: ragt2ridges
## ----------------------------------------------------------------------------
## Usage: function(p = 25, n = 15, data = "sim5", topologyA = "clique", 
##                 topologyP = "banded")
##            p: dimension of the matrix autoregressive coefficient matrix A
##            n: number of smaples 
##         data: indicate level of sparsity of A
##    topologyA: structure of autoregressive coeficient matrix A: "clique", 
##               "hub", "chain", or "random"                 
##    topologyP: structure of precision matrix P: "complete" and "banded"
## ----------------------------------------------------------------------------
## Value: autoregression coefficient matrix A and precision matrix P
## ----------------------------------------------------------------------------
    
    if (data == "sim5") {
        trueA <- createA(p, 
                         topology=topologyA,
                         nonzeroA=0.3, 
                         nCliques=8, 
                         nHubs=p%/%10,
                         percZeros=0.8,
                         stationary=F
        )
        trueP <- createS(n, 
                         p,
                         topology=topologyP, 
                         nonzero=0.5,
                         banded.n=3, 
                         precision=T
        )
    }
    if (data == "sim25") {
        trueA <- createA(p,
                         topology=topologyA, 
                         nonzeroA=0.3,
                         nCliques=2, 
                         nHubs=p%/%2,
                         percZeros=0.55,
                         stationary=F
        )
        trueP <- createS(n,
                         p,
                         topology=topologyP,
                         nonzero=0.5,
                         banded.n=3, 
                         precision=T
        )
    }
    if (data == "real") {
        data(hpvP53)
        Y <- longitudinal2array(t(exprs(hpvP53)))[1:p, , ]
        VAR1hat <- ridgeVAR1(Y, 0.3, 0.1)  #100,0.0024
        trueA = VAR1hat$A
        trueP = VAR1hat$P
    }
    return(list(A=trueA, P=trueP))
}


createA <- function (p, 
                     topology,
                     nonzeroA=0, 
                     nCliques=1,
                     nHubs=1,
                     nBands=1, 
                     percZeros=0.9,
                     stationary=TRUE){
## ----------------------------------------------------------------------------
## Title: createA
## ----------------------------------------------------------------------------
## Authors:  Wessel van Wieringen
## ----------------------------------------------------------------------------
## Description: Generation of the VAR(1) autoregressive coeficient matrix  
## with various types of topologies
## ----------------------------------------------------------------------------
## Required Packages: ragt2ridges
## ----------------------------------------------------------------------------
## Usage: function(p = 25, n = 15, data = "sim5", topologyA = "clique", 
##                 topologyP = "banded")
##            p: dimension of the matrix autoregressive coefficient matrix A
##    topologyA: structure of autoregressive coeficient matrix A: "clique", 
##               "hub", "chain", or "random"                 
##     nonzeroA: value of nonzero elements of A
##     nCliques: number of cliques when topology="clique"
##        nHubs: number of hubs when topology="hub"
##       nBands: number of bands when topology="chain"
##    precZeros: probability with shich zero elements are samples when 
##               topology="random"
##   stationary: should the generated A be stationary
## ----------------------------------------------------------------------------
## Value: Matrix with auto-regression coefficient matrix A of the VAR(1) model
## ----------------------------------------------------------------------------

    if (as.character(class(p)) != "numeric") {
        stop("Input (p) is of wrong class.")
    }
    if (length(p) != 1) {
        stop("Input (p) is of wrong length.")
    }
    if (is.na(p)) {
        stop("Input (p) is not a positive integer.")
    }
    if (p < 0) {
        stop("Input (p) is not a positive integer.")
    }
    if (as.character(class(topology)) != "character") {
        stop("Input (topology) is of wrong class.")
    }
    if (as.character(class(topology)) == "character") {
        if (!(topology %in% c("clique", "chain", "hub", "random"))) {
            stop("Input (topology) ill-specified.")
        }
    }
    if (as.character(class(nonzeroA)) != "numeric") {
        stop("Input (nonzeroA) is of wrong class.")
    }
    if (length(nonzeroA) != 1) {
        stop("Input (nonzeroA) is of wrong length.")
    }
    if (is.na(nonzeroA)) {
        stop("Input (nonzeroA) is not a non-negative number.")
    }
    if (as.character(class(nCliques)) != "numeric") {
        stop("Input (nCliques) is of wrong class.")
    }
    if (length(nCliques) != 1) {
        stop("Input (nCliques) is of wrong length.")
    }
    if (is.na(nCliques)) {
        stop("Input (nCliques) is not a positive integer.")
    }
    if (nCliques < 0) {
        stop("Input (nCliques) is not a positive integer.")
    }
    if (nCliques > p) {
        stop("Input (nCliques) is not smaller than (or equal to) p.")
    }
    if (as.character(class(nHubs)) != "numeric") {
        stop("Input (nHubs) is of wrong class.")
    }
    if (length(nHubs) != 1) {
        stop("Input (nHubs) is of wrong length.")
    }
    if (is.na(nHubs)) {
        stop("Input (nHubs) is not a positive integer.")
    }
    if (nHubs < 0) {
        stop("Input (nHubs) is not a positive integer.")
    }
    if (nHubs > p) {
        stop("Input (nHubs) is not smaller than (or equal to) p.")
    }
    if (as.character(class(nBands)) != "numeric") {
        stop("Input (nBands) is of wrong class.")
    }
    if (length(nBands) != 1) {
        stop("Input (nBands) is of wrong length.")
    }
    if (is.na(nBands)) {
        stop("Input (nBands) is not a positive integer.")
    }
    if (nBands < 0) {
        stop("Input (nBands) is not a positive integer.")
    }
    if (nBands > p) {
        stop("Input (nBands) is not smaller than (or equal to) p.")
    }
    if (as.character(class(percZeros)) != "numeric") {
        stop("Input (percZeros) is of wrong class.")
    }
    if (length(percZeros) != 1) {
        stop("Input (percZeros) is of wrong length.")
    }
    if (is.na(percZeros)) {
        stop("Input (percZeros) is not a non-negative number.")
    }
    if (percZeros <= 0) {
        stop("Input (percZeros) is not a positive number.")
    }
    if (percZeros >= 1) {
        stop("Input (percZeros) is not smaller than one.")
    }
    if (as.character(class(stationary)) != "logical") {
        stop("Input (stationary) is of wrong class.")
    }
    again = TRUE
    while (again) {
        if (nonzeroA == 0) {
            nonzeroA = runif(1, -1, 1)
        }
        if (topology == "chain") {
            diags = list()
            for (d in 0:nBands) {
                diags[[d + 1]] = rep(nonzeroA, p - d)
            }
            A = as.matrix(bandSparse(p, k = -c(0:nBands), diagonals=diags,
                symmetric=FALSE))
        }
        if (topology == "hub") {
            if (p%%nHubs == 0) {
                hubIDs = 1 + p/nHubs * c(0:(nHubs - 1))
            }
            else {
                hubIDs = 1 + floor(p/nHubs) * c(0:(nHubs - 1))
            }
            A = matrix(0, p, p)
            A[, hubIDs] = nonzeroA
            A[upper.tri(A)] = 0
        }
        if (topology == "clique") {
            if (p%%nCliques == 0) {
                cliqueSizes = rep(p/nCliques, nCliques)
            }
            else {
                cliqueSizes = rep(floor(p/nCliques), nCliques)
                cliqueSizes[nCliques] = cliqueSizes[nCliques] +
                  p%%nCliques
            }
            A = as.matrix(bdiag(lapply(cliqueSizes, function(x) {
                matrix(nonzeroA, x, x)
            })))
            A[upper.tri(A)] = 0
        }
        if (topology == "random") {
            A = matrix(sample(0:1,
                               p * p, 
                               replace=TRUE,
                               prob=c(percZeros,
                                      1 - percZeros)
                        ), 
                        nrow=p,
                        ncol=p
            )
            A[A != 0] = nonzeroA
            A[(row(A)<col(A))] <- 0
        }
        evs = abs(eigen(A, only.values=TRUE)$values)
        if (max(evs) < 1 || !stationary) {
            again = FALSE
        }
        else {
            print("non-stationary A generated: trying again")
        }
    }
    return(A)
}