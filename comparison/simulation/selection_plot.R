################################################################################
#                                                                              #
#                                                                              #       
#   Filename  :	  selection_plot.R											                       #
#                                                                              #       
#   Project   :   BiomJ article "Ridge estimation of VAR(1) model and its time #
#                 series chain graph from multivariate time-course omics data" #
#   Date      :   11.08.2016                                                   #
#   Purpose   :   As described in BiomJ article                                #
#																				                                       #
#   R Version :   R-3.2.2                                                      #          
#                                                                              #
#                                                                              #
#   Input data files  :    ROCbandedHub25, ROCbandedClique25,                  #
#                          ROCbandedRandom25, ROCcompleteHub25                 #
#                          ROCcompleteClique25, ROCcompleteRandom25,           #
#                          ROCrealData25                                       #                           
#   Output data files :    Figure_10a_SM-VII, Figure_10b_SM-VII                #
#                                                                              #
#   Required R packages :  ---                                                 #
#                                                                              #
#                                                                              #
################################################################################

# set the directory for imput        
setwd("./simulation/intermediate_results/Comparison_output_5/T20N15ROC")

load("ROCbandedHub25.RData")
load("ROCbandedClique25.RData")
load("ROCbandedRandom25.RData")
load("ROCcompleteHub25.RData")
load("ROCcompleteClique25.RData")
load("ROCcompleteRandom25.RData")

# set the directory for output
setwd("./results")

h1 <- ROCbandedHub25
c1 <- ROCbandedClique25
r1 <- ROCbandedRandom25
h2 <- ROCcompleteHub25
c2 <- ROCcompleteClique25
r2 <- ROCcompleteRandom25


# comparison of false positive rate between ridge and SCAD estimatros
# Figure 24a SM-VII
postscript(file="Figure_24a_SM-VII.eps")
rng <- range(h1$scadF, h1$ridgeFa, c1$scadF, c1$ridgeFa, r1$scadF, r1$ridgeFa,
    h2$scadF, h2$ridgeFa, c2$scadF, c2$ridgeFa, r2$scadF, r2$ridgeFa)

boxplot(as.matrix(h1$scadF)[, 1], at = 1, xlim = c(0.5, 17), ylim = rng,
    xaxt = "n", col = "pink", border = "red", ylab = "False positive rate",
    main = expression(paste("False positive rate : p=25, T=20, n=15")),
    pch = 20)
boxplot(as.matrix(h1$ridgeFa)[, 1], at = 2, xaxt = "n", add = TRUE, col = "blue",
    border = "lightblue", pch = 20)
boxplot(as.matrix(h2$scadF)[, 1], at = 4, xaxt = "n", add = TRUE, col = "pink",
    border = "red", pch = 20)
boxplot(as.matrix(h2$ridgeFa)[, 1], at = 5, xaxt = "n", add = TRUE, col = "blue",
    border = "lightblue", pch = 20)

boxplot(as.matrix(c1$scadF)[, 1], at = 7, xaxt = "n", add = TRUE, col = "pink",
    border = "red", pch = 20)
boxplot(as.matrix(c1$ridgeFa)[, 1], at = 8, xaxt = "n", add = TRUE, col = "blue",
    border = "lightblue", pch = 20)
boxplot(as.matrix(c2$scadF)[, 1], at = 10, xaxt = "n", add = TRUE, col = "pink",
    border = "red", pch = 20)
boxplot(as.matrix(c2$ridgeFa)[, 1], at = 11, xaxt = "n", add = TRUE, col = "blue",
    border = "lightblue", pch = 20)

boxplot(as.matrix(r1$scadF)[, 1], at = 13, xaxt = "n", add = TRUE, col = "pink",
    border = "red", pch = 20)
boxplot(as.matrix(r1$ridgeFa)[, 1], at = 14, xaxt = "n", add = TRUE, col = "blue",
    border = "lightblue", pch = 20)
boxplot(as.matrix(r2$scadF)[, 1], at = 16, xaxt = "n", add = TRUE, col = "pink",
    border = "red", pch = 20)
boxplot(as.matrix(r2$ridgeFa)[, 1], at = 17, xaxt = "n", add = TRUE, col = "blue",
    border = "lightblue", pch = 20)

mtext(1, at = c(1.5, 3, 4.5, 6, 7.5, 9, 10.5, 12, 13.5, 15, 16.5),
    text = c("A(hub),", "", "A(hub),", "", "A(clique),", "", "A(clique),", "",
    "A(random),", "", "A(random),"), line = 1)
mtext(1, at = c(1.5, 3, 4.5, 6, 7.5, 9, 10.5, 12, 13.5, 15, 16.5),
    text = c(expression(paste(Omega,"(band)")), "",
    expression(paste(Omega, "(complete)")), "", expression(paste(Omega,
    "(band)")), "", expression(paste(Omega, "(complete)")), "",
    expression(paste(Omega, "(band)")), "",
    expression(paste(Omega, "(complete)"))), line = 2)
legend("topright", c("SCAD", "ridge"), fill = c("pink", "blue"),
    border = c("red", "lightblue"), box.lwd = 0)
dev.off()

# comparison of true positive rate between ridge and SCAD estimatros
# Figure_24b_SM-VII
postscript(file="Figure_24b_SM-VII.eps")
rng <- range(h1$scadT, h1$ridgeTa, c1$scadT, c1$ridgeTa, r1$scadT,
             r1$ridgeTa, h2$scadT, h2$ridgeTa, c2$scadT, c2$ridgeTa,
             r2$scadT, r2$ridgeTa
)
boxplot(as.matrix(h1$scadT)[, 1],
        at = 1,
        xlim = c(0.5, 17),
        ylim = rng,
        xaxt = "n",
        col = "pink", 
        border = "red", 
        ylab = "True positive rate",
        main = expression(paste("True positive rate : p=25, T=20, n=15")),
        pch = 20
)
boxplot(as.matrix(h1$ridgeTa)[, 1],
        at = 2, 
        xaxt = "n",
        add = TRUE,
        col = "blue",
        border = "lightblue",
        pch = 20
)
boxplot(as.matrix(h2$scadT)[, 1], 
        at = 4, 
        xaxt = "n", 
        add = TRUE,
        col = "pink",
        border = "red",
        pch = 20
)
boxplot(as.matrix(h2$ridgeTa)[, 1],
        at = 5,
        xaxt = "n", 
        add = TRUE, 
        col = "blue",
    border = "lightblue",
    pch = 20
)
boxplot(as.matrix(c1$scadT)[, 1], 
        at = 7, 
        xaxt = "n",
        add = TRUE,
        col = "pink",
        border = "red",
        pch = 20
)
boxplot(as.matrix(c1$ridgeTa)[, 1], 
        at = 8,
        xaxt = "n", 
        add = TRUE,
        col = "blue",
        border = "lightblue", 
        pch = 20
)
boxplot(as.matrix(c2$scadT)[, 1], 
        at = 10,
        xaxt = "n",
        add = TRUE,
        col = "pink",
        border = "red", 
        pch = 20
)
boxplot(as.matrix(c2$ridgeTa)[, 1], 
        at = 11, 
        xaxt = "n",
        add = TRUE, 
        col = "blue",
        border = "lightblue",
        pch = 20
)
boxplot(as.matrix(r1$scadT)[, 1], 
        at = 13, 
        xaxt = "n",
        add = TRUE,
        col = "pink",
        border = "red",
        pch = 20
)
boxplot(as.matrix(r1$ridgeTa)[, 1],
        at = 14, 
        xaxt = "n",
        add = TRUE,
        col = "blue",
        border = "lightblue",
        pch = 20
)
boxplot(as.matrix(r2$scadT)[, 1], 
        at = 16,
        xaxt = "n", 
        add = TRUE,
        col = "pink",
        border = "red",
        pch = 20
)
boxplot(as.matrix(r2$ridgeTa)[, 1], 
        at = 17, 
        xaxt = "n", 
        add = TRUE, 
        col = "blue",
        border = "lightblue",
        pch = 20
)


mtext(1, 
      at=c(1.5, 3, 4.5, 6, 7.5, 9, 10.5, 12, 13.5, 15, 16.5),
      text=c("A(hub),","",
             "A(hub),", "",
             "A(clique),", "", 
             "A(clique),", "",
             "A(random),", "",
             "A(random),"),
      line=1
)
mtext(1, 
      at=c(1.5, 3, 4.5, 6, 7.5, 9, 10.5, 12, 13.5, 15, 16.5),
      text=c(expression(paste(Omega, "(band)")), "",
               expression(paste(Omega, "(complete)")), "",
               expression(paste(Omega, "(band)")), "", 
               expression(paste(Omega, "(complete)")), "",
               expression(paste(Omega, "(band)")), "",
               expression(paste(Omega, "(complete)"))),
      line=2
)
legend("topright",
       c("SCAD", "ridge"), 
       fill=c("pink", "blue"),
       border=c("red", "lightblue"),
       box.lwd=0
)
dev.off()
