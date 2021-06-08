################################################################################
#                                                                              #
#                                                                              #       
#   Filename  :	  loss_plot.R   												                       #
#                                                                              #       
#   Project   :   BiomJ article "Ridge estimation of VAR(1) model and its time #
#                 series chain graph from multivariate time-course omics data" #
#   Date      :   11.08.2016                                                   #
#   Purpose   :   As described in BiomJ article                                #
#																				                                       #
#   R Version :   R-3.2.2                                                      #          
#                                                                              #
#                                                                              #
#   Input data files  :    LOSSbandedHub25, LOSSbandedClique25,                #
#                          LOSSbandedRandom25, LOSScompleteHub25               #
#                          LOSScompleteClique25, LOSScompleteRandom25,         #
#                          LOSSrealData25                                      #                           
#   Output data files :    Figure_10a_SM-VII, Figure_10b_SM-VII                #
#                                                                              #
#   Required R packages :  ---                                                 #
#                                                                              #
#                                                                              #
################################################################################


# set the directory for imput   
setwd("./simulation/intermediate_results/Comparison_output_5/T20N5LOSS")

load("LOSSbandedHub25.Rdata")
load("LOSSbandedClique25.Rdata")
load("LOSSbandedRandom25.Rdata")
load("LOSScompleteHub25.Rdata")
load("LOSScompleteClique25.Rdata")
load("LOSScompleteRandom25.Rdata")
load("LOSSrealData25.Rdata")

# set the directory for output
setwd("./results")

h1 <- LOSSbandedHub25
c1 <- LOSSbandedClique25
r1 <- LOSSbandedRandom25
t1 <- LOSSrealData25
h2 <- LOSScompleteHub25
c2 <- LOSScompleteClique25
r2 <- LOSScompleteRandom25


# frobenius loss comparison between SCAD and ridge estimatiors for precision
# Figure_8a_SM-VII
postscript(file="Figure_8a_SM-VII.eps")
rng <- range(h1$lossAs, h1$lossAr, c1$lossAs, c1$lossAr, r1$lossAs, r1$lossAr,
    t1$lossAs, t1$lossAr, h2$lossAs, h2$lossAr, c2$lossAs, c2$lossAr, r2$lossAs,
    r2$lossAr)
boxplot(as.matrix(h1$lossAs)[1, ], 
        at=1,
        xlim=c(0.5, 20),
        ylim=rng,
        xaxt="n",
        col="pink",
        border="red",
        ylab="Frobenius loss",
        main=expression(paste("Loss of A : p=25, T=20, n=5")),
        pch=20
)
boxplot(as.matrix(h1$lossAr)[1, ],
        at=2,
        xaxt="n",
        add=TRUE,
        col="blue",
        border="lightblue",
        pch=20
)
boxplot(as.matrix(h2$lossAs)[1, ], 
        at=4,
        xaxt="n",
        add=TRUE,
        col="pink",
        border="red",
        pch=20
)
boxplot(as.matrix(h2$lossAr)[1, ],
        at = 5,
        xaxt = "n",
        add = TRUE,
        col = "blue",
        border = "lightblue",
        pch = 20
)
boxplot(as.matrix(c1$lossAs)[1, ],
        at=7, 
        xaxt="n",
        add=TRUE,
        col="pink",
        border="red", 
        pch=20
)
boxplot(as.matrix(c1$lossAr)[1, ],
        at=8, 
        xaxt="n",
        add=TRUE, 
        col="blue",
        border="lightblue", 
        pch=20
)
boxplot(as.matrix(c2$lossAs)[1, ], 
        at=10,
        xaxt="n",
        add=TRUE,
        col="pink",
        border="red",
        pch=20
)
boxplot(as.matrix(c2$lossAr)[1, ], 
        at=11,
        xaxt="n", 
        add=TRUE,
        col="blue",
        border="lightblue",
        pch=20
)
boxplot(as.matrix(r1$lossAs)[1, ], 
        at=13,
        xaxt="n", 
        add=TRUE, 
        col="pink",
        border="red", 
        pch=20
)
boxplot(as.matrix(r1$lossAr)[1, ], 
        at=14,
        xaxt="n", 
        add=TRUE, 
        col="blue",
        border="lightblue",
        pch=20
)
boxplot(as.matrix(r2$lossAs)[1, ], 
        at=16, 
        xaxt="n",
        add=TRUE,
        col="pink",
        border="red",
        pch=20
)
boxplot(as.matrix(r2$lossAr)[1, ], 
        at=17,
        xaxt="n", 
        add=TRUE, 
        col="blue",
        border="lightblue",
        pch=20
)
boxplot(as.matrix(t1$lossAs)[1, ],
        at=19,
        xaxt="n",
        add=TRUE,
        col="pink",
        border="red",
        pch=20
)
boxplot(as.matrix(t1$lossAr)[1, ],
        at=20,
        xaxt="n",
        add=TRUE,
        col="blue",
        border="lightblue",
        pch=20
)

mtext(1, 
      at=c(1.5, 3, 4.5, 6, 7.5, 9, 10.5, 12, 13.5, 15, 16.5, 18, 19.5),
      text=c("A(hub),", "",
             "A(hub),", "", 
             "A(clique),", "",
             "A(clique),", "",
             "A(random),", "",
             "A(random),", "",
             "A(data)"), 
      line=1
)
mtext(1, 
      at=c(1.5, 3, 4.5, 6, 7.5, 9, 10.5, 12, 13.5, 15, 16.5, 18, 19.5),
      text=c(expression(paste(Omega, "(band)")), "",
             expression(paste(Omega, "(complete)")), "",
             expression(paste(Omega, "(band)")), "", 
             expression(paste(Omega, "(complete)")), "",
             expression(paste(Omega, "(band)")), "",
             expression(paste(Omega, "(complete)")), "",
             expression(paste(Omega, "(data)"))), 
      line=2
)
legend("topright", 
       c("SCAD", "ridge"), 
       fill=c("pink", "blue"), 
       border=c("red", "lightblue"),
       box.lwd=0
)
dev.off()



# Frobenius loss comparison between SCAD and ridge estimatiors for precision
# Figure_8b_SM-VII
postscript(file="Figure_8b_SM-VII.eps")
rng <- range(h1$lossPs, h1$lossPr, c1$lossPs, c1$lossPr, r1$lossPs, r1$lossPr,
             t1$lossPs, t1$lossPr, h2$lossPs, h2$lossPr, c2$lossPs, c2$lossPr,
             r2$lossPs, r2$lossPr)
boxplot(as.matrix(h1$lossPs)[1, ],
        at=1,
        xlim=c(0.5, 20), 
        ylim=rng,
        xaxt="n", 
        col="pink", 
        border="red",
        ylab="Frobenius loss",
        main=expression(paste("Loss of ", Omega, " : p=25, T=20, n=5")),
        pch=20
)
boxplot(as.matrix(h1$lossPr)[1, ],
        at=2, 
        xaxt="n",
        add=TRUE,
        col="blue",
        border="lightblue",
        pch=20
)
boxplot(as.matrix(h2$lossPs)[1, ],
        at=4,
        xaxt="n", 
        add=TRUE,
        col="pink",
        border="red",
        pch=20
)
boxplot(as.matrix(h2$lossPr)[1, ],
        at=5, 
        xaxt="n", 
        add=TRUE, 
        col="blue",
        border="lightblue", 
        pch=20
)
boxplot(as.matrix(c1$lossPs)[1, ], 
        at=7,
        xaxt="n",
        add=TRUE, 
        col="pink",
        border="red", 
        pch=20
)
boxplot(as.matrix(c1$lossPr)[1, ],
        at=8,
        xaxt="n",
        add=TRUE,
        col="blue",
        border="lightblue",
        pch=20
)
boxplot(as.matrix(c2$lossPs)[1, ],
        at=10, 
        xaxt="n",
        add=TRUE,
        col="pink",
        border="red",
        pch=20
)
boxplot(as.matrix(c2$lossPr)[1, ],
        at=11, 
        xaxt="n", 
        add=TRUE,
        col="blue",
        border="lightblue",
        pch=20
)

boxplot(as.matrix(r1$lossPs)[1, ],
        at=13,
        xaxt="n",
        add=TRUE,
        col="pink",
        border="red",
        pch=20
)
boxplot(as.matrix(r1$lossPr)[1, ],
        at=14,
        xaxt="n",
        add=TRUE,
        col="blue",
        border="lightblue",
        pch=20
)
boxplot(as.matrix(r2$lossPs)[1, ],
        at=16,
        xaxt="n",
        add=TRUE,
        col="pink",
        border="red",
        pch=20
)
boxplot(as.matrix(r2$lossPr)[1, ],
        at=17,
        xaxt="n",
        add=TRUE,
        col="blue",
        border="lightblue",
        pch=20
)
boxplot(as.matrix(t1$lossPs)[1, ],
        at=19,
        xaxt="n", 
        add=TRUE,
        col="pink",
        border="red",
        pch=20
)
boxplot(as.matrix(t1$lossPr)[1, ], 
        at=20, 
        xaxt="n",
        add=TRUE,
        col="blue",
        border="lightblue",
        pch=20
)

mtext(1, 
      at=c(1.5, 3, 4.5, 6, 7.5, 9, 10.5, 12, 13.5, 15, 16.5, 18, 19.5),
      text=c("A(hub),", "",
               "A(hub),", "", 
               "A(clique),", "", 
               "A(clique),", "",
               "A(random),", "", 
               "A(random),", "",
               "A(data)"), 
      line=1
)
mtext(1, 
      at=c(1.5, 3, 4.5, 6, 7.5, 9, 10.5, 12, 13.5, 15, 16.5, 18, 19.5),
      text=c(expression(paste(Omega, "(band)")), "", 
               expression(paste(Omega, "(complete)")), "",
               expression(paste(Omega, "(band)")), "",
               expression(paste(Omega, "(complete)")), "",
               expression(paste(Omega, "(band)")), "", 
               expression(paste(Omega, "(complete)")), "",
               expression(paste(Omega, "(data)"))),
      line=2
)
legend("topright", 
       c("SCAD", "ridge"), 
       fill=c("pink", "blue"), 
       border=c("red", "lightblue"), 
       box.lwd=0
)
dev.off()


