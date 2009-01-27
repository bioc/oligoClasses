par(mfrow=c(1, 1))
plot(0:1, xlim=c(-0.3, 1), ylim=c(-0.5,1), xaxt="n", yaxt="n", bty="n", type="n",
     xaxs="i", yaxs="i", xlab="", ylab="")
eset.coords <- xy.coords(c(0.40, 0.6, 0.6, 0.4),
                         c(0.8, 0.8, 0.99, 0.99))
polygon(eset.coords, col="grey70")
text(mean(eset.coords$x), mean(eset.coords$y), "eSet")  

sls.coords <- xy.coords(eset.coords$x,
                        c(0.5, 0.5, 0.69, 0.69))
polygon(sls.coords, col="grey70")
text(mean(sls.coords$x), mean(sls.coords$y), "SnpLevelSet")

arrows(x0=mean(eset.coords$x), y0=0.80,
       x1=mean(eset.coords$x), y1=0.69, length=0.1)


y <- c(0, 0, 0.15, 0.15)
callset.coords <- xy.coords(c(0, 0.25, 0.25, 0),
                            y)
polygon(callset.coords)
text(mean(callset.coords$x), mean(callset.coords$y), "SnpCallSet")  
arrows(x0=mean(sls.coords$x), y0=0.5,
       x1=mean(callset.coords$x), y1=y[3], length=0.1)

yPlus <- c(-0.4, -0.4,-0.25, -0.25)
callsetPlus.coords <- xy.coords(c(-0.20, .1, .1, -0.20), yPlus)
polygon(callsetPlus.coords)
text(mean(callsetPlus.coords$x), mean(callsetPlus.coords$y), "SnpCallSetPlus")  
arrows(x0=mean(callset.coords$x), y0=y[1],
       x1=mean(callsetPlus.coords$x), y1=-0.25, length=0.1)


cnset.coords <- xy.coords(c(0.30, 0.65, 0.65, 0.30),
                          y)
polygon(cnset.coords)
text(mean(cnset.coords$x), mean(cnset.coords$y), "SnpCopyNumberSet")
arrows(x0=mean(sls.coords$x), y0=0.5,
       x1=mean(cnset.coords$x), y1=y[3], length=0.1)


snpset.coords <- xy.coords(c(0.7, 0.95, 0.95, 0.7),
                           y)
polygon(snpset.coords)
text(mean(snpset.coords$x), mean(snpset.coords$y), "oligoSnpSet")
arrows(x0=mean(sls.coords$x), y0=0.5,
       x1=mean(snpset.coords$x), y1=y[3], length=0.1)

exprset.coords <- xy.coords(c(0.7, 0.95, 0.95, 0.7),
                            c(0.6, 0.6, 0.69, 0.69))
polygon(exprset.coords)
text(mean(exprset.coords$x), mean(exprset.coords$y), "ExpressionSet")
arrows(x0=mean(eset.coords$x), 0.8,
       x1=mean(exprset.coords$x), 0.69, length=.1)

rnalevel <- xy.coords(c(0.65, 1, 1, 0.65),
                      c(0.5, 0.5, 0.79, 0.79))
polygon(rnalevel, lty=2)
text(mean(rnalevel$x), 0.75, "Gene Expression")



