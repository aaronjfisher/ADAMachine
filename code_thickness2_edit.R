# 756 final 
# functional regression of thickness of CC
# part one simulation

# generate "5" shape function
# I <- 300
# N <- 500
# set.seed(999)
# y.out <- rpois(I, 3)
# datx <- NULL
# daty <- NULL
# set.seed(999)
# err.base <- runif(300,0.03,0.10)
# for (i in 1:I){
#   err <- rnorm(3000,0,err.base[i])
#   err[901:1350] <- runif(450, -(y.out[i])/20, (y.out[i])/20)
#   x1 <- runif(900, 0, 1)
#   x2 <- rep(0, 450)
#   x3 <- runif(450, 0, 0.5)
#   theta <- runif(750, -pi/2, pi/2)
#   x4 <- 0.5+(0.5+err[1801:2550])*cos(theta)
#   x5 <- runif(450, 0, 0.5)
#   z1 <- rep(0, 900)
#   z2 <- runif(450, -1, 0)
#   z3 <- rep(-1, 450)
#   z4 <- -1.5+(0.5+err[1801:2550])*sin(theta)
#   z5 <- rep(-2, 450)
#   x <- c(x1,x2+err[901:1350],x3,x4,x5)
#   y <- c(z1+err[1:900],z2,z3+err[1351:1800],z4,z5+err[2551:3000])
#   datx <- cbind(datx, x)
#   daty <- cbind(daty, y)
#   print(i)
# }
# path1 <- 'C:/Users/Chen/Documents/Coursework/Advanced Statistcal Methods/756 final'
# setwd(path1)
# save(datx, daty, file='simul_2.RData')


# step 1 fit the principal curve
I <- 300
N <- 500
library(mgcv)
path2 <- 'C:/Users/Chen/Documents/Research/Modelling the Corpus Callosum/R/110512Fast/principal curve'
setwd(path2)
source('principalcurve.R')
source('lmstep.R')
source('pjstep.R')
thickness <- NULL
for (i in 1:I){
  curr.dat <- cbind(datx[,i],daty[,i])
  pc <- PrinCurve(curr.dat, sub.size=500, grid.n=100, 
                  max.iter=30, err.tol=0.02,  
                  neighbor=0.05, kn=15)
  proj <- pc$proj
  curr.dist <- pc$height
  curr.dist[proj$X1%in%c(0,1)] <- mean(curr.dist[!(proj$X1%in%c(0,1))])
  grid.x <- seq(0,1, length.out=101)
  y <- sapply(grid.x, function(a){
    quantile(curr.dist[abs(a-proj)<0.025],prob=0.95)
  })
  y <- smooth.spline(grid.x, y, spar=.6)$y
  thickness <- cbind(thickness, y)
  cat('Subject ',i,'\ ')
}

setwd(path1)
save(thickness, file='thickness2.RData')

# Try plot some of the fitted curves
# plot(curr.dat[,1]-mean(curr.dat[,1]), curr.dat[,2]-mean(curr.dat[,2]))
# points(pc$dat.curve, col='blue')










###########
#START HERE!!
mydir <- "/Users/aaronfisher/Documents/Aaron/Chen's Code/projectmaterials"
setwd(mydir)
load('thickness2.Rdata')
load('simul_2.RData')
source('principalcurve.R')
source('lmstep.R')
source('pjstep.R')
I <- 300
N <- 500
set.seed(999)
y.out <- rpois(I, 3)


if(dim(thickness)[1]<300) thickness <- t(thickness)
library(refund)
result <- pfr(y.out, funcs=thickness, kz=30, kb=30, family='poisson', method=REML, smooth.cov=T)

matplot(cbind(result$BetaHat[[1]], result$Bounds[[1]]),
  type = 'l', lty = c(1,2,2), col = c(1,2,2), ylab = "BetaHat", 
  main = "refund plot")

p.value <- 2*(1-pnorm(abs(result$BetaHat[[1]]/((result$BetaHat[[1]]-result$Bounds[[1]][,1])/qnorm(0.975)))))


# plot session

# original data with fitted curves
set.seed(9)
samp <- runif(4,1,301)
samp <- floor(samp)
png('principalcurve2.png', width=600, height=600, pointsize=12)
layout(matrix(1:4,2,2))
par(mar=c(2,2,1,1))
for (i in samp){
  curr.dat <- cbind(datx[,i], daty[,i])
  pc <- PrinCurve(curr.dat, sub.size=500, grid.n=100, 
                  max.iter=30, err.tol=0.02,  
                  neighbor=0.05, kn=15)
  plot(curr.dat[,1]-mean(curr.dat[,1]), curr.dat[,2]-mean(curr.dat[,2]),
       pch=19, cex=.7, col=rainbow(3000), xlab='', ylab='', xaxt='n', yaxt='n', bty='n')
  points(pc$dat.curve, col='blue', pch=19)
}
dev.off()

# plot of the thickness function
png('thickness2.png', width=600, height=600, pointsize=12)
plot(thickness[1,], ylim=c(0,.8),  type='l', 
     xlab='Along the curve', ylab='Thickness', col='green', cex=.3, main='Thickness Function')
for (i in 2:300){
  lines(thickness[i,], ylim=c(0,0.7),  type='l', cex=.3,col=rgb(i/300,1-i/300,0))
}
dev.off()

source('colmixer.r')
# plot of the beta(t), its confidence bound and its corresponding p-value
png('betaest2.png',width=700, height=400, pointsize=12)
cl.id <- apply(result$Bounds[[1]], 1, min)<0
c.l <- rep('blue',101)
c.l[cl.id] <- 'black'
c.l.b <- rep('red',101)
c.l.b[cl.id] <- 'green'
layout(matrix(1:2,1,2))
par(mar=c(2,2,1,1))
plot(result$BetaHat[[1]], type='b', col=c.l, pch=19, xlab='Along the curve', 
     ylab=expression(beta(t)), ylim=c(-0.1, .5), cex=.6, main='Functional estimates')
points(apply(result$Bounds[[1]], 1, min),pch=19, col=c.l.b,cex=.5)
points(apply(result$Bounds[[1]] o, 1, max), pch=19, col=c.l.b, cex=.5)

plot(p.value, col=colmixer(p.value), pch=19, xlab='Along the curve', type='b',ylab='P-value', cex=.6, main='Functional P-Value')

lines(1:101, rep(0.05, 101), col='red', lwd=1)
dev.off()

# Backwards mapping

set.seed(9)
samp <- runif(4,1,301)
samp <- floor(samp)
#png('backmapping2.png', width=600, height=600, pointsize=12)
layout(matrix(1:4,2,2))
par(mar=c(2,2,1,1))
for (i in samp){
  curr.dat <- cbind(datx[,i], daty[,i])
  pc <- PrinCurve(curr.dat, sub.size=500, grid.n=100, 
                  max.iter=30, err.tol=0.02,  
                  neighbor=0.08, kn=15)
  c.ll <- (pc$proj$X1*100+1)%in%which(cl.id==1)
  c.l.r <- rep('blue',3000)
  c.l.r[c.ll] <- 'yellow'
  plot(curr.dat[,1]-mean(curr.dat[,1]), curr.dat[,2]-mean(curr.dat[,2]),pch='.', cex=2, col=rainbow(3000), xlab='', ylab='', xaxt='n', yaxt='n', bty='n')
  points(pc$dat.curve, col=c.l.r, pch=19)
}
dev.off()




