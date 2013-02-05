# 756 final 
# functional regression of thickness of CC
# part one simulation

# generate "C" shape function
I <- 300
N <- 500
set.seed(999)
y.out <- rpois(I, 3)
fix <- runif(I,.1,.3)
datx <- NULL
daty <- NULL
for (i in 1:I){
  theta <- runif(N, pi/3, 7*pi/4)
  thick <- rep(N, 0)
  thick[theta<(3*pi/4)] <- 4*log(y.out[i]+1)/10
  thick[(theta<(5*pi/4))&(theta>(3*pi/4))] <- fix[i]+rnorm(sum((theta<(5*pi/4))&(theta>3*pi/4)),0,0.01)
  thick[theta>(5*pi/4)] <- 2*log(y.out[i]+1)/10
  r <- 1+runif(N, -thick/2, thick/2)
  x <- r*cos(theta)
  y <- r*sin(theta)
  plot(x,y)
  datx <- cbind(datx, x)
  daty <- cbind(daty, y)
}
path1 <- 'C:/Users/Chen/Documents/Coursework/Advanced Statistcal Methods/756 final'
setwd(path1)
save(datx, daty, file='simul_1.RData')
load('simul_1.RData')

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
                  neighbor=0.08, kn=15)
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

for (i in 1:I){
  if (mean(thickness[c(1:10),i])<mean(thickness[c(92:101),i])){
    thickness[,i] <- thickness[c(101:1),i]
  }
}
setwd(path1)
save(thickness, file='thickness.RData')

# Try plot some of the fitted curves
# plot(curr.dat[,1]-mean(curr.dat[,1]), curr.dat[,2]-mean(curr.dat[,2]))
# points(pc$dat.curve, col='blue')

thickness <- t(thickness)
library(refund)
result <- pfr(y.out, funcs=thickness, kz=30, kb=30, family='poisson', method=REML, smooth.cov=T)
p.value <- 2*(1-pnorm(abs(result$BetaHat[[1]]/((result$BetaHat[[1]]-result$Bounds[[1]][,1])/qnorm(0.975)))))


# plot session

# original data with fitted curves
set.seed(9)
samp <- runif(4,1,301)
samp <- floor(samp)
png('principalcurve.png', width=600, height=600, pointsize=12)
layout(matrix(1:4,2,2))
par(mar=c(2,2,1,1))
for (i in samp){
  curr.dat <- cbind(datx[,i], daty[,i])
  pc <- PrinCurve(curr.dat, sub.size=500, grid.n=100, 
                  max.iter=30, err.tol=0.02,  
                  neighbor=0.08, kn=15)
  plot(curr.dat[,1]-mean(curr.dat[,1]), curr.dat[,2]-mean(curr.dat[,2]),
       pch='.', cex=6, col=rainbow(500), xlab='', ylab='', xaxt='n', yaxt='n', bty='n')
  points(pc$dat.curve, col='blue', pch=19)
}
dev.off()

# plot of the thickness function
png('thickness.png', width=600, height=600, pointsize=12)
plot(thickness[1,], ylim=c(0,.8),  type='l', 
     xlab='Along the curve', ylab='Thickness', col='green', cex=.3, main='Thickness Function')
for (i in 2:300){
  lines(thickness[i,], ylim=c(0,0.7),  type='l', cex=.3,col=rgb(i/300,1-i/300,0))
}
dev.off()

# plot of the beta(t), its confidence bound and its corresponding p-value
png('betaest.png',width=700, height=400, pointsize=12)
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
points(apply(result$Bounds[[1]], 1, max), pch=19, col=c.l.b, cex=.5)
plot(p.value, col=c.l, pch=19, xlab='Along the curve', type='b',ylab='P-value', cex=.6, main='Functional P-Value')
lines(1:101, rep(0.05, 101), col='red', lwd=1)
dev.off()

# Backwards mapping

set.seed(9)
samp <- runif(4,1,301)
samp <- floor(samp)
png('backmapping.png', width=600, height=600, pointsize=12)
layout(matrix(1:4,2,2))
par(mar=c(2,2,1,1))
for (i in samp){
  curr.dat <- cbind(datx[,i], daty[,i])
  pc <- PrinCurve(curr.dat, sub.size=500, grid.n=100, 
                  max.iter=30, err.tol=0.02,  
                  neighbor=0.08, kn=15)
  c.ll <- (pc$proj$X1*100+1)%in%which(cl.id==1)
  c.l.r <- rep('blue',500)
  c.l.r[c.ll] <- 'yellow'
  plot(curr.dat[,1]-mean(curr.dat[,1]), curr.dat[,2]-mean(curr.dat[,2]),
       pch='.', cex=2, col=rainbow(500), xlab='', ylab='', xaxt='n', yaxt='n', bty='n')
  points(pc$dat.curve, col=c.l.r, pch=19)
}
dev.off()





