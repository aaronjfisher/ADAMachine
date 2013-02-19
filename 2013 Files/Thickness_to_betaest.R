# Code name PyschAPATH
# Analysis of simulated data 
# 2013 Feb 08 edit - Chen Yue

# be sure to change date tag in files saved

# set working directory
setwd('C:/Users/Chen/Documents/GitHub/research/ADAMachine/2013 Files')

# load thickness data for "c"
rm(list=ls())
load('temp_output/thickness_est_C.RData')

thick_func <- thickness
load('Simulate_C_all_data_feb4.RData')

# check the dimension
dim(thickness)

# plot all the thicknesses
matplot(1:101, thickness, type='l')

# compare the thickness with the original thickness
set.seed(998)
samp <- sample(1:200, 9)
layout(matrix(1:9,3,3))
for (i in 1:9){
  plot((1:100)/100,2*thicknesses[samp[i],], type='l', ylim=range(thicknesses[samp[i],]*2, thickness[,samp[i]]),
       xlab='Coor along the curve', ylab='thickness', col='red')
  lines((101:1)/100,thickness[,samp[i]])
}
# Similar pattern, different value, works
#We use the refund package
library(refund)

PFR.fit <- pfr(Y=y, funcs=t(thickness), kz=35, kb=35, method='REML')
beta.fit <- as.vector(PFR.fit$BetaHat[[1]])
matplot(cbind(PFR.fit$BetaHat[[1]][101:1], PFR.fit$Bounds[[1]][101:1,]),
        type = 'l', lty = c(1,2,2), col = c(1,2,2), ylab = "BetaHat", xlab='Along PC',
        main = "beta function", ylim=range(PFR.fit$Bounds, betat))
lines(1:101,rep(0,101), col='blue', lty='dashed')
lines(betat, col='blue', lty='dotdash')
null.beta.true <- which(betat==0)
null.beta.test <- 101-which(PFR.fit$Bounds[[1]][,2]<0)

# The beta function looks similar but not the same.
# however, we can compare the significant part on the line
# following code does this
plot((1:100)/100, rep(.4,100), cex=2, pch=19,col=rgb(betat==0, betat!=0,1), type='p', ylim=c(-.2,1.2))
lines((101:1)/101, rep(.6,101), cex=2, pch=19, col=rgb(PFR.fit$Bounds[[1]][,2]<0, PFR.fit$Bounds[[1]][,2]>=0, 1), type='p')
text(.2,.7, 'True zero-nonzero Beta')
text(.2,.3, 'Estimated zero-nonzero Beta')

# we see good prediction, good job

# how about use the true thickness? this is gonna be a test of pfr function
PFR.fit <- pfr(Y=y, funcs=thicknesses, kz=35, kb=35, method='REML')
beta.fit <- as.vector(PFR.fit$BetaHat[[1]])
matplot(cbind(PFR.fit$BetaHat[[1]], PFR.fit$Bounds[[1]]),
        type = 'l', lty = c(1,2,2), col = c(1,2,2), ylab = "BetaHat", xlab='Along PC',
        main = "beta function", ylim=range(PFR.fit$Bounds, betat))
lines(1:101,rep(0,101), col='blue', lty='dashed')
lines(betat, col='blue', lty='dotdash')
# not perfect, but good on the edge. 


# Now we use bayesian method to obtain the pointwise credible interval and joint credible band
# winbugs will be used
# We use pca for the thickness decomposition
thick.demean <- t(scale(t(thickness)))
matplot(1:101, thick.demean, type='l')

# demeaned function works the same as the original thickness function
# but the estimate for beta 0 will be different
# Now we obtain the pca scores for the thickness function
dim(thick.demean)
SVD <- svd(thick.demean)

# the top J principal curves, J is the number of PC's that we choose
J <- 35
p.comp.J <- SVD$u[,1:J]
dim(p.comp.J)
p.compscore.J <- (diag(SVD$d)%*%t(SVD$v))[1:J,]
dim(p.compscore.J)

# Now we get the b-spline basis for the beta coefficient
library(splines)
K <- 35
bspl.beta.K <- bs(1:101,df=K)
dim(bspl.beta.K)
matplot(1:101, bspl.beta.K, type='l')

# Now we construct the functional regression form
basis.mat <- t(p.comp.J)%*%bspl.beta.K
design <- t(p.compscore.J)%*%basis.mat

# continue working on 021213
# using winbugs to check
MCMC.result <- bugs(data=list(CJ=C%*%J, y=y.out, K_b=K_b, I=I), inits=list(
  list(B0=0, beta=rep(0,15), taubeta=.1)), 
                    parameters.to.save=c('beta', 'B0', 'taubeta'),
                    model.file="bugcode.txt",n.chains=1, n.burnin=10000, n.thin=1,
                    n.iter=20000, bugs.seed=999, bugs.directory="C:/Users/Chen/Documents/R/winbugs14/",
                    working.directory = path1, debug=F)
save(MCMC.result, file='C_shape.RData')







# Was fitting PCA's by hand here, but gave up
# eigenList<-eigen(t(thickZero)%*%thickZero) 

# svd(thickZero)
# svdt<-svd(t(thickZero)%*%thickZero)
# svdt<-svd(thickZero)

# svdt$d^2/cumsum(svdt$d^2)

# eigenList$values
# eigenList$values^2/cumsum(eigenList$values^2)

#Just use Pre-sets

pcx<-princomp(x=thickZero)

summary(pcx)
names(pcx)

ups<-pcx$load #Matrix of eigen vectors (by column) denoted 
image(ups)
image(t(ups)%*%ups)
t(ups[,1])%*%ups[,1] #all ¦×_k'¦×_k terms will cancel to 1

#First 5 PCs
plot(ups[,1])
plot(ups[,2])
plot(ups[,3])
plot(ups[,4])
plot(ups[,5])

###################################
####### PICKING K_x
#How much did each subj weight on these 5 pcs?
K_x<-15 #How many PC's to use? Once we hit 15 it looks super close
xScoreKx<-pcx$score[,1:K_x]

plot(1:101,ylim=range(thickZero))
for(l in 1:101){
  lines(thickZero[,l])
}

plot(1:101,ylim=range(thickness))
for(l in 1:101){
  lines(thickness[,l])
}


#Est people based on PCs
estCurve<-array(dim=dim(thickZero))
for(i in 1:300){
  estCurvei<-rep(0,p)
  for(k in 1:K_x){
    estCurvei<-estCurvei+xScoreKx[i,k]*ups[,k]
  }
  estCurve[i,]<-estCurvei
}

#Check how well the fit is on 4 random Xi's
par(mfrow=c(2,2))
for(i in sample(1:I,4)){
  plot(thickZero[i,],main=paste('subject ',i),pch=4)
  lines(estCurve[i,])
}

######## END OF PICKING K_x
###################################

mpca<-glm(y.out~xScoreKx,family='poisson') #model PCA

summary(mpca)

#recover beta(t)
betaPCA<-rep(0,p)
for(kb in 1:K_x) betaPCA<-betaPCA+mpca$coef[kb+1]*ups[,kb] #don't include intercept
plot(betaPCA)

####################################################
####################################################
####################################################
# Adding B Splines to Beta representation

library(splines)

#Get splines
K_b<-15 # # of slides
bsplines<-bs(1:p,df=K_b,degree=3,intercept=T) #degree 3 = cubic, could also do quadratic splines and it would look almost identical.
plot(bsplines[,1],type='l',ylim=range(bsplines), main='bspline basis')
for(i in 1:K_b) lines(bsplines[,i])

t(bsplines)%*%bsplines
image(t(bsplines)%*%bsplines) #OK, not exactly orthogonal?? But still a basis!

#Get J matrix (integrated thing)
J<-matrix(nrow=K_x,ncol=K_b)
for(kx in 1:K_x){
  for (kb in 1:K_b){
    J[kx,kb]<-sum(ups[,kx] * bsplines[,kb])/p #*1/p to integrate from 0 to 1
  }  
}
C<-xScoreKx

# Using winbugs do the penalized functional regression
MCMC.result <- bugs(data=list(CJ=C%*%J, y=y.out, K_b=K_b, I=I), inits=list(
  list(B0=0, beta=rep(0,15), taubeta=.1)), 
                    parameters.to.save=c('beta', 'B0', 'taubeta'),
                    model.file="bugcode.txt",n.chains=1, n.burnin=10000, n.thin=1,
                    n.iter=20000, bugs.seed=999, bugs.directory="C:/Users/Chen/Documents/R/winbugs14/",
                    working.directory = path1, debug=F)
save(MCMC.result, file='C_shape.RData')

# Check the chains for beta
name <- sapply(1:15, function(x){paste('beta',as.character(x))})
layout(matrix(1:4,2,2))
for (i in 1:15){
  plot(MCMC.result$sims.list$beta[,i], type='l', ylab='', main=name[i])
}
# beta mixed very well !

# Check the chains for B0 and taubeta
plot(MCMC.result$sims.list$B0, type='l', ylab='', main='B0')
plot(MCMC.result$sims.list$taubeta, type='l', ylab='', main=expression(tau[beta]))

# B0 mixed well while taubeta is not quite good

# Reconstruct all the curves
beta.hist <- MCMC.result$sims.list$beta
dim(beta.hist)
dim(bsplines)

height.hist <- beta.hist%*%t(bsplines)
sigma <- var(height.hist)
# check the correlation coefficient of all the curve grid points
image(cov2cor(sigma))