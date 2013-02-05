# # 756 final 
# # functional regression of thickness of CC
# # part one simulation
# 

#2013 Jan 29 edit - AF
#Improve data simulation


#MODEL
#Y_i=β_0 + ∫X_i(t)β(t)

# # generate "C" shape function
#NOTE THIS IS NOT THE USUAL NOTATION FOR I????? CHANGE????
I <- 100 #I is the length of the curve
N <- 200 #N is the number of curves
set.seed(999)
y.out <- rpois(I, 3)
fix <- runif(I,.1,.3)
datx <- NULL
daty <- NULL

library(RColorBrewer)
darkCol<-brewer.pal(10,'Paired')
lightCol<-brewer.pal(10,'Set3')


######
#Generate $k$ signals to mix into the simulated thickness

k<-7 #number of signals to mix for X
signals<-array(dim=c(k,I))
signals[1,]<-sin((1:I)/I*2*pi)
signals[2,]<-seq(-1,3,length=I)

library(splines)
bsSig<-t(bs(1:I,df=7,degree=3,intercept=TRUE)) # "b spline signals"
matplot(t(bsplineSignals),type='l')

signals[3,]<-bsSig[1,]+bsSig[3,]
signals[4,]<-bsSig[2,]-bsSig[6,]
signals[5,]<-bsSig[3,]+bsSig[6,]+bsSig[7,]
signals[6,]<-bsSig[4,]
signals[7,]<-bsSig[7,]

for(i in 1:k) signals[i,]<-signals[i,]/max(abs(signals[i,]))
signals<-signals-rowMeans(signals)
matplot(t(signals),type='l',lty=1,col=darkCol)
##############


############
#Generate true B(t) & functional model parameters
bsBeta<-t(bs(1:I,df=15,degree=3,intercept=TRUE))/2 # "b spline signals"
betat<-matrix(colSums(bsBeta[c(1:3,12:15),]))
dim(betat)<-c(I,1)
plot(betat)
beta0<-3
sd_e<-1
###########

###################
#Make X_i(t)'s (thicknesses)
#defining rwalk, a random walk generator.
#to be used to add SMOOTH random noise to the covariate function, beyond our defined principle components
rwalk<-function(I,amplitude=1){
  x<-rep(0,I)
  x[1]<-rnorm(1)
  for(i in 2:I){
    x[i]<-rnorm(1,x[i-1],1)
  }
  outpre1<-lowess(x,f=1/8)$y
  outpre2<-(outpre1-mean(outpre1))
  out<-outpre2*amplitude/max(abs(outpre2))
  return(out)
}
#plot(rwalk(I,amplitude=.2)) #test

#write thicknesses, generate outcomes (y)
y<-rep(NA,N)
e_y<-rnorm(N,0,sd=sd_e) #noise in the model for y, see above for sd_e
thicknesses<-matrix(nrow=N,ncol=I)
for(i in 1:N){
  #Get thicknesses:
  #noise<-rnorm(I) ???
  int<-runif(1,.9,1.1) #some variation in the overall height.
  sigWeightsi<-runif(k,-.2,.2)
  noisei<-rwalk(I,amplitude=.15)
  thicki<-int+matrix(sigWeightsi,nrow=1)%*%signals + noisei
  thicknesses[i,]<-thicki
  
  #get outcome
  y[i]<-beta0+thicknesses[i,]%*%betat+e_y[i]
}
matplot(t(thicknesses[1:40,]),lty=1,type='l',lwd=2.5,col=lightCol)
hist(y)
#sanity check
plot(rowMeans(thicknesses),y)
plot(rowMeans(thicknesses[,betat>mean(betat)]),y) #rough way to check

#####generate 2-d shapes######
cShape.theta<-seq(pi/4,7*pi/4,length=I)
R<-5 #Distance from the center of the circle


#Get X & Y Data
obsPerRow<-8
dat.x<-matrix(nrow=N,ncol=I*obsPerRow)
dat.y<-matrix(nrow=N,ncol=I*obsPerRow)
for(i in 1:N){ #going across subjects
  i.ticker<-0
  for(t in 1:I){ #going across angles
    for(o in 1:obsPerRow){ #going across observations per angle
      i.ticker<-i.ticker+1
      obsThickness.ito<-runif(1,R-thicknesses[i,t],R+thicknesses[i,t])
      dat.x[i,i.ticker]<-cos(cShape.theta[t])*obsThickness.ito
      dat.y[i,i.ticker]<-sin(cShape.theta[t])*obsThickness.ito
    }
  }
}
plot(dat.x[1,],dat.y[1,])
lines(rowMeans(dat.x[1,]),dat.y[1,])

save(file='C_xy_data_feb4.RData',list=)

#END OF GENERATING DATA
########################################################
########################################################
########################################################
########################################################
########################################################



########################################################
             ######### WORK SPACE ##########
             ## IGNORE IF YOU WANT TO :) ###
             ###############################
########################################################

cBaseShape.x<-cos(cShape.theta)*R
cBaseShape.y<-sin(cShape.theta)*R
plot(cBaseShape.x,cBaseShape.y)

#Plot some full shapes
#Note, this isn't how we'll generate data. We'll instead
#generate random radiuses for each theta.
r.i.out<- 5+thicknesses[1,]
c.i.x<-cos(cShape.theta)*r.i.out
c.i.y<-sin(cShape.theta)*r.i.out
plot(c.i.x,c.i.y,type='l',col=darkCol[1])
for(i in 1:60){
  for(j in 1:2){
  r.i.out<- 5+thicknesses[i,]*(j==1)-thicknesses[i,]*(j==2)
  c.i.x<-cos(cShape.theta)*r.i.out
  c.i.y<-sin(cShape.theta)*r.i.out
  lines(c.i.x,c.i.y,type='l',col=lightCol[i])
  }
}



########################################################
########################################################
########################################################
########################################################
########################################################
#OLD STUFF FROM LAST VERSION OF FILE

# 
# path1 <- "/Users/aaronfisher/Documents/Aaron/Chen's Code/projectmaterials"
# setwd(path1)
# save(datx, daty, file='simul_1.RData')
# load('simul_1.RData')
# 
# 
# # step 1 fit the principal curve
# I <- 300
# N <- 500
# library(mgcv)
# path2 <- 'C:/Users/Chen/Documents/Research/Modelling the Corpus Callosum/R/110512Fast/principal curve'
# setwd(path2)
# source('principalcurve.R')
# source('lmstep.R')
# source('pjstep.R')
# thickness <- NULL
# for (i in 1:I){
#   curr.dat <- cbind(datx[,i],daty[,i])
#   pc <- PrinCurve(curr.dat, sub.size=500, grid.n=100, 
#                   max.iter=30, err.tol=0.02,  
#                   neighbor=0.08, kn=15)
#   proj <- pc$proj
#   curr.dist <- pc$height
#   curr.dist[proj$X1%in%c(0,1)] <- mean(curr.dist[!(proj$X1%in%c(0,1))])
#   grid.x <- seq(0,1, length.out=101)
#   y <- sapply(grid.x, function(a){
#     quantile(curr.dist[abs(a-proj)<0.025],prob=0.95)
#   })
#   y <- smooth.spline(grid.x, y, spar=.6)$y
#   thickness <- cbind(thickness, y)
#   cat('Subject ',i,'\ ')
# }
# 
# for (i in 1:I){
#   if (mean(thickness[c(1:10),i])<mean(thickness[c(92:101),i])){
#     thickness[,i] <- thickness[c(101:1),i]
#   }
# }
# setwd(path1)
# save(thickness, file='thickness.RData')

# Try plot some of the fitted curves
# plot(curr.dat[,1]-mean(curr.dat[,1]), curr.dat[,2]-mean(curr.dat[,2]))
# points(pc$dat.curve, col='blue')



mydir <- "/Users/aaronfisher/Documents/Aaron/Chen's Code/projectmaterials"
setwd(mydir)
load('thickness.Rdata') #wrong size?

 #copied from above
  I <- 300
 N <- 500
  set.seed(999)
 y.out <- rpois(I, 3)
#
 

if(dim(thickness)[1]>50)thickness <- t(thickness)
library(refund)
result <- pfr(y.out, funcs=thickness, kz=30, kb=30, family='poisson', method=REML, smooth.cov=TRUE)
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





