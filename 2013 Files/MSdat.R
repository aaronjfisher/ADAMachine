# Code for analyzing the MS data with EDSS score
# Thickness analysis of the MS data
# Chen Yue on Feb. 22nd, 2013
rm(list=ls())
setwd('C:/Users/Chen/Documents/Research/Modeling HDLFPCA/Masked')
files <- dir(, pattern='.nii')
# let's try one subject
i <- 3
dat.dti <- f.read.nifti.volume(files[i])[,,,1]
mask <- which(dat.dti!=0, arr.ind=T, useNames=F)
dim(mask)
mean(range(mask[,1]))
range(mask[,1])

# get the mid sagittal slice
slice.id <- mask[,1]%in%c(floor(mean(range(mask[,1]))), floor(mean(range(mask[,1])))+1)
plot(mask[slice.id,2], mask[slice.id,3])
setwd('C:/Users/Chen/Documents/GitHub/Research/ADAMachine/2013 Files/core_code_for_alg')
source('principalcurve.R')
source('lmstep.R')
source('pjstep.R')

setwd('C:/Users/Chen/Documents/GitHub/Research/ADAMachine/2013 Files/temp_output')
# do the PC fitting
dat.x <- mask[slice.id,2]
dat.y <- mask[slice.id,3]

# it is too flat, which ruins the curve fitting algorithm, we adjust the y range: four times of the original scale
scale.fac <- 4
curr.dat <- cbind(dat.x, scale.fac*dat.y)

# use angle to determine the intial parameter
initial.para <- acos((dat.x-mean(dat.x))/((dat.x-mean(dat.x))^2+(dat.y-mean(dat.y))^2)^.5)
initial.para <- ((dat.y-mean(dat.y))>0)*initial.para + ((dat.y-mean(dat.y))<0)*((dat.x-mean(dat.x))<0)*(2*pi-initial.para)+
  ((dat.y-mean(dat.y))<0)*((dat.x-mean(dat.x))>0)*(-initial.para)
initial.para <- (initial.para-min(initial.para))/diff(range(initial.para))


pc <- PrinCurve(curr.dat, sub.size=dim(curr.dat)[1], grid.n=100, 
                max.iter=30, err.tol=0.005,  
                neighbor=0.05, kn=12, initial.para=initial.para)
id <- !(pc$proj$X1)%in%c(0,1)
# get rid of the ending point
proj <- pc$proj[id,1]
# since we incorporated a scale factor, we need to transform it back and obtain the height

cdat <- scale(curr.dat, scale=F)
dat.thick <- 2*apply(cbind(cdat[,1], cdat[,2]/scale.fac)-cbind(pc$dat.curve[,1], pc$dat.curve[,2]/scale.fac), 1, function(x){(sum(x^2))^.5})


curr.dist <- dat.thick[id]
grid.x <- seq(0,1, length.out=101)
y <- sapply(grid.x, function(a){
  quantile(curr.dist[abs(a-proj)<0.04],prob=0.95)
})

# 95% quantile is defined for the thickness, PC fitting is biased
ysm <- smooth.spline(grid.x, y, spar=.6)$y


# check the fitting result in new scale
ip.col <- initial.para
plot(curr.dat[,1]-mean(curr.dat[,1]), (curr.dat[,2]-mean(curr.dat[,2]))/scale.fac,
     col=rgb(ip.col,1-ip.col,0), cex=.5, pch=19, ylim=range(curr.dat[,2]-mean(curr.dat[,2]))/scale.fac, asp=4)
lines(pc$result$prin.curve[,1], pc$result$prin.curve[,2]/scale.fac,lty='dashed', pch=19, cex=1, col='black', lwd=2)

# check the fitting result in orginal scale
ip.col <- initial.para
plot(curr.dat[,1]-mean(curr.dat[,1]), (curr.dat[,2]-mean(curr.dat[,2])),
     col=rgb(ip.col,1-ip.col,0), cex=.5, pch=19, asp=1, ylim=range(curr.dat[,2]-mean(curr.dat[,2]))/scale.fac)
lines(pc$result$prin.curve[,1], pc$result$prin.curve[,2], lty='dashed', pch=19, cex=1, col='black', lwd=2)


plot((101:1)/101,y, type='b', lty='solid', pch=19, cex=1, col='black', lwd=2)
points((101:1)/101, ysm, pch=19, col='blue')


# do four subjects
set.seed(1000)
samp <- sample(1:466,4)
layout(matrix(1:8,4,2))
setwd('C:/Users/Chen/Documents/Research/Modeling HDLFPCA/Masked')
for (i in samp){
  dat.dti <- f.read.nifti.volume(files[i])[,,,1]
  mask <- which(dat.dti!=0, arr.ind=T, useNames=F)
#   slice.id <- mask[,1]%in%c(floor(mean(range(mask[,1]))), floor(mean(range(mask[,1])))+1)
  slice.id <- mask[,1]==floor(mean(range(mask[,1])))
  dat.x <- mask[slice.id,2]
  dat.y <- mask[slice.id,3]
  scale.fac <- 4
  curr.dat <- cbind(dat.x, scale.fac*dat.y)
  initial.para <- acos((dat.x-mean(dat.x))/((dat.x-mean(dat.x))^2+(dat.y-mean(dat.y))^2)^.5)
  initial.para <- ((dat.y-mean(dat.y))>0)*initial.para + ((dat.y-mean(dat.y))<0)*((dat.x-mean(dat.x))<0)*(2*pi-initial.para)+
    ((dat.y-mean(dat.y))<0)*((dat.x-mean(dat.x))>0)*(-initial.para)
  initial.para <- (initial.para-min(initial.para))/diff(range(initial.para))
  pc <- PrinCurve(curr.dat, sub.size=dim(curr.dat)[1], grid.n=100, 
                  max.iter=30, err.tol=0.005,  
                  neighbor=0.05, kn=12, initial.para=initial.para)
  id <- !(pc$proj$X1)%in%c(0,1)
  proj <- pc$proj[id,1]
  cdat <- scale(curr.dat, scale=F)
  dat.thick <- 2*apply(cbind(cdat[,1], cdat[,2]/scale.fac)-cbind(pc$dat.curve[,1], pc$dat.curve[,2]/scale.fac), 1, function(x){(sum(x^2))^.5})
  curr.dist <- dat.thick[id]
  grid.x <- seq(0,1, length.out=101)
  y <- sapply(grid.x, function(a){
    quantile(curr.dist[abs(a-proj)<0.05],prob=0.95)
  })
  ysm <- smooth.spline(grid.x, y, spar=.6)$y
  ip.col <- initial.para
  plot(curr.dat[,1]-mean(curr.dat[,1]), (curr.dat[,2]-mean(curr.dat[,2]))/scale.fac,
       col=rgb(ip.col,1-ip.col,0), cex=.5, pch=19, ylim=range(curr.dat[,2]-mean(curr.dat[,2]))/scale.fac, asp=4,
       xaxt='n',yaxt='n', main=paste('Subject',i,sep=''), xlab='', ylab='')
  lines(pc$result$prin.curve[,1], pc$result$prin.curve[,2]/scale.fac,lty='dashed', pch=19, cex=1, col='black', lwd=2)
  plot((101:1)/101,y, type='b', lty='solid', pch=19, cex=1, col='black', lwd=2, xaxt='n', yaxt='n', 
       main='Thickness function', xlab='n', ylab='n')
  points((101:1)/101, ysm, pch=19, col='blue')  
}
