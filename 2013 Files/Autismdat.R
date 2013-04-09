# Code name PyschAPATH
# Analysis of real data 
# 2013 Feb 14 edit - Chen Yue
# get the names of all mask of CC
files <- dir(path='C:/Users/Chen/Documents/Research/Thickness analysis/Data/cleaned/mask', pattern='.nii')
length(files)


# let's do one subject, later we will do the loop.
library(AnalyzeFMRI)
path='C:/Users/Chen/Documents/Research/Thickness analysis/Data/cleaned/mask'
setwd(path)
i <- 20
dat.mri <- f.read.nifti.volume(files[i])[,,,1]
dim(dat.mri)

# check how many are
sum(dat.mri!=0)
# get the mask location and plot
mask <- which(dat.mri!=0, arr.ind=T, useNames=F)
dim(mask)

# now get the mid-sagittal slice
mean(range(mask[,1]))
range(mask[,1])
slice.id <- mask[,1]==90
plot(mask[slice.id,2], mask[slice.id,3])


setwd('C:/Users/Chen/Documents/GitHub/Research/ADAMachine/2013 Files/core_code_for_alg')
source('principalcurve.R')
source('lmstep.R')
source('pjstep.R')

setwd('C:/Users/Chen/Documents/GitHub/Research/ADAMachine/2013 Files/temp_output')
# do the PC fitting
dat.x <- mask[slice.id,2]
dat.y <- mask[slice.id,3]
curr.dat <- cbind(dat.x, dat.y)

# use angle to determine the intial parameter
initial.para <- acos((dat.x-mean(dat.x))/((dat.x-mean(dat.x))^2+(dat.y-mean(dat.y))^2)^.5)
initial.para <- ((dat.y-mean(dat.y))>0)*initial.para + ((dat.y-mean(dat.y))<0)*((dat.x-mean(dat.x))<0)*(2*pi-initial.para)+
  ((dat.y-mean(dat.y))<0)*((dat.x-mean(dat.x))>0)*(-initial.para)
initial.para <- (initial.para-min(initial.para))/diff(range(initial.para))


pc <- PrinCurve(curr.dat, sub.size=dim(curr.dat)[1], grid.n=100, 
                max.iter=30, err.tol=0.01,  
                neighbor=0.04, kn=10, initial.para=initial.para)
id <- !(pc$proj$X1)%in%c(0,1)
# get rid of the ending point
proj <- pc$proj[id,1]
curr.dist <- pc$height[id]
grid.x <- seq(0,1, length.out=101)
y <- sapply(grid.x, function(a){
  quantile(curr.dist[abs(a-proj)<0.025],prob=0.95)
})

# 95% quantile is defined for the thickness, PC fitting is biased
ysm <- smooth.spline(grid.x, y, spar=.6)$y


# check the fitting result
ip.col <- initial.para
plot(curr.dat[,1]-mean(curr.dat[,1]), curr.dat[,2]-mean(curr.dat[,2]),
     col=rgb(ip.col,1-ip.col,0), cex=.5, pch=19, asp=1, ylim=range(curr.dat[,2]-mean(curr.dat[,2])))
lines(pc$result$prin.curve[,], lty='dashed', pch=19, cex=1, col='black', lwd=2)


plot((101:1)/101,y, type='b', lty='solid', pch=19, cex=1, col='black', lwd=2)
points((101:1)/101, ysm, pch=19, col='blue')


# Now we do the loop, save the thickness
I <- length(files)
path='C:/Users/Chen/Documents/Research/Thickness analysis/Data/cleaned/mask'
setwd(path)
thick.dat <- NULL
for (i in 1:I){
  dat.mri <- f.read.nifti.volume(files[i])[,,,1]
  # dim(dat.mri)
  # check how many are
  sum(dat.mri!=0)
  # get the mask location and plot
  mask <- which(dat.mri!=0, arr.ind=T, useNames=F)
  # dim(mask)
  # ww get the mid-sagittal slice
  slice.id <- mask[,1]==floor(mean(range(mask[,1])))
  dat.x <- mask[slice.id,2]
  dat.y <- mask[slice.id,3]
  curr.dat <- cbind(dat.x, dat.y)
  # use angle to determine the intial parameter
  initial.para <- acos((dat.x-mean(dat.x))/((dat.x-mean(dat.x))^2+(dat.y-mean(dat.y))^2)^.5)
  initial.para <- ((dat.y-mean(dat.y))>0)*initial.para + ((dat.y-mean(dat.y))<0)*((dat.x-mean(dat.x))<0)*(2*pi-initial.para)+
    ((dat.y-mean(dat.y))<0)*((dat.x-mean(dat.x))>0)*(-initial.para)
  initial.para <- (initial.para-min(initial.para))/diff(range(initial.para))
  pc <- PrinCurve(curr.dat, sub.size=dim(curr.dat)[1], grid.n=100, 
                  max.iter=30, err.tol=0.01,  
                  neighbor=0.04, kn=10, initial.para=initial.para)
  id <- !(pc$proj$X1)%in%c(0,1)
  # get rid of the ending point
  proj <- pc$proj[id,1]
  curr.dist <- pc$height[id]
  grid.x <- seq(0,1, length.out=101)
  y <- sapply(grid.x, function(a){
    quantile(curr.dist[abs(a-proj)<0.025],prob=0.95)
  })
  
  # 95% quantile is defined for the thickness, PC fitting is biased
  ysm <- smooth.spline(grid.x, y, spar=.6)$y
  thick.dat <- cbind(thick.dat, ysm)
}

setwd('C:/Users/Chen/Documents/GitHub/Research/ADAMachine/2013 Files/temp_output')
save(thick.dat, file='thickness_est_dat.RData')
load('thickness_est_dat.RData')

## Now thickness analysis

# load the response Y into R
setwd('C:/Users/Chen/Documents/Research/Thickness analysis/Data')
dat.tab <- read.csv('subject_info.csv', header=T)
names(dat.tab)

# a quick look at the data, we are interested in ados and paness
ados <- dat.tab$ADOS.Total.Score
paness <- dat.tab$Total.PANESS
id <- (ados!=0)&(paness!=0)
plot(ados[id], paness[id])

# first let's do ados
length(ados)
ados.id <- which(!is.na(ados))
thick.ados <- thick.dat[,ados.id]
dim(thick.ados)
gender <- dat.tab$Gender
dx <- dat.tab$Primary.Dx
# match the name
dat.tab$MRI.[ados.id]
files[ados.id]
#check!

#
matplot(thick.ados, type='l')
# now let's do functional regression
library(refund)

PFR.fit <- pfr(Y=ados[ados.id], funcs=t(thick.ados), kz=30, kb=30, method='REML')
beta.fit <- as.vector(PFR.fit$BetaHat[[1]])
matplot(cbind(PFR.fit$BetaHat[[1]][101:1], PFR.fit$Bounds[[1]][101:1,]),
        type = 'l', lty = c(1,2,2), col = c(1,2,2), ylab = "BetaHat", xlab='Along PC',
        main = "beta function", ylim=range(PFR.fit$Bounds))
lines(1:101,rep(0,101), col='blue', lty='dashed')

# not good!!!!!

# do analysis for paness

length(paness)
paness.id <- which(!is.na(paness))
thick.paness <- thick.dat[,paness.id]
matplot(thick.paness, type='l')
dim(thick.paness)
PFR.fit <- pfr(Y=paness[paness.id], funcs=t(thick.paness), kz=10, kb=10, method='REML')
beta.fit <- as.vector(PFR.fit$BetaHat[[1]])
matplot(cbind(PFR.fit$BetaHat[[1]][101:1], PFR.fit$Bounds[[1]][101:1,]),
        type = 'l', lty = c(1,2,2), col = c(1,2,2), ylab = "BetaHat", xlab='Along PC',
        main = "beta function", ylim=range(PFR.fit$Bounds))
lines(1:101,rep(0,101), col='blue', lty='dashed')

# not good either....

# added on 02/27/13
# pointwise analysis and bin analysis
# first we do pointwise


p.value.ados <- apply(thick.ados, 1, function(x){ summary(lm(ados[ados.id]~x))$coefficients[2,4] })
plot(p.value.ados, type='p', pch=19, cex=.5)
p.value.paness <- apply(thick.paness, 1, function(x){ summary(lm(paness[paness.id]~x+dx[paness.id]))$coefficients[2,4] })
plot(p.value.paness, type='p', pch=19, cex=.5)

# calculate the average thickness
ave.thick.ados <- colMeans(thick.ados)
plot(ave.thick.ados, ados[ados.id])
summary(lm(ados[ados.id]~ave.thick.ados))
ave.thick.paness <- colMeans(thick.paness)
plot(ave.thick.paness, paness[paness.id])
summary(lm(paness[paness.id]~ave.thick.paness+dx[paness.id]))

# try bin analysis
Bin <- 5
group.fac <- floor(seq(1,Bin+.9999,length.out=101))
# unlist(by(thick.ados[,1], group.fac, mean))
# ados
thick.group <- apply(thick.ados, 2, function(x){by(data=data.frame(x=x, levl=as.factor(group.fac)), INDICES=group.fac, function(a){mean(a[,1])})})
p.value.ados.bin <- apply(thick.group, 1, function(x){ summary(lm(ados[ados.id]~x))$coefficients[2,4] })
plot(p.value.ados.bin, type='b', pch=19, cex=.5)

# try bin analysis for paness
Bin <- 50
group.fac <- floor(seq(1,Bin+.9999,length.out=101))
# unlist(by(thick.ados[,1], group.fac, mean))
# ados
thick.group <- apply(thick.paness, 2, function(x){by(data=data.frame(x=x, levl=as.factor(group.fac)), INDICES=group.fac, function(a){mean(a[,1])})})
p.value.ados.bin <- apply(thick.group, 1, function(x){ summary(lm(paness[paness.id]~x))$coefficients[2,4] })
plot(p.value.ados.bin, type='b', pch=19, cex=.5)