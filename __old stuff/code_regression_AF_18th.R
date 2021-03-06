###########
#START HERE!!

# C DATA

path1 <- '/Users/aaronfisher/Documents/JH/Biostatistics V & VI/Bst VI Final Project with Chen'
setwd(path1)
load('thickness.Rdata')
load('simul_1.RData')
# source('principalcurve.R')
# source('lmstep.R')
# source('pjstep.R')
I <- 300
N <- 500
set.seed(999)
y.out <- rpois(I, 3)
if(diff(dim(thickness))>0) thickness <- t(thickness)


#NOW ROWS ARE SUBJ
#COLS ARE COORDS

#Look at five different subjects
dim(thickness)
plot(thickness[1,]) 
plot(thickness[2,])
plot(thickness[3,])
plot(thickness[4,])
plot(thickness[5,])

########################################################################
########################################################################
########################################################################
#GLMS


#TRY BIN METHOD FOR DIFFERENT BIN SIZES
#Show below the instability of the coef's at higher bin size

p<-dim(thickness)[2]

binmodel<-function(nbins){
  binMat<-matrix(nrow=dim(thickness)[1],ncol=nbins)
  thickbreaks<-seq.int(1,p,length=nbins+1)
  for(i in 1:300){
    for(j in 1:nbins){
      binMat[i,j]<- mean(thickness[i,thickbreaks[j]:thickbreaks[j+1]])
    }
  }
  
  
  mpois<-glm(y.out~binMat,family='poisson') #poisson Model
  if(nbins>1) plot(mpois$coef[-1],type='l',main=paste('β fits,  nbins = ',nbins),xlab='bin',ylab='beta') #not the intercept
  
  if(nbins<10) print(mpois$coef)
  
  #SEE HOW CLOSELY IT'S FITTING
  #par(mfrow=c(2,1))
  if(nbins>1){
    plot(thickness[1,],type='l',main=paste('X[1,] fits,  nbins = ',nbins))
    points(thickbreaks[1:nbins]+p/(2*nbins),binMat[1,],pch=19)
  }
  #plot(thickness[2,])
  #points(thickbreaks[1:nbins]+p/(2*nbins),binMat[2,],pch=19)
  return(mpois)	
}

par(mfrow=c(4,2))
for(numberbins in c(1,5,10,20,40)) {binmodel(numberbins)}

mBIN<-binmodel(5)
betaBIN<-mBIN$coef[-1]
dev.off()
plot(betaBIN)


########################################################################
########################################################################
########################################################################
#PCA


# plot(colMeans(thickness))
# 
# colMeanMat<-matrix(colMeans(thickness),byrow=T,nrow=dim(thickness)[1],ncol=dim(thickness)[2])
# 
# thickZero<-thickness-colMeanMat

thickZero <- scale (thickness, center=T, scale=F)
matplot(t(matrix(rep(1:101, each=300), 300, 101)), t(thickZero), type='l')
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
t(ups[,1])%*%ups[,1] #all ψ_k'ψ_k terms will cancel to 1

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
  lines(thickZero[l,])
}

plot(1:101,ylim=range(thickness))
for(l in 1:101){
  lines(thickness[l,])
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

##########################################################
##########################################################
##########################################################
# AARON'S VERSION, COMMENT THIS OUT

# Using winbugs do the penalized functional regression
# MCMC.result <- bugs(data=list(CJ=C%*%J, y=y.out, K_b=K_b, I=I), inits=list(
  # list(B0=0, beta=rep(0,15), taubeta=.1)), 
                    # parameters.to.save=c('beta', 'B0', 'taubeta'),
                    # model.file="bugcode.txt",n.chains=1, n.burnin=10000, n.thin=1,
                    # n.iter=20000, bugs.seed=999, bugs.directory="C:/Users/Chen/Documents/R/winbugs14/",
                    # working.directory = path1, debug=F)
# save(MCMC.result, file='C_shape.RData')


##########################################################
##########################################################
##########################################################
load( file='C_shape.RData')
H<-MCMC.result$n.keep
p<-101

# Check the chains for beta
#just functions for plotting here.
plotInd<-floor(seq(1,H,length=min(10000,H)))
p<-101 #length of the curve t=1,...p
name <- sapply(1:15, function(x){paste('beta',as.character(x))})
layout(matrix(1:4,2,2))
for (i in 1:15){
  plot(MCMC.result$sims.list$beta[,i][plotInd], type='l', ylab='', main=name[i])
}
# beta mixed very well !

# Check the chains for B0 and taubeta
plot(MCMC.result$sims.list$B0[plotInd], type='l', ylab='', main='B0')
plot(MCMC.result$sims.list$taubeta[plotInd], type='l', ylab='', main=expression(tau[beta]))
plot(1/MCMC.result$sims.list$taubeta[plotInd], type='l', ylab='', main=expression(1/tau[beta]))

# B0 mixed well while taubeta is not quite good


# Reconstruct all the curves
beta.hist <- MCMC.result$sims.list$beta
dim(beta.hist)
dim(bsplines)

height.hist <- beta.hist%*%t(bsplines)

sigma <- var(height.hist)
sd_t<-sqrt(diag(sigma))
# check the correlation coefficient of all the curve grid points
image(cov2cor(sigma))

#plot some histories
par(mfrow=c(1,2))
matplot(t(height.hist[1:500,]),type='l',main='raw') #plots the columns

BtEst<-colMeans(height.hist)
Bt.means<-array(rep(BtEst,each=H),dim=dim(height.hist))
Bt.hist0<-height.hist-Bt.means
matplot(t(Bt.hist0[1:500,]),type='l',main='zeroed') 

#Uncomment these to do Ciprian's Way
#library(MASS)
#Bt.hist0<-mvrnorm(50000,rep(0,101),sigma)



Bt.sds<-array(rep(sd_t,each=H),dim=dim(height.hist))
d<- abs(Bt.hist0/Bt.sds)

par(mfrow=c(1,1))
bigd<-apply(d,1,max)
band95dval<-quantile(bigd,.95)
band95width_t<-band95dval*sd_t
plot(band95width_t)

matplot(t(d[1:500,]),type='l',main='standardized & zeroed') 
abline(h=band95dval,lwd=2)

matplot(t(Bt.hist0[1:500,]),type='l',main='zeroed') 
lines(band95width_t,lwd=2)
lines(-band95width_t,lwd=2)

##########################
######## COOL PLOT!!!
matplot(t(height.hist[1:800,]),lwd=.3 ,type='l',main='Target B draws, with bands & pointwise CIs',ylim=c(-20,40),xlab='t',ylab='B(t)',) #plots the 
lines(BtEst,lwd=2,col='black')
lines(BtEst+band95width_t,lty=2,col='red',lwd=2,)
lines(BtEst-band95width_t,lty=2,col='red',lwd=2)
abline(h=0)
d.pt.L<- apply(height.hist,2,function(x) quantile(x,probs=.025)) #pointwise lower
d.pt.U<- apply(height.hist,2,function(x) quantile(x,probs=.975)) 
#lines(BtEst+2*sd_t,lwd=2)
#lines(BtEst-2*sd_t,lwd=2)
lines(d.pt.L,col='black',lwd=2,lty=2)
lines(d.pt.U,col='black',lwd=2,lty=2)
legend('topleft',c('Estimate','Pointwise CI','Band'),lty=c(1,2,2),col=c('black','black','red'),lwd=2)
#########################

#dev.copy2pdf(file='Bands_C_Shape_12-18-12_1.pdf')

band.dval.i<-quantile(bigd,c(.5,.6,.7,.8,.9,.95))
band.i.width.t<-matrix(nrow=6,ncol=p)
for(i in 1:6) band.i.width.t[i,]<-band.dval.i[i]*sd_t
row.names(band.i.width.t)<-c("50%","60%","70%","80%","90%","95%")
band.i.upper<-band.i.width.t #just to capture names and dimension
band.i.lower<-band.i.width.t #just to capture names and dimension
for(i in 1:6){
	band.i.upper[i,]<-	BtEst+band.i.width.t[i,]
	band.i.lower[i,]<-	BtEst-band.i.width.t[i,]
}
matlines(cbind(t(band.i.upper),t(band.i.lower)),type='l')

dev.off()

#dev.copy2pdf(file='Chart to show bands -- 5 shape.pdf')
save(list=c('band.i.upper','band.i.lower','BtEst'),file='bands_&_Beta_Estimate_C_Shape.RData')
#save.image('Final_Results_Rdata_Five')







#Fit unpenalized model
# mbs<-glm(y.out~1+C %*% J,family='poisson') #bspline model
# summary(mbs)
# 
# #recover beta(t)
# betaBS<-rep(0,p)
# for(kb in 1:K_b) betaBS<-betaBS+mbs$coef[kb+1]*bsplines[,kb]
# 
# par(mfrow=c(2,1))
# plot(rep(betaBIN,each=floor(p/length(betaBIN))),type='l',main='Five bin model')
# plot(betaPCA,type='l',col='red',main='PCA β basis (red) and B-Spline β basis (blue)')
# lines(betaBS,type='l',col='blue')
# #Totally different scale??? What's up with that???
# #NOTE!!!!! When K_x \neq K_b, these graphs will look much different
# 
# 
# #JEFF'S FUNCTION
# library(refund)
# result <- pfr(y.out, funcs=thickness, kz=30, kb=30, family='poisson', method=REML, smooth.cov=T)
# 
# 
# ######################################################
# ######################################################
# ######################################################
# #WINBUGS
# library(R2WinBUGS)
# 
# 
# #NO ERROR VERSION
# dat4bugs_noError<-list(y.out,C%*%J,I,K_b)
# names(dat4bugs_noError)<-c('y',"CJ",'N','K_b')
# 
# save.image('forChen.Rdata')
# load('forChen.Rdata')
# 
# #bugs(data=dat4bugs_noError, inits=list(B0=0, beta=c(0,0,0,0,0,0,0,0,0,0), taubeta=1) , ...)
# 
# 
# 
# #WITH MEASUREMENT ERROR
# #this code still has some errors in it
# dat4winbugs<-list(y.out,betaPCA,0,1,C,C,J,I,K_b,K_x,ups[,1:K_x])
# names(dat4winbugs)<-c('y','beta','B0','taubeta','WC',"C","J",'N_subj','K_b','K_x','ups') #if we do measure error, this is just WC, not C
# bugs.data(data=dat4winbugs,data.file='r2wintest.txt')

######################################################
######################################################
######################################################
# SPLINES WITH BASIC PENALTY
#HOW DID JEFF DO THIS WITH NMLE?

# library(nmle)

# CJ<-C %*% J
# id<- 1:I
# mlme<-glmer(y.out~1+(CJ[,]|id),family='poisson') #bspline model, where 
# summary(mbs)
