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
#fix <- runif(I,.1,.3) #Not sure what this is doing here
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
matplot(t(bsSig),type='l')

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
  int.i<-runif(1,.9,1.1) #some variation in the overall height.
  sigWeights.i<-runif(k,-.2,.2)
  noise.i<-rwalk(I,amplitude=.15)
  thick.i<-int.i+matrix(sigWeights.i,nrow=1)%*%signals + noise.i
  thicknesses[i,]<-thick.i
  
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
plot(dat.x[1,],dat.y[1,],pch=18,cex=.7,col=darkCol[2])
points(dat.x[2,],dat.y[2,],pch=18,cex=.7,col=darkCol[4])
plot(dat.x[3,],dat.y[3,],pch=18,cex=.7,col=darkCol[2])
points(dat.x[4,],dat.y[4,],pch=18,cex=.7,col=darkCol[4])
plot(dat.x[5,],dat.y[5,],pch=18,cex=.7,col=darkCol[2])
points(dat.x[6,],dat.y[6,],pch=18,cex=.7,col=darkCol[4])
plot(dat.x[7,],dat.y[7,],pch=18,cex=.7,col=darkCol[2])
points(dat.x[8,],dat.y[8,],pch=18,cex=.7,col=darkCol[4])


setwd("/Users/aaronfisher/Documents/JH/Biostatistics V & VI/Bst VI Final Project with Chen/2013 Files/")
save.image(file='Simulate_C_all_data_feb4.RData')
save(file='Simulate_C_only_XY_data_feb4.RData',list=list(dat.x,dat.y))

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




