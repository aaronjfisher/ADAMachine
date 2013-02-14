#Code name PyschAPATH
#Data simulation/generation
#2013 Jan 29 edit - AF

#### Notes 2/13 ###
#I'm no longer dating the output files
#"I" is greater here just to ensure a better resolution truth. Since the output is a grid, it shouldn't matter what I we generate with. We'll need to scale the fitted thickness indeces up a bit though, if we want to compare against the true thicknesses.
#resSquare sets the resolution of the grid. resSquare<-100 gives a 100x100 resolution grid.
####################

setwd("/Users/aaronfisher/Documents/JH/ADAMachine/2013 Files/")
#setwd('C:/Users/Chen/Documents/GitHub/ADAMachine/2013 Files')
#MODEL
#Y_i=β_0 + ∫X_i(t)β(t) +e_i

#shapeType<-'c'; set.seed(34871)
shapeType<-'snail'; set.seed(40871)

library(RColorBrewer)
darkCol<-brewer.pal(10,'Paired')
lightCol<-brewer.pal(10,'Set3')

#NOTE THIS IS NOT THE USUAL NOTATION FOR I????? CHANGE????
N <- 200 #N is the number of curves
I <- 500 #I is the length of the curve (final true resolution will be determined by resSquare)
resSquare<-100 #20 takes<5sec; 60 takes 13sec; 100 takes 36sec (for snail)
#This is the grid resolution of both axes in the final plots




########################
#Generate $k$ signals to mix into the simulated thickness

k <- 7 #number of signals to mix for X
signals <- array(dim=c(k,I))
signals[1,]<-sin((1:I)/I*2*pi)
signals[2,]<-seq(-1,3,length=I)

# plot(signals[2,])

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
######################



########################
#Generate true B(t) & functional model parameters
bsBeta<-t(bs(1:I,df=15,degree=3,intercept=TRUE))/2 # "b spline signals"
betat<-matrix(colSums(bsBeta[c(1:3,12:15),]))
dim(betat)<-c(I,1)
plot(betat)
beta0<-3
sd_e<-1
#######################




########################
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
plot(rwalk(I,amplitude=.2)) #test


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
############################




####################################
#####generate 2-d shapes######
#First define the core shape
#then get a function to go from a grid point to the closest theta angle we've defined.
#The check radius at that theta angle to see if it's consistent with the simulated shape (i.e core shape's spine + or - thickness at that point)

#######
#Core Shape
if(shapeType=='c') {
  shapeTheta<-seq(pi/4,7*pi/4,length=I)
  rSpine<-rep(5,I)
}
if(shapeType=='snail') {
  shapeTheta<-seq(0,2*pi,length=I)
  rSpine<-3.5+4*(1:I/I)
}
maxRange<-max(rSpine)+max(thicknesses) 


########
#Function to recover binned θ for a given x, y coord 
recoverThetaInd<-function(x,y){
  theta.xy<-atan2(x=x,y=y)
  #make sure we go from (0,2pi), not (-pi,pi)
  if(theta.xy<0) theta.xy<- (2*pi + theta.xy)
  binInd.xy<-cut(theta.xy,breaks=shapeTheta,labels=F)
  if(is.na(binInd.xy)) return('out of range')
  thetaBin.xy<-shapeTheta[binInd.xy]
  #we have the floor on the in, but check to see if the ceiling fits better
  if(binInd.xy==I) return(binInd.xy)
  moveUp<-abs(theta.xy-thetaBin.xy) > abs(theta.xy-shapeTheta[binInd.xy+1])
  if(moveUp){
    binInd.xy<-binInd.xy+1
    thetaBin.xy<-shapeTheta[binInd.xy]
  }
  return(binInd.xy)
}
########

########
#Generate Shapes
#recall, resSquare is defined at the top of the file, as the grid resolution on for both axes.
resRange<-seq(-maxRange,maxRange,length=resSquare) #same range for all shape grids
dat.x<-matrix(nrow=N,ncol=resSquare^2)
dat.y<-dat.x
maxITicker<-1

system.time({## time in
print(paste("PROGRESS --  resSqure =",resSquare,'  N =',N))
pb<-txtProgressBar(min = 1, max = N, initial = 1, char = "=", style = 3)
for(i in 1:N){
  iTicker<-1
  #to help weed out grid points:
  maxR.i<-max(rSpine)+max(thicknesses[i,]) #biggest radius
  minR.i<-min(rSpine)-max(thicknesses[i,]) #biggest radius
  for(x in resRange){
    for(y in resRange){
      #First pass
      r.xy<-sqrt(x^2+y^2) #radius
      if(r.xy>maxR.i | r.xy<minR.i) next
      #check if it's in the shape
      binInd.xy<-recoverThetaInd(x=x,y=y)
      if(binInd.xy=='out of range') next
      xyInUp<-r.xy <= rSpine[binInd.xy]+thicknesses[i,binInd.xy]
      xyInDown<- r.xy >= rSpine[binInd.xy]-thicknesses[i,binInd.xy]
      xyIn<-xyInUp&xyInDown
      if(xyIn){
        dat.x[i,iTicker]<-x
        dat.y[i,iTicker]<-y
        iTicker<-iTicker+1
      }
    }
  }
  maxITicker<-max(maxITicker,iTicker)
  setTxtProgressBar(pb,i)
}
})## time out
close(pb)

#trim some fat
dat.x<-dat.x[,1:maxITicker]
dat.y<-dat.y[,1:maxITicker]
#######################################




######################################
#Check your work!
plot(dat.x[1,],dat.y[1,],pch=11,cex=.7,col=darkCol[2])
points(dat.x[2,],dat.y[2,],pch=1,cex=.7,col=darkCol[4])
plot(dat.x[3,],dat.y[3,],pch=18,cex=.7,col=darkCol[2])
points(dat.x[4,],dat.y[4,],pch=11,cex=.7,col=darkCol[4])
plot(dat.x[5,],dat.y[5,],pch=18,cex=.7,col=darkCol[2])
points(dat.x[6,],dat.y[6,],pch=18,cex=.7,col=darkCol[4])
plot(dat.x[7,],dat.y[7,],pch=18,cex=.7,col=darkCol[2])
points(dat.x[8,],dat.y[8,],pch=18,cex=.7,col=darkCol[4])

plotsubj<-function(i){
  par(mfrow=c(1,2))
  plot(thicknesses[i,],type='l',cex=.7,col=darkCol[2],main=paste0('Thickness i=',i))
  plot(dat.x[i,],dat.y[i,],pch=18,cex=.7,col=darkCol[2],main=paste0('C Shape i=',i))
}
plotsubj(1)
plotsubj(2)
plotsubj(7)

save.image(file=paste0('Simulate_',shapeType,'_all_data.RData')) 
save(file=paste0('Simulate_',shapeType,'_only_xyOutcome_data.RData'),list=c("dat.x","dat.y","y"))

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
workshop<-FALSE

if(workshop)
{ 

  
####### OLD NON GRID DATA GEN ###########
####### OLD NON GRID DATA GEN ###########
####### OLD NON GRID DATA GEN ###########
####### OLD NON GRID DATA GEN ###########



resSquare<-50
maxObsPerRow<-max(30)#,sqrt(2)*resSquare)
maxThick<-max(thicknesses) #could change this a bit later
raw.x<-matrix(nrow=N,ncol=I*maxObsPerRow)
raw.y<-matrix(nrow=N,ncol=I*maxObsPerRow)
for(i in 1:N){ #going across subjects
  iTicker<-0
  for(t in 1:I){ #going across angles
    obsPerRow.it<-floor(thicknesses[i,t]/maxThick * maxObsPerRow)
    obsR.it<-seq(rSpine[t]-thicknesses[i,t],rSpine[t]+thicknesses[i,t],length=obsPerRow.it)
    for(o in 1:obsPerRow.it){ #going across observations per angle
      iTicker<-iTicker+1
      raw.x[i,iTicker]<-cos(shapeTheta[t])*obsR.it[o]
      raw.y[i,iTicker]<-sin(shapeTheta[t])*obsR.it[o]
    }
  }
  if(maxObsPerRow*N>5000) print(i)
}

#Remove most of the NAs
all(is.na(raw.y)==is.na(raw.x))
naDatCols<-apply(raw.x,2,function(x) {all(is.na(x))} )
raw.x2<-raw.x[,!naDatCols]
raw.y2<-raw.y[,!naDatCols]
min.xy<-min(c(raw.x2,raw.y2),na.rm=TRUE)
max.xy<-max(c(raw.x2,raw.y2),na.rm=TRUE)

#convert to grid
gridBreaks<-seq(min.xy,max.xy,length=resSquare)
dat.x<-array(dim=dim(raw.x2))
dat.y<-array(dim=dim(raw.y2))
for(i in 1:N){ #cut only works on vectors, not matrices
  dat.x[i,]<-cut(raw.x2[i,],breaks=gridBreaks)
  dat.y[i,]<-cut(raw.y2[i,],breaks=gridBreaks)
  #remove duplicates from lower rez
  dupes.xy<-duplicated(  cbind(dat.x[i,],dat.y[i,]), margin=1)
  dat.x[i,dupes.xy]<-NA
  dat.y[i,dupes.xy]<-NA
  #re-order
  dat.x[i,]<-dat.x[i,order(is.na(dat.x[i,]))]
  dat.y[i,]<-dat.y[i,order(is.na(dat.y[i,]))]
  if(maxObsPerRow*N>5000) print(i)
}
rbind(dat.x[1,],dat.y[1,])
par(mfrow=c(1,2))
plot(raw.x[1,],raw.y[1,],col=darkCol[2])
plot(dat.x[1,],dat.y[1,],col=darkCol[2])
all(is.na(dat.y[1,])==is.na(dat.x[1,]))
naDatCols<-apply(raw.x,2,function(x) {all(is.na(x))} )




############### OTHER STUFF #################
############### OTHER STUFF #################
############### OTHER STUFF #################



cBaseShape.x<-cos(shapeTheta)*rSpine[t]
cBaseShape.y<-sin(shapeTheta)*rSpine[t]
plot(cBaseShape.x,cBaseShape.y)

#Plot some full shapes
#Note, this isn't how we'll generate data. We'll instead
#generate random radiuses for each theta.
r.i.out<- 5+thicknesses[1,]
c.i.x<-cos(shapeTheta)*r.i.out
c.i.y<-sin(shapeTheta)*r.i.out
plot(c.i.x,c.i.y,type='l',col=darkCol[1])
for(i in 1:60){
  for(j in 1:2){
  r.i.out<- 5+thicknesses[i,]*(j==1)-thicknesses[i,]*(j==2)
  c.i.x<-cos(shapeTheta)*r.i.out
  c.i.y<-sin(shapeTheta)*r.i.out
  lines(c.i.x,c.i.y,type='l',col=lightCol[i])
  }
}

}


