#NOTES & MAJOR STUFF


####
#Errors: Something's wrong with permutation test? The cutoff value is off the charts tiny
#The beta estimates from redoing voxelwise regressions on X data reconstructed from the principle components are coming up shifted down (centered?) why??? Shouldn't the beta's be the same, even if the projected stuff is lower?
	#tag-aosdfoao

####
#Options:
useRandInt = F  #returns a standard lm analysis (basic iid errors)
#Whether or not to use tInd in the permutation test section.

#OTHER NOTES ON CHANGES
#5/27/13 - Chen was standardizing FA values, which, biologically, is a little crazy.
#I think the reason why the beta estimates from the basis (functional) models are centered around zero is because they're all adjusting for eachother.. while in the voxelwise regression they aren't.
#He had code: 











########################################################
##############       OLD   FUNCTIONS      
########################################################


# Load Aaron's little functions for plotting with colorspace
library(colorspace)

disp<-function(){ #make new window, cross mac v pc
  if(.Platform$OS.type=="windows") windows()
  else quartz()
}

pal <- function(col, border = "light gray", ...)
#Copy pasted from HCL-Based Color Palettes in R
#http://cran.r-project.org/web/packages/colorspace/vignettes/hcl-colors.pdf
 {
 n <- length(col)
 plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
 axes = FALSE, xlab = "", ylab = "", ...)
 rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)
 }
#rainbow_hcl(10)
#pal(rainbow_hcl(50,c=100,l=80))
#pal(sequential_hcl(50))

c.r.hcl<-function(x,n=1000, ...){ #cut rainbow_hcl
	xCut<-cut(x,breaks=n)
	colors<-rainbow_hcl(n=n, ...)
	out<-colors[xCut]
	return(out)
}
c.s.hcl<-function(x,n=1000,only.up.to=n, ...){ #cut sequantial_hcl
	xCut<-1 #these two lines guard against errors caused when this function is piped into c.d.hcl
	if(only.up.to>1) xCut <-cut(x,breaks=only.up.to)
	colors<-sequential_hcl(n=n, ...)[n:(n-only.up.to)]#reverse the order
	#if you don't want the full sequence, this will chop it short!
	out<-colors[xCut]
	return(out)
}

c.d.hcl<-function(x,n=1000, h=c(0,260),c=80,...){ #cut divergent_hcl
#c is the max chroma, must be zero where they meet in the middle
#hues: 0 is low, 260 is the higher.
	xNegInd <- which(x<0)
	xPosInd <- which(x>=0)
	
	nNeg <- length(xNegInd) #so we only get more faded parts on the side of the spectrum that doesn't have the same magnitude
	nPos <- length(xPosInd)
	biggerN<-max(nPos,nNeg)
	
	out<-rep(0,length(x))	
	if(!length(xNegInd)==0){
		xNegCol<-c.s.hcl(abs(x[xNegInd]),n=biggerN,only.up.to=nNeg,h=h[1],c=c(c,0),...)
		out[xNegInd]<-xNegCol
	}
	if(!length(xPosInd)==0){
		xPosCol<-c.s.hcl(x[xPosInd],n=biggerN,only.up.to=nPos,h=h[2],c=c(c,0),...)
		out[xPosInd]<-xPosCol
	}

	return(out)
}

fastSvdWide<-function(A,talk=FALSE){ 
#I future, this function could be adjusted to do matrix multiplication more efficiently! Needed for super massive !A!
#takes a wide centered Nxp matrix with N<<p, and does an svd for it
#Let A=UDV'; where A=(Nxp); D=(NxN) diagonal; V'=(Nxp); V=(pxN)
	N<-dim(A)[1]
	p<-dim(A)[2]
	if(p < N){ print('error - must be a wide matrix'); return()}
	if(talk) print("Let's get ready to party!")
	#get AA'
	if(talk){
		maxIndex<-N
		pb<-txtProgressBar(min = 1, max = maxIndex,char = "=",  style = 3)
	}
	AAt<- A %*% t(A)
	if(talk) print("Got AA^T!!")
	#get svd for AA' = UDV'VDU'=U D^2 U'
	svdAAt<-svd(AAt)
	#A=UDV' => (D^-1)(U^-1)A=V' => (D^-1)(U')A=V'
	#svd for AA'= XYZ'
	#with UD^2 U' = XYZ', so X=D^2, and Z=X=U
	#Technically the U above might be not unique up to sign (it may be equal to Uδ, where δ is an identity matrix with some negative ones that make some columns neg. However, in this case, 
	#(Uδ')D(δU')=UδDδU' b/c δ'=δ
	#=UDU' b/c δ'δ=I
	#=So our method would give us:
	#t(Uδ)DInv A = δV, so we'd end up with a totally valid V, just with some different signed vectors also.
	D<-sqrt(svdAAt$d)  #note, D is a vector here
	if(talk) print('Got svd(AA^T)!!!')
	DInv<-1/D
	U<- svdAAt$u
	DInvUt<- diag(DInv) %*% t(U)
	Vt<- DInvUt %*% A
	V<-t(Vt)
	return(list(v=V,u=U,d=D))
}


getBetaSigma<-function(Sigma.phi,basis,diag.only=TRUE){
  p<-length(basis)
  V<-matrix(nrow=p,ncol=p)
  if(diag.only){
    V<-rep(NA,p)
    for(i in 1:p) V[i]<-t(basis[i,]) %*% Sigma.phi %*% basis[i,]
    return(V)
  }
  for(i in 1:p){
    for(j in 1:p){
      V[i,j]<- t(basis[i,]) %*% Sigma.phi %*% basis[j,]  
    }
  }
  return(V)
}

#A<-matrix(rnorm(3*4),ncol=4)


imageVec<-function(v,color='other', col=c.s.hcl(1:100,h=150),...){
	p<-length(v)
	vMat<-matrix(v, sqrt(p),sqrt(p))
	if(color=='d') col<-c.d.hcl(  seq(min(v,na.rm=T),max(v,na.rm=T),length=100),...  )
	if(color=='s') col<-c.s.hcl(  seq(min(v,na.rm=T),max(v,na.rm=T),length=100),...  )
	if(color =='other') image(vMat,  asp=1, bty='n', xlab='', ylab='', xaxt='n', yaxt='n',col=col)
	if(color!='other') image(vMat,col=vCol,  asp=1, bty='n', xlab='', ylab='', xaxt='n', yaxt='n',col=col)
}









###################################################
###########          READ DATA         ############
###################################################

# Reading data for FA analysis of CC
# written by Chen on Feb. 27th, 2013, for Aaron. Thanks Chen!
# A few edits from Aaron


# load the image data
# change your working directory
# setwd('C:/Users/Chen/Dropbox/CC_surface/paper/Paper version 020813')
setwd('/Users/aaronfisher/Documents/JH/ADAMachine/FA analysis/data')
load('PSCC.RData')
ls()
# we are interested in the fa.map
# let try to plot FA image of several scans
# Each row of fa.map is a subject
# it's a carpet of 101x101=10201 grid points. 
# The numbers go down the rows, and the over to the next column 
# and down again, the same way R orders it [1,1], [1,2]...[1,101], [2,1].... [101,101]

#Check the ordering of fa.map
#HERE SWITCH fa.map so that SUBJECTS ARE ROWS AND COLS ARE VOXELS
if(dim(fa.map)[1]>dim(fa.map)[2]) fa.map<-t(fa.map)
dim(fa.map) #again, each row is an entry
mean(is.na(fa.map)) #no NA's
N<-dim(fa.map)[1] #N=466 
p<-dim(fa.map)[2] #p=10201
orderFaMap<-1:(101^2)
matrix(fa.map[1,], nrow=101,ncol=101)[1:4,1:4]  #this is how we map back to matrixes, or "carpets"
matrix(orderFaMap, 101,101)[1:4,1:4]


#lets visualize some of these carpets
library(colorspace)
library(plotrix)
set.seed(999)
# x <- seq(0, 1, length.out=101) #these lines turn out to be redundant w/ image() defaults
# y <- seq(0, 1, length.out=101)
layout(matrix(1:4,2,2))
samp <- sample(c(1:466),4) #pick 4 people from the data for a demo
for (i in samp){
  aint <- matrix(fa.map[i,], 101,101) #removed a t() from Chen's code here
  image(aint, asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.s.hcl(1:100),main=paste('Subj',i))
  color.legend(1.01, 0,1.05,1,legend=paste0('FA=',signif(range(aint),2)), rect.col=c.s.hcl(1:100), align='rb',gradient='y', cex=.5)
}
par(mfrow=c(1,1))

## now load the data of EDSS & identity data
edss <- c(t(read.table('edss_data.txt'))) #length of N = 466
load(file='id.RData') #contains a vector "identity" with the ID of patients (repeated measures), and ALSO the time of visit
#first 5 digits indicate subject id, last 2 digits indicate visits. I forgot to send this to you last time.
head(identity) #strings where first 5 digits are id, last 2 are visit number
id<-sapply(identity,FUN=function(x) {substr(x,1,5)}  )
names(id)<-c()
time<-sapply(identity,FUN=function(x) {substr(x,6,7)}  )
time<-as.numeric(time)
head(cbind(id,time),20)

n.id<-length(unique(id)) # #of people
length(edss) #=N

mean(is.na(edss))    # 3.4% missing
length(unique(id))   #176 subjects
hist(time)           #distribution of times
hist(c(edss))
summary(edss)

idObs<-id[!is.na(edss)]
edssObs<-edss[!is.na(edss)]
fa.map.obs<-fa.map[!is.na(edss),]




################################################
#Look at reliability of surfaces across time
#compare image[subj=i,visit=1] to image[i,visit=2:3] and also to image[j≠i, 2:3]

#find three images from three subjects, compare subtraction images against 1st subject.
id[1:50]
#compare col 2:4 against 1, 11:12 against 10, 21:23 against 24, 29:31 against 32, and 36:38 against 39
compareVector<-c(2:4,11:13,21:23,29:31,36:38)
matrix(id[compareVector],nrow=3)
refInd<-32
id[refInd]
par(mfcol=c(3,5),mar=c(2,1,2,5))
#just plain stuff
for (i in compareVector){
  aint <- matrix(fa.map[i,], 101,101) #removed a t() from Chen's code here
  tempcol<-c.s.hcl(seq(min(aint),max(aint),length=100))
  image(aint, asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=tempcol,main=paste('FA matrix row',i,'\n id',id[i]))
  color.legend(1.01, 0,1.05,1,legend=paste0('FA=',signif(range(aint),2)), rect.col=tempcol, align='rb',gradient='y', cex=.5)
}
for (i in compareVector){
  aint <- matrix(fa.map[i,]-fa.map[refInd,], 101,101) #removed a t() from Chen's code here
  tempcol<-c.d.hcl(seq(min(aint),max(aint),length=100))
  image(aint, asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=tempcol,main=paste('FA matrix row',i,'\n id',id[i]))
  color.legend(1.01, 0,1.05,1,legend=paste0('FA=',signif(range(aint),2)), rect.col=tempcol, align='rb',gradient='y', cex=.5)
}
for (i in compareVector){
  hist(fa.map[i,]-fa.map[refInd,],main=paste('FA matrix row',i,'absdiff',round(sum(abs(fa.map[i,]-fa.map[refInd,])))),xlim=c(-.4,.4),breaks=30)
}
abline(0,1)
par(mfrow=c(1,1))
#within subj seems normal?
qqnorm(scale(fa.map[1,]-fa.map[2,]))
abline(0,1)
qqnorm(scale(fa.map[11,]-fa.map[12,]))
abline(0,1)
qqnorm(scale(fa.map[19,]-fa.map[20,]))
abline(0,1)

#across subj also seems normal?
qqnorm(scale(fa.map[1,]-fa.map[20,]))
abline(0,1)





########################################################
################          BASIC ANALYSIS      
########################################################



library(lme4)
##############################


#Function for voxel-wise regression w/ optional random intercept
#useRandInt = F  returns a standard lm analysis (basic iid errors)
#useRandInt = T  uses lme4 to add a random intercept
#Need to add "if(skip]v]) next " line so that you can use it in the telescope function
vlm<-function(x,y,useRandInt=F,id=NA,talk=TRUE,skip='none'){ #voxelwise lm
  N<-dim(x)[1] #N=466 
  p<-dim(x)[2] #p=10201
  if(skip[1]=='none') skip <- rep(FALSE, p)
  
  #set up vectors to record results
  ests.voxelwise<-rep(NA,times=p)
  stdErrs.voxelwise<-rep(NA,times=p)
  tvals.voxelwise<-rep(NA,times=p)
  pvals.voxelwise<-rep(NA,times=p)

  if(talk){
  	print(paste("PROGRESS --  p =",p,'  N =',N))
  	pb<-txtProgressBar(min = 1, max = p,  char = "=", style = 3)
  }

  for(v in 1:p){
  	if(skip[v]) next 
  	if(useRandInt){
  		m.v<-lmer(y~ (1|id) + x[,v]) #model for i.th voxel
  		m.v.summary<- attributes(summary(m.v))$coef[2,]
  	}
  	if(!useRandInt){
  		m.v 				<-lm(y~x[,v])
  		m.v.summary			<- summary(m.v)$coef[2,]
  		pvals.voxelwise[v]  <- m.v.summary[4]
  	}
  	ests.voxelwise[v]    <- m.v.summary[1]
  	stdErrs.voxelwise[v] <- m.v.summary[2]
  	tvals.voxelwise[v]   <- m.v.summary[3]
  	
  	if(talk) setTxtProgressBar(pb,v)
  }
  output<-list("ests.voxelwise"=ests.voxelwise, "stdErrs.voxelwise"=stdErrs.voxelwise, "tvals.voxelwise"=tvals.voxelwise,"pvals.voxelwise"=pvals.voxelwise)
  return(output)
}

#Run voxelwise regression, store results.
vlm.basic<-vlm(x=fa.map,y=edss,useRandInt=useRandInt,id=NA,talk=TRUE)

ests.v<-vlm.basic$ests.voxelwise
stdErrs.v<-vlm.basic$stdErrs.voxe 
tvals.v<-vlm.basic$tvals.voxelw
pvals.v<-vlm.basic$pvals.vox


#a quick function to permute data randomly, just once
#Permute outcomes function, maintaining ids. (a function for it)
#by default assume all id's are different
#nApprox = how many times should we permute it?
#When possible, we will NOT use replacement when reassigning. Replacement will
#essentially serve as interpolation
#id won't change, but the values will switch around.
permuteY<-function(y,id=1:length(y)){
	#potential errors
	if(sum(is.na(y) | is.na(id))  > 0  ){
		print('Please do not enter any NA values in "y" or "id"')
		return()
	}
	if(length(y)!=length(id) ){
		print('"y" and "id" must have same length')
		return()
	}

	#permute
	N  		<-length(y)
	uId  	<-unique(id)
	nGroups <-length(uId)
	newOrder  <-sample(nGroups)
	yMixed    <-rep(NA,N)

	for(i in 1:nGroups){ #for each id, take their values and assign them to the new person
		#sample with replacement if nessecary, sample without replacement if possible.
		oldPersonInd <- which(id==uId[i])
		newPersonInd <- which(id==uId[newOrder[i]])

		newPersonLength <- length(newPersonInd)
		oldPersonLength <- length(oldPersonInd)

		oldPersonVals   <- y[oldPersonInd]

		if(newPersonLength<=oldPersonLength) newPersonVals<-sample(oldPersonVals, size=newPersonLength,replace=FALSE)
		if(newPersonLength>oldPersonLength) newPersonVals<-sample(oldPersonVals, size=newPersonLength,replace=TRUE)
		#adjust for the case when oldPersonVals has length 1, and sample() treats it like an integer rather than a vector
		if(length(oldPersonVals)==1) newPersonVals<- rep(oldPersonVals,newPersonLength)

		yMixed[newPersonInd]<- newPersonVals
	}
	
	return(yMixed)	
}





#small scale simulation to see if permutation tests work well
#they seem to, not sure why it's so low power with voxelwise.
x<-matrix(rnorm(200*300),nrow=200)
y<-x[,1]+rnorm(200)

nPerms<-200
minPvalPerm<-rep(NA,nPerms)

pb<-txtProgressBar(min = 1, max = nPerms,  char = "=", style = 3)
for(r in 1:nPerms){
	setTxtProgressBar(pb,r)

	y.r<-permuteY(y)#,id=idObs)
	vlm.basic.r<-vlm(x=x,y=y.r,useRandInt=F,id=NA,talk=FALSE)
	minPvalPerm[r]<-min(vlm.basic.r$pvals.vox)
}
quantile(log(minPvalPerm),c(.05,.1,.2,.3))
min(log(vlm(x=x,y=y,useRandInt=F,id=NA,talk=FALSE)$pvals))
log(.05/200)
#OK, here the permutation test is more powerful







##############################
#get dist of max pvals
setwd('/Users/aaronfisher/Documents/JH/ADAMachine/FA analysis')

#test index, to make sure it works OK
tInd<- (101*50+1):(101*51)
matplot(t(fa.map[2:20,tInd]),type='l')
tInd<-(101*50+1):(101*52)
tInd<-1:10201


######  This will take a long time, so we'll just deactivate this chunk for now:
if(FALSE){

nPerms<-20
minPvalPerm<-rep(NA,nPerms)

print(paste("PROGRESS --  nPerms =",nPerms))
pb<-txtProgressBar(min = 1, max = nPerms,  char = "=", style = 3)
for(r in 1:nPerms){
	setTxtProgressBar(pb,r)

	y.r<-permuteY(edssObs)#,id=idObs)
	vlm.basic.r<-vlm(x=fa.map.obs[,tInd],y=y.r,useRandInt=F,id=NA,talk=FALSE)
	min(vlm.basic.r$pvals.vox)
	minPvalPerm[r]<-min(vlm.basic.r$pvals.vox)
}
quantile(log(minPvalPerm,.05))
min(log(pvals.v[tInd]))
setwd('/Users/aaronfisher/Documents/JH/ADAMachine/FA analysis')
#save(list='minPvalPerm',file='Permutation_test_p-values_voxelwiselm_2013-05-27.RData')
}
######


#save(list='minPvalPerm',file='Permutation_test_p-values_voxelwiselm_2013-05-27.RData')
load(file='Permutation_test_p-values_voxelwiselm_2013-05-27.RData')



##########################
#Plot results


setwd('/Users/aaronfisher/Documents/JH/ADAMachine/FA analysis')

hist(pvals.v) # Huh, *lots* of signal?
hist(log(pvals.v))
abline(v=log(.05),lwd=2,col='blue')   #signif cutoff line
abline(v=log(.05/p),lwd=2,col='blue') #bonferonni correction
abline(v=log(quantile(minPvalPerm,.05)),lwd=2,col='blue') #permutation test correction
mean(pvals.v<.05) #43 are pointwise sig? Wow!
#dev.copy2pdf(file='Hist of p-values.pdf')
#dev.off()

hist(ests.v) #almost all negative

pal(c.s.hcl(100:1))


png(file='plots/voxelwise pvals -- tight layout -- 2013-5-27.png',height=420,width=1750)
#layout(matrix(c(1,3,2,4,6,5),2,3))
par(mfrow=c(1,5),mar=c(2,1.9,4.3,4.8))

maincex<-2.3

image(matrix(fa.map[1,],101,101), asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.s.hcl(1:100), main='Sample Subject FA Map \n(subject 1)',cex.main=maincex)
color.legend(1.01, 0,1.05,1,legend=signif(range(fa.map[1,]),2), rect.col=c.s.hcl(1:100), align='rb',gradient='y', cex=1.15)

estsMat.v<-matrix(ests.v,101,101)
estsRange<-seq(min(ests.v),max(ests.v),length=100)
image(estsMat.v, asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.d.hcl(estsRange,h=c(295,40)) ,main='Beta Estimates',cex.main=maincex)
color.legend(1.01, 0,1.05,1,legend=c(round(range(ests.v),digits=2)), rect.col=c.d.hcl(estsRange,h=c(295,40)), align='rb',gradient='y', cex=1.15)

pvalsMat.v<-matrix(pvals.v,101,101)
image(pvalsMat.v, asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.s.hcl(100:1,h=150),main='p-values',cex.main=maincex)
color.legend(1.01, 0,1.05,1,legend=c('0','1'), rect.col=c.s.hcl(100:1,h=150), align='rb',gradient='y', cex=1.15)

image(pvalsMat.v>=.05, asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.s.hcl(100:1,h=150),main='p-values<.05',cex.main=maincex)
color.legend(1.01, 0,1.05,1,legend=c('<.05','≥.05'), rect.col=c.s.hcl(2:1,h=150), align='rb',gradient='y', cex=1.15)

image(log(pvalsMat.v), asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.s.hcl((100:1),h=150),main='log(p-values)',cex.main=maincex)
color.legend(1.01, 0,1.05,1,legend=round(range(log(pvals.v))), rect.col=c.s.hcl(100:1,h=150), align='rb',gradient='y', cex=1.15)

dev.off()
#dev.copy2pdf(file='voxelwise pvals - better labels 2013-4-4.pdf')
#dev.copy2pdf(file='voxelwise pvals 2013-3-8.pdf')
#dev.copy2pdf(file='voxelwise pvals 2013-2-27.pdf')



#######################
#abstract verion
makeAbstract<-FALSE
if(makeAbstract){
	png(file='plots/voxelwise pvals -- abstract 2013-05-30.png',height=400,width=1750)
		par(mfrow=c(1,5),mar=c(.4,.4,4,.4))

		image(matrix(fa.map[1,],101,101), asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.s.hcl(1:100), main='')

		estsMat.v<-matrix(ests.v,101,101)
		estsRange<-seq(min(ests.v),max(ests.v),length=100)
		image(estsMat.v, asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.d.hcl(estsRange,h=c(295,40)) ,main='')

		pvalsMat.v<-matrix(pvals.v,101,101)
		image(pvalsMat.v, asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.s.hcl(100:1,h=150),main='')

		image(pvalsMat.v>.05, asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.s.hcl(100:1,h=150),main='')

		image(log(pvalsMat.v), asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.s.hcl((100:1),h=150),main='')

	dev.off()
#	dev.copy2pdf(file='voxelwise pvals - abstract -- 2013-4-4.pdf')
}















####################################################################
##########     Multi-resolution stuff
####################################################################




#funciton to pixelize a design matrix (no analysis yet)
#takes a WIDE form n by p^2 matrix
#resolution is a vector of resolution to form the sqrt of the total pixels in the desired image
#outputs a list of pixelized wide design matrices
#resolution of 1 or p (returning the full matrix) are fine to enter here
#returnMat will have the function just return the FIRST pixelizing matrix. Then totally stop.
pixelize<-function(x,resolution=c(2),talk=TRUE,returnMat=FALSE) {
	xRes<-list()
	p<-dim(x)[2]
	N<-dim(x)[1]
	for(r in 1:length(resolution)){
		res<-resolution[r]
		if(talk) print(paste0('working on res=',res))
		charRes<-as.character(res)
		xRes[[charRes]]<-matrix(NA,nrow=N,ncol=res^2)
		if(res^2>=p) {xRes[[charRes]]<-x; next}

		#build a vector, then a mat, then a vector, then a mat to multiply x by.
		rIndVec1<-round(seq(1,res,length=sqrt(p) ))
		rIndMat<- outer(rIndVec1, rIndVec1, FUN=paste)
		rInd<-c(rIndMat)

		pixelizeMat<- matrix(nrow=p,ncol=res^2)
		for(m in 1:(res^2) ){
			thisSquare<-as.numeric(  rInd== unique(rInd)[m]  )
			thisSquare<-thisSquare/sum(thisSquare)
			pixelizeMat[,m]<-thisSquare
		}

		if(returnMat) return(pixelizeMat)
		xRes[[charRes]] <- x %*% pixelizeMat

	}
	return(xRes)
}

resolution<-c(1,2,4,8,16,32,64,101) #sqrt root of total # pixels per data version
xRes<-pixelize(x=fa.map,resolution=resolution, talk=TRUE)
names(xRes)

#plot em!
par(mfrow=c(1,length(resolution)-1),mar=c(1,1,4,5))
for(r in 2:length(resolution)){
	tempMat<-  matrix(xRes[[r]][9,],resolution[r],resolution[r])  
	image(tempMat, asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.s.hcl(1:100), main=paste0('res=',resolution[r],'x',resolution[r]))
	xPosition<-1+1/resolution[r]+.01
	color.legend(xPosition, 0,xPosition+.05,1,legend=signif(range(tempMat[1,]),2), rect.col=c.s.hcl(1:100), align='rb',gradient='y', cex=.7)
}

#Do new voxelwise regressions for lower res images
ests.v.multires<-list() #v for voxelwise, one entry in the list level for each resolution
pvals.v.multires<-list()
for(r in 1:length(resolution)){
	res<-resolution[r]
	if(res==1) {
		resOut<-summary(lm(edss~xRes[[r]]))$coef
		ests.v.multires[[r]]<-resOut[2,1]
		pvals.v.multires[[r]]<-resOut[2,4]
		next
	}
	if(res>1){
		resOut<-vlm(x=xRes[[r]],y=edss,useRandInt=useRandInt,id=NA,talk=TRUE)
		ests.v.multires[[r]]<-resOut$ests.voxelwise
		pvals.v.multires[[r]]<-resOut$pvals.voxelwise
	}
}



##################################
#Simplest case -- just regress against mean across image?
ests.v.multires[[1]]
pvals.v.multires[[1]]
plot(xRes[[1]],jitter(edss))
abline(coefficients(lm(edss~xRes[[1]])))

ests.v.multires[[2]]
pvals.v.multires[[2]]

##################################


##################################
##########     REDO COOL PLOTS
##################################




setwd('/Users/aaronfisher/Documents/JH/ADAMachine/FA analysis')
png(file='plots/voxelwise multires pvals -- 2013-5-29.png',height=2000,width=2000)
par(mfrow=c(5,5),mar=c(1,1,4,7))
maincex<-2.55

for(r in c(2:5,8)){
	res<-resolution[r]
	xPosition<-1+1.01/resolution[r]

	image(matrix(xRes[[r]][1,],res,res), asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.s.hcl(1:100), main='Sample Subject FA Map \n(subject 1)',cex.main=maincex)
	color.legend(xPosition, 0,xPosition+.05,1,legend=signif(range(xRes[[r]][1,]),2), rect.col=c.s.hcl(1:100), align='rb',gradient='y', cex=maincex/2)

	estsMat.v<-matrix(ests.v.multires[[r]],res,res)
	estsRange<-seq(min(ests.v.multires[[r]]),max(ests.v.multires[[r]]),length=100)
	image(estsMat.v, asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.d.hcl(estsRange,h=c(295,40)) ,main='Estimates',cex.main=maincex)
	color.legend(xPosition, 0,xPosition+.05,1,legend=c(round(range(ests.v.multires[[r]]),digits=2)), rect.col=c.d.hcl(estsRange,h=c(295,40)), align='rb',gradient='y', cex=maincex/2)

	pvalsMat.v<-matrix(pvals.v.multires[[r]],res,res)
	image(pvalsMat.v, asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.s.hcl(100:1,h=150),main='p-values',cex.main=maincex)
	color.legend(xPosition, 0,xPosition+.05,1,legend=c('0','1'), rect.col=c.s.hcl(100:1,h=150), align='rb',gradient='y', cex=maincex/2)

	image(pvalsMat.v>.05, asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.s.hcl(2:1,h=150),main='p-values<.05',cex.main=maincex)
	color.legend(xPosition, 0,xPosition+.05,1,legend=c('<.05','>.05'), rect.col=c.s.hcl(2:1,h=150), align='rb',gradient='y', cex=maincex/2)

	image(log(pvalsMat.v), asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.s.hcl((100:1),h=150),main='log(p-values)',cex.main=maincex)
	color.legend(xPosition, 0,xPosition+.05,1,legend=round(range(log(pvals.v.multires[[r]]))), rect.col=c.s.hcl(100:1,h=150), align='rb',gradient='y', cex=maincex/2)
}

dev.off()





##################################
##########     Looking at why that lower left square is so significant
##################################
#in image(), the lower left is the 1,1 index. (near zero/origin)


idNum1<-!duplicated(id)
idNum2<-id
idNum2[!idNum1]<-NA
idNum2[idNum1]<-1:sum(idNum1)
idNum2<-as.numeric(idNum2)
for(i in 2:length(idNum2)) if(is.na(idNum2[i])) idNum2[i]<-idNum2[i-1]

#In lowest resolution
#RES = 1
plot(xRes[['1']],jitter(edss),pch=19,col=colors()[idNum2])
#RES = 2
which(pvals.v.multires[[2]]==min(pvals.v.multires[[2]]))
plot(xRes[['2']][,1],jitter(edss),pch=19,col=colors()[idNum2])
#RES = 4
which(pvals.v.multires[[3]]==min(pvals.v.multires[[3]]))
plot(xRes[['4']][,1],jitter(edss),pch=19,col=colors()[idNum2])
#RES = 8
which(pvals.v.multires[[4]]==min(pvals.v.multires[[4]]))
plot(xRes[['8']][,2],jitter(edss),pch=19,col=colors()[idNum2])
#RES = 16
which(pvals.v.multires[[5]]==min(pvals.v.multires[[5]]))
plot(xRes[['16']][,3],jitter(edss),pch=19,col=colors()[idNum2])
#RES = 32
which(pvals.v.multires[[6]]==min(pvals.v.multires[[6]]))
plot(xRes[['32']][,38],jitter(edss),pch=19,col=colors()[idNum2])
#RES = 64
which(pvals.v.multires[[7]]==min(pvals.v.multires[[7]]))
plot(xRes[['64']][,138],jitter(edss),pch=19,col=colors()[idNum2])
#RES = 101
which(pvals.v.multires[[8]]==min(pvals.v.multires[[8]]))
plot(xRes[['101']][,621],jitter(edss),pch=19,col=colors()[idNum2])



#Do a telescoping voxelwise regression
#resolution should be sqrt of total # pixels
#x is the unpixelized matrix
#resolution are the different resolution levels you want to look at
#percentile is the proportion of the image you explore further at each stage
	#can be a vector of length = resolution, or a scalar
#Threshold vector can also be specified. If percentile is also specified, we use an "either or" selection criteria, not a "both" criterea
telescope<-function(y,x,resolution=2,percentile=.5,threshold=0,talk=TRUE,currList=list(),skip='none', ...){
	#just some book-keeping/ordering
	lenRes<-length(resolution)
	resOrd<-order(resolution,decreasing=FALSE)
	resolution<-unique(resolution[resOrd])
	if(length(percentile==1)) percentile<-rep(percentile,lenRes)
	if(length(threshold==1)) threshold<-rep(threshold,lenRes)
	percentile<-percentile[resOrd]
	threshold<-threshold[resOrd]
	p<-dim(x)[2]
	if(skip[1]=='none') skip<-rep(FALSE,p)
	######

	##### 
	#Actual stuff
	if(talk) print(paste('Working on res =',resolution[1]))
	pixelizeMat<-pixelize(x,resolution[1],talk=FALSE,returnMat=TRUE)
	if(resolution[1]>=p) pixelizeMat <- diag(rep(1,p)) #as the above function actually returns a list xRes
	xRes<- x %*% pixelizeMat
	xLabels<-pixelizeMat>0 #each column is now a label for that column of xRes
	#xRes[,j] corresponds to x[,xLabels[,j]]
	skipRes<- (skip %*% xLabels)>0
	vlmRes<-vlm(x=xRes,y=y,skip=skipRes,useRandInt=F,id=NA,talk=FALSE)
	percThres<-quantile(vlmRes$pvals,percentile[1], na.rm=TRUE)
	pick<- vlmRes$pvals<=percThres | vlmRes$pvals<=threshold[1] 
	pickLabs<-xLabels[,(!pick)&(!skipRes)]
	skip<-apply(pickLabs,1,any) | skip
	skip[is.na(skip)]<-TRUE #at this point, some of the NA pvals will bleed into the skip vector. Skip em!
	currList[[as.character(resolution[1])]]<-vlmRes

	#if we're done:
	if(lenRes==1) {if(talk) print('last'); return(currList)}
	#if we've still got layers:
	if(lenRes>1){
		currList<-telescope(y=y,x=x,skip=skip, resolution=resolution[-1], percentile=percentile[-1], threshold=threshold[-1], talk=talk,currList=currList)
		if(talk) print('next')
		return(currList)
	}
}

xtest<-pixelize(fa.map,64)[[1]]
tel<-telescope(y=edss,x=xtest,resolution=c(2,4,8,16,32),percentile=.5,threshold=0,talk=TRUE,currList=list())
imageVec(tel[[1]]$pval, col=c.s.hcl(100:1,h=150))
imageVec(tel[[2]]$pval,col= c.s.hcl(100:1,h=150))
imageVec(tel[[3]]$pval,col=c.s.hcl(100:1,h=150))
imageVec(tel[[4]]$pval,col=c.s.hcl(100:1,h=150))
imageVec(tel[[5]]$pval,col=c.s.hcl(100:1,h=150))

tel[[4]]$pval[!is.na(tel[[4]]$pval)]
(tel[[5]]$pval[!is.na(tel[[5]]$pval)])
min(tel[[5]]$pval,na.rm=TRUE)



##########################################
#    LOOOONG STEP!!!!
##########################################
#get cutoff from permutation test:
nPerms<-800
minPvalPerm<-rep(NA,nPerms)

pb<-txtProgressBar(min = 1, max = nPerms,  char = "=", style = 3)
for(i in 1:nPerms){
	setTxtProgressBar(pb,i)

	y.i<-permuteY(edss[!is.na(edss)])#,id=idObs)
	tel.i<-telescope(y=y.i,x=xtest[!is.na(edss),],resolution=c(2,4,8,16,32),percentile=.5,threshold=0,talk=FALSE,currList=list())
	minPvalPerm[i]<-min(tel.i[[5]]$pval,na.rm=TRUE)
	imageVec(tel.i[[5]]$pval)
}
#save(list='minPvalPerm',file='data/telescope_permutation_null_dist_2013-05-30.RData')
quantile(minPvalPerm,c(.05,.1,.2,.3))
permTelThreshold<-quantile(minPvalPerm,.05)
permTelThreshold
numTests<-0
for(i in 1:5) numTests<-numTests + sum(!is.na(tel[[i]]$pval))
.05/numTests
#Permutation test is less powerfuL? just a little bit?
mean(tel[[5]]$pval<permTelThreshold,na.rm=TRUE)
hist(log(tel[[5]]$pval))
abline(v=log(permTelThreshold),col='blue',lwd=2)
legend('topright','Permutation Cut-off',lwd=2,col='blue')
#dev.copy2pdf(file='plots/telescope_pval_hist_2013-05-30.png',height=1400,width=600)
hist(log(minPvalPerm))
abline(v=log(median(tel[[5]]$pval,na.rm=TRUE)))
imageVec(tel[[5]]$pval,col=c.s.hcl(100:1,h=150))
imageVec(tel[[5]]$pval[tel[[5]]$pval<permTelThreshold],col=c.s.hcl(100:1,h=150))

png(file='plots/telescope_2013-05-30.png',height=1400,width=600)
par(mfrow=c(5,2),mar=c(2,2,3,5))
for(r in 2:6){
	res<-resolution[r]
	xPosition<-1+1.01/resolution[r]

	pvalsMat.v<-matrix(pvals.v.multires[[r]],res,res)
	image(log(pvalsMat.v), asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.s.hcl((100:1),h=150),main='log(p-values)',cex.main=maincex)
	color.legend(xPosition, 0,xPosition+.05,1,legend=round(range(log(pvals.v.multires[[r]]))), rect.col=c.s.hcl(100:1,h=150), align='rb',gradient='y', cex=maincex/2)

	telescopeMat.plot<-matrix(tel[[as.character(res)]]$pval,res,res)
	image(log(telescopeMat.plot), asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.s.hcl(100:1,h=150),main='log(p-values) Telescope',cex.main=maincex)
	color.legend(xPosition, 0,xPosition+.05,1,legend=round(range(log(telescopeMat.plot),na.rm=T)), rect.col=c.s.hcl(100:1,h=150), align='rb',gradient='y', cex=maincex/2)
}
dev.off()


image(pc.i, asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=pc.i.col,main=paste('PC',i,'; Var Explained =',round(varexp[i],digits=3)),cex.main=1.6)
imageVec(x[1,])
imageVec(x2[1,])
imageVec(xRes[2,])






####################################################################
##########     More Complex Stuff               ####################
####################################################################

library(refund)
if(FALSE){
	faPfr <- pfr(edss, funcs=fa.map, kz=35, kb=35, family='poisson', method=REML, smooth.cov=F)
	faPfr <- pfr(edss, funcs=fa.map, kz=35, kb=35, family='poisson', method=REML, smooth.cov=T)
	faPfr <- pfr(edss, funcs=fa.map, kz=35, kb=35, family='gaussian', method=REML, smooth.cov=F)
	faPfr <- pfr(edss, funcs=projFa, kz=35, kb=35, family='gaussian', method=REML, smooth.cov=F)
}


faMapZero<-scale(fa.map, center=T, scale=F)
#faMapZero.t1<-faMapZero[time==1,]

system.time({
faSvd<-fastSvdWide(faMapZero,talk=TRUE)
})# Awesome! 5 seconds & done with full SVD!!!!!!

sum(faSvd$d^2)-sum(faMapZero^2)
varexp<-cumsum(faSvd$d^2)/sum(faMapZero^2)
plot(varexp)

#How much is explained?
#See work check, total sum of the diagonal squared eigenvalues will be:
# = ΣX[,j] by taking Z_j to the δ_j vector

#disp()
png(file='plots/telescope_2013-05-30.png')(file='plots/PC plot 2013-05-29.png',height=750,width=750)
par(mfrow=c(3,3), mar=c(1,1,4,3.6))
#plot some principle components!
for(i in 1:9){
	pc.i <- matrix(faSvd$v[,i], 101,101)
	pc.i.col<-c.d.hcl(  seq(min(pc.i),max(pc.i),length=100)  )
	image(pc.i, asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=pc.i.col,main=paste('PC',i,'; Var Explained =',round(varexp[i],digits=3)),cex.main=1.6)
	color.legend(1.01, 0,1.05,1,legend=round(range(pc.i),3), rect.col=pc.i.col, align='rb',gradient='y', cex=.6)
}
dev.off()
#dev.copy2pdf(file='plots/PCA_test.pdf')
##Smoothed version of data
library(fields)
library(ks)

#####trying to use kde for bases
#v<-faMapZero[3,]
#vmat<-matrix(v,nrow=101,ncol=101)
#image(vmat)
#sv<-kde(vmat, H=1, h=1, gridsize=101, xmin=1, xmax=101, supp=3.7, eval.points=c(10,20,30,40,50))


#################################
#Get Scores for PCs
#Build Basis for Beta

k.b.sqrt<-6 # # of splines; some squared to get both dimensions of the final basis
k.b<-k.b.sqrt^2
k.x<-k.b #num of PC's to use.

getScores<-function(basis,data){
  #inputs: basis is pXk.b; data is nXp
  #outputs: scores are nXk.b
  n<-dim(data)[1]
  q<-dim(basis)[2]
  scores<-matrix(nrow=n,ncol=q) #rows are the score vectors for each person
  print(paste("PROGRESS --  q =",q,'  N =',n))
  pb<-txtProgressBar(min = 1, max = n,  char = "=", style = 3)
  for(i in 1:n){
    scores[i,]<-lm(data[i,]~basis-1)$coef
    setTxtProgressBar(pb,i)
  }
  return(scores)
}

faPcs<-faSvd$v[,1:k.x]
C<-getScores(basis=faPcs,data=fa.map)
projFa<-C%*%t(faPcs) #in wide form
projFaZero<-scale(projFa, center=T, scale=F) #not actually used right now

#how well do the projections look?
par(mfrow=c(3,2),mar=c(5,4,4,4))
for(subj in c(64,200,343)){  
  subj.i.col<-c.d.hcl(  seq(min(fa.map[subj,]),max(fa.map[subj,]),length=100) )
  image(matrix(fa.map[subj,],101,101),col=subj.i.col, main=paste('Raw Data Subj',subj))
  color.legend(1.01, 0,1.05,1,legend=round(range(fa.map[subj,]),2), rect.col=c.s.hcl(1:100), align='rb',gradient='y', cex=.7)
  proj.i.col<-c.d.hcl(  seq(min(projFa[subj,]),max(projFa[subj,]),length=100)  )
  image(matrix(projFa[subj,],101,101),col=proj.i.col,main=paste('Projected FA Subj',subj))
  color.legend(1.01, 0,1.05,1,legend=round(range(projFa[subj,]),2), rect.col=c.s.hcl(1:100), align='rb',gradient='y', cex=.7)
}
par(mfrow=c(1,1))





#recreating "resolution=1" analysis with projected FA
plot(apply(fa.map,1,mean),apply(projFa,1,mean))
plot(apply(fa.map,1,mean),edss)
plot(apply(projFa,1,mean),edss)
plot(fa.map[,621],edss)
plot(projFa[,621],edss)
plot(ests.v,vlmProj$est)
#the estimates are all moved UP a bit. This sort of makes sense? Smoothed?

#Does voxelwise on projected data match up with voxelwise above?
vlmProj<-vlm(x=projFa,y=edss,useRandInt=F,id=NA)
disp()
par(mfrow=c(2,2))

#tag-aosdfoao

estsMat.v<-matrix(ests.v,101,101)
estsRange<-seq(min(ests.v),max(ests.v),length=100)
image(estsMat.v,bty='n',asp=1, xlab='', ylab='', xaxt='n', yaxt='n', col=c.d.hcl(estsRange,h=c(295,40)) ,main='Estimates')
color.legend(1.01, 0,1.05,1,legend=c('low FA','high FA'), rect.col=c.d.hcl(estsRange,h=c(295,40)), align='rb',gradient='y', cex=.5)

image(matrix(vlmProj$est,nrow=101), asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.d.hcl(estsRange,h=c(295,40)) ,main='Estimates')
color.legend(1.01, 0,1.05,1,legend=c('low FA','high FA'), rect.col=c.d.hcl(estsRange,h=c(295,40)), align='rb',gradient='y', cex=.5)

image(log(pvalsMat.v), asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.s.hcl((100:1),h=150),main='log(p-values)')
color.legend(1.01, 0,1.05,1,legend=c(paste('log(p)=',round(min(log(pvals.v)))),'log(p)=0'), rect.col=c.s.hcl(100:1,h=150), align='rb',gradient='y', cex=.5)

image(log(matrix(vlmProj$p,nrow=101)), asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.s.hcl((100:1),h=150),main='log(p-values)')
color.legend(1.01, 0,1.05,1,legend=c(paste('log(p)=',round(min(log(pvals.v)))),'log(p)=0'), rect.col=c.s.hcl(100:1,h=150), align='rb',gradient='y', cex=.5)

plot(quantile(pvalsMat.v))
lines(quantile(vlmProj$p))
#pvalues line up well
par(mfrow=c(1,1))
#################################
#Multidimensional Basis

library(mgcv)
library(splines)

#Get splines
bsplines<-bs(1:sqrt(p),df=k.b.sqrt,degree=2,intercept=T) #degree 3 = cubic, could also do quadratic splines and it would look almost identical.
plot(bsplines[,1],type='l',ylim=range(bsplines), main='bspline basis')
for(i in 1:k.b.sqrt) lines(bsplines[,i],col=c.r.hcl(1:k.b.sqrt)[i],lwd=2)

bsplines2<-matrix(nrow=p,ncol=k.b)

kb2Ticker<-0
for(i in 1:k.b.sqrt){
  for(j in 1:k.b.sqrt){
    kb2Ticker<-kb2Ticker+1
    squareBasis.ij<-bsplines[,i]%*%t(bsplines[,j])
    bsplines2[,kb2Ticker]<-c(squareBasis.ij)
  }
}
image(matrix(bsplines2[,1],101,101), col=c.s.hcl(1:100))
image(matrix(bsplines2[,2],101,101), col=c.s.hcl(1:100))
image(matrix(bsplines2[,8],101,101), col=c.s.hcl(1:100))
image(matrix(bsplines2[,9],101,101), col=c.s.hcl(1:100))
image(matrix(bsplines2[,10],101,101), col=c.s.hcl(1:100))
image(matrix(bsplines2[,16],101,101), col=c.s.hcl(1:100))
image(matrix(bsplines2[,22],101,101), col=c.s.hcl(1:100))
image(matrix(bsplines2[,28],101,101), col=c.s.hcl(1:100))
image(matrix(bsplines2[,34],101,101), col=c.s.hcl(1:100))
image(matrix(bsplines2[,k.b],101,101),col=c.s.hcl(1:100))
image(squareBasis.ij,col=c.s.hcl(1:100))

image(t(bsplines2)%*%bsplines2) #OK, not exactly orthogonal?? But still a basis!

#Get J matrix (integrated thing)
J<-matrix(nrow=k.x,ncol=k.b)
for(kxi in 1:k.x){
  for (kbi in 1:k.b){
    J[kxi,kbi]<-sum(faPcs[,kxi] * bsplines2[,kbi])/p #*1/p to integrate from 0 to 1
  }  
}
test<-t(faPcs)%*%bsplines2
image(t(faPcs)%*%faPcs)
image(faPcs%*%t(faPcs))
image(t(faPcs)%*%bsplines2)
image(J)

disp()

#NEW -- What if we don't bother to approximate X and just use the original thing?
CJEmp<-fa.map%*% bsplines2 #J "empirical" with just 
CJEmpPcBasis<-fa.map%*% faPcs #J also empirical, using pc's as basis functions

par(mfcol=c(1,3))


########## SUPER SUPER Basic Model
#t stat estimations are still totally off:
#WORKS FOR T & PVALUES ONLY VERY ROUGHLY!!!! Uses Z approximation



#############
#############
#############
#Some legwork before comparing models:

plotModelVector<-function(phi,se=NA,beta.basis,title='plot',pval=FALSE,col, ...){
	vec<-beta.basis%*%phi
	if(!is.na(se[1])) vec<-beta.basis%*%(phi/se)
	z<-dnorm(1-abs(vec))
	if(pval)vec<-dnorm(1-z)
	vmat<-matrix(c(vec),nrow=101,ncol=101)
	image(vmat,asp=1,main=title,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=col,...)
	color.legend(1.01, 0,1.05,1,legend=round(range(vmat),digits=4), rect.col=col, align='rb',gradient='y', cex=.5)
}
#A quick function for testing estimates on empirical X(t) matrix for beta coefficeints, from different beta functions.
#plots a beta matrix.
plotEmpModelVector<-function(beta.basis,title='plot', ...){
	CJEmp<-fa.map %*% beta.basis
	model<-lm(edss~CJEmp)
	phi<-model$coefficients[-1]
	vec<-beta.basis%*%phi
	vmat<-matrix(c(vec),nrow=101,ncol=101)
	estsRange<-seq(min(vmat),max(vmat),length=20)
	col<-c.d.hcl(estsRange,h=c(295,40))
	image(vmat,asp=1,main=title,bty='n', xlab='', ylab='', xaxt='n', yaxt='n',col=col,...)
	color.legend(1.01, 0,1.05,1,legend=round(range(vmat),digits=4), rect.col=col, align='rb',gradient='y', cex=.5)
}

## Build up a "square basis" to see how it compares
squareBasisGen<-matrix(round(seq(.5,sqrt(k.b)+.5,length=sqrt(p))),nrow=sqrt(p),ncol=sqrt(k.b))
for(i in 1:sqrt(k.b)){
	thiscol<-squareBasisGen[,i]
	squareBasisGen[thiscol!=i,i]<-0
	squareBasisGen[thiscol==i,i]<-1
}
squareBasis<-matrix(ncol=k.b,nrow=p)
kb2Ticker<-0
for(i in 1:k.b.sqrt){
  for(j in 1:k.b.sqrt){
    kb2Ticker<-kb2Ticker+1
    squareBasis.ij<-squareBasisGen[,i]%*%t(squareBasisGen[,j])
    squareBasis[,kb2Ticker]<-c(squareBasis.ij)
  }
}
image(matrix(squareBasis[,8],101,101), col=c.s.hcl(1:100))
image(matrix(squareBasis[,20],101,101), col=c.s.hcl(1:100))
image(matrix(squareBasis[,k.b],101,101),col=c.s.hcl(1:100) )

#compare to the cubic
par(mfrow=c(1,2))
image(matrix(bsplines2[,1],101,101), col=c.s.hcl(1:100));image(matrix(squareBasis[,1],101,101), col=c.s.hcl(1:100))
image(matrix(bsplines2[,2],101,101), col=c.s.hcl(1:100));image(matrix(squareBasis[,2],101,101), col=c.s.hcl(1:100))
image(matrix(bsplines2[,8],101,101), col=c.s.hcl(1:100));image(matrix(squareBasis[,8],101,101), col=c.s.hcl(1:100))
image(matrix(bsplines2[,9],101,101), col=c.s.hcl(1:100));image(matrix(squareBasis[,9],101,101), col=c.s.hcl(1:100))
image(matrix(bsplines2[,10],101,101), col=c.s.hcl(1:100));image(matrix(squareBasis[,10],101,101), col=c.s.hcl(1:100))
image(matrix(bsplines2[,16],101,101), col=c.s.hcl(1:100));image(matrix(squareBasis[,16],101,101), col=c.s.hcl(1:100))
image(matrix(bsplines2[,22],101,101), col=c.s.hcl(1:100));image(matrix(squareBasis[,22],101,101), col=c.s.hcl(1:100))
image(matrix(bsplines2[,28],101,101), col=c.s.hcl(1:100));image(matrix(squareBasis[,28],101,101), col=c.s.hcl(1:100))
image(matrix(bsplines2[,34],101,101), col=c.s.hcl(1:100));image(matrix(squareBasis[,34],101,101), col=c.s.hcl(1:100))
image(matrix(bsplines2[,k.b],101,101),col=c.s.hcl(1:100));image(matrix(squareBasis[,k.b],101,101), col=c.s.hcl(1:100))

#############
#############
#############
#Now actually make some basic models and compare beta estimates

png(file='plots/compare_beta_estimates_2013-05-30.png',height=250,width=1000)
par(mfrow=c(1,4))
plotEmpModelVector(beta.basis=bsplines2,title='b splines²')
plotEmpModelVector(beta.basis=squareBasis,title='Square basis')
plotEmpModelVector(beta.basis=faPcs,title='PCs')
plotEmpModelVector(beta.basis=faPcs[,1:20],title='PCs')
estsMat.v<-matrix(ests.v,101,101)
estsRange<-seq(min(ests.v),max(ests.v),length=20)
image(estsMat.v, asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.d.hcl(estsRange,h=c(295,40)) ,main='Voxelwise Estimates')
dev.off()


CJ<-C%*%J
nieveLm<-lm(edss~CJ)
nieveEmpLm<-lm(edss~CJEmp)
summary(nieveLm)$coef
summary(nieveEmpLm)$coef
names(summary(nieveEmpLm))

plotModelVector(phi=summary(nieveLm)$coef[-1,1],beta.basis=bsplines2,col=c.d.hcl(estsRange,h=c(295,40)),title='nieve ests')
plotModelVector(phi=summary(nieveEmpLm)$coef[-1,1],beta.basis=bsplines2,col=c.d.hcl(estsRange,h=c(295,40)),title='nieve ests')
plotModelVector(phi=summary(nieveLm)$coef[-1,4],beta.basis=bsplines2,col=c.d.hcl(estsRange,h=c(295,40)),title='nieve ests')
plotModelVector(phi=summary(nieveEmpLm)$coef[-1,4],beta.basis=bsplines2,col=c.d.hcl(estsRange,h=c(295,40)),title='nieve ests',pval=TRUE)
plotModelVector(phi=summary(nieveLm)$coef[-1,1],se=summary(nieveLm)$coef[-1,2],beta.basis=bsplines2,col=c.d.hcl(estsRange,h=c(190,60)),title='nieve tstats')






######## Regress outcome on the scores?
# Y = scores %*% PCs %*% t[phi %*% PCs]
# Y = scores %*% PCs %*% t[PCs] %*% t[phi]
# Y = scores %*% t[phi]
image(cov(C))
scoreLm<-lm(edss~C)
scoreEmpLm<-lm(edss~CJEmpPcBasis) #ends up withs same results!!
cbind(summary(scoreLm)$coef,summary(scoreEmpLm)$coef)
#only the first PC is significant, then nonsig for a while, then sig again at maybe the 19th PC.
phiVecScore<-faPcs%*%summary(nieveLm)$coef[-1,1]


estsRange<-seq(min(ests.v),max(ests.v),length=20)
plotModelVector(phi=summary(scoreLm)$coef[-1,1],beta.basis=faPcs,col=c.d.hcl(estsRange,h=c(295,40)),title='score ests')
plotModelVector(phi=summary(scoreEmpLm)$coef[-1,1],beta.basis=faPcs,col=c.d.hcl(estsRange,h=c(295,40)),title='score ests')

seVec<-sqrt(getBetaSigma(basis=faPcs,Sigma.phi=summary(scoreLm)$cov.un[-1,-1]))
plotModelVector(phi=summary(scoreLm)$coef[-1,3],se=seVec,beta.basis=faPcs,col=c.d.hcl(estsRange,h=c(190,60)),title='score tstats')


estsMat.v<-matrix(ests.v,101,101)
estsRange<-seq(min(ests.v),max(ests.v),length=20)
image(estsMat.v, asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.d.hcl(estsRange,h=c(295,40)) ,main='Voxelwise Estimates')
color.legend(1.01, 0,1.05,1,legend=c(round(range(ests.v),digits=2)), rect.col=c.d.hcl(estsRange,h=c(295,40)), align='rb',gradient='y', cex=.5)

#t.v<-matrix(tvals.v,101,101)
#tRange<-seq(min(t.v),max(t.v),length=20)
#image(t.v, asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.d.hcl(estsRange,h=c(190,60)) ,main='Voxelwise t stats')
#color.legend(1.01, 0,1.05,1,legend=c(round(range(t.v),digits=2)), rect.col=c.d.hcl(estsRange,h=c(190,60)), align='rb',gradient='y', cex=.5)

dev.copy2pdf(file='plots/mismatchingBetas.pdf')
#way off????????????


############# WINBUGS

##########################################################
##########################################################
##########################################################
### Still buggy connection

# Using winbugs do the penalized functional regression
MCMC.result <- bugs(data=list(CJ=C%*%J, y=edss, K_b=K.b, I=p), inits=list(
  list(B0=0, beta=rep(0,15), taubeta=.1)), 
                    parameters.to.save=c('beta', 'B0', 'taubeta'),
                    model.file="bugcode.txt",n.chains=1, n.burnin=10000, n.thin=1,
                    n.iter=20000, bugs.seed=999, bugs.directory="/Users/aaronfisher/Documents/JH/ADAMachine/FA analysis/BUGS"
,working.directory = "/Users/aaronfisher/Documents/JH/ADAMachine/FA analysis/BUGS", debug=F)
save(MCMC.result, file='FA_Surface_MCMC.RData')

##########################################################
##########################################################
##########################################################






########################################
########################################
########################################
############## WORK CHECK ##############
#Show that the sume squared diagonals svd$d^2
#= sum of the squared f terms (f are the scores for each component, a vector with one element for each i)
#f = (U_1'U_1)U_1'X from projecting X onto the 1st PC
# = U_1'X by orthonormality
# ΣU_j'X = ΣZ_j'X where Z_j's are some other basis for X (j=1,..n)
# = ΣX[,j] by taking Z_j to the δ_j vector


## Generate X with a 3 signals, rows are people
nr<-60
nc<-90
x<-matrix(floor(rnorm(nr*nc)),nrow=60,ncol=nc)
sig0<-rep(1,nc)
sig1<-sin(seq(0,2*pi,length=nc))
sig2<-rep(c(1,0),each=nc/2)
mag0<-runif(nr,-4,4)
mag1<-runif(nr,-5,5)
mag2<-runif(nr,-5,5)
x<-x+mag0%*%t(sig0)+mag1%*%t(sig1)+mag2%*%t(sig2)

##analyze it
svdx<-svd(x)
sum(svdx$d^2) - sum(x^2) #get same result! summing the squared entries works.
plot(svdx$v[,1],type='l')
plot(svdx$v[,2],type='l')

plot(svdx$u[1,],type='l')
plot(svdx$u[,1],type='l')

#yea, it's the v's we want

