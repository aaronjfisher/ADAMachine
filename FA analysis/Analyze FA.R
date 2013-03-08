########################################################
##############          FUNCTIONS          #############
########################################################


# Load Aaron's little functions for plotting with colorspace
library(colorspace)

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
c.s.hcl<-function(x,n=1000, ...){ #cut sequantial_hcl
	xCut<-cut(x,breaks=n)
	colors<-sequential_hcl(n=n, ...)[n:1]#reverse the order
	out<-colors[xCut]
	return(out)
}
c.d.hcl<-function(x,n=1000, ...){ #cut divergent_hcl
	xCut<-cut(x,breaks=n)
	colors<-diverge_hcl(n=n, ...)
	out<-colors[xCut]
	return(out)
}

fastSvdWide<-function(A){ #takes a wide centered Nxp matrix with N<<p, and does an svd for it
	#Let A=UDV'; where A=(Nxp); D=(NxN) diagonal; V'=(Nxp); V=(pxN)
	N<-dim(A)[1]
	p<-dim(A)[2]
	if(p<=N) return('error - must be a wide matrix')
	print("Let's get ready to party!")
	#get AA'
	AAt<- matrix(nrow=N,ncol=N) #AA'
	for(i in 1:N){
		for(j in 1:N){
			AAt[i,j]<-sum(A[i,]*A[j,])
		}
	}
	print("Got AA'!!")
	#get svd for AA' = UDV'VDU'=U D^2 U'
	svdAAt<-svd(AAt)
	#A=UDV' => (D^-1)(U^-1)A=V' => (D^-1)(U')A=V'
	#svd for AA'= XYZ'
	#with UD^2 U' = XYZ', so X=D^2, and Z=X=U
	#adjust for sign of U? could be -U? doesn't matter, some of the vectors can be mult by -1 and we still get fine PCs
	D<-sqrt(svdAAt$d)  #note, D is a vector here
	DInv<-1/D
	U<- svdAAt$v
	Vt<-matrix(nrow=N,ncol=p)
	for(i in 1:N){
		for(j in 1:p){
			Vt[i,j]<- D[i]*sum(t(U)[i,]%*%A[,j])
		}
	}
	V<-t(Vt)
	return(list
}


#A-U%*%diag(D)%*%t(V)
#A-svd(A)$u%*%diag(svd(A)$d)%*%svd(A)
#
#A<-matrix(rnorm(4*6),ncol=6,nrow=4)
#t(U)%*%solve(diag(D))%*%A
#



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
matrix(fa.map[1,], 101,101)[1:4,1:4]  #this is how we map back to matrixes, or "carpets"
matrix(orderFaMap, 101,101)[1:4,1:4]


#lets visualize some of these carpets
library(colorspace)
library(plotrix)
set.seed(999)
x <- seq(0, 1, length.out=101)
y <- x
layout(matrix(1:4,2,2))
samp <- sample(c(1:466),4) #pick 4 people from the data for a demo
for (i in samp){
  aint <- matrix(fa.map[i,], 101,101) #removed a t()
  aint <- (aint-min(aint))/diff(range(aint)) #standardize
  image(x,y, aint, asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.s.hcl(1:100),main=paste('Subj',i))
  color.legend(1.01, 0,1.05,1,legend=c('low FA','high FA'), rect.col=c.s.hcl(1:100), align='rb',gradient='y', cex=.5)
}
par(mfrow=c(1,1))

## now load the data of EDSS & identity data
edss <- c(t(read.table('edss_data.txt'))) #length of N = 466
load(file='id.RData') #contains a vector "identity" with the ID of patients (repeated measures), and ALSO the time of visit
#first 5 digits indicate subject id, last 2 digits indicate visits. I forgot to send this to you last time.
head(identity)
id<-sapply(identity,FUN=function(x) {substr(x,1,5)}  )
names(id)<-c()
time<-sapply(identity,FUN=function(x) {substr(x,6,7)}  )
time<-as.numeric(time)

n.id<-length(unique(id)) # #of people
length(edss) #=N

mean(is.na(edss)) # 3.4% missing
sum(duplicated(id)) # num w/ repeated visits
sum(!duplicated(id)) #num w/o repeated visits
hist(c(edss))
summary(edss)
hist(time)






########################################################
################          ANALYSIS      ################
########################################################


library(lme4)
##############################
# Voxel-wise regression w/ random intercept
useRandInt<- TRUE  #lme4 model
#useRandInt<- FAlSE #basic iid errors model


#set up vectors to record results
ests.voxelwise<-rep(NA,times=p)
stdErrs.voxelwise<-rep(NA,times=p)
tvals.voxelwise<-rep(NA,times=p)
pvals.voxelwise<-rep(NA,times=p)

#v=1 Try this to test the loop

print(paste("PROGRESS --  p =",p,'  N =',N))
pb<-txtProgressBar(min = 1, max = p,  char = "=", style = 3)

for(v in 1:p){
	if(useRandInt==T){
		m.v<-lmer(edss~ (1|id) + fa.map[,v]) #model for i.th voxel
		m.v.summary<- attributes(summary(m.v))$coef[2,]
	}
	if(useRandInt==F){
		m.v<-lm(edss~fa.map[,v])
		m.v.summary<- summary(m.v)$coef[2,]
		pvals.voxelwise[v]   <- m.v.summary[4]
	}
	ests.voxelwise[v]    <- m.v.summary[1]
	stdErrs.voxelwise[v] <- m.v.summary[2]
	tvals.voxelwise[v]   <- m.v.summary[3]
	
	setTxtProgressBar(pb,v)
}

hist(unlist(Fvals.voxelwise))



##########################
#Plot results

hist(tvals.voxelwise) # Huh, *lots* of signal?
hist(log(tvals.voxelwise))
abline(v=log(.05),lwd=2,col='blue')
mean(tvals.voxelwise<.05) #43 are pointwise sig? Wow!

pal(c.s.hcl(100:1))




setwd('/Users/aaronfisher/Documents/JH/ADAMachine/FA analysis/plots')
#jpeg(file='voxelwise pvals 2013-2-27.jpg',height=750,width=1200)
#layout(matrix(c(1,3,2,4,6,5),2,3))
par(mfrow=c(1,5))

image(matrix(fa.map[1,],101,101), asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.s.hcl(1:100), main='Sample Subject FA Map (subject 1)')
color.legend(1.01, 0,1.05,1,legend=c('low FA','high FA'), rect.col=c.s.hcl(1:100), align='rb',gradient='y', cex=.5)

estsMat.voxelwise<-matrix(ests.voxelwise,101,101)
image(estsMat.voxelwise, asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.s.hcl(1:100),main='Estimates')
color.legend(1.01, 0,1.05,1,legend=c('low FA','high FA'), rect.col=c.s.hcl(1:100), align='rb',gradient='y', cex=.5)

pvalsMat.voxelwise<-matrix(tvals.voxelwise,101,101)
image(pvalsMat.voxelwise, asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.s.hcl(100:1,h=150),main='p-values')
color.legend(1.01, 0,1.05,1,legend=c('p=0','p=1'), rect.col=c.s.hcl(100:1,h=150), align='rb',gradient='y', cex=.5)


image(pvalsMat.voxelwise>.05, asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.s.hcl(100:1,h=150),main='p-values<.05')
color.legend(1.01, 0,1.05,1,legend=c('p=0','p=1'), rect.col=c.s.hcl(100:1,h=150), align='rb',gradient='y', cex=.5)

image(log(pvalsMat.voxelwise), asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.s.hcl((100:1),h=150),main='log(p-values)')
color.legend(1.01, 0,1.05,1,legend=c(paste('log(p)=',round(min(log(tvals.voxelwise)))),'log(p)=0'), rect.col=c.s.hcl(100:1,h=150), align='rb',gradient='y', cex=.5)

#dev.off()
#dev.copy2pdf(file='voxelwise pvals 2013-2-27.pdf')


####################################################################
##########     More Complex Stuff               ####################
####################################################################

library(refund)
faPfr <- pfr(edss, funcs=fa.map, kz=35, kb=35, family='poisson', method=REML, smooth.cov=F)
faPfr <- pfr(edss, funcs=fa.map, kz=35, kb=35, family='poisson', method=REML, smooth.cov=T)


library(irlba) # For fast SVD

system.time({ #30 seconds
#faMapVar<-cov(fa.map)
})
dim(faMapVar)
testk<-500
faMapVarTest<-faMapVar[1:testk,1:testk]
system.time({ 
svdTest<-svd(fa.map,nu=20,nv=0)
})


faMapZero<-scale(fa.map, center=T, scale=F)
faMapZero.t1<-faMapZero[time==1,]

system.time({
faFastSvd<-irlba(faMapZero, nu=0, nv=220)
})# Awesome! 3 seconds for kb=35; 12sec/100kb; 172 sec for 350kb

system.time({
faFastSvd.t1<-irlba(faMapZero.t1, nu=0, nv=150)
})# Awesome! 3 seconds for kb=35; 12sec/100kb; 172 sec for 350kb


cumsum(faFastSvd.t1$d^2)/sum(faMapZero.t1^2)
varexp<-cumsum(faFastSvd.t1$d^2)/sum(faMapZero.t1^2)


cumsum(faFastSvd$d^2)/sum(faFastSvd$d^2)
cumsum(faFastSvd$d^2)/sum(faMapZero^2)
varexp<-cumsum(faFastSvd$d^2)/sum(faMapZero^2)

plot(faFastSvd$v[,1])

#How much is explained?
#See work check, total sum of the diagonal squared eigenvalues will be:
# = ΣX[,j] by taking Z_j to the δ_j vector


par(mfrow=c(3,3), mar=c(1,1,4,2.5))
#plot some principle components!
x<-seq(0,1,length=101)
y<-seq(0,1,length=101)
for(i in 1:9){
	pc.i <- matrix(faFastSvd$v[,i], 101,101)
	#aint <- (aint-min(aint))/diff(range(aint)) #standardize
	image(x,y, pc.i, asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.s.hcl(1:100),main=paste('PC',i,'; Var Explained =',round(varexp[i],digits=3)))
	color.legend(1.01, 0,1.05,1,legend=c('low FA','high FA'), rect.col=c.s.hcl(1:100), align='rb',gradient='y', cex=.5)
}

pal(c.s.hcl(1:100,h=360))
##Smoothed version of data
library(fields)
library(ks)

v<-faMapZero[3,]
vmat<-matrix(v,nrow=101,ncol=101)
image(vmat)
sv<-kde(vmat, H=1, h=1, gridsize=101, xmin=1, xmax=101, supp=3.7, eval.points=c(10,20,30,40,50))

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


#################################


