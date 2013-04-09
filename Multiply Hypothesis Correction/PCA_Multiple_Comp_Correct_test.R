#COVERAGE MEANS: What % of the time is the true beta fully contained in our joint ci

########################################################
##############        LIBRARIES            #############
########################################################

library(MASS)


########################################################
##############          FUNCTIONS          #############
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
	if((p < N )& (N>1000)){ print('error - must be a wide matrix'); return()}
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


#A<-matrix(rnorm(3*4),ncol=4)


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

δ<-function(i,length){
  x<-rep(0,length)
  x[i]<-1
  return(x)
}

############################################
############################################
# Functions to simulate X & Y

#given n, p, #smooth x inputs, #smooth y inputs, form of β 


#defining rwalk, a random walk generator.
#to be used to add SMOOTH random noise to the covariate function, beyond our defined principle components
rWalk<-function(p,amplitude=1){
  x<-rep(0,p)
  x[1]<-rnorm(1)
  for(i in 2:p){
    x[i]<-rnorm(1,x[i-1],1)
  }
  outpre1<-lowess(x,f=1/8)$y
  outpre2<-(outpre1-mean(outpre1))
  out<-outpre2*amplitude/max(abs(outpre2))
  return(out)
}
plot(rWalk(p=100,amplitude=.2)) #test

genX<-function(n,p,noise.sd=1,k.x.s=10,k.x.r=0){ #k.x.r is # of rough x signals; k.x.s = #smooth signals
  #get signals
  k.x<-k.x.s+k.x.r
  xSmoothSig<-c()
  xRoughSig<-c()
  if(k.x.s>0) xSmoothSig<-apply(matrix(rep(p,k.x.s),ncol=1),1,rWalk) #p X k.x.s
  if(k.x.r>0) xRoughSig<-matrix(rnorm(k.x.r*p),ncol=k.x.r)
  signals<-cbind(xSmoothSig,xRoughSig)
  for(i in 1:k.x) signals[i,]<-signals[i,]/max(abs(signals[i,])) #somewhat redundant but OK, (smooth sig is already scaled)
  signals<-scale(signals,scale=FALSE,center=TRUE)
  
  #weight signals
  intercepts<-runif(n,-1,1) #some variation in the overall height.
  sigWeights<-matrix(runif(k.x*n,-1,1),nrow=n)
  noise<-matrix(rnorm(100,mean=0,sd=noise.sd),ncol=p,nrow=n)
  X<- sigWeights %*% t(signals) + noise
  return(X)
}


n<-500
p<-500
noise.sd<-0
k.x.s<-0
k.x.r<-10
sd.e<-1

x<-genX(n=n,p=p,noise.sd=noise.sd, k.x.s=k.x.s, k.x.r=k.x.r)
x0<-scale(x,center=TRUE,scale=FALSE)
system.time({
  xsvd<-fastSvdWide(x0,talk=TRUE)
})
varExp<-cumsum(xsvd$d^2)/sum(xsvd$d^2)
plot(varExp)
k<-20
pcs<-xsvd$v[,1:k]
plot(pcs[,1])

betat<-c(rep(1,20),rep(0,p-20)) #<-pcs[,1]
#betat<-pcs[,1]
beta0<-runif(1,0,3)
intVec<-matrix(rep(1,n),nrow=n)
e<-rnorm(n,sd=sd.e)
y<-intVec*beta0+x%*%betat+e


#####################################################
#####################################################
#####################################################
#LOAD X & Y (and true beta)
#Y=B0+∫X(t)B(t)+ε 

#load('/Users/aaronfisher/Documents/JH/ADAMachine/Multiply Hypothesis Correction/Test_Data_For_PCA_MultCorrection.RData')

dim(x)
dim(y)
length(y)

#p<-dim(x)[2]
#N<-dim(x)[1]

#PROBLEM? DO I SUBTRACT MEAN???
#IS THIS THE PROBLEM I HAD BEFORE??? WHEN YOU SHOWED CIPRIAN THAT JUNK?

if(k==1) pcs<-matrix(pcs,ncol=1)
pcsTpcs<-t(pcs) %*% pcs
pcspcsT<-pcs %*% t(pcs)

par(mfrow=c(1,2))
dim(pcsTpcs)
image(pcsTpcs)
diag(pcsTpcs)

#dim(pcspcsT)
#image(pcspcsT)
#diag(pcspcsT)
par(mfrow=c(1,1))

############################################################
#" After regression:
# Y = ∫ X(t)B(t) + e
# Y ≅ X'B + e
# Y ≅ ξψ'ψ ϕ+e
# = ξϕ+e
#hat(B)=ψhat(ϕ)
#Σ_β = ψΣ_ϕψ'
#svd(Σ_ϕ) = D'diag(L)D


############################
#Get Σ_ϕ

#Option 1: test with no corr on hat(ϕ)
#Sigma.phi<-diag(rep(1,k))
#m<-NA

#Option 2: test with some small corr
# set.seed(4)
# junk1<-mvrnorm(k,mu=rep(0,k),Sigma=diag(rep(1,k)))
# Sigma.phi<-cov(junk1)
# image(Sigma.phi)
# m<-NA

#Option 3: get Sigma.phi from an actual regression:
xi<-getScores(data=x,basis=pcs)
m<-lm(y~1+xi)
phi<-m$coef[-1]
beta<- c(pcs %*% phi)
Sigma.phi<-summary(m)$cov.unscaled[-1,-1]

plot(beta,type='l')
lines(betat,col='blue')
image(Sigma.phi)
############################

#Write a function to finish it from here

svd.Sigma.phi<-svd(Sigma.phi)
D<-svd.Sigma.phi$v
L<-svd.Sigma.phi$d
Q<-pcs %*% D %*% diag(sqrt(L))

#if p is small enough that we can compute this, let's try it!
#V<-pcs %*% Sigma.phi %*% t(pcs)
#hist(V-Q%*%t(Q)) #nice!

sds<-rep(NA,p)
for(i in 1:p){
  varB.i<-t(pcs[i,]) %*% Sigma.phi %*% pcs[i,]  
  sds[i]<-sqrt(varB.i)
}

Q.scaled<- diag(1/sds) %*% Q
#corV<-cov2cor(V)
#hist(corV-Q.scaled%*%t(Q.scaled)) #nice!

nreps<-2000
#each col of Z will be a draw from a std normal length k dist
Z<-matrix(rnorm(k*nreps,mean=0,sd=1) ,ncol=nreps)
#make a function to get the scaled Bdraw from each Z, and the max_t [scaled B(t)] from that draw. Then apply it to all the reps
Zi2MaxBeta<-function(Zi){
    Bdraw<-abs(Q.scaled %*% Zi)
    return(max(Bdraw))
}
system.time({
  maxBdraws<-apply(Z,2,Zi2MaxBeta)
})
hist(maxBdraws)

qm<-quantile(maxBdraws,.95)
qnorm(1-.05/(2*p))
qnorm(1-.05/(2*k))
qm #almost equal to bonf @ level k, but MUCH CLOSER IF you omit the absolute value requirement in Zi2MaxBeta

plotlim<-max(sds)*max(c(qm,qnorm(1-.05/(2*p))))

pointwise95band<-sds*qnorm(1-.05/2)
joint95band<-sds*qm
#bonf95band<-sds*qnorm(.975/p) #almost equal to joing95band????
#????????
bonf95band.p<-sds*qnorm(1-.05/(2*p))
#our joint is exactly equal to the bonferroni with correction level k????
bonf95band.k<-sds*qnorm(1-.05/(2*k))

colpal<-c.r.hcl(1:10,l=60)
pal(colpal[1:5*2])
par(bg=rgb(.95,.95,.95))

#CI
#pdf(file='/Users/aaronfisher/Documents/JH/ADAMachine/Multiply Hypothesis Correction/TestBands.pdf',height=7,width=7)
plot(pointwise95band,ylim=c(-plotlim*1.5,plotlim),type='l',lwd=2,xlab='Beta Index',ylab='Band Value')
lines(-pointwise95band,ylim=c(-plotlim,plotlim),type='l',lwd=2)
lines(joint95band,lwd=2,col=colpal[8])
lines(-joint95band,lwd=2,col=colpal[8])
lines(bonf95band.p,lwd=2,col=colpal[10])
lines(-bonf95band.p,lwd=2,col=colpal[10])
lines(bonf95band.k,lwd=2,col=colpal[4],lty=2)
lines(-bonf95band.k,lwd=2,col=colpal[4],lty=2)
#legend('bottomright',c('pointwise','q_m band','bonferonni-p','bonferonni-k'),lwd=2,col=c('black',colpal[c(8,10,4)]),lty=c(1,1,1,2))
#dev.off()
lines(betat)

#coverage
plot(betat,type='l',lwd=2)
lines(beta+pointwise95band,ylim=c(-plotlim*1.5,plotlim),type='l',lwd=1,xlab='Beta Index',ylab='Band Value')
lines(beta-pointwise95band,ylim=c(-plotlim,plotlim),type='l',lwd=1)
lines(beta+joint95band,lwd=1,col=colpal[8])
lines(beta-joint95band,lwd=1,col=colpal[8])
lines(beta+bonf95band.p,lwd=1,col=colpal[10])
lines(beta-bonf95band.p,lwd=1,col=colpal[10])
lines(beta+bonf95band.k,lwd=1,col=colpal[4],lty=2)
lines(beta-bonf95band.k,lwd=1,col=colpal[4],lty=2)

if(!is.na(m)[1]){
  lines(betat)
  lines(trueBeta)
}






#WORKSPACE
i<-22;j<-36
(V[i,j] == t(pcs[i,]) %*% Sigma.phi %*% pcs[j,]  )
(V[i,j] - t((pcs%*%D)[i,]) %*% diag(L) %*% (pcs%*%D)[j,] ) 
#so the formula for V[i,j] is 

########################
#generate β from b splines
#Generate true B(t) & functional model parameters
bsBeta<-t(bs(1:I,df=15,degree=3,intercept=TRUE))/2 # "b spline signals"
betat<-matrix(colSums(bsBeta[c(1:3,12:15),]))
dim(betat)<-c(I,1)
plot(betat)
beta0<-3
sd_e<-1
        
#test genX
X<-genX(n=20,p=40,k.x.s=1,k.x.r=1,noise.sd=.01)
dim(X)
matplot(t(X),type='l')
#######################

getBetaSigma<-function(Sigma.phi,basis,diag.only=TRUE){
  p<-min(c(dim(basis)[1],length(basis)))
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


#image(V,col=c.d.hcl(-10:10))
vc<-cov2cor(V)
#image(vc,col=c.d.hcl(-10:10))





