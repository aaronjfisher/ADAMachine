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
layout(matrix(1:4,2,2))
set.seed(999)
x <- seq(0, 1, length.out=101)
y <- x
samp <- sample(c(1:466),4) #pick 4 people from the data for a demo
for (i in samp){
  aint <- matrix(fa.map[i,], 101,101) #removed a t()
  aint <- (aint-min(aint))/diff(range(aint)) #standardize
  image(x,y, aint, asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.s.hcl(1:100))
  color.legend(1.01, 0,1.05,1,legend=c('low FA','high FA'), rect.col=c.s.hcl(1:100), align='rb',gradient='y', cex=.5)
}
par(mfrow=c(1,1))

## now load the data of EDSS
# change folder here
edss <- c(t(read.table('edss_data.txt'))) #length of N = 466
length(edss) #=N
mean(is.na(edss)) # 3.4% missing
hist(c(edss))
summary(edss)






########################################################
################          ANALYSIS      ################
########################################################



##############################
# Voxel-wise regression

#set up vectors to record results
ests.voxelwise<-rep(NA,times=p)
stdErrs.voxelwise<-rep(NA,times=p)
pvals.voxelwise<-rep(NA,times=p)

#v=1 Try this to test the loop

print(paste("PROGRESS --  p =",p,'  N =',N))
pb<-txtProgressBar(min = 1, max = p,  char = "=", style = 3)

for(v in 1:p){
	m.v<-lm(edss~fa.map[,v]) #model for i.th voxel
	ests.voxelwise[v]   <- summary(m.v)$coef[2,1]
	stdErrs.voxelwise[v] <- summary(m.v)$coef[2,2]
	pvals.voxelwise[v] <- summary(m.v)$coef[2,4]
	
	setTxtProgressBar(pb,v)
}


##########################
#Plot results

hist(pvals.voxelwise) # Huh, *lots* of signal?
hist(log(pvals.voxelwise))
abline(v=log(.05),lwd=2,col='blue')
mean(pvals.voxelwise<.05) #43 are pointwise sig? Wow!

pal(c.s.hcl(100:1))




setwd('/Users/aaronfisher/Documents/JH/ADAMachine/FA analysis/plots')
#jpeg(file='voxelwise pvals 2013-2-27.jpg',height=750,width=1200)
layout(matrix(c(1,3,2,4,6,5),2,3))


image(matrix(fa.map[1,],101,101), asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.s.hcl(1:100), main='Sample Subject FA Map (subject 1)')
color.legend(1.01, 0,1.05,1,legend=c('low FA','high FA'), rect.col=c.s.hcl(1:100), align='rb',gradient='y', cex=.5)

estsMat.voxelwise<-matrix(ests.voxelwise,101,101)
image(estsMat.voxelwise, asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.s.hcl(1:100),main='Estimates')
color.legend(1.01, 0,1.05,1,legend=c('low FA','high FA'), rect.col=c.s.hcl(1:100), align='rb',gradient='y', cex=.5)

pvalsMat.voxelwise<-matrix(pvals.voxelwise,101,101)
image(pvalsMat.voxelwise, asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.s.hcl(100:1,h=150),main='p-values')
color.legend(1.01, 0,1.05,1,legend=c('p=0','p=1'), rect.col=c.s.hcl(100:1,h=150), align='rb',gradient='y', cex=.5)


image(pvalsMat.voxelwise>.05, asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.s.hcl(100:1,h=150),main='p-values<.05')
color.legend(1.01, 0,1.05,1,legend=c('p=0','p=1'), rect.col=c.s.hcl(100:1,h=150), align='rb',gradient='y', cex=.5)

image(log(pvalsMat.voxelwise), asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.s.hcl((100:1),h=150),main='log(p-values)')
color.legend(1.01, 0,1.05,1,legend=c(paste('log(p)=',round(min(log(pvals.voxelwise)))),'log(p)=0'), rect.col=c.s.hcl(100:1,h=150), align='rb',gradient='y', cex=.5)

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

system.time({
faFastSvd<-irlba(faMapZero, nu=0, nv=35)
})# Awesome! 3 seconds!
faFastSvd$d^2/sum(faFastSvd$d^2)

plot(faFastSvd$v[,1])

#How much is explained?
#See work check, total sum of the diagonal squared eigenvalues will be:
# = ΣX[,j] by taking Z_j to the δ_j vector


par(mfrow=c(3,3), mar=c(1,1,4,2.5))
#plot some principle components!
x<-seq(0,1,length=101)
y<-seq(0,1,length=101)
for(i in 20:28){
	pc.i <- matrix(faFastSvd$v[,i], 101,101)
	#aint <- (aint-min(aint))/diff(range(aint)) #standardize
	image(x,y, pc.i, asp=1,bty='n', xlab='', ylab='', xaxt='n', yaxt='n', col=c.s.hcl(1:100),main=paste('PC',i))
	color.legend(1.01, 0,1.05,1,legend=c('low FA','high FA'), rect.col=c.s.hcl(1:100), align='rb',gradient='y', cex=.5)
}

sum(faMapZero^2)

cumsum(faFastSvd$d^2)/sum(faMapZero^2) #actual var explained


############## WORK CHECK ##############
#Show that the sume squared diagonals svd$d^2
#= sum of the squared f terms (f are the scores for each component, a vector with one element for each i)
#f = (U_1'U_1)U_1'X from projecting X onto the 1st PC
# = U_1'X by orthonormality
# ΣU_j'X = ΣZ_j'X where Z_j's are some other basis for X (j=1,..n)
# = ΣX[,j] by taking Z_j to the δ_j vector

x<-matrix(rnorm(200),nrow=10,ncol=20)
varx<-var(x)
svdx<-svd(x)
head(svdx$d^2/sum(svdx$d^2))
head(svdx$d^2/sum(varx))


sum((t(svdx$u[,1])%*%x)^2)
svdx$d[1]^2
sum((t(svdx$u[,2])%*%x)^2)
svdx$d[2]^2
sum((t(svdx$u[,3])%*%x)^2)
svdx$d[3]^2

sum(svdx$d^2)
sum(x^2)

totalvar<-0
for(i in 1:d){
	totalvar<-totalvar+sum((t(svdx$u[,i])%*%x)^2)
}
totalvar

Z.3<-rep(0,d)
Z.3[3]<-1
t(Z.3)%*%x

#what about a different basis?
totalvarδ<-0
for(i in 1:d){
	totalvarδ<-totalvarδ+sum(x[,i]^2) #note however that x is a symmetric matrix
}
totalvarδ

sum(svdx$d^2)
#################################







 svdt<-svd(t(faMapZero)%*%faMapZero)
 svdt<-svd(faMapZero)

 svdt$d^2/cumsum(svdt$d^2)

 eigenList$values
 eigenList$values^2/cumsum(eigenList$values^2)

#Just use Pre-sets

pcx<-princomp(x=(faMapZero))

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