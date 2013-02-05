
###########
#START HERE!!
mydir <- "/Users/aaronfisher/Documents/Aaron/Chen's Code/projectmaterials"
setwd(mydir)
load('thickness2.Rdata')
load('simul_2.RData')
source('principalcurve.R')
source('lmstep.R')
source('pjstep.R')
I <- 300
N <- 500
set.seed(999)
y.out <- rpois(I, 3)
if(dim(thickness)[1]<300) thickness <- t(thickness)


#NOW ROWS ARE SUBJ
#COLS ARE COORDS

dim(thickness)
plot(thickness[1,]) #Five different subjects
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

par(mfcol=c(2,2))
for(nbins in c(1,5,10,20,40)){
	binMat<-matrix(nrow=dim(thickness)[1],ncol=nbins)
	p<-dim(thickness)[2]
	thickbreaks<-seq.int(1,p,length=nbins+1)
	for(i in 1:300){
		for(j in 1:nbins){
			binMat[i,j]<- mean(thickness[i,thickbreaks[j]:thickbreaks[j+1]])
		}
	}
	
	#SEE HOW CLOSELY IT'S FITTING
	#par(mfrow=c(2,1))
	#plot(thickness[1,])
	#points(thickbreaks[1:nbins]+p/(2*nbins),binMat[1,],pch=19)
	#plot(thickness[2,])
	#points(thickbreaks[1:nbins]+p/(2*nbins),binMat[2,],pch=19)
	
	mpois<-glm(y.out~binMat,family='poisson') #poisson Model
	if(nbins>1) plot(mpois$coef[-1],type='l',main=paste('nbins = ',nbins),xlab='bin',ylab='beta') #not the intercept

	if(nbins<10) print(mpois$coef)
	
}


########################################################################
########################################################################
########################################################################
#PCA


# plot(colMeans(thickness))

# colMeanMat<-matrix(colMeans(thickness),byrow=T,nrow=dim(thickness)[1],ncol=dim(thickness)[2])

# thickZero<-thickness-colMeanMat

# plot(thickZero[1,]) #Five different subjects
# plot(thickZero[2,])
# plot(thickZero[3,])
# plot(thickZero[4,])
# plot(thickZero[5,])

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

#First 5 PCs
plot(pcx$load[,1])
plot(pcx$load[,2])
plot(pcx$load[,3])
plot(pcx$load[,4])
plot(pcx$load[,5])

#How much did each subj weight on these 5 pcs?
pcaScore5<-pcx$score[,1:5] #using just these 5 components

#Est people based on PCs
estCurve<-array(dim=dim(thickZero))
for(i in 1:300){
	estCurvei<-rep(0,p)
	for(k in 1:5){
		estCurvei<-estCurvei+pcaScore5[i,k]*pcx$load[,k]
	}
	estCurve[i,]<-estCurvei
}

plot(thickZero[1,])
lines(estCurve[1,])
plot(thickZero[2,])
lines(estCurve[2,])
plot(thickZero[3,])
lines(estCurve[3,])
plot(thickZero[4,])
lines(estCurve[4,])
plot(thickZero[5,])
lines(estCurve[5,])


mpcaScore5<-glm(y.out~pcaScore5,family='poisson')

summary(mpcaScore5)

matrix(rep(mpcaScore5$coef,times=101),byrow=T,nrow=101,ncol=6)

