# Calculate the I2C2 of FA maps
# Chen Yue 03/25/2013

# Source the I2C2 function
setwd('C:/Users/Chen/Documents/Research/Thickness analysis/FAmap/I2C2_software')
source('I2C2_inference_comment.R')

# load data
setwd('C:/Users/Chen/Documents/Research/Modelling the Corpus Callosum/corpus callosum paper/code for plot')
load('PSCC.RData')
dim(fa.map)
fa <- t(fa.map)

# get the file names
path <- "C:/Users/Chen/Documents/Research/Modeling HDLFPCA/Masked"
files.masked <- dir(path, glob2rx("*.nii")) # read in all scans

# get the I, J, id and visit from the file names

raw.id<-sapply(files.masked, FUN=function(x) {substr(x,1,5)}  )
names(raw.id)<-c()
visit<-sapply(files.masked, FUN=function(x) {substr(x,6,7)}  )
visit<-as.numeric(visit)

I <- length(unique(raw.id))
J <- as.vector(table(raw.id, dnn=NULL))
id <- NULL
for (i in 1:I){
  id <- c(id, rep(i, times=J[i]))
}
system.time( i2c2.est <- I2C2(fa, I = I, J = J, id = id, visit = visit, demean = TRUE) )
i2c2.est$lambda

# exclude the data with only 1 visit
# each variable name is added by ".2" 
# indicating getting rid of subjest less than 2 visit. (does not contribute to ICC)
valid <- which(as.vector(table(raw.id))!=1)
I.2 <- length(valid)
J.2 <- (as.vector(table(raw.id)))[valid]
visit.2 <- NULL
id.2 <- NULL
for (i in 1:I.2){
  visit.2 <- c(visit.2, 1:J.2[i])
  id.2 <- c(id.2, rep(i, J.2[i]))
}
valid.id <- names(which(table(raw.id)!=1))
valid.scan <- raw.id%in%valid.id
fa.2 <- fa[valid.scan,]

# calculate the I2C2
system.time( i2c2.est.2 <- I2C2(fa.2, I = I.2, J = J.2, id = id.2, visit = visit.2, demean = TRUE) )
i2c2.est.2$lambda

# calculate the confidence interval, need some time (10 minutes in my laptop)
system.time( i2c2.int <- I2C2.mcCI(fa.2, I = NULL, J = NULL, id = id.2, visit = visit.2, R = 100, rseed = 1, ncores = 1, demean = TRUE, ci = 0.95 ) )
hist(unlist(i2c2.int))

# calculate the null distribution
system.time( i2c2.null <- I2C2.mcNulldist( fa.2, I = I.2, J = J.2, id = id.2, visit = visit.2, R = 100, rseed = 1, ncores = 1, demean = FALSE ) )  
hist(unlist(i2c2.null))

# do the plot for I2C2 analysis
library(beanplot)
beanplot( data.frame(VN = i2c2.null), border = 8, ylim = c(-0.2,.8), cex.main = 1, cex.axis = .8, main = "I2C2-MAmaps", ll = 0.001, col = c(8, 8, "#B2DF8A") )
lines( rep(1, 2), i2c2.int$CI, col = 1, lwd = 3, lty = 1)
lines( c(0.9, 1.1), rep(i2c2.int$CI[1], 2), col = 1, lty = 1, lwd = 2)
lines( c(0.9, 1.1), rep(i2c2.int$CI[2], 2), col = 1, lty = 1, lwd = 2)
lines(c(0.95, 1.05), rep(i2c2.est.2$lambda,2), col = 2, lwd = 3, lty = 1)
legend('topright',c("I2C2","95% CI","Null"),lwd=2,cex=.6,col=c(2,1,8))

# save the result
setwd('C:/Users/Chen/Documents/Research/Thickness analysis/FAmap')
save(i2c2.est.2, i2c2.int, i2c2.null, file='i2c2_result.RData')