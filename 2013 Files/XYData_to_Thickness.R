# Code name PyschAPATH
# Analysis of simulated data 
# 2013 Feb 08 edit - Chen Yue

# be sure to change date tag in files saved

# set working directory
setwd('C:/Users/Chen/Documents/GitHub/research/ADAMachine/2013 Files')

# load the xy data AF generated
# choose c or snail to load
rm(list=ls())
load('Simulate_c_only_xyOutcome_data_feb4.RData')
load('Simulate_C_all_data_feb4.RData')

# load('Simulate_snail_only_xyOutcome_data_feb4.RData')
# load('Simulate_snail_all_data_feb4.RData')

# a quick look at the data, four random plot
par(mar=c(2,2,2,0.5))
layout(matrix(1:4,2,2))
set.seed(999)
sam.ple <- sample(1:200, 4)
for (i in sam.ple){
  plot(dat.x[i,], dat.y[i,], pch=19, col='blue', cex=.5)  
}

## curving fitting and height function estimation

# source the core code of curve fitting
setwd('C:/Users/Chen/Documents/GitHub/Research/ADAMachine/2013 Files/core_code_for_alg')
source('principalcurve.R')
source('lmstep.R')
source('pjstep.R')

setwd('C:/Users/Chen/Documents/GitHub/Research/ADAMachine/2013 Files/temp_output')
thickness <- NULL

## Let's try for one subject now, then we will do loop
i <- 1

## if the data does not have NA's, we can ignore this line
indices <- !is.na(dat.x[i,])
curr.dat <- cbind(dat.x[i,indices],dat.y[i,indices])
# delete all the NA's
pc <- PrinCurve(curr.dat, sub.size=sum(indices), grid.n=100, 
                max.iter=30, err.tol=0.02,  
                neighbor=0.05, kn=10)
# some notes in here, we use all the none NA data, the grid is also 100, the number of knots is 10.

# plot the fitting result
plot(curr.dat[,1]-mean(curr.dat[,1]), curr.dat[,2]-mean(curr.dat[,2]),
     col='blue', cex=.5, pch=19)
lines(pc$result$prin.curve[,], lty='dashed', pch=19, cex=1, col='black', lwd=2)

# obtain the thickness function
proj <- pc$proj
curr.dist <- pc$height
grid.x <- seq(0,1, length.out=101)
y <- sapply(grid.x, function(a){
  quantile(curr.dist[abs(a-proj)<0.025],prob=.95)
})
y <- smooth.spline(grid.x, y, spar=0.6)$y
plot((101:1)/101,y, type='b', lty='solid', pch=19, cex=1, col='black', lwd=2, ylim=range(thicknesses[i,]*2, y))
lines((1:100)/100,thicknesses[i,]*2)
# plot the thickness funtion
# compare this thickness to the original thickness


# The problem comes from the edge at the beginning, want to deal with it later on
I <- dim(dat.x)[1]
# Now we do loop and create thickness
thickness <- NULL
for (i in 1:I){
  indices <- !is.na(dat.x[i,])
  curr.dat <- cbind(dat.x[i,indices],dat.y[i,indices])
  # delete all the NA's
  pc <- PrinCurve(curr.dat, sub.size=sum(indices), grid.n=100, 
                  max.iter=30, err.tol=0.02,  
                  neighbor=0.05, kn=10)
  proj <- pc$proj
  curr.dist <- pc$height
  grid.x <- seq(0,1, length.out=101)
  y <- sapply(grid.x, function(a){
    quantile(curr.dist[abs(a-proj)<0.025],prob=0.95)
  })
  # 90% quantile is defined for the thickness, PC fitting is biased
  y <- smooth.spline(grid.x, y, spar=.6)$y
  thickness <- cbind(thickness, y)
  cat('Subject ',i,'\ ')
}
save(thickness, file='thickness_est_C.RData')
# save(thickness, file='thickness_est_snail.RData')

# Now do snail
# Repeat the procedure. 
# A special line, which is pretty ~~~~~~~~~~~~ :)
rm(list=ls())
setwd('C:/Users/Chen/Documents/GitHub/research/ADAMachine/2013 Files')

rm(list=ls())
# load the xy data AF generated

load('Simulate_snail_only_xyOutcome_data_feb4.RData')
load('Simulate_snail_all_data_feb4.RData')

# a quick look at the data, four random plot
par(mar=c(2,2,2,0.5))
layout(matrix(1:4,2,2))
set.seed(999)
sam.ple <- sample(1:200, 4)
for (i in sam.ple){
  plot(dat.x[i,], dat.y[i,], pch=19, col='blue', cex=.5)  
}

## curving fitting and height function estimation

# source the core code of curve fitting
setwd('C:/Users/Chen/Documents/GitHub/Research/ADAMachine/2013 Files/core_code_for_alg')
source('principalcurve.R')
source('lmstep.R')
source('pjstep.R')

setwd('C:/Users/Chen/Documents/GitHub/Research/ADAMachine/2013 Files/temp_output')
thickness <- NULL

## Let's try for one subject now, then we will do loop
i <- 1

## if the data does not have NA's, we can ignore this line
indices <- !is.na(dat.x[i,])
curr.dat <- cbind(dat.x[i,indices],dat.y[i,indices])

## we define the initial parametrization of the curve
initial.para <- seq(1,0, length.out=sum(indices))
plot(initial.para)
# delete all the NA's
pc <- PrinCurve(curr.dat, sub.size=sum(indices), grid.n=100, 
                max.iter=30, err.tol=0.01,  
                neighbor=0.03, kn=10, initial.para=initial.para)
# some notes in here, we use all the none NA data, the grid is also 100, the number of knots is 10.

# plot the fitting result
plot(curr.dat[,1]-mean(curr.dat[,1]), curr.dat[,2]-mean(curr.dat[,2]),
     col='blue', cex=.5, pch=19)
lines(pc$result$prin.curve, lty='dashed', pch=19, cex=1, col='black', lwd=2)

# obtain the thickness function
proj <- pc$proj
curr.dist <- pc$height
grid.x <- seq(0,1, length.out=101)
y <- sapply(grid.x, function(a){
  quantile(curr.dist[abs(a-proj)<0.025],prob=0.90)
})
y <- smooth.spline(grid.x, y, spar=.6)$y

# plot the thickness funtion
# compare this thickness to the original thickness
plot(y, type='l', pch=19, cex=1, col='black', lwd=.5, ylim=range(y, thicknesses[i,]*2))
lines(thicknesses[i,100:1]*2, col='blue', lwd=.5, lty='dashed')


# The problem comes from the edge at the beginning, want to deal with it later on
I <- dim(dat.x)[1]
# Now we do loop and create thickness
thickness <- NULL
for (i in 1:I){
  indices <- !is.na(dat.x[i,])
  curr.dat <- cbind(dat.x[i,indices],dat.y[i,indices])
  # delete all the NA's
  # delete all the NA's
  initial.para <- seq(1,0, length.out=sum(indices))
  pc <- PrinCurve(curr.dat, sub.size=sum(indices), grid.n=100, 
                  max.iter=30, err.tol=0.01,  
                  neighbor=0.03, kn=10, initial.para=initial.para)
  proj <- pc$proj
  curr.dist <- pc$height
  grid.x <- seq(0,1, length.out=101)
  y <- sapply(grid.x, function(a){
    quantile(curr.dist[abs(a-proj)<0.025],prob=0.95)
  })
  # 90% quantile is defined for the thickness, PC fitting is biased
  y <- smooth.spline(grid.x, y, spar=.6)$y
  thickness <- cbind(thickness, y)
  #   plot(thickness, type='b', lty='solid', pch=19, cex=1, col='black', lwd=2)
  cat('Subject ',i,'\ ')
}
save(thickness, file='thickness_est_snail.RData')
# save(thickness, file='thickness_est_snail.RData')
