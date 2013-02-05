## This is the main function for calculating the principal surface of a given 3D dataset

PrinCurve <- function(dat, sub.size=100, grid.n=100, max.iter=20, err.tol=0.05, neighbor=0.03, kn=10){
  set.seed(0)
  samp <- sample(1:(dim(dat))[1], sub.size)
  samp.dat <- dat[samp,]
  center <- apply(samp.dat, 2, mean)
  csamp  <- cbind( (samp.dat[,1]-center[1]),  (samp.dat[,2]-center[2]))
  cdat  <- cbind( (dat[,1]-center[1]),  (dat[,2]-center[2]))
  VarCov <- crossprod(csamp)
  eivec <- eigen(VarCov)$vectors
  score1 <- csamp%*%eivec[,1]
  score1 <- (score1-min(score1))/(max(score1)-min(score1))
  proj <- score1
  j <- 1
  diff <- 1
  library(mgcv)
  while((j<=max.iter)&&(diff>err.tol)){
    l.mean <- lm.step(dat=csamp, proj=proj, neighbor=neighbor)
    result <- pj.step(csamp, l.mean, proj, n.search=grid.n+1, kn=kn)
    proj.new <- result$proj
    diff <- sum((proj.new-proj)^2)
    proj <- proj.new
    j <- j+1
    print(j)
    print(diff)
  }
  grid.fit <- result$prin.curve
  new <- apply(cdat, 1, function(x){which.min(colSums((t(grid.fit)-x)^2))})
  n.search <- grid.n+1
  x.p <- (new-1)/grid.n
  grid.x <- seq(0,1, length.out=grid.n+1)
  projection <- data.frame(x.p)
  names(projection) <- c('X1')
  x.gam <- result$curve.x
  y.gam <- result$curve.y
  curve.dat.x <- predict(x.gam, projection)
  curve.dat.y <- predict(y.gam, projection)
  dat.curve <- cbind(curve.dat.x, curve.dat.y)
  height <- 2*apply(cdat-dat.curve, 1, function(x){(sum(x^2))^.5})
  PrinSurf <- list(result = result, height = height, proj = projection, dat.curve = dat.curve, iter = j-1)
}