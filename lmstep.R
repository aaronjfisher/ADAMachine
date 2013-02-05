lm.step <- function(dat, proj, neighbor=0.03){
  n <- nrow(dat)
  p <- ncol(dat)
  lmean <- sapply(proj, function(x){
    c(mean(dat[abs(proj-x)<=neighbor,1]), mean(dat[abs(proj-x)<=neighbor,2]))
  }
  )
  lm.step <- t(lmean)
}