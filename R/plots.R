#' Functions to be tested
#' @import sads
#' @export
diagPlots <- function(which=1:4) {
  if (length(which) > 2)
    par(mfrow=c(2,2))
  else if (length(which) > 1)
    par(mfrow=c(1,2))
  if(1 %in% which) {
    S <- function(r) sum(r>0)
    my.S <- apply(history(), 1, S)
    plot(my.S, type='l', main="Species richness", xlab="Time", ylab="")
    cat("Species:\n")
    print(summary(my.S))
  }
  if(2 %in% which) {
    my.N <- apply(history(), 1, sum)
    plot(my.N, type='l', main="Total individuals", xlab="Time", ylab="")
    abline(h=sum(K()), lty=2, lwd=0.8) #Can we estimate the expected number of individuals from K and alpha??
    cat("Individuals:\n")
    print(summary(my.N))
  }
  if(3 %in% which) {
    my.alpha <- apply(history(), 1, get.alpha)
    plot(my.alpha, type='l', main="Fisher's alpha", xlab="Time", ylab="")
    cat("Fisher's alpha:\n")
    print(summary(my.alpha))
  }
  if(4 %in% which) {
    my.sdlog <- apply(history(), 1, get.sdlog)
    plot(my.sdlog, type='l', main="Log-normal sd", xlab="Time", ylab="")
    cat("Log-normal sd:\n")
    print(summary(my.sdlog))
  }
}

get.alpha <- function (r) {
  r <- as.numeric(r[r>0])
  a <- tryCatch({f <- sads::fitls(r); return(bbmle::coef(f)[2]);}, error=function(x) return(NA));
  return(a)
}

get.sdlog <- function (r) {
  r <- as.numeric(r[r>0])
  a <- tryCatch({f <- sads::fitlnorm(r); return(bbmle::coef(f)[2]);}, error=function(x) return(NA));
  return(a)
}
