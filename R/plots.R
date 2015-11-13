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

#' @import grDevices
# TODO: documentation
# TODO: input parameters
radOverTime <- function(steps) {
  if (missing(steps)) steps <- time() 
  palette <- grDevices::colorRampPalette(c("gray90", "gray10"))(steps)
  palette <- grDevices::adjustcolor(palette, alpha.f=0.5)
  h <- history()
  now <- dim(h)[1]
  J <- dim(h)[2]
  if(now < steps) stop("Not enough simulated data for this number of steps, check history()")
  tinc <- floor(now / steps)
  ab <- as.numeric(abundance())
  plot(0, type='n', log="y", xlab="Species Rank", ylab="Species Abundance", xlim=c(0, J), ylim=c(1, max(h)))
  for (i in 1:steps) {
    lines(rad(h[i * tinc, ]), col=palette[i])
  }
  lines(rad(ab), type='l', col='blue4', lwd=2)
}


octavOverTime <- function(steps, prop=FALSE) {
  if (missing(steps)) steps <- time() 
  par.axis <- list() # TODO
  palette <- grDevices::colorRampPalette(c("gray90", "gray10"))(steps)
  palette <- grDevices::adjustcolor(palette, alpha.f=0.5)
  h <- history()
  now <- dim(h)[1]
  J <- dim(h)[2]
  if(now < steps) stop("Not enough simulated data for this number of steps, check history()")
  plot.new()
  plot(octav(as.numeric(abundance())), col="blue4", prop=prop)
  tinc <- floor(now / steps)
  ab <- as.numeric(abundance())
  for (i in 1:steps) {
    lines(octav(h[i * tinc, ]), col=palette[i], prop=prop)
  }
  lines(octav(as.numeric(abundance())), col="blue4", prop=prop, lwd=1.5)
}

myplot <- function(x, prop=FALSE, x.oct=FALSE, par.axis=list(), ...){
            dots <- list(...)
            x.hist <- rep(as.integer(as.character(x$octave)), as.integer(as.character(x$Freq)))
            h1 <- hist(x=x.hist,
                           breaks = c((min(as.integer(as.character(x$octave)))-1),as.integer(as.character(x$octave))),
                           plot=FALSE)
            if(prop) h1$counts <- h1$counts/sum(h1$counts)
            if(x.oct) xlab <- x[seq(1,length(x[,1]),2),1]
            if(!x.oct) xlab <- x[seq(1,length(x[,1]),2),2]
            if(!"col" %in% names(dots)) dots$col = "gray"
            if(!"main" %in% names(dots)) dots$main = ""
            if(!"ylab" %in% names(dots) & !prop) dots$ylab = "N of species"
            if(!"ylab" %in% names(dots) & prop) dots$ylab = "Proportion of species"
            if(!"xlab" %in% names(dots) & !x.oct) dots$xlab = "Abundance class"
            if(!"xlab" %in% names(dots) & x.oct) dots$xlab = "Abundance class (log2)"
            if(!"axes" %in% names(dots)){
                do.call(plot, c(list(x=h1, axes=FALSE, type='l'),dots, add=TRUE))
                n <- as.numeric(as.character(x[,1]))
                do.call(axis, c(list(side=2), par.axis))
                do.call(axis, c(list(side=1,at=n[seq(1,length(x[,1]),2)],
                     labels=xlab),par.axis))
            }
            else
              do.call(plot, c(list(x=h1,dots)))
          }

