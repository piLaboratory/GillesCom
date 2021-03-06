#' Plotting functions
#'
#' These functions plot relevant information about a simulation
#' 
#' The function \code{diagPlots} plots up to four different diagnostic plots (controlled by the
#' parameter "which":
#' 1 - Species richness over time
#' 2 - Total individuals over time
#' 3 - Number of species equivalents calculated from Shannon information index
#' 4 - Kolmogorov-Smirnoff statistic between SADs at each time compared to the SAD at the last time
#' 5 - Estimated Fisher's alpha over time (from a log-series fit)
#' 6 - Estimated sdlog over time (from a lognormal fit)
#'
#' Notice that the default invocation of this function displays only the plots 1 to 4.
#'
#' The functions \code{radOverTime} and \code{octavOverTime} provide a superimposing plot
#' with the rad and octav, respectively, at distinct points in time.
#' @import sads
#' @rdname plots
#' @export
#' @param which numerical, may be a single number or a vector: which plots should be done (see Details)
diagPlots <- function(which=1:4) {
    opar <- par(no.readonly=TRUE)
    on.exit(par(opar))
  if (length(which) > 4)
    par(mfrow=c(2,3))
  else if (length(which) > 2)
    par(mfrow=c(2,2))
  else if (length(which) > 1)
    par(mfrow=c(1,2))
  if(1 %in% which) {
    S <- function(r) sum(r>0)
    my.S <- apply(trajectories(), 1, S)
    plot(my.S, type='l', main="Species richness", xlab="Time", ylab="")
    cat("Species:\n")
    print(summary(my.S))
  }
  if(2 %in% which) {
    my.N <- apply(trajectories(), 1, sum)
    plot(my.N, type='l', main="Total individuals", xlab="Time", ylab="")
    abline(h=sum(K()), lty=2, lwd=0.8) #Can we estimate the expected number of individuals from K and alpha??
    cat("Individuals:\n")
    print(summary(my.N))
  }
  if(3 %in% which) {
    my.H <- apply(trajectories(), 1, get.Heq)
    plot(my.H, type='l', main="Shannon's species equivalent", xlab="Time", ylab="")
    cat("Shannon's species equivalent:\n")
    print(summary(my.H))
  }
  if(4 %in% which) {
    my.D <- get.ks(trajectories())
    plot(ks ~ tempo, data=my.D , type='l', main="KS distance to final SAD", xlab="Time", ylab="")
    cat("Ks distance to final SAD:\n")
    print(summary(my.D$ks))
  }
  if(5 %in% which) {
    my.alpha <- apply(trajectories(), 1, get.alpha)
    plot(my.alpha, type='l', main="Fisher's alpha", xlab="Time", ylab="")
    cat("Fisher's alpha:\n")
    print(summary(my.alpha))
  }
  if(6 %in% which) {
    my.sdlog <- apply(trajectories(), 1, get.sdlog)
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
## Shannon's species equivalents
get.Heq <- function(r){
    r <- as.numeric(r[r>0])
    rp <- r/sum(r)
    a <- sum(-rp*log(rp))
    return(exp(a))
}
## Kolmogorov-Smirnoff distance from the sad at the highest time
get.ks <- function(m, lag=1){
    tempo <- seq(1,nrow(m), by=lag)
    kst <- c()
    ksp <- c()
    for(i in 1:(length(tempo)-1)){
        teste <- ks.test(m[tempo[i],], m[nrow(m),])
        kst[i] <- teste$statistic
        ksp[i] <- teste$p.value
    }
    data.frame(tempo=tempo[1:length(tempo)-1], ks=kst, ksp=ksp)
}

#' @import grDevices
#' @rdname plots
#' @param steps number of intervals in which to cut the trajectories
#' @param col For both \code{radOverTime} and \code{octavOverTime}, the colors are chosen by three values: [1] the color of the most ancestral rad/octav, [2] the color of the latest rad/octav and [3] a highlight color for the current rad/octav. The colors for the remaining rad/octavs to be plotted are generated with a ramp palette from col[1] to col[2].
#' @param par.axis Additional graphical parameters for handling the axes of the plot.
#' @param \dots Additional graphical parameters for the plots, such as main title, labels, font size, etc; EXCEPT for the axis.
radOverTime <- function(steps, col=c("gray90", "gray10", "blue4"), par.axis=list(), ...) {
    dots <- list(...)
    if (length(col) != 3) stop ("The col argument must have exactly three elements; see help")
    if (missing(steps)) steps <- elapsed_time() 
    palette <- grDevices::colorRampPalette(c(col[1], col[2]))(steps)
  palette <- grDevices::adjustcolor(palette, alpha.f=0.5)
  h <- trajectories()
  now <- dim(h)[1]
  J <- dim(h)[2]
  Jmax <- max(apply(h, 1, function(x) sum(x>0)))
  if(now < steps) stop("Not enough simulated data for this number of steps, check trajectories()")
  tinc <- floor(now / steps)
  ab <- as.numeric(abundance())
  if(!"main" %in% names(dots)) dots$main = "Simulated rank abundances over time"
  if(!"ylab" %in% names(dots)) dots$ylab = "Species Abundance"
  if(!"xlab" %in% names(dots)) dots$xlab = "Species Rank"
  do.call(plot, c(list(x=1, type='n', axes=FALSE, log="y", xlim=c(0,Jmax), ylim=c(1, max(h))), dots))
  do.call(axis, c(list(1), par.axis))
  do.call(axis, c(list(2), par.axis))
  for (i in 1:steps) {
    if(sum(h[i*tinc,]) >0) lines(rad(h[i * tinc, ]), col=palette[i])
  }
  lines(rad(ab), type='l', col=col[3], lwd=2)
}

#' @rdname plots
#' @param prop Logical. Should the octav be plotted using proportions (as opposed to absolute numbers)?
octavOverTime <- function(steps, prop=TRUE, col=c("gray90", "gray10", "blue4"), par.axis=list(), ...) {
  dots <- list(...) 
  if (length(col) != 3) stop ("The col argument must have exactly three elements; see help")
  if (missing(steps)) steps <- elapsed_time() 
  palette <- grDevices::colorRampPalette(c(col[1], col[2]))(steps)
  palette <- grDevices::adjustcolor(palette, alpha.f=0.5)
  h <- trajectories()
  now <- dim(h)[1]
  maxO <- ceiling(max(log2(trajectories()))+1)
  tinc <- floor(now / steps)
  # Finds the maximum scale for y
  maxY = 0
  for (i in 1:steps) {
    if(sum(h[i*tinc,])>0) {
        o <- octav(h[i * tinc, ])
        if (prop) newy = max(o$Freq)/sum(o$Freq)
        else newy = max(o$Freq)
        if (newy > maxY) maxY = newy
    }
  }
  o <- octav(as.numeric(abundance()))
  if (prop) newy = max(o$Freq)/sum(o$Freq)
  else newy = max(o$Freq)
  if (newy > maxY) maxY = newy

  J <- dim(h)[2]
  if(now < steps) stop("Not enough simulated data for this number of steps, check trajectories()")
  if(!"main" %in% names(dots)) dots$main = "Simulated octaves over time"
  if(!"ylab" %in% names(dots) & prop) dots$ylab = "Proportion of species"
  if(!"ylab" %in% names(dots) & !prop) dots$ylab = "Number of species"
  if(!"xlab" %in% names(dots)) dots$xlab = "Abundance class"
  do.call(plot, c(list(x=0, axes=FALSE, type='n', xlim=c(-0.5, maxO), ylim=c(0,maxY)),dots))
  ab <- as.numeric(abundance())
  for (i in 1:steps) {
    if(sum(h[i*tinc,])>0) lines(octav(h[i * tinc, ]), col=palette[i], prop=prop, type='l')
  }
  lines(octav(as.numeric(abundance())), col=col[3], prop=prop, lwd=1.5)
  x <- octav(as.numeric(abundance()))
  xlab <- x[seq(1,length(x[,1]),2),2]
    n <- as.numeric(as.character(x[,1]))
    do.call(axis, c(list(side=2), par.axis))
    do.call(axis, c(list(side=1,at=n[seq(1,length(x[,1]),2)],
                         labels=xlab),par.axis))
}

