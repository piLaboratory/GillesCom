#' Plotting functions
#'
#' These functions plot relevant information about a simulation
#' 
#' The function \code{diagPlots} plots up to four different diagnostic plots (controlled by the
#' parameter "which":
#' 1 - Total individuals over time
#' 2 - Species- equivalents over time: number of species, number of species equivalents calculated by Shannon's and Simpson's diversity indexes.
#' 3 - Mean abundance (mean number of individual per species)
#' 4 - Abundances standard deviation (per species)
#'
#' Notice that the default invocation of this function displays plots 1 to 4.
#'
#' The functions \code{radOverTime} and \code{octavOverTime} provide a superimposing plot
#' with the rad and octav, respectively, at distinct points in time.
#' @param community An object of class \code{\link{Community}}
#' @param which numerical, may be a single number or a vector: which plots should be done (see Details)
#' @examples
#' Com = Init_Community(50)
#' bdm(Com, 1e4, progress="none")
#' diagPlots(Com)
#' radOverTime(Com)
#' octavOverTime(Com)
#' @rdname plots
#' @export
diagPlots <- function(community, which=1:4) {
    opar <- par(no.readonly=TRUE)
    on.exit(par(opar))
    if (length(which) > 2)
        par(mfrow=c(2,2))
    else if (length(which) > 1)
        par(mfrow=c(1,2))
    if(1 %in% which) {
        my.N <- apply(community$trajectories[[1]], 1, sum)
        plot(my.N ~ community$trajectories[[2]], type='l', main="Total individuals", xlab="Time", ylab="")
        cat("Individuals:\n")
        print(summary(my.N))
    }
    if(2 %in% which) {
        S <- function(r) sum(r>0)
        my.S <- apply(community$trajectories[[1]], 1, S)
        plot(my.S ~ community$trajectories[[2]], type='l', main="Species-equivalent \n Richness, sHannon, Simpson", xlab="Time", ylab="", col=1)
        my.H <- apply(community$trajectories[[1]], 1, get.Heq)
        lines(my.H ~ community$trajectories[[2]], type='l', col=2)
        my.D <- apply(community$trajectories[[1]], 1, get.Deq)
        lines(my.D ~ community$trajectories[[2]], type='l', col=3)
        legend("topleft", c("R", "H", "S"), lty=1, col=1:3, bty="n")
        cat("Species:\n")
        print(summary(my.S))
        cat("Shannon's Species-equivalent:\n")
        print(summary(my.H))
        cat("Simpson's Species-equivalent:\n")
        print(summary(my.D))
    }
    if(3 %in% which) {
        my.mean <- apply(community$trajectories[[1]], 1, mean)
        plot(my.mean ~ community$trajectories[[2]], type='l', main="Abundance mean", xlab="Time", ylab="")
        cat("Mean abundance / species:\n")
        print(summary(my.mean))
    }
    if(4 %in% which) {
        my.sd <- apply(community$trajectories[[1]], 1, sd)
        plot(my.sd ~ community$trajectories[[2]], type='l', main="Abundance standard dev", xlab="Time", ylab="")
        cat("Mean sd / species:\n")
        print(summary(my.sd))
    }
}

## Shannon's species equivalents
get.Heq <- function(r){
    r <- as.numeric(r[r>0])
    rp <- r/sum(r)
    a <- sum(-rp*log(rp))
    return(exp(a))
}
## Simpson's species equivalents
get.Deq <- function(r){
    r <- as.numeric(r[r>0])
    rp <- r/sum(r)
    a <- sum(rp^2)
    return(1/a)
}

#' @import grDevices
#' @param steps number of intervals in which to cut the trajectories
#' @param col For both \code{radOverTime} and \code{octavOverTime}, the colors are chosen by three values: [1] the color of the most ancestral rad/octav, [2] the color of the latest rad/octav and [3] a highlight color for the current rad/octav. The colors for the remaining rad/octavs to be plotted are generated with a ramp palette from col[1] to col[2].
#' @param par.axis Additional graphical parameters for handling the axes of the plot.
#' @param \dots Additional graphical parameters for the plots, such as main title, labels, font size, etc; EXCEPT for the axis.
#' @rdname plots
#' @export
radOverTime <- function(community, steps, col=c("gray90", "gray10", "blue4"), par.axis=list(), ...) {
    dots <- list(...)
    if (length(col) != 3) stop ("The col argument must have exactly three elements; see help")
    if (missing(steps)) steps <- community$time
    palette <- grDevices::colorRampPalette(c(col[1], col[2]))(steps)
    palette <- grDevices::adjustcolor(palette, alpha.f=0.5)
    h <- community$trajectories[[1]]
    now <- dim(h)[1]
    J <- dim(h)[2]
    Jmax <- max(apply(h, 1, function(x) sum(x>0)))
    if(now < steps) stop("Not enough simulated data for this number of steps, check trajectories()")
    tinc <- floor(now / steps)
    ab <- as.numeric(community$abundance)
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

#' @param prop Logical. Should the octav be plotted using proportions (as opposed to absolute numbers)?
#' @rdname plots
#' @export
octavOverTime <- function(community, steps, prop=TRUE, col=c("gray90", "gray10", "blue4"), par.axis=list(), ...) {
    dots <- list(...) 
    if (length(col) != 3) stop ("The col argument must have exactly three elements; see help")
    if (missing(steps)) steps <- community$time
    palette <- grDevices::colorRampPalette(c(col[1], col[2]))(steps)
    palette <- grDevices::adjustcolor(palette, alpha.f=0.5)
    h <- community$trajectories[[1]]
    now <- dim(h)[1]
    maxO <- ceiling(max(log2(h))+1)
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
    o <- octav(as.numeric(community$abundance))
    if (prop) newy = max(o$Freq)/sum(o$Freq)
    else newy = max(o$Freq)
    if (newy > maxY) maxY = newy

    J <- dim(h)[2]
    if(now < steps) stop("Not enough simulated data for this number of steps, check your community trajectories")
    if(!"main" %in% names(dots)) dots$main = "Simulated octaves over time"
    if(!"ylab" %in% names(dots) & prop) dots$ylab = "Proportion of species"
    if(!"ylab" %in% names(dots) & !prop) dots$ylab = "Number of species"
    if(!"xlab" %in% names(dots)) dots$xlab = "Abundance class"
    do.call(plot, c(list(x=0, axes=FALSE, type='n', xlim=c(-0.5, maxO), ylim=c(0,maxY)),dots))
    ab <- as.numeric(community$abundance)
    for (i in 1:steps) {
        if(sum(h[i*tinc,])>0) lines(octav(h[i * tinc, ]), col=palette[i], prop=prop, type='l')
    }
    lines(octav(ab), col=col[3], prop=prop, lwd=1.5)
    x <- octav(ab)
    xlab <- x[seq(1,length(x[,1]),2),2]
    n <- as.numeric(as.character(x[,1]))
    do.call(axis, c(list(side=2), par.axis))
    do.call(axis, c(list(side=1,at=n[seq(1,length(x[,1]),2)],
                         labels=xlab),par.axis))
}

