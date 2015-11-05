#' Community
#' 
#' Functions for generating and altering the simulated Community. Function \code{Init_Community}
#' initializes a simulation.
#' 
#' Because of the way the interaction with the underlying c++ code is implemented, only one
#' community may be simulated at a time. Calling \code{Init_Community} more than once will
#' overwrite the previous simulation objects!
#' @param abundance vector of initial abundances of species in the community (set species not present to zero). For convenience, a single number is expanded as \code{rep(1,N)} (this also happens with parameters K, d0, b and m).
#' @param interaction matrix of interaction coefficients, see \code{\link{Interaction}}.
#' @param K carrying capacities of each species
#' @param d0 death rate when N=0
#' @param b birth rates (constant)
#' @param m per capita migration rate in the metacommunity
#' @param save.int History saving interval (in simulated time units)
#' @examples
#' # Initializes the community (in a global object)
#' Init_Community(100)
#' # Runs 1e6 iteractions of the birth-death-migration process
#' bdm(1e6)
#' # Gets and analyzes the abundance vector
#' (ab <- as.numeric(abundance()))
#' f <- sads::fitlnorm(ab[ab>0])
#' plot(f, which=1)
#' # Simulation internal time elapsed
#' time()
#' # History saves a line for each time period elapsed (starting with 0):
#' dim(history())
#' @export
#' @import graphics
#' @useDynLib GillesCom
#' @rdname Community
Init_Community <- function(abundance, interaction, K = 1000, b = 1, m = 0.1, d0 = 0, save.int = 1) {
  # Error checking, etc
  if (length(abundance)==1) abundance <- rep(1, abundance)
  J <- length(abundance)
  if(missing(interaction)) interaction <- Interaction(J)
  if(J > 100000) stop("Maximum number of species reached!")
  if (class(interaction) != "matrix") stop("Interaction must be a matrix!")
  if (length(K)==1) K <- rep(K, J)
  if (length(d0)==1) d0 <- rep(d0, J)
  if (length(b)==1) b <- rep(b, J)
  if (length(m)==1) m <- rep(m, J)
  if (any(abundance < 0)) stop ("Abundances must be positive integers or zero")
  if (length(K) != J || length(d0) != J || length(b) != J || length(m) != J || dim(interaction) != c(J,J))
     stop("All objects must have the same dimension as the abundance vector")
  create_community(abundance, interaction, K, d0, b, m, save.int)
}

#' Interaction matrix
#' 
#' \code{Interaction} generates a Lotka-Volterra interaction matrix, following May. For a 
#' Caswell matrix, use \code{diag(J)}, and for a Hubbell matrix, use \code{ones(J)}.
#'
#' @param J size of metacommunity
#' @param con connectance of the interaction matrix. Defaults to 1 (totally connected)
#' @param stren strength of interaction matrix, which is the standard deviation of the Gaussin from which the values are drawn
#' @param comp Logical. Use \code{TRUE} for a competition only matrix (all entries are positive); \code{FALSE} for otherwise
#' @export
#' @import stats
Interaction <- function(J, stren = 0.1, con = 1, comp = TRUE) {
  alphas <- matrix(rnorm(J*J,sd=stren), ncol=J)
  if(comp) alphas <- abs(alphas)
  diag(alphas) <- 1
  if(con<1){
    indexes <- expand.grid(1:J,1:J)
    indexes <- indexes[indexes[,1]!=indexes[,2],]
    ind.i <- sample(c(TRUE,FALSE), nrow(indexes), replace=TRUE, prob=c(1-con, con))
    if(sum(ind.i)>0){
      indexes <- indexes[ind.i]
      alphas[indexes[,1], indexes[,2]] <- 0
    }
  }
  alphas
}

#' @export
#' @rdname Interaction
ones <- function(J) matrix(rep(1, J*J), ncol=J)

##' Function \code{bdm} runs one interaction of a Gillespie Algorithm of birth death and migration process in 
##' a system of generalized Lotka-Volterra system of competing species
##' @rdname Community
##' @param count Number of cycles to be simulated
##' @export
#"bdm"
#
##' Function \code{abundance} returns the current abundance vector for the community.
##' @rdname Community
##' @export
#"abundance"
#
##' Function \code{time} returns the current simulation time for the community.
##' @rdname Community
##' @export
#"time"
