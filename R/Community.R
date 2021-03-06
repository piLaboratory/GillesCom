#' Community
#' 
#' Functions for generating and altering the simulated Community. Function \code{Init_Community}
#' initializes a simulation.
#' 
#' Because of the way the interaction with the underlying c++ code is implemented, only one
#' community may be simulated at a time. Calling \code{Init_Community} more than once will
#' overwrite the previous simulation objects!
#' @param abundance vector of initial abundances of species in the community (set species not present to zero). For convenience, a single number is expanded as \code{rep(1,N)} (this also happens with parameters K, d0, b and m).
#' @param interaction matrix of interaction coefficients, see \code{\link{interaction}}.
#' @param K carrying capacities of each species
#' @param d0 death rate when N=0
#' @param b birth rates (constant)
#' @param m per capita migration rate in the metacommunity. May be the given as the 
#' resulting list of the \code{\link{ls_migration}} function
#' @param save.int History saving interval (in simulated time units)
#' @examples
#' # Initializes the community (in a global object)
#' Init_Community(100)
#' # Runs 1e6 iteractions of the birth-death-migration process
#' bdm(1e6)
#' # Gets and analyzes the abundance vector
#' (ab <- as.numeric(abundance()))
#' require(sads)
#' f <- sads::fitlnorm(ab[ab>0])
#' plot(f, which=1)
#' # Simulation internal time elapsed
#' elapsed_time()
#' # History saves a line for each time period elapsed (starting with 0):
#' dim(trajectories())
#' @export
#' @import graphics
#' @useDynLib GillesCom
#' @rdname Community
Init_Community <- function(abundance, interaction, K = 1000, b = 1, m = 0.1, d0 = 0, save.int = 1) {
  # Error checking, etc
  if (length(abundance)==1) abundance <- rep(0, abundance)
  if (length(abundance) == 0) stop ("Please provide an abundance vector or a positive number of species")
  J <- length(abundance)
  if(missing(interaction)) interaction <- interaction_matrix(J)
  if(save.int <= 0) stop("Save interval must be strictly positive")
  if(J > 100000) stop("Maximum number of species reached!")
  if (class(interaction) != "matrix") stop("interaction must be a matrix!")
  # Helper for using results from ls_migration
  if (class(m) == "list" & "m" %in% names(m)) m = m$m 
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
#' The function \code{interaction_matrix} generates a Lotka-Volterra interaction matrix, following May. For a 
#' Caswell matrix, use \code{diag(J)}, and for a Hubbell matrix, use \code{ones(J)}.
#'
#' @param J size of metacommunity
#' @param con connectance of the interaction matrix. Defaults to 1 (totally connected)
#' @param stren strength of interaction matrix, which is the standard deviation of the Gaussian from which the values are drawn
#' @param comp Logical. Use \code{TRUE} for a competition only matrix (all entries are positive); \code{FALSE} for otherwise
#' @export
#' @import stats
interaction_matrix <- function(J, stren = 0.1, con = 1, comp = TRUE) {
  alphas <- matrix(rnorm(J*J,sd=stren), ncol=J)
  if(comp) alphas <- abs(alphas)
  diag(alphas) <- 1
  if(con<1){
    indexes <- expand.grid(1:J,1:J)
    indexes <- indexes[indexes[,1]!=indexes[,2],]
    ind.i <- sample(c(TRUE,FALSE), nrow(indexes), replace=TRUE, prob=c(1-con, con))
    if(sum(ind.i)>0){
      indexes <- indexes[ind.i, ]
      alphas[indexes[,1], indexes[,2]] <- 0
    }
  }
  alphas
}

#' @export
#' @rdname interaction_matrix
ones <- function(J) matrix(rep(1, J*J), ncol=J)

#' Function \code{bdm} runs one interaction of a Gillespie Algorithm of birth death and migration process in 
#' a system of generalized Lotka-Volterra system of competing species
#' @rdname Community
#' @param count Number of cycles to be simulated
#' @param progress Should a text bar be used? Currently, "text" will produce a text based bar, and \code{NULL} will produce none.
#' @export
#' @import utils

bdm <- function(count=1, progress="text") {
  if (count < 100)
    return(Cbdm(count))
  step <- count / 100
  if(progress=="text") pb <- utils::txtProgressBar(style=3)
  for (i in 1:100) {
    Cbdm(step)
    if(progress=="text") setTxtProgressBar(pb, i/100)
  }
  if(progress=="text") cat ("\n") 
}

#' Function \code{abundance} returns the current abundance vector for the community.
#' @rdname Community
#' @export
"abundance"

#' Function \code{elapsed_time} returns the current simulation time for the community.
#' @rdname Community
#' @export
"elapsed_time"

#' Function \code{trajectories} returns a data frame in which each line corresponds to the species abundance
#' distribution at a different time.
#' @rdname Community
#' @export
"trajectories"

#' Helper functions
#' 
#' Generates migration rates from a log-series metacommunity. The user is expected
#' to provide either S or alpha, but not both.
#'
#' @param J expected size of metacommunity (total number of individuals)
#' @param S Expected number of species in the metacommunity
#' @param alpha Fisher's alpha of the metacommunity
#' @param m per species migration rate
#' @export
#' @import sads
ls_migration <- function(J, S, alpha, m){
    if(missing(alpha)&missing(S))
        stop("Please provide alpha or S")
    ## Fisher's alpha
    if(!missing(S)){
        f1 <- function(a){ S + a * log((a/(a + J))) }
        sol <- uniroot(f1, interval = c(1/J, J))
        alpha <- sol$root
        warning(paste("alpha = ", alpha, " calculated from S and J"))
    }
    ## expected number of species
    else
        S <- ceiling(alpha*log(1+J/alpha))
    ## Sampling S abundances from a logseries
    N <- qls(runif(S), N=J, alpha=alpha) #change to rls when sads 0.3 is released
    ## Migration rates are wheighted by the abundances in the community
    return(list(m=m*N/sum(N), N = N, alpha = alpha))
}

