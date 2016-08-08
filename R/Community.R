#' Rcpp module: Community
#' 
#'    The Community class (or Rcpp_Community) implements a Rcpp interface to simulate the dynamics of an ecological
#'    community. It should be created with \code{\link{Init_Community}}, and the simulation should be carried out
#'    using the \code{\link{bdm}} function, but the data may be retrieved directly using the accessors described
#'    below.
#' @details
#' The class Community contains the following slots:
#' \itemize{
#'  \item{\code{$abundance}}{contains the current abundance vector for the community.}
#'
#'  \item{\code{$time}} {contains the current simulation time for the community.}
#'
#'  \item{\code{$cycles} }{contains the number of elapsed simulation cycles.}
#'
#'  \item{\code{$trajectories}}{ contains a list with the relevant data to reconstruct the system history. The first 
#' element in this list is a data frame in which each line corresponds to the species abundance
#' distribution at a different time, and the second and third elements contains vectors with the time and cycle elapsed
#' at each "snapshot".}
#'
#'  \item{  \code{$K}}{ contains a vector with the carrying capacity of each species.}
#'
#'  \item{  \code{$b}}{ contains a vector with the birth rate of each species.}
#'
#'  \item{  \code{$d0}}{ contains a vector with the base death rate of each species (the actual death rate is affected by the
#'relation between total populations and carrying capacities).}
#'
#'  \item{  \code{$interaction}}{ contains an interaction matrix detailing how each species affects one another. }
#'
#'  \item{  \code{$m}}{ contains a vector with the migration rates of each species. See \code{\link{ls_migration}}.}
#'
#'  \item{  \code{$save_int}}{ contains the saving interval, which is the simulation time that must go by between "snapshots"
#'    that are saved in the slot \code{$trajectories}.}
#'
#' \item{   \code{environmental}}{ contains a \code{data.frame} related to how demographic environmentality or
#' environmental changes affect the community. See the vignettes for details.}
#' }
#' @examples
#' show(Community)
#' @name Community
#' @export
NULL
# Loads the Rcpp Module Community
loadModule("Community", TRUE)

#' Community generating functions
#' 
#' Functions for generating and altering the simulated Community. Function \code{Init_Community}
#' initializes a simulation.
#' 
#' Because of the way the interaction with the underlying c++ code is implemented, saving a Community
#' object using save() or exiting R will not work! See \code{\link{GillesComToFile}} for an alternative.
#' @param abundance vector of initial abundances of species in the community (set species not present to zero). 
#' For convenience, a single number is expanded as \code{rep(1,N)} (this also happens with parameters K, d0, b and m).
#' @param interaction matrix of interaction coefficients, see \code{\link{interaction_matrix}}.
#' @param K vector with carrying capacities of each species
#' @param d0 vector with death rate when N=0
#' @param b vector with birth rates (constant)
#' @param m vector with per capita migration rate in the metacommunity. May be the given as the 
#' resulting list of the \code{\link{ls_migration}} function
#' @param save.int History saving interval (in simulated time units)
#' @param environmental Optional; only use if you want to include demographic stocasticity or environmental
#' fluctuations. A data.frame consisting of 
#' two columns. The first represents the time intervals of the environmental variation, and the second represents the
#' multiplier to carrying capacities at those time intervals. See the vignettes for a more in-depth explanation.
#' @examples
#' # Initializes the community
#' Com = Init_Community(100)
#' # Runs some iteractions of the birth-death-migration process
#' bdm(Com, 5e5, progress = "none")
#' # Gets and analyzes the abundance vector
#' print(ab <- as.numeric(Com$abundance))
#' require(sads)
#' f <- sads::fitlnorm(ab[ab>0])
#' plot(f, which=1)
#' # Simulation internal time elapsed
#' print(Com$time)
#' # Trajectories saves a line for each time period elapsed (starting with 0), along with 
#' # the exact time at each snapshot:
#' print(length(Com$trajectories[[2]]))
#' @return \code{Init_Community} returns an object of class \code{\link{Community}}. 
#' \code{bdm} returns nothing. This function is called by its side effects only.
#' @export
#' @import stats sads Rcpp RcppArmadillo methods graphics
#' @useDynLib GillesCom
Init_Community <- function(abundance, interaction, K = 1000, b = 1, m = 0.1, d0 = 0, save.int = 1, 
                           environmental = data.frame()) {
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
#  create_community(abundance, interaction, K, d0, b, m, save.int, as.matrix(environmental))
  # TODO: Move around some of the parameters to the constructor?
  a = new (Community, abundance, interaction, save.int)
  a$K = K; a$b = b; a$m = m; a$d0 = d0; a$environmental = as.matrix(environmental)
  return(a)
}

#' Interaction matrix
#' 
#' The function \code{interaction_matrix} generates a Lotka-Volterra interaction matrix, following May
#' with interaction coefficients drawn from a Normal or a Gamma distribution. For a 
#' Caswell matrix, use \code{diag(J)}, and for a Hubbell matrix, use \code{ones(J)}.
#'
#' @param J size of metacommunity
#' @param con connectance of the interaction matrix. Defaults to 1 (totally connected)
#' @param stren strength of interaction matrix, which is the standard deviation of the distribution from which the values are drawn
#' @param comp Logical. Use \code{TRUE} for a competition only matrix (all entries are positive); \code{FALSE} for otherwise. Ignored with Gamma distribution.
#' @param distr character, the names of the distribution to draw the interaction coefficients.
#' @param meanInter the mean value of interspecific competition coefficients for Gamma distribution.
#' @note Coefficients drawn from the Gamma distribution are always positive numbers (implying only competion).
#' In this case the mean value of the interspecific competition coefficients can be set with the argument 'meanInter'.
#' The intra-species coefficients are always set to one.
#' @export
#' @import stats
interaction_matrix <- function(J, stren = 0.1, con = 1, comp = TRUE, distr = c("normal", "gamma"), meanInter = 1) {
    distr <- match.arg(distr)
    if(distr=="normal"){
        alphas <- matrix(rnorm(J*J, sd=stren), ncol=J)
        if(comp) alphas <- abs(alphas)
        }
    if(distr=="gamma")
        alphas <- matrix(rgamma(J*J, scale = stren^2, shape = stren^(-2)), ncol=J) * meanInter
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
#' a system of generalized Lotka-Volterra system of competing species. You may provide the number of simulated
#' cycles or the end time for the simulation, but not both.
#' @rdname Init_Community
#' @param community An object of class \code{\link{Community}}
#' @param count Number of cycles to be simulated
#' @param time Total time to be simulated (i.e., if the simulation clock is in t=10, bdm(time=10) will simulate up to t=20)
#' @param progress Should a text bar be used? Currently, "text" will produce a text based bar, and \code{NULL} will produce none.
#' @export
#' @import utils

bdm <- function(community, count, time, progress=c("text", "none")) {
    if (!missing(count) & !missing(time)) stop("Provide count or time, but not both")
    if (missing(count) & missing(time)) count = 1
    progress <- match.arg(progress)
    if(missing(count)) {  ## This section is for specified time
        now = community$time
        step <- time / 100
        if(progress=="text") pb <- utils::txtProgressBar(style=3)
        for (i in 1:100) {
            community$Tbdm(now + i*step)
            if(progress=="text") setTxtProgressBar(pb, i/100)
        }
        if(progress=="text") cat ("\n") 
    } else {
        if (count < 100)
            return(community$Cbdm(count))
        step <- count / 100
        if(progress=="text") pb <- utils::txtProgressBar(style=3)
        for (i in 1:100) {
            community$Cbdm(step)
            if(progress=="text") setTxtProgressBar(pb, i/100)
        }
        if(progress=="text") cat ("\n") 
    }
}


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
#' @import stats sads
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

