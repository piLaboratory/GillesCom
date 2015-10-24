#' Community
#' 
#' Functions for generating and altering \code{Community} objects
#' 
#' The \code{Community} object is one of the fundamental objects in the GillesCom package.
#' @export
#' @useDynLib GillesCom
Init_Community <- function(abundance, interaction, K = 100, b = 1, m = 0.1, d0 = 0) {
  # Error checking, etc
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
  create_community(abundance, interaction, K, d0, b, m)
}

#' @param J size of metacommunity
#' @param con connectance of the interaction matrix. Defaults to 1 (totally connected)
#' @param stren strength of interaction matrix, which is the standard deviation of the Gaussin from which the values are drawn
#' @param comp Logical. Use \code{TRUE} for a competition only matrix (all entries are positive); \code{FALSE} for otherwise
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
