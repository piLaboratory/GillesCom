#' Community
#' 
#' Functions for generating and altering \code{Community} objects
#' 
#' The \code{Community} object is one of the fundamental objects in the GillesCom package.
#' @export
#' @useDynLib GillesCom
Init_Community <- function(abundance, interaction, K, b, m, d0 = 0) {
  # Error checking, etc
  J <- length(abundance)
  if(J > 100000) stop("Maximum number of species reached!")
  if (class(interaction) != "matrix") stop("Interaction must be a matrix!")
  if (length(K)==1) K <- rep(K, J)
  if (length(d0)==1) d0 <- rep(d0, J)
  if (length(b)==1) b <- rep(b, J)
  if (length(m)==1) m <- rep(m, J)
  create_community(abundance, interaction, K, d0, b, m)
}
