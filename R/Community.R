#' Community
#' 
#' Functions for generating and altering \code{Community} objects
#' 
#' The \code{Community} object is one of the fundamental objects in the GillesCom package.
#' @export
#' @useDynLib GillesCom
Init_Community <- function(abundance, interaction) {
  if(length(abundance) > 100000) stop("Maximum number of species reached!")
  # Error checking, etc
  create_community(abundance, interaction)
}
