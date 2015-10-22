#' Community
#' 
#' Functions for generating and altering \code{Community} objects
#' 
#' The \code{Community} object is one of the fundamental objects in the GillesCom package.

setClass("Community", representation("list"))

Community <- function() {
  new("Community", list())
}
