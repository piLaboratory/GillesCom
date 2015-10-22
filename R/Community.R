#' Community
#' 
#' Functions for generating and altering \code{Community} objects
#' 
#' The \code{Community} object is one of the fundamental objects in the GillesCom package.
#' @export
#' @useDynLib GillesCom
Community <- function() {
 print(test_basic()) 
  new("Community", list())
}

setClass("Community", representation("list"))
