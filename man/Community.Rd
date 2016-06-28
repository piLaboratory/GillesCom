\name{Rcpp_Community}
\alias{Rcpp_Community}
\alias{Community}
\docType{data}
\title{
    Community class
}
\description{
    The Community class (or Rcpp_Community) implements a Rcpp interface to simulate the dynamics of an ecological
    community. It should be created with \code{\link{Init_Community}}, and the simulation should be carried out
    using the \code{\link{bdm}} function, but the data may be retrieved directly using the accessors described
    below.

}
\usage{Com = new(Community)}
\format{
#' Function \code{abundance} returns the current abundance vector for the community.
#' @rdname Community
#' @export
"abundance"

#' Function \code{elapsed_time} returns the current simulation time for the community.
#' @rdname Community
#' @export
"elapsed_time"

#' Function \code{elapsed_cycles} returns the number of elapsed simulation cycles.
#' @rdname Community
#' @export
"elapsed_cycles"

#' Function \code{trajectories} returns a data frame in which each line corresponds to the species abundance
#' distribution at a different time.
#' @rdname Community
#' @export
"trajectories"
  The format is:
Formal class 'C++Class' [package "Rcpp"] with 11 slots
  ..@ .Data       : chr "Rcpp_Community"
  ..@ pointer     :<externalptr> 
  ..@ module      :<externalptr> 
  ..@ fields      :List of 11
  .. ..$ K           :Reference class 'C++Field' [package "Rcpp"] with 5 fields
  .. .. ..$ pointer      :<externalptr> 
  .. .. ..$ cpp_class    : chr "arma::Col<double>"
  .. .. ..$ read_only    : logi FALSE
  .. .. ..$ class_pointer:<externalptr> 
  .. .. ..$ docstring    : chr ""
  .. .. ..and 14 methods.
  .. ..$ abundance   :Reference class 'C++Field' [package "Rcpp"] with 5 fields
  .. .. ..$ pointer      :<externalptr> 
  .. .. ..$ cpp_class    : chr "arma::Col<double>"
  .. .. ..$ read_only    : logi FALSE
  .. .. ..$ class_pointer:<externalptr> 
  .. .. ..$ docstring    : chr ""
  .. .. ..and 14 methods.
  .. ..$ b           :Reference class 'C++Field' [package "Rcpp"] with 5 fields
  .. .. ..$ pointer      :<externalptr> 
  .. .. ..$ cpp_class    : chr "arma::Col<double>"
  .. .. ..$ read_only    : logi FALSE
  .. .. ..$ class_pointer:<externalptr> 
  .. .. ..$ docstring    : chr ""
  .. .. ..and 14 methods.
  .. ..$ cycles      :Reference class 'C++Field' [package "Rcpp"] with 5 fields
  .. .. ..$ pointer      :<externalptr> 
  .. .. ..$ cpp_class    : chr "int"
  .. .. ..$ read_only    : logi FALSE
  .. .. ..$ class_pointer:<externalptr> 
  .. .. ..$ docstring    : chr ""
  .. .. ..and 14 methods.
  .. ..$ d0          :Reference class 'C++Field' [package "Rcpp"] with 5 fields
  .. .. ..$ pointer      :<externalptr> 
  .. .. ..$ cpp_class    : chr "arma::Col<double>"
  .. .. ..$ read_only    : logi FALSE
  .. .. ..$ class_pointer:<externalptr> 
  .. .. ..$ docstring    : chr ""
  .. .. ..and 14 methods.
  .. ..$ interaction :Reference class 'C++Field' [package "Rcpp"] with 5 fields
  .. .. ..$ pointer      :<externalptr> 
  .. .. ..$ cpp_class    : chr "arma::Mat<double>"
  .. .. ..$ read_only    : logi FALSE
  .. .. ..$ class_pointer:<externalptr> 
  .. .. ..$ docstring    : chr ""
  .. .. ..and 14 methods.
  .. ..$ m           :Reference class 'C++Field' [package "Rcpp"] with 5 fields
  .. .. ..$ pointer      :<externalptr> 
  .. .. ..$ cpp_class    : chr "arma::Col<double>"
  .. .. ..$ read_only    : logi FALSE
  .. .. ..$ class_pointer:<externalptr> 
  .. .. ..$ docstring    : chr ""
  .. .. ..and 14 methods.
  .. ..$ save_int    :Reference class 'C++Field' [package "Rcpp"] with 5 fields
  .. .. ..$ pointer      :<externalptr> 
  .. .. ..$ cpp_class    : chr "double"
  .. .. ..$ read_only    : logi FALSE
  .. .. ..$ class_pointer:<externalptr> 
  .. .. ..$ docstring    : chr ""
  .. .. ..and 14 methods.
  .. ..$ stochastic  :Reference class 'C++Field' [package "Rcpp"] with 5 fields
  .. .. ..$ pointer      :<externalptr> 
  .. .. ..$ cpp_class    : chr "arma::Mat<double>"
  .. .. ..$ read_only    : logi FALSE
  .. .. ..$ class_pointer:<externalptr> 
  .. .. ..$ docstring    : chr ""
  .. .. ..and 14 methods.
  .. ..$ time        :Reference class 'C++Field' [package "Rcpp"] with 5 fields
  .. .. ..$ pointer      :<externalptr> 
  .. .. ..$ cpp_class    : chr "double"
  .. .. ..$ read_only    : logi FALSE
  .. .. ..$ class_pointer:<externalptr> 
  .. .. ..$ docstring    : chr ""
  .. .. ..and 14 methods.
  .. ..$ trajectories:Reference class 'C++Field' [package "Rcpp"] with 5 fields
  .. .. ..$ pointer      :<externalptr> 
  .. .. ..$ cpp_class    : chr "arma::Mat<double>"
  .. .. ..$ read_only    : logi FALSE
  .. .. ..$ class_pointer:<externalptr> 
  .. .. ..$ docstring    : chr ""
  .. .. ..and 14 methods.
  ..@ methods     :List of 2
  .. ..$ Cbdm:Reference class 'C++OverloadedMethods' [package "Rcpp"] with 8 fields
  .. .. ..$ pointer      :<externalptr> 
  .. .. ..$ class_pointer:<externalptr> 
  .. .. ..$ size         : int 1
  .. .. ..$ void         : logi TRUE
  .. .. ..$ const        : logi FALSE
  .. .. ..$ docstrings   : chr ""
  .. .. ..$ signatures   : chr "void Cbdm(int)"
  .. .. ..$ nargs        : int 1
  .. .. ..and 15 methods, of which 1 is  possibly relevant:
  .. .. ..  info
  .. ..$ Tbdm:Reference class 'C++OverloadedMethods' [package "Rcpp"] with 8 fields
  .. .. ..$ pointer      :<externalptr> 
  .. .. ..$ class_pointer:<externalptr> 
  .. .. ..$ size         : int 1
  .. .. ..$ void         : logi TRUE
  .. .. ..$ const        : logi FALSE
  .. .. ..$ docstrings   : chr ""
  .. .. ..$ signatures   : chr "void Tbdm(double)"
  .. .. ..$ nargs        : int 1
  .. .. ..and 15 methods, of which 1 is  possibly relevant:
  .. .. ..  info
  ..@ constructors:List of 1
  .. ..$ :Reference class 'C++Constructor' [package "Rcpp"] with 5 fields
  .. .. ..$ pointer      :<externalptr> 
  .. .. ..$ class_pointer:<externalptr> 
  .. .. ..$ nargs        : int 3
  .. .. ..$ signature    : chr "Community(arma::Col<double>, arma::Mat<double>, double)"
  .. .. ..$ docstring    : chr ""
  .. .. ..and 14 methods.
  ..@ generator   :Formal class 'refObjectGenerator' [package "methods"] with 4 slots
  .. .. ..@ .Data    :function (...)  
  .. .. ..@ generator:Reference class 'refGeneratorSlot' [package "methods"] with 2 fields
  .. .. .. ..$ def      :Formal class 'refClassRepresentation' [package "methods"] with 15 slots
  .. .. .. .. .. ..@ fieldClasses   :List of 11
  .. .. .. .. .. .. ..$ K           : chr "activeBindingFunction"
  .. .. .. .. .. .. ..$ abundance   : chr "activeBindingFunction"
  .. .. .. .. .. .. ..$ b           : chr "activeBindingFunction"
  .. .. .. .. .. .. ..$ cycles      : chr "activeBindingFunction"
  .. .. .. .. .. .. ..$ d0          : chr "activeBindingFunction"
  .. .. .. .. .. .. ..$ interaction : chr "activeBindingFunction"
  .. .. .. .. .. .. ..$ m           : chr "activeBindingFunction"
  .. .. .. .. .. .. ..$ save_int    : chr "activeBindingFunction"
  .. .. .. .. .. .. ..$ stochastic  : chr "activeBindingFunction"
  .. .. .. .. .. .. ..$ time        : chr "activeBindingFunction"
  .. .. .. .. .. .. ..$ trajectories: chr "activeBindingFunction"
  .. .. .. .. .. ..@ fieldPrototypes:<environment: 0x458b940> 
  .. .. .. .. .. ..@ refMethods     :<environment: 0x458c620> 
  .. .. .. .. .. ..@ refSuperClasses: chr "envRefClass"
  .. .. .. .. .. ..@ slots          :List of 1
  .. .. .. .. .. .. ..$ .xData: atomic [1:1] environment
  .. .. .. .. .. .. .. ..- attr(*, "package")= chr "methods"
  .. .. .. .. .. ..@ contains       :List of 6
  .. .. .. .. .. .. ..$ C++Object   :Formal class 'SClassExtension' [package "methods"] with 10 slots
  .. .. .. .. .. .. .. .. ..@ subClass  : chr "Rcpp_Community"
  .. .. .. .. .. .. .. .. ..@ superClass: chr "C++Object"
  .. .. .. .. .. .. .. .. ..@ package   : chr "GillesCom"
  .. .. .. .. .. .. .. .. ..@ coerce    :function (from, strict = TRUE)  
  .. .. .. .. .. .. .. .. ..@ test      :function (object)  
  .. .. .. .. .. .. .. .. ..@ replace   :function (from, to, value)  
  .. .. .. .. .. .. .. .. ..@ simple    : logi TRUE
  .. .. .. .. .. .. .. .. ..@ by        : chr(0) 
  .. .. .. .. .. .. .. .. ..@ dataPart  : logi FALSE
  .. .. .. .. .. .. .. .. ..@ distance  : num 1
  .. .. .. .. .. .. ..$ envRefClass :Formal class 'SClassExtension' [package "methods"] with 10 slots
  .. .. .. .. .. .. .. .. ..@ subClass  : atomic [1:1] Rcpp_Community
  .. .. .. .. .. .. .. .. .. ..- attr(*, "package")= chr "GillesCom"
  .. .. .. .. .. .. .. .. ..@ superClass: atomic [1:1] envRefClass
  .. .. .. .. .. .. .. .. .. ..- attr(*, "package")= chr "methods"
  .. .. .. .. .. .. .. .. ..@ package   : chr "GillesCom"
  .. .. .. .. .. .. .. .. ..@ coerce    :function (from, strict = TRUE)  
  .. .. .. .. .. .. .. .. ..@ test      :function (object)  
  .. .. .. .. .. .. .. .. ..@ replace   :function (from, to, value)  
  .. .. .. .. .. .. .. .. ..@ simple    : logi TRUE
  .. .. .. .. .. .. .. .. ..@ by        : chr(0) 
  .. .. .. .. .. .. .. .. ..@ dataPart  : logi FALSE
  .. .. .. .. .. .. .. .. ..@ distance  : num 1
  .. .. .. .. .. .. ..$ .environment:Formal class 'SClassExtension' [package "methods"] with 10 slots
  .. .. .. .. .. .. .. .. ..@ subClass  : atomic [1:1] Rcpp_Community
  .. .. .. .. .. .. .. .. .. ..- attr(*, "package")= chr "GillesCom"
  .. .. .. .. .. .. .. .. ..@ superClass: chr ".environment"
  .. .. .. .. .. .. .. .. ..@ package   : chr "methods"
  .. .. .. .. .. .. .. .. ..@ coerce    :function (from, strict = TRUE)  
  .. .. .. .. .. .. .. .. ..@ test      :function (object)  
  .. .. .. .. .. .. .. .. ..@ replace   :function (from, to, value)  
  .. .. .. .. .. .. .. .. ..@ simple    : logi TRUE
  .. .. .. .. .. .. .. .. ..@ by        : chr "envRefClass"
  .. .. .. .. .. .. .. .. ..@ dataPart  : logi FALSE
  .. .. .. .. .. .. .. .. ..@ distance  : num 2
  .. .. .. .. .. .. ..$ refClass    :Formal class 'SClassExtension' [package "methods"] with 10 slots
  .. .. .. .. .. .. .. .. ..@ subClass  : atomic [1:1] Rcpp_Community
  .. .. .. .. .. .. .. .. .. ..- attr(*, "package")= chr "GillesCom"
  .. .. .. .. .. .. .. .. ..@ superClass: chr "refClass"
  .. .. .. .. .. .. .. .. ..@ package   : chr "methods"
  .. .. .. .. .. .. .. .. ..@ coerce    :function (from, strict = TRUE)  
  .. .. .. .. .. .. .. .. ..@ test      :function (object)  
  .. .. .. .. .. .. .. .. ..@ replace   :function (from, to, value)  
  .. .. .. .. .. .. .. .. ..@ simple    : logi TRUE
  .. .. .. .. .. .. .. .. ..@ by        : chr "envRefClass"
  .. .. .. .. .. .. .. .. ..@ dataPart  : logi FALSE
  .. .. .. .. .. .. .. .. ..@ distance  : num 2
  .. .. .. .. .. .. ..$ environment :Formal class 'SClassExtension' [package "methods"] with 10 slots
  .. .. .. .. .. .. .. .. ..@ subClass  : atomic [1:1] Rcpp_Community
  .. .. .. .. .. .. .. .. .. ..- attr(*, "package")= chr "GillesCom"
  .. .. .. .. .. .. .. .. ..@ superClass: chr "environment"
  .. .. .. .. .. .. .. .. ..@ package   : chr "methods"
  .. .. .. .. .. .. .. .. ..@ coerce    :function (from, strict = TRUE)  
  .. .. .. .. .. .. .. .. ..@ test      :function (object)  
  .. .. .. .. .. .. .. .. ..@ replace   :function (from, to, value)  
  .. .. .. .. .. .. .. .. ..@ simple    : logi FALSE
  .. .. .. .. .. .. .. .. ..@ by        : chr "envRefClass"
  .. .. .. .. .. .. .. .. ..@ dataPart  : logi FALSE
  .. .. .. .. .. .. .. .. ..@ distance  : num 3
  .. .. .. .. .. .. ..$ refObject   :Formal class 'SClassExtension' [package "methods"] with 10 slots
  .. .. .. .. .. .. .. .. ..@ subClass  : atomic [1:1] Rcpp_Community
  .. .. .. .. .. .. .. .. .. ..- attr(*, "package")= chr "GillesCom"
  .. .. .. .. .. .. .. .. ..@ superClass: chr "refObject"
  .. .. .. .. .. .. .. .. ..@ package   : chr "methods"
  .. .. .. .. .. .. .. .. ..@ coerce    :function (from, strict = TRUE)  
  .. .. .. .. .. .. .. .. ..@ test      :function (object)  
  .. .. .. .. .. .. .. .. ..@ replace   :function (from, to, value)  
  .. .. .. .. .. .. .. .. ..@ simple    : logi TRUE
  .. .. .. .. .. .. .. .. ..@ by        : atomic [1:1] envRefClass
  .. .. .. .. .. .. .. .. .. ..- attr(*, "package")= chr "methods"
  .. .. .. .. .. .. .. .. ..@ dataPart  : logi FALSE
  .. .. .. .. .. .. .. .. ..@ distance  : num 3
  .. .. .. .. .. ..@ virtual        : logi FALSE
  .. .. .. .. .. ..@ prototype      :Formal class 'S4' [package ""] with 0 slots
 list()
  .. .. .. .. .. ..@ validity       : NULL
  .. .. .. .. .. ..@ access         : list()
  .. .. .. .. .. ..@ className      : atomic [1:1] Rcpp_Community
  .. .. .. .. .. .. ..- attr(*, "package")= chr "GillesCom"
  .. .. .. .. .. ..@ package        : chr "GillesCom"
  .. .. .. .. .. ..@ subclasses     : list()
  .. .. .. .. .. ..@ versionKey     :<externalptr> 
  .. .. .. .. .. ..@ sealed         : logi FALSE
  .. .. .. ..$ className: chr "Rcpp_Community"
  .. .. .. ..and 20 methods, of which 6 are  possibly relevant:
  .. .. .. ..  accessors, fields, help, lock, methods, new
  .. .. ..@ className: atomic [1:1] Rcpp_Community
  .. .. .. ..- attr(*, "package")= chr "GillesCom"
  .. .. ..@ package  : chr "GillesCom"
  ..@ docstring   : chr ""
  ..@ typeid      : chr "9Community"
  ..@ enums       : Named list()
  ..@ parents     : chr(0) 
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(Rcpp_Community)
## maybe str(Rcpp_Community) ; plot(Rcpp_Community) ...
}
\keyword{datasets}
