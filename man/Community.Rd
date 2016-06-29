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
\usage{new(Community)}
\format{
    The class Community contains the following slots:

    \code{$abundance} contains the current abundance vector for the community.

    \code{$time} contains the current simulation time for the community.

    \code{$cycles} contains the number of elapsed simulation cycles.

    \code{$trajectories} contains a data frame in which each line corresponds to the species abundance
distribution at a different time.

    \code{$K} contains a vector with the carrying capacity of each species.

    \code{$b} contains a vector with the birth rate of each species.

    \code{$d0} contains a vector with the base death rate of each species (the actual death rate is affected by the
relation between total populations and carrying capacities).

    \code{$interaction} contains an interaction matrix detailing how each species affects one another. 

    \code{$m} contains a vector with the migration rates of each species. See \code{\link{ls_migration}}.

    \code{$save_int} contains the saving interval, which is the simulation time that must go by between "snapshots"
    that are saved in the slot \code{$trajectories}.

    \code{stochastic} contains a \code{data.frame} related to how demographic stochasticity affects the
    community. See the vignettes for details.
}
\details{
    Objects of the class \code{Community} should _only_ be created with the \code{\link{Init_Community}} function
    or loaded from a saved file with the \code{\link{GillesComFromFile}} function.
}
\examples{
    com = Init_Community(10)
    bdm(com, 99)
    com$K
    com$abundance
    com$time
}
\keyword{datasets}
