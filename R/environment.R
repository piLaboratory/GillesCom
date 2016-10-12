#' Random environment generator
#'
#' This function is a simple generator for a stochastic data frame
#' that can be used to represent a fluctuating environment
#' that performs a random walk
#' for use in \code{\link{Init_Community}}.
#'
#' @export
#' @rdname environment
#' @param time Final time
#' @param sd Standard deviation of the random walk process
#' @param intervals Number of steps of the random walk
stoch_environment <- function(time, sd, intervals = 10000) {
    K <- exp(cumsum(rnorm(intervals, sd=sd)))
    data.frame(time = seq(0, time, length=intervals), K)
}
