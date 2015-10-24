
# Sample data
#cs <- 100
#K <- 100
#d0=0
#b=rep(1, cs)
#m=rep(0.1,cs)
#N0 <- rep(10,cs)
#alphas <- matrix(rnorm(cs*cs)^2, ncol=cs)
#alpha <- c()
#for (j in 1:100) {
#  for (i in 1:1000) N0 <- bdm(N0, alphas, K, d0, b, m)
#  alpha <- c(alpha, coef(fitls(N0[N0>0]))[2])
#}
#plot(alpha)

#' Run one interaction of a Gillespie Algorithm of birth death and migration process in a system of generalized Lotka-Volterra system of competing species
#' @param N0 vector of initial abundances of species in the community (set species not present to zero)
#' @param alphas matrix of interaction coefficients
#' @param K carrying capacities of species
#' @param d0 death rate when N=0
#' @param b birth rates (constant)
#' @param m per capita migration rate in the metacommunity
# TBI: include continuous time record (sampling from an exponential)
Rbdm <- function(N0, alphas, K, d0=0, b, m){
    d <- (b-d0)/K # slope of the density-dependent linear relation of death rate to N
    N <- N0
    N[N0>0] <- N0[N0>0] %*% alphas[N0>0,N0>0]
    dt <- d0+d*N # death rates of each species
    w <- N0*(b+dt) + m  # Gillespie weights for each specie, which are the sum of their rates
    i <- sample((1:length(N)), size=1, prob=w) ## sampling which species will suffer the next action, proportionaly to their weights
    N0[i] <- N0[i] + sample(c(1, -1), size=1, prob= c(b[i]*N0[i]+m[i], dt[i]*N0[i])) ## Sampling if the selected species will gain or loss an individual
    return(N0)
}

#' Generates migration rates from a log-series metacommunity
#' @param J size of metacommunity
#' @param alpha Fisher's alpha
#' @param m per species migration rate
ls.m <- function(J, alpha, m){
    ## expected number of species
    S <- ceiling(alpha*log(1+J/alpha))
    ## Sampling S abundances from a logseries
    N <- qls(runif(S), N=J, alpha=alpha)
    ## Migration rates are wheighted by the abundances in the community
    return(m*N/sum(N))
}


## Running many times
#' @param con connectance of the interaction matrix
#' @param stren strength of interaction matrix, which is the standard deviation of the Gaussin from which the values are drawn
run.bdm <- function(alphas, N0, K, d0=0, b, m, con, stren=0.1, comp=TRUE,  nrep, rec.step, file="teste.dat", return.df=FALSE){
    if(return.df) file=gsub(" ", "_", date()) ## um jeito tosco de fazer arquivos de saida
    J <- length(N0)
    if(length(K)==1) K <- rep(K,J)
    if(length(d0)==1) d0 <- rep(d0,J)
    if(length(b)==1) b <- rep(b,J)
    step <- 0
    ## Assembling the interaction matrix (should be an auxiliary function)
    ##
    df <- data.frame(time=rep(step,J), sp=1:J, N=N0)
    write.table(df, file=file, row.names=FALSE)
    s1 <- bdm(N0=N0, alphas=alphas, K, d0, b, m=m, con=con, stren=stren, comp=comp)
    for(i in 1:nrep){
        step <- step+1
        if(step%%rec.step==0){
            df <- data.frame(time=rep(step,J), sp=1:J, N=s1)
            write.table(df, file=file, append=TRUE, col.names=FALSE, row.names=FALSE) ## 10 times faster than rbinding a dataframe
        }
        s1 <- bdm(N0=s1, alphas=alphas, K=K, d0=d0, b=b, m=m, con=con, stren=stren, comp=comp)
    }
    if(return.df) read.table(file)
}




