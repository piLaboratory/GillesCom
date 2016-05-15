#source("bdm.R")
## Migration rates for each specie in the metacommunity
#m.test <- ls_migration(1e6, alpha=50, 1e-5)
#run.bdm(N0=rep(c(1,0),c(1,length(m.test)-1)),
#        K=1e4, b=0.1, m=m.test, con=1, stren=0.1, nrep=1e6, rec.step=1e4, file="teste.dat")
## Hubbell's neutral model
#a2 <- matrix(1, ncol=length(m.test), nrow=length(m.test))
#run.bdm(N0=rep(c(1,0),c(1,length(m.test)-1)), alphas=a2,
#        K=1e4, b=0.1, m=m.test, con=1, stren=0.1, nrep=1e6, rec.step=1e4, file="teste2.dat")
## Caswell's Neutral model
#a3 <- matrix(0, ncol=length(m.test), nrow=length(m.test))
#diag(a3) <- 1
#run.bdm(N0=rep(c(1,0),c(1,length(m.test)-1)), alphas=a3,
#        K=1e4, b=0.1, m=m.test, con=1, stren=0.1, nrep=1e6, rec.step=1e4, file="teste3.dat")
#library(parallel) ## parece que deu algum erro
#cl <- makeCluster(3)
#clusterExport(cl, c("a2", "a3", "m.test", "bdm"))
#teste.cl <- parLapply(cl, list(NULL, a2, a3), run.bdm, N0=rep(c(1,0),c(1,length(m.test)-1)),
#          K=1e4, b=0.1, m=m.test, con=1, stren=0.1, nrep=1e6, rec.step=1e4, return.df=TRUE) ## retornou erro
#stopCluster(cl)
## Checking if Number of species stabilized
##results <- read.table("teste.dat", header=TRUE)
##results <- read.table("teste2.dat", header=TRUE)
##results <- read.table("teste3.dat", header=TRUE)
#results <- read.table("Sat_Aug_15_05:21:09_2015", header=TRUE)
#results <- read.table("Sat_Aug_15_05:22:45_2015", header=TRUE)
#results <- read.table("Sat_Aug_15_05:25:18_2015", header=TRUE)
## Checking if S and N stabilized    
#results.S <- aggregate(results$N, list(time=results$time), function(x)sum(x>0)) ## algum problema
#plot(results.S, type="b")
#results.N <- aggregate(results$N, list(time=results$time), sum)
#plot(results.N, type="b")
## sads
#plot(rad(results$N[results$time==max(results$time)]))
## Trajectories
#ggplot(results, aes(x=time, y=N, group=factor(sp))) + geom_line()
