## Como calcular tempo até um evento, quando a taxa do evento varia no tempo
## Necessário para inclusão de estocasticidade demográfica no modelo
## eq 9 Gibson J. Phys. Chem. A 2000 (ver pasta stats/stochast_simula)

## Criando uma taxa de qua varia no tempo de acordo com Wiener process
## supondo uma taxa de b=0.1, fazemos o processo por tempo para quantil de 99.99% (este numero tem que ser ajustado em funcao de t0)
(N <- ceiling(qexp(.9999, 0.1)))
b <- 0.1 + cumsum(rnorm(N*10000, sd=0)) ## para verificar deixe sd=0 , tem que bater com pexp
tempo <- seq(0,N, length=N*10000)
##plot(b~tempo, type="l")
f1 <- approxfun(tempo, b)
## Funcao aproximada para eq 9
f9 <- function(tau, t0, div=1e6){
    index <- seq(t0, tau, length=div)
    x <- f1(index) #a_mu(S,tau)
    cx.dt <- cumsum((tau-t0)*x/div) #integrais de lambda(t) de t=t0 a t=tau
    ecx <- exp(-cx.dt) # exponenciais das taxas até o tempo acumulado
    x.ecx <- ecx*x ## Produto do valor da taxa em cada tempo e a integral acima
    approxfun(y=cumsum((tau-t0)*x.ecx/div), x=index)
    }
f2 <- f9(N,0)
f2(1)
f2(10)
f2(N)
## deve dar valores aproximados
pexp(1,rate=mean(b))
pexp(10,rate=mean(b))
pexp(N,rate=mean(b))
curve(f2, 0,N)
## Agora é generalizar e incluir isso em uma funçao
## ANTES: verificar se está funcionando corretamente para t0>0 (Acho que tem que corrigir o tau, e talvez definir que funcao retorna zero para valores abaixo de t0 e acima de tau)
f2 <- f9(N,80)
f2(N)
plot(f2, 0,N)
