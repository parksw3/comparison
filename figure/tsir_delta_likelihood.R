library(Deriv)

load("../sim/gillespie_sim.rda")
load("../data/gillespie_data.rda")

N <- 1e5

i <- 1

dd <- datalist[[i]][1:20,]
rr <- reslist[[i]]

I <- dd$incidence/0.7
Inew <- tail(I, -1)
Iprev <- head(I, -1)

## assume that we know S exactly
S <- approx(x=rr$time, y=N-rr$incidence-10, xout=1:19)$y

lfit <- lm(log(Inew) ~ 1 + log(Iprev) + offset(log(S/N)))

hatbeta <- fixfun(coef(lfit)[[1]], coef(lfit)[[2]], N, I)

dl <- fixfun_deriv(coef(lfit)[[1]], coef(lfit)[[2]], N, I) 

hatbeta_sd <- sqrt(t(dl) %*% vcov(lfit) %*% dl)[[1]]

pp <- profile_likelihood(betavec=seq(1.5, 2.5, by=0.01), N=N, I=I, S=S, lfit=lfit)

plot(pp$beta, pp$nll + logLik(lfit), type="l")
curve((x-hatbeta)^2/(2 * hatbeta_sd^2), add=T, col=2)
