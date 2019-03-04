library(Deriv)

load("../sim/gillespie_sim.rda")
load("../data/gillespie_data.rda")

fixfun <- function(logbeta, alpha, N, I) {
	mI <- mean(I)
	beta <- exp(logbeta)
	
	-log(1 - beta * mI^alpha/N) * N/mI
}

fixfun_deriv <- Deriv(fixfun, c("logbeta", "alpha"))

reversefun <- function(beta, alpha, N, I) {
	mI <- mean(I)
	
	(1 - exp(-beta*mI/N)) * N/mI^alpha
}

reversenll <- function(beta, alpha, N, I, S, lfit) {
	bb <- reversefun(beta, alpha, N, I)
	Iprev <- head(I, -1)
	Inew <- tail(I, -1)
	
	lfit2 <- lm(log(Inew) ~ -1 + offset(log(bb) + alpha * log(Iprev) + log(S/N)))
	
	-logLik(lfit2)[[1]]
}

profile_likelihood <- function (betavec=seq(1, 3, by=0.1),
								N,
								I,
								S,
								lfit) {
	reslist <- vector('list', length(betavec))
	
	for (i in 1:length(betavec)) {
		beta <- betavec[i]
		
		oo <- optim(par=c(alpha=0.95), reversenll, beta=beta, N=N, I=I, S=S, lfit=lfit,
			  method="Brent", lower=0.9, upper=1.1)
		
		reslist[[i]] <- data.frame(
			beta=beta,
			nll=oo$value
		)
	}
	
	do.call("rbind", reslist)
}

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
