fixfun <- function(logbeta, alpha, N, I) {
	mI <- mean(I)
	beta <- exp(logbeta)
	
	-log(1 - beta * mI^alpha/N) * N/mI
}

fixfun_deriv <- Deriv::Deriv(fixfun, c("logbeta", "alpha"))

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
