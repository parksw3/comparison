library(mgcv)
load("../data/gillespie_sinusoidal_data.rda")

N <- 5e6
nsim <- length(datalist)

fitlist <- vector('list', nsim)
translist <- vector('list', nsim)

for (i in 1:nsim) {
	print(i)
	dd <- datalist[[i]]
	
	regdata <- data.frame(
		x=cumsum(dd$incidence),
		y=cumsum(dd$birth)
	)
	
	regfit <- lm(y~x, data=regdata)
	
	Z <- regfit$residuals
	rho <- 1/coef(regfit)[[2]]
	
	I <- dd$incidence/rho

	S0vec <- seq(0.03, 0.08, by=0.001)
	llvec <- rep(NA, length(S0vec))
	
	for (j in 1:length(S0vec)) {
		
		fitdata <- data.frame(
			logInew=tail(log(I), -1),
			logIprev=head(log(I), -1),
			S=head(S0vec[j]*N + Z, -1),
			N=N,
			biweek=head(dd$biweek, -1)
		)
		
		fitdata$offterm <- log(fitdata$S) - log(fitdata$N)
		
		lfit <- gam(logInew ~  s(biweek, bs="cc") + logIprev + offset(offterm), data=fitdata)
		
		llvec[j] <- logLik(lfit)
	}
	
	j <- which.max(llvec)
	
	fitdata <- data.frame(
		logInew=tail(log(I), -1),
		logIprev=head(log(I), -1),
		S=head(S0vec[j]*N + Z, -1),
		N=N,
		biweek=head(dd$biweek, -1)
	)
	
	fitdata$offterm <- log(fitdata$S) - log(fitdata$N)
	
	lfit <- gam(logInew ~ s(biweek, bs="cc") + logIprev + offset(offterm), data=fitdata)
	
	logR0 <- coef(lfit)[[1]]
	
	cdata <- data.frame(
		param=c("R0", "rprob", "alpha"),
		mean=c(exp(logR0), 1/coef(regfit)[[2]], coef(lfit)[[2]]),
		lwr=c(exp(confint.default(lfit, 1)[[1]]) , 1/confint(regfit, 2)[[2]],  confint.default(lfit, 2)[1]),
		upr=c(exp(confint.default(lfit, 1)[[2]]), 1/confint(regfit, 2)[[1]], confint.default(lfit, 2)[2])
	)
	
	cdata$coverage <- c(
		cdata$lwr[1] < 500/26 && 500/26 < cdata$upr[1],
		cdata$lwr[2] < 0.7 && 0.7 < cdata$upr[2],
		NA
	)
	
	fitlist[[i]] <- cdata
	translist[[i]] <- data.frame(
		time=seq(1, 26, by=0.01),
		beta=exp(predict(lfit, newdata=data.frame(biweek=seq(1, 26, by=0.01), logIprev=0, offterm=0)))
	)
}

save("fitlist", "translist", file="sinusoidal_tsir_fit_cyclic.rda")
