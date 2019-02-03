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
	
	logI <- log(dd$incidence)
	
	gfit <- MASS::glm.nb(incidence~ns(time, knots=seq(1, 260, by=6)), data=dd)
	
	S0vec <- seq(0.03, 0.07, by=0.001)
	llvec <- rep(NA, length(S0vec))
	
	for (j in 1:length(S0vec)) {
		fitdata <- data.frame(
			grad=diff(predict(gfit, newdata = data.frame(time=1:261)))+1,
			logS=log(S0vec[j]*N + Z),
			biweek=dd$biweek
		)
		
		fitdata$offterm <- fitdata$logS - log(N)
		
		lfit <- gam(grad ~ s(biweek, bs="cc") + offset(offterm), data=fitdata,
					family = gaussian("log"))
		
		llvec[j] <- logLik(lfit)
	}
	
	j <- which.max(llvec)
	
	fitdata <- data.frame(
		grad=diff(predict(gfit, newdata = data.frame(time=1:261)))+1,
		logS=log(S0vec[j]*N + Z),
		biweek=dd$biweek
	)
	
	fitdata$offterm <- fitdata$logS - log(N)
	
	lfit <- gam(grad ~ s(biweek, bs="cc") + offset(offterm), data=fitdata,
				family = gaussian("log"))
	
	logR0 <- coef(lfit)[[1]]
	
	cdata <- data.frame(
		param=c("R0", "rprob"),
		mean=c(exp(logR0), 1/coef(regfit)[[2]]),
		lwr=c(exp(confint.default(lfit, 1)[[1]]) , 1/confint(regfit, 2)[[2]]),
		upr=c(exp(confint.default(lfit, 1)[[2]]), 1/confint(regfit, 2)[[1]])
	)
	
	cdata$coverage <- c(
		cdata$lwr[1] < 500/26 && 500/26 < cdata$upr[1],
		cdata$lwr[2] < 0.7 && 0.7 < cdata$upr[2]
	)
	
	fitlist[[i]] <- cdata
	translist[[i]] <- data.frame(
		time=seq(1, 26, by=0.01),
		beta=exp(predict(lfit, newdata=data.frame(biweek=seq(1, 26, by=0.01), offterm=0)))
	)
}

save("fitlist", "translist", file="sinusoidal_gradient_fit.rda")
