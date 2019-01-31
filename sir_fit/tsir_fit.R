load("../data/gillespie_data.rda")

N <- 1e5
nsim <- length(datalist)

fitlist <- vector('list', nsim)

for (i in 1:nsim) {
	dd <- datalist[[i]][1:20,]
	
	I <- dd$incidence/0.7
	Inew <- tail(I, -1)
	Iprev <- head(I, -1)
	
	S <- N - cumsum(Iprev) - 10
	
	lfit <- lm(log(Inew) ~ 1 + log(Iprev) + offset(log(S)))

	cdata <- data.frame(
		param=c("beta", "alpha"),
		mean=c(exp(coef(lfit)[[1]])*N, coef(lfit)[[2]]),
		lwr=c(exp(confint(lfit)[1,1])*N, confint(lfit)[2,1]),
		upr=c(exp(confint(lfit)[1,2])*N, confint(lfit)[2,2])
	)
	
	cdata$coverage <- c(
		cdata$lwr[1] < 2 && 2 < cdata$upr[1],
		NA
	)
	
	fitlist[[i]] <- cdata
}

save("fitlist", file="tsir_fit.rda")
