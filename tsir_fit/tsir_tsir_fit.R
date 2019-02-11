load("../data/tsir_data.rda")

N <- 1e5
nsim <- length(datalist)

fitlist <- vector('list', nsim)

for (i in 1:nsim) {
	dd <- datalist[[i]][1:20,]
	
	I <- dd$incidence/0.7
	Inew <- tail(I, -1)
	Iprev <- head(I, -1)
	
	## assume that we know S exactly
	S <- head(dd$S, -1)
	
	lfit <- lm(log(Inew) ~ 1 + log(Iprev) + offset(log(S/N)))

	cdata <- data.frame(
		param=c("beta", "alpha"),
		mean=c(exp(coef(lfit)[[1]]), coef(lfit)[[2]]),
		lwr=c(exp(confint(lfit)[1,1]), confint(lfit)[2,1]),
		upr=c(exp(confint(lfit)[1,2]), confint(lfit)[2,2])
	)
	
	cdata$coverage <- c(
		cdata$lwr[1] < 2 && 2 < cdata$upr[1],
		cdata$lwr[2] < 0.97 && 0.97 < cdata$upr[2]
	)
	
	fitlist[[i]] <- cdata
}

save("fitlist", file="tsir_tsir_fit.rda")
