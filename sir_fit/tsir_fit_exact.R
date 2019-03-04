library(Deriv)

load("../sim/gillespie_sim.rda")
load("../data/gillespie_data.rda")

fixfun <- function(logbeta, alpha, N, I) {
	mI <- mean(I)
	beta <- exp(logbeta)
	
	-log(1 - beta * mI^alpha/N) * N/mI
}

fixfun_deriv <- Deriv(fixfun, c("logbeta", "alpha"))

N <- 1e5
nsim <- length(datalist)

fitlist <- vector('list', nsim)

for (i in 1:nsim) {
	dd <- datalist[[i]][1:20,]
	rr <- reslist[[i]]
	
	ii <- c(rr$incidence[1], diff(rr$incidence))
	tt <- rr$time
	
	dd <- as.data.frame(table(ceiling(tt[ii==1])))
	
	I <- dd$Freq[1:20]
	
	Inew <- tail(I, -1)
	Iprev <- head(I, -1)
	
	S <- head(N - cumsum(I) - 10, -1)
	
	lfit <- lm(log(Inew) ~ 1 + log(Iprev) + offset(log(S/N)))
	
	hatbeta <- fixfun(coef(lfit)[[1]], coef(lfit)[[2]], N, I)
	
	dl <- fixfun_deriv(coef(lfit)[[1]], coef(lfit)[[2]], N, I) 
	
	hatbeta_sd <- sqrt(t(dl) %*% vcov(lfit) %*% dl)[[1]]
	
	cdata <- data.frame(
		param=c("beta", "alpha"),
		mean=c(hatbeta, coef(lfit)[[2]]),
		lwr=c(hatbeta - 1.96 * hatbeta_sd, confint(lfit)[2,1]),
		upr=c(hatbeta + 1.96 * hatbeta_sd, confint(lfit)[2,2])
	)
	
	cdata$coverage <- c(
		cdata$lwr[1] < 2 && 2 < cdata$upr[1],
		NA
	)
	
	fitlist[[i]] <- cdata
}

save("fitlist", file="tsir_fit_exact.rda")
