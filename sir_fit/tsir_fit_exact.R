library(Deriv)
source("../R/tsir_util.R")

load("../sim/gillespie_sim.rda")
load("../data/gillespie_data.rda")

N <- 1e5
nsim <- length(datalist)

fitlist <- vector('list', nsim)

for (i in 1:nsim) {
	print(i)
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
	
	pp <- profile_likelihood(betavec=seq(1.5, 3, by=0.05), N=N, S=S, I=I)
	pp$nll <- pp$nll + logLik(lfit)
	pp$nll[pp$beta < hatbeta] <- - pp$nll[pp$beta < hatbeta]
	
	ci <- sort(approx(x=pp$nll, y=pp$beta, xout=c(-qchisq(0.95, 1)/2, qchisq(0.95, 1)/2))$y)
	
	cdata <- data.frame(
		param=c("beta", "alpha"),
		mean=c(hatbeta, coef(lfit)[[2]]),
		lwr=c(ci[1], confint(lfit)[2,1]),
		upr=c(ci[2], confint(lfit)[2,2])
	)
	
	cdata$coverage <- c(
		cdata$lwr[1] < 2 && 2 < cdata$upr[1],
		NA
	)
	
	fitlist[[i]] <- cdata
}

save("fitlist", file="tsir_fit_exact.rda")
