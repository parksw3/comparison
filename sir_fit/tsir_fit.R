library(Deriv)

load("../sim/gillespie_sim.rda")
load("../data/gillespie_data.rda")

fixfun <- function(beta, alpha, N, I) {
	mI <- mean(I)
	
	-log(1 - beta * mI^alpha/N) * N/mI
}

fixfun_deriv <- Deriv(fixfun, c("beta", "alpha"))

N <- 1e5
nsim <- length(datalist)

fitlist <- vector('list', nsim)

for (i in 1:nsim) {
	dd <- datalist[[i]][1:20,]
	rr <- reslist[[i]]
	
	I <- dd$incidence/0.7
	Inew <- tail(I, -1)
	Iprev <- head(I, -1)
	
	## assume that we know S exactly
	S <- approx(x=rr$time, y=N-rr$incidence-10, xout=1:19)$y
	
	lfit <- lm(log(Inew) ~ 1 + log(Iprev) + offset(log(S/N)))

	hatbeta <- fixfun(exp(coef(lfit))[[1]], coef(lfit)[[2]], N, I)
	
	dl <- fixfun_deriv(exp(coef(lfit))[[1]], coef(lfit)[[2]], N, I) 
	
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

save("fitlist", file="tsir_fit.rda")
