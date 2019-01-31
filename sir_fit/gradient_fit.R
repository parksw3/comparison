library(mgcv)
library(deSolve)
library(bbmle)
load("../data/gillespie_data.rda")

sir <- function(t, y, pars, transfun) {
	with(as.list(c(y, pars)), {
		I <- exp(logI)
		inf <- transfun(t)
		
		dS <- - inf
		dlogI <- inf/I - gamma
		dCI <- inf
		
		list(c(dS, dlogI, dCI))
	})
}

objfun <- function(param, 
				   data,
				   N=1e5,
				   gfit) {
	beta <- exp(param[["log.beta"]])
	rprob <- plogis(param[["logit.rprob"]])
	size <- exp(param[["log.size"]])
	
	transfun <- function(t) predict(gfit, type="response", newdata=data.frame(time=t))/rprob
	
	times <- seq(0, nrow(data), by=1)
	
	out <- with(as.list(param), {
		par <- c(gamma=1, N=N)
		I0 <- plogis(logit.I0)
		y <- c(S=N*(1-I0), logI=log(N*I0), CI=0)
		
		out <- as.data.frame(rk4(y, times, sir, par, transfun=transfun))
		
		out
	})
	
	ind <- seq(1, nrow(data)*1+1, by=1)
	
	nll <- -sum(dnbinom(data$incidence, mu=diff(out$CI[ind])*rprob, size=size, log=TRUE)) # +
#		sum((transfun(times)*rprob - beta * out$S * exp(out$logI)/N*rprob)^2)
	
	nll
}

parnames(objfun) <- c("log.beta", "logit.rprob", "logit.I0", "log.size")

N <- 1e5
nsim <- length(datalist)

for (i in 1:nsim) {
	print(i)
	dd <- datalist[[i]][1:20,]
	
	gfit <- gam(incidence~s(time, bs="cr"), data=dd, family=nb)
	
	start <- c(log.beta=log(2), logit.rprob=qlogis(0.7), logit.I0=qlogis(10/1e5), log.size=log(1))
	
	m <- mle2(objfun, start, data=list(data=dd, gfit=gfit), vecpar=TRUE)
	
	cc <- confint(m)
	
	cdata <- data.frame(
		param=c("rprob", "I0", "size"),
		mean=c(plogis(coef(m))[[1]], plogis(coef(m))[[2]], exp(coef(m))[[3]]),
		lwr=c(plogis(cc[1,1]), plogis(cc[2,1]), exp(cc[3,1])),
		upr=c(plogis(cc[1,2]), plogis(cc[2,2]), exp(cc[3,2]))
	)
	
	fitlist[[i]] <- cdata
}

save("fitlist", file="gradient_fit.rda")
