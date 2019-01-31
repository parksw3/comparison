library(deSolve)
library(bbmle)
load("../data/gillespie_data.rda")

if (fe <- file.exists("trajectory_fit.rda")) {
	load("trajectory_fit.rda")
}

sir <- function(t, y, pars) {
	with(as.list(c(y, pars)), {
		inf <- beta * S/ N
		I <- exp(logI)
		dS <- - inf * I
		dlogI <- inf - gamma
		dCI <- inf * I
		
		list(c(dS, dlogI, dCI))
	})
}

objfun <- function(param, 
				   data,
				   N=1e5) {
	out <- with(as.list(param), {
		par <- c(beta=exp(log.beta), gamma=1, N=N)
		I0 <- plogis(logit.I0)
		y <- c(S=N*(1-I0), logI=log(N*I0), CI=0)
		
		times <- c(0:nrow(data))
		
		out <- as.data.frame(rk4(y, times, sir, par))
		
		out
	})
	
	rprob <- plogis(param[["logit.rprob"]])
	size <- exp(param[["log.size"]])

	nll <- -sum(dnbinom(data$incidence, mu=diff(out$CI)*rprob, size=size, log=TRUE))
	
	nll
}

parnames(objfun) <- c("log.beta", "logit.rprob", "logit.I0", "log.size")

N <- 1e5
nsim <- length(datalist)

if (!fe) {
	fitlist <- vector('list', nsim)
}

init <- sum(sapply(fitlist, length)!=0)

for (i in (init+1):nsim) {
	print(i)
	dd <- datalist[[i]][1:20,]
	
	start <- c(log.beta=log(2), logit.rprob=qlogis(0.7), logit.I0=qlogis(10/1e5), log.size=log(1))
	
	m <- mle2(objfun, start, data=list(data=dd), vecpar=TRUE)
	
	cc <- confint(m)
	
	cdata <- data.frame(
		param=c("beta", "rprob", "I0", "size"),
		mean=c(exp(coef(m))[[1]], plogis(coef(m))[[2]], plogis(coef(m))[[3]], exp(coef(m))[[4]]),
		lwr=c(exp(cc[1,1]), plogis(cc[2,1]), plogis(cc[3,1]), exp(cc[4,1])),
		upr=c(exp(cc[1,2]), plogis(cc[2,2]), plogis(cc[3,2]), exp(cc[4,2]))
	)
	
	cdata$coverage <- c(
		cdata$lwr[1] < 2 && 2 < cdata$upr[1],
		cdata$lwr[2] < 0.7 && 0.7 < cdata$upr[2],
		cdata$lwr[3] < 1e-4 && 1e-4 < cdata$upr[3],
		NA
	)
	
	fitlist[[i]] <- cdata
	
	save("fitlist", file="trajectory_fit.rda")
}

save("fitlist", file="trajectory_fit.rda")
