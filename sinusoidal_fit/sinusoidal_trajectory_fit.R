library(Matrix)
library(deSolve)
library(bbmle)
load("../data/gillespie_sinusoidal_data.rda")

sir <- function(t, state, parameters) {
	with(as.list(c(state, parameters)), {
		beta <- b0 * (1 + b1 * cos(2 * pi * t/ 26))
		I <- exp(logI)
		dS <- mu * N - beta * S * I/N
		dlogI <- beta * S/N - gamma
		dR <- gamma * I
		dC <- beta * S * I/N
		list(c(dS, dlogI, dR, dC))
	})
}

objfun <- function(param, 
				   data,
				   N=5e6) {
	out <- with(as.list(param), {
		par <- c(b0=exp(log.b0), b1=plogis(logit.b1), mu=1/(50*26), gamma=1, N=N)
		I0 <- plogis(logit.I0)
		S0 <- plogis(logit.S0)
		y <- c(S=N*S0, logI=log(N*I0), R=0, C=0)
		
		times <- c(0:nrow(data))
		
		out <- as.data.frame(rk4(y, times, sir, par))
		
		out
	})
	
	rprob <- plogis(param[["logit.rprob"]])
	size <- exp(param[["log.size"]])

	nll <- -sum(dnbinom(data$incidence, mu=diff(out$C)*rprob, size=size, log=TRUE))
	
	nll
}

parnames(objfun) <- c("log.b0", "logit.b1", "logit.rprob", 
					  "logit.S0", "logit.I0", "log.size")

N <- 5e6
nsim <- length(datalist)

reslist <- fitlist <- translist <- vector('list', nsim)

for (i in 1:nsim) {
	print(i)
	dd <- datalist[[i]]
	
	start <- c(log.b0=log(500/26), logit.b1=qlogis(0.15), logit.rprob=qlogis(0.7),
			   logit.S0=qlogis(0.05), logit.I0=qlogis(1e-4), log.size=log(10))
	
	m <- mle2(objfun, start, data=list(data=dd), vecpar=TRUE, method="BFGS")
	
	m@details$hessian <- as.matrix(nearPD(m@details$hessian)$mat)
	
	pp <- profile(m, which=c(1, 3), trace=TRUE, zmax=2, continuation="naive")
	
	cc <- confint(pp)
	
	cdata <- data.frame(
		param=c("R0", "rprob"),
		mean=c(exp(coef(m))[[1]], plogis(coef(m))[[3]]),
		lwr=c(exp(cc[1,1]), plogis(cc[2,1])),
		upr=c(exp(cc[1,2]), plogis(cc[2,2]))
	)
	
	cdata$coverage <- c(
		cdata$lwr[1] < 500/26 && 500/26 < cdata$upr[1],
		cdata$lwr[2] < 0.7 && 0.7 < cdata$upr[2]
	)
	
	fitlist[[i]] <- cdata
	translist[[i]] <- data.frame(
		time=seq(1, 26, by=0.01),
		beta=exp(coef(m)[[1]]) * (1 + plogis(coef(m)[[2]]) * cos(2 * pi * seq(1, 26, by=0.01)/26))
	)
	
	reslist[[i]] <- m
	
	save("fitlist", "translist", "reslist", file="sinusoidal_trajectory_fit.rda")
}

save("fitlist", "translist", "reslist", file="sinusoidal_trajectory_fit.rda")
