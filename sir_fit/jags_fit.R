library(rjags)
load("../data/gillespie_data.rda")

nsim <- length(datalist)

fitlist <- vector('list', nsim)

for (i in 1:nsim) {
	print(i)
	
	dd <- datalist[[i]][1:20,]
	
	jagsdata <- list(
		ell=40,
		N0=1e5,
		tmax=20,
		nstep=10,
		incidence=dd$incidence
	)
	
	jags <- jags.model('../bugscode/sir.bug',
					   data = jagsdata,
					   n.chains = 1,
					   n.adapt = 2000)
	
	update(jags, 2000)
	
	jj <- jags.samples(jags,c("I", "mu", "R0", "reporting", "Gmean", "Gvar", "r",
							  "khat", "k", "I0"), 2000)
	
	
	cdata <- data.frame(
		param=c("R0", "rprob", "I0", "size"),
		mean=c(mean(jj$R0), mean(jj$reporting), mean(jj$I0), mean(jj$r)),
		lwr=c(quantile(jj$R0, 0.025), quantile(jj$reporting, 0.025), NA, NA),
		upr=c(quantile(jj$R0, 0.975), quantile(jj$reporting, 0.975), NA, NA)
	)
	
	
	cdata$coverage <- c(
		cdata$lwr[1] < 2 && 2 < cdata$upr[1],
		cdata$lwr[2] < 0.7 && 0.7 < cdata$upr[2],
		cdata$lwr[3] < 1e-4 && 1e-4 < cdata$upr[3],
		NA
	)
	
	fitlist[[i]] <- cdata
	
	save("fitlist", file="jags_fit.rda")
}

save("fitlist", file="jags_fit.rda")
