library(rjags)
load("../data/gillespie_data.rda")

argvals <- commandArgs(trailingOnly=TRUE)
batch_num <- as.numeric(argvals[1])

fn <- paste0("jags_fit_", batch_num, ".rda")

nsim <- 10

jagslist <- fitlist <- vector('list', nsim)

set.seed(101)
for (i in 1:nsim) {
	print(i)
	
	j <- batch_num * 10 + i
	
	dd <- datalist[[j]][1:20,]
	
	jagsdata <- list(
		ell=40,
		N0=1e5,
		tmax=20,
		nstep=10,
		incidence=dd$incidence
	)
	
	jags <- jags.model('../bugscode/sir.bug',
					   data = jagsdata,
					   n.chains = 4,
					   n.adapt = 2000)
	
	update(jags, 2000)
	
	jj <- jags.samples(jags,c("R0", "r", "reporting", "I0"), 2000, thin=10)
	
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
	jagslist[[i]] <- jj
	
	save("fitlist", "jagslist", file=fn)
}

save("fitlist", "jagslist", file=fn)
