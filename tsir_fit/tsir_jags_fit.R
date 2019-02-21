library(rjags)
load("../data/tsir_data.rda")

argvals <- commandArgs(trailingOnly=TRUE)
batch_num <- as.numeric(argvals[1])

fn <- paste0("tsir_jags_fit_", batch_num, ".rda")

nsim <- 10

jagslist <- fitlist <- vector('list', nsim)

set.seed(101)
for (i in 1:nsim) {
	print(i)
	
	j <- batch_num * 10 + i
	
	dd <- datalist[[j]]
	
	jagsdata <- list(
		ell=40,
		Gmean=1,
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
	
	jj <- jags.samples(jags,c("R0", "r", "reporting", "I0", "Gvar"), 2000, thin=10)
	
	cdata <- data.frame(
		param=c("R0", "rprob", "I0", "size"),
		mean=c(mean(jj$R0), mean(jj$reporting), mean(jj$I0), mean(jj$r)),
		lwr=c(quantile(jj$R0, 0.025), NA, NA, NA),
		upr=c(quantile(jj$R0, 0.975), NA, NA, NA)
	)
	
	cdata$coverage <- c(
		cdata$lwr[1] < 2 && 2 < cdata$upr[1],
		NA,
		NA,
		NA
	)
	
	fitlist[[i]] <- cdata
	jagslist[[i]] <- jj
	
	save("fitlist", "jagslist", file=fn)
}

save("fitlist", "jagslist", file=fn)
