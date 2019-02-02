library(dplyr)
library(pomp)
source("../R/fitfun_pomp_sir.R")
source("../R/confint_pomp.R")

load("../data/gillespie_data.rda")

argvals <- commandArgs(trailingOnly=TRUE)
batch_num <- as.numeric(argvals[1])

fn <- paste0("pomp_0_5_fit_", batch_num, ".rda")

delta_t <- 0.5
gamma <- -1/delta_t * log(1 - delta_t/1)

globals <- Csnippet(paste0("double N=100000; double Gamma=",gamma,";"))
pomp_arg$rprocess <- pomp::euler.sim(rprocess, delta.t=5/10)

nsim <- 10

fitlist <- vector('list', nsim)

set.seed(101)
for (i in 1:nsim) {
	print(i)
	
	j <- batch_num * 10 + i
	
	dd <- datalist[[j]][1:20,]
	
	pomp_arg2 <- append(pomp_arg, list(data=dd, globals=globals, t0=0))
	
	pomp_model <- do.call(pomp, pomp_arg2)
	
	start <- c(beta=2, rho=0.7, I0=1e-4, disp=40)
	
	rwsd_arg <- list(beta=0.1, rho=0.1, I0=0.1, disp=0.1)
	
	m <- mif2(
		pomp_model,
		Nmif=300,
		start=start,
		Np=1000,
		cooling.fraction.50=0.95,
		rw.sd=do.call(rw.sd, rwsd_arg),
		transform=TRUE)
	
	cc <- confint_pomp(m, rwsd_arg=rwsd_arg)
	
	cdata <- data.frame(
		param=c("beta", "rprob", "I0", "size"),
		mean=coef(m),
		lwr=c(cc$lwr, NA, NA),
		upr=c(cc$upr, NA, NA)
	)
	
	cdata$coverage <- c(
		cdata$lwr[1] < 2 && 2 < cdata$upr[1],
		cdata$lwr[2] < 0.7 && 0.7 < cdata$upr[2],
		NA,
		NA
	)
	
	fitlist[[i]] <- cdata
	
	save("fitlist", file=fn)
}

save("fitlist", file=fn)
