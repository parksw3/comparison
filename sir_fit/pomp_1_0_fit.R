library(dplyr)
library(pomp)
source("../R/fitfun_pomp_sir.R")
source("../R/confint_pomp.R")

load("../data/gillespie_data.rda")

argvals <- commandArgs(trailingOnly=TRUE)
batch_num <- as.numeric(argvals[1])

fn <- paste0("pomp_1_0_fit_", batch_num, ".rda")

delta_t <- 1
gamma <- 1000

globals <- Csnippet(paste0("double N=100000; double Gamma=",gamma,";"))
pomp_arg$rprocess <- pomp::euler.sim(rprocess, delta.t=1)

nsim <- 10

fitlist <- vector('list', nsim)

set.seed(101)
for (i in 1:nsim) {
	print(i)
	
	j <- batch_num * 10 + i
	
	dd <- datalist[[j]][1:20,]
	
	pomp_arg2 <- append(pomp_arg, list(data=dd, globals=globals, t0=0))
	
	pomp_model <- do.call(pomp, pomp_arg2)
	
	start <- c(beta=2, rho=0.7, I0=10e-5, disp=40)
	
	rwsd_arg <- list(beta=0.01, rho=0.01, I0=0.01, disp=0.01)
	
	m <- mif2(
		pomp_model,
		Nmif=100,
		start=start,
		Np=1000,
		cooling.fraction.50=0.95,
		rw.sd=do.call(rw.sd, rwsd_arg),
		transform=TRUE) %>%
		continue(Nmif=100, cooling.fraction=0.8) %>%
		continue(Nmif=100, cooling.fraction=0.6)
	
	cc <- confint_pomp2(m, rwsd_arg=rwsd_arg, par=1)
	
	cdata <- data.frame(
		param=c("beta", "rprob", "I0", "size"),
		mean=coef(m),
		lwr=c(cc$lwr, NA, NA, NA),
		upr=c(cc$upr, NA, NA, NA)
	)
	
	cdata$coverage <- c(
		cdata$lwr[1] < 2 && 2 < cdata$upr[1],
		NA,
		NA,
		NA
	)
	
	fitlist[[i]] <- cdata
	
	save("fitlist", file=fn)
}

save("fitlist", file=fn)
