library(dplyr)
library(pomp)
source("../R/fitfun_pomp_renewal.R")
source("../R/confint_pomp.R")

load("../data/gillespie_data.rda")

argvals <- commandArgs(trailingOnly=TRUE)
batch_num <- as.numeric(argvals[1])

fn <- paste0("pomp_renewal_fit_", batch_num, ".rda")

globals <- Csnippet(paste0("double N0=100000; double Gmean=1;"))

nsim <- 10

fitlist <- vector('list', nsim)

set.seed(101)
for (i in 1:nsim) {
	print(i)
	
	j <- batch_num * 10 + i
	
	dd <- datalist[[j]][1:20,]
	
	pomp_arg <- make_pomp_renewal()
	
	pomp_arg2 <- append(pomp_arg, list(data=dd, globals=globals, t0=0))
	
	pomp_model <- do.call(pomp, pomp_arg2)
	
	start <- c(R0=2, rho=0.7, I0=1e-4, disp=10, Gvar=1)
	
	rwsd_arg <- list(R0=0.05, rho=0.05, I0=0.05, disp=0.05, Gvar=0.01)
	
	m <- mif2(
		pomp_model,
		Nmif=50,
		start=start,
		Np=1000,
		cooling.fraction.50=0.95,
		rw.sd=do.call(rw.sd, rwsd_arg),
		transform=TRUE) %>%
		continue(Nmif=50, cooling.fraction=0.8) %>%
		continue(Nmif=50, cooling.fraction=0.6) %>%
		continue(Nmif=50, cooling.fraction=0.2) %>%
		continue(Nmif=50, cooling.fraction=0.1)
	
	cc <- confint_pomp(m, rwsd_arg=rwsd_arg, par=1, trace=TRUE)
	
	cdata <- data.frame(
		param=names(start),
		mean=coef(m),
		lwr=c(cc$lwr, NA, NA, NA, NA),
		upr=c(cc$upr, NA, NA, NA, NA)
	)
	
	cdata$coverage <- c(
		cdata$lwr[1] < 2 && 2 < cdata$upr[1],
		NA,
		NA,
		NA,
		NA
	)
	
	fitlist[[i]] <- cdata
	
	save("fitlist", file=fn)
}

save("fitlist", file=fn)
