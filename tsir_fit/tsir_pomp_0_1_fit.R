library(dplyr)
library(pomp)
source("../R/fitfun_pomp_sir.R")
source("../R/confint_pomp.R")

load("../data/tsir_data.rda")

argvals <- commandArgs(trailingOnly=TRUE)
batch_num <- as.numeric(argvals[1])

fn <- paste0("tsir_pomp_0_1_fit_", batch_num, ".rda")

delta_t <- 0.1
gamma <- -1/delta_t * log(1 - delta_t/1)

globals <- Csnippet(paste0("double N=100000; double Gamma=",gamma,";"))

nsim <- 10

fitlist <- vector('list', nsim)

set.seed(101)
for (i in 1:nsim) {
	print(i)
	
	j <- batch_num * 10 + i
	
	dd <- datalist[[j]][1:20,]
	
	dd <- data.frame(
		time=dd$t,
		incidence=dd$incidence
	)
	
	pomp_arg2 <- append(pomp_arg, list(data=dd, globals=globals, t0=0))
	
	pomp_model <- do.call(pomp, pomp_arg2)
	
	start <- c(beta=1.5, rho=0.7, I0=1e-4, disp=40)
	
	rwsd_arg <- list(beta=0.05, rho=0.05, I0=0.05, disp=0.05)
	
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
	
	cc <- confint_pomp(m, rwsd_arg=rwsd_arg, par=1)
	
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
