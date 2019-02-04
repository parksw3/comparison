library(dplyr)
library(pomp)
source("../R/fitfun_pomp_sinusoidal.R")

load("../data/gillespie_sinusoidal_data.rda")

argvals <- commandArgs(trailingOnly=TRUE)
batch_num <- as.numeric(argvals[1])

fn <- paste0("sinusoidal_pomp_0_1_fit_", batch_num, ".rda")

delta_t <- 0.1
gamma <- -1/delta_t * log(1 - delta_t/1)

globals <- Csnippet(paste0("double N=5000000; double Gamma=",gamma,"; double mu=0.0007692308; double pi=3.14159265;"))

nsim <- 10

fitlist <- vector('list', nsim)
translist <- vector('list', nsim)
reslist <- vector('list', nsim)

set.seed(101)
for (i in 1:nsim) {
	print(i)
	
	j <- batch_num * 10 + i
	
	dd <- datalist[[j]]
	
	pomp_arg2 <- append(pomp_arg, list(data=dd[,c("time", "incidence")], globals=globals, t0=0))
	
	pomp_model <- do.call(pomp, pomp_arg2)
	
	start <- c(b0=500/26, b1=0.15, rho=0.7, S0=0.05, I0=1e-4, disp=5)
	
	rwsd_arg <- list(b0=0.01, b1=0.01, rho=0.01, S0=0.01, I0=0.01, disp=0.01)
	
	m <- mif2(
		pomp_model,
		Nmif=100,
		start=start,
		Np=2000,
		cooling.fraction.50=0.95,
		rw.sd=do.call(rw.sd, rwsd_arg),
		transform=TRUE) %>%
		continue(Nmif=100, cooling.fraction=0.8) %>%
		continue(Nmif=100, cooling.fraction=0.6) %>%
		continue(Nmif=100, cooling.fraction=0.2) %>%
		continue(Nmif=100, cooling.fraction=0.1)
	
	cdata <- data.frame(
		param=c("R0", "rprob"),
		mean=coef(m)[[1]], coef(m)[[3]]
	)

	fitlist[[i]] <- cdata
	translist[[i]] <- data.frame(
		time=seq(1, 26, by=0.01),
		beta=coef(m)[[1]] * (1 + coef(m)[[2]] * cos(seq(1, 26, by=0.01) * 2 * pi/26))
	)
	reslist[[i]] <- m
	
	save("fitlist", "translist", "reslist", file=fn)
}

save("fitlist", "translist", "reslist", file=fn)
