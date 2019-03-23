library(dplyr)
library(mgcv)
library(pomp)
source("../R/fitfun_pomp_renewal_basis.R")
source("../R/basis.R")

load("../data/gillespie_sinusoidal_data.rda")

argvals <- commandArgs(trailingOnly=TRUE)
batch_num <- as.numeric(argvals[1])

fn <- paste0("sinusoidal_pomp_renewal_fit_", batch_num, ".rda")

globals <- Csnippet(paste0("double N0=5000000; double Gmean=1; double mu=0.0007692308;"))

k <- seq(0, 26, by=2)

BX <- mgcv::cSplineDes(seq(0, 260, by=0.1)%%26,k)
basis_covar <- as.data.frame(BX)
colnames(basis_covar) <- paste0("B", 1:ncol(basis_covar))
basis_covar$time <- seq(0, 260, by=0.1)

bb <- estimate_basis(500/26 *(1 + 0.15 * cos(1:26 * 2 * pi/26)))

nsim <- 10

fitlist <- vector('list', nsim)
translist <- vector('list', nsim)
reslist <- vector('list', nsim)

set.seed(101)
for (i in 1:nsim) {
	print(i)
	
	j <- batch_num * 10 + i
	
	dd <- datalist[[j]][1:(9*26),]
	
	pomp_arg <- make_pomp_renewal_basis()
	
	pomp_arg2 <- append(pomp_arg, list(data=dd[,c("time", "incidence")], globals=globals, t0=0, covar=basis_covar, tcovar="time"))
	
	pomp_model <- do.call(pomp, pomp_arg2)
	
	start <- c(bb, rho=0.7, S0=0.05, I0=1e-4, disp=5, Gvar=1)
	
	rwsd_bb <- as.list(rep(0.01, length(bb)))
	names(rwsd_bb) <- names(bb)
	
	rwsd_arg <- c(rwsd_bb, rho=0.01, S0=0.01, I0=0.01, disp=0.01, Gvar=0.01)
	
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
	
	cdata <- data.frame(
		param=names(start),
		mean=coef(m)
	)
	
	fitlist[[i]] <- cdata
	translist[[i]] <- data.frame(
		time=seq(1, 26, by=0.1),
		beta=(BX %*% coef(m)[1:13])[1:251]
	)
	
	reslist[[i]] <- m
	
	save("fitlist", "translist", "reslist", file=fn)
}

save("fitlist", "translist", "reslist", file=fn)
