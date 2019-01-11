library(subplex)
library(dplyr)
library(tsiR)
library(pomp)

source("../R/fitfun_pomp_exp_fixed_det.R")
source("../R/basis.R")

measles_data <- read.csv("../data/measlesUKUS.csv")

measles_UK <- measles_data %>% 
	filter(country=="UK") %>%
	mutate(cases=ifelse(is.na(cases), 0, cases))

measles_list <- measles_UK %>%
	split(as.character(.$loc))

ld <- measles_list$LONDON

ld$time <- ld$index <- 1:nrow(ld)

fit_tsir <- runtsir(
	data=data.frame(
		time=ld$time,
		cases=ld$cases,
		pop=ld$pop,
		births=ld$rec
	),
	alpha=0.97,
	sbar=0.035,
	inits.fit = FALSE,
	method="det",
	nsim=1
)

ld_obs <- select(ld, time, cases)
ld_covar <- select(ld, pop, rec, index, biweek)

ld_covar <- append_basis(ld_covar)

ld_pomp_arg <- append(pomp_arg, list(data=ld_obs, covar=ld_covar, t0=min(ld$time)))

ld_pomp <- do.call(pomp, ld_pomp_arg)

trans.param <- estimate_basis(fit_tsir$contact$beta * mean(ld_covar$pop))
names(trans.param) <- paste0("b", 1:length(trans.param))

params <- c(trans.param, alpha=unname(fit_tsir$alpha), m=10, 
		   S0=fit_tsir$inits[1],
		   I0=fit_tsir$inits[2],
		   rho=mean(1/fit_tsir$rho),
		   disp=5)

london_pomp_exp_det_fixed_frequency_list <- vector('list', 100)

set.seed(101)
for (i in 1:100) {
	print(i)
	start <- params
	start <- rlnorm(length(params), meanlog=log(params), sdlog=0.01)
	names(start) <- names(params)
	
	m <- traj.match(ld_pomp, start=start, est=names(params),
					transform=TRUE,
					maxit=1e5)
	
	ll <- logLik(m)
	
	london_pomp_exp_det_fixed_frequency_list[[i]] <- list(mle=m, ll=ll)
}

save("trans.param", "london_pomp_exp_det_fixed_frequency_list", "ld_pomp_arg", 
	 file="london_pomp_exp_det_fixed_frequency.rda")
