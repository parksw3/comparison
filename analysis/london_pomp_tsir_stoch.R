library(dplyr)
library(tsiR)
library(pomp)

source("../R/fitfun_pomp_tsir.R")

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
	inits.fit = TRUE,
	method="negbin"
)

ld_obs <- select(ld, time, cases)
ld_covar <- select(ld, pop, rec, index, biweek)

ld_pomp_arg <- apapend(pomp_arg, list(data=ld_obs, covar=ld_covar, t0=min(ld$time)))

ld_pomp <- do.call(pomp, ld_pomp_arg)

trans.param <- fit_tsir$contact$beta * mean(ld$pop)
names(trans.param) <- paste0("b", 1:26)

params <- c(trans.param, alpha=unname(fit_tsir$alpha), m=10, 
		   S0=fit_tsir$inits[1],
		   I0=fit_tsir$inits[2],
		   rho=mean(1/fit_tsir$rho),
		   disp=20)

ld_rwsd_arg <- as.list(abs(log(params)) * 0.001)

mflist_tsir <- vector('list', 10)

set.seed(101)
for (i in 1:10) {
	print(i)
	start <- params
	start <- rlnorm(length(params), meanlog=log(params), sdlog=0.05)
	names(start) <- names(params)
	
	m <- mif2(
		ld_pomp,
		Nmif=30,
		start=start,
		Np=2000,
		cooling.fraction.50=0.95,
		rw.sd=do.call(rw.sd, ld_rwsd_arg),
		transform=TRUE) %>%
		continue(Nmif=30, cooling.fraction=0.8) %>%
		continue(Nmif=30, cooling.fraction=0.6) %>%
		continue(Nmif=30, cooling.fraction=0.2)
	
	ll <- replicate(n=10,logLik(pfilter(m,Np=10000)))
	
	mflist_tsir[[i]] <- list(mif=m, ll=logmeanexp(ll,se=TRUE))
}

save("trans.param", "mflist_tsir", "ld_pomp_arg", file="london_pomp_tsir.rda")
