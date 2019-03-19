library(tidyr)
library(dplyr)
library(tsiR)
library(pomp)
library(ggplot2); theme_set(theme_bw())
source("../R/fitfun_pomp_tsir_lognormal_process.R")

rr <- read.csv("../data/measlesUKUS.csv")

boston <- rr %>%
	filter(loc=="BOSTON") %>%
	filter(year >= 1920, year < 1940) %>%
	rename(
		time=decimalYear,
		cases=cases,
		pop=pop,
		births=rec
	)

bs_fit <- runtsir(boston, alpha=0.975, sbar=0.035, inits.fit=TRUE)

globals <- Csnippet(paste0("double S0=", bs_fit$inits[1], "; double I0=", bs_fit$inits[2], ";"))

pomp_arg <- make_pomp_tsir_lognormal_process()

fitdata <- data.frame(
	time=1:nrow(boston),
	cases=boston$cases
)

fitdata$cases[fitdata$cases==0] <- 1

alphavec <- seq(0.9, 1, by=0.01)
sigmavec <- seq(0, 0.9, by=0.1)

liklist <- tsirlist <- simlist <- vector('list', length(alphavec))

set.seed(101)
for (i in 1:length(alphavec)) {
	print(i)
	alpha <- alphavec[i]
	
	tsir_fit <- runtsir(boston, alpha=alpha, sbar=0.035, userYhat=bs_fit$Yhat, regtype="user", nsim=1)
	
	pomp_covar <- data.frame(
		ctime=1:nrow(boston),
		B=round(boston$births),
		N=boston$pop,
		Beta=rep(tsir_fit$beta, 100)[1:nrow(boston)],
		rho=1/bs_fit$rho
	)
	
	pomp_model <- do.call(pomp, append(
		pomp_arg, list(
			tcovar="ctime",
			data=fitdata,
			t0=1,
			covar=pomp_covar,
			globals=globals
		))) 
	
	tsir_var <- var(residuals(tsir_fit$glmfit))
	
	fitlist <- lapply(sigmavec, function(x) {
		pp <- logmeanexp(replicate(10, logLik(pfilter(pomp_model, param=c(alpha=alpha, sigma_p=sqrt(x * tsir_var), sigma_o=sqrt(1-x * tsir_var)), Np=10000, tol=1e-301))), se=TRUE)
		
		data.frame(
			logLik=pp[1],
			se=pp[2],
			amount=x,
			alpha=alpha
		)
	})
	
	liklist[[i]] <- fitlist %>% bind_rows
	
	tsirlist[[i]] <- tsir_fit
	
	save("liklist", "tsirlist", file="compare_likelihood_lognormal_process.rda")
}

save("liklist", "tsirlist", file="compare_likelihood_lognormal_process.rda")
