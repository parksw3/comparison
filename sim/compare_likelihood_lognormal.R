library(tidyr)
library(dplyr)
library(tsiR)
library(pomp)
library(ggplot2); theme_set(theme_bw())
source("../R/fitfun_pomp_tsir_lognormal.R")

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

globals <- Csnippet(paste0("double I0=", bs_fit$inits[2], ";"))

pomp_arg <- make_pomp_tsir_lognormal()

fitdata <- data.frame(
	time=1:nrow(boston),
	cases=boston$cases
)

fitdata$cases[fitdata$cases==0] <- 1

alphavec <- seq(0.9, 0.98, by=0.002)

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
		rho=1/bs_fit$rho,
		S=0.051 * mean(rr$pop) + bs_fit$Z
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
	
	p0 <- logmeanexp(replicate(10, logLik(pfilter(pomp_model, param=c(alpha=alpha, sigma_p=sqrt(0 * tsir_var), sigma_o=sqrt(1 * tsir_var)), Np=10000, tol=1e-301))), se=TRUE)
	p1 <- logmeanexp(replicate(10, logLik(pfilter(pomp_model, param=c(alpha=alpha, sigma_p=sqrt(0.1 * tsir_var), sigma_o=sqrt(0.9 * tsir_var)), Np=10000, tol=1e-301))), se=TRUE)
	p2 <- logmeanexp(replicate(10, logLik(pfilter(pomp_model, param=c(alpha=alpha, sigma_p=sqrt(0.5 * tsir_var), sigma_o=sqrt(0.5 * tsir_var)), Np=10000, tol=1e-301))), se=TRUE)
	p3 <- logmeanexp(replicate(10, logLik(pfilter(pomp_model, param=c(alpha=alpha, sigma_p=sqrt(0.8 * tsir_var), sigma_o=sqrt(0.2 * tsir_var)), Np=10000, tol=1e-301))), se=TRUE)
	p4 <- logmeanexp(replicate(10, logLik(pfilter(pomp_model, param=c(alpha=alpha, sigma_p=sqrt(0.9 * tsir_var), sigma_o=sqrt(0.1 * tsir_var)), Np=10000, tol=1e-301))), se=TRUE)
	p5 <- logmeanexp(replicate(10, logLik(pfilter(pomp_model, param=c(alpha=alpha, sigma_p=sqrt(0.95 * tsir_var), sigma_o=sqrt(0.05 * tsir_var)), Np=10000, tol=1e-301))), se=TRUE)
	
	liklist[[i]] <- data.frame(
		logLik=c(p0[1], p1[1], p2[1], p3[1], p4[1], p5[1]),
		se=c(p0[2], p1[2], p2[2], p3[2], p4[2], p5[2]),
		amount=c(0, 0.1, 0.5, 0.8, 0.9, 0.95),
		alpha=alpha
	)
	
	tsirlist[[i]] <- tsir_fit
}

save("liklist", "tsirlist", file="compare_likelihood_lognormal.rda")
