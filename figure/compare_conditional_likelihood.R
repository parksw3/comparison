library(dplyr)
library(ggplot2); theme_set(theme_bw())

likelihood_list <- list()

load("../sim/compare_likelihood_lognormal.rda")

likelihood_list$conditional <- liklist

load("../sim/compare_likelihood_lognormal_process.rda")

likelihood_list$unconditional <- liklist

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

bs_fit <- runtsir(boston, alpha=0.975, sbar=0.051, inits.fit=TRUE)

fitdata <- data.frame(
	time=1:nrow(boston),
	cases=boston$cases
)

fitdata$cases[fitdata$cases==0] <- 1

alphavec <- seq(0.9, 1, by=0.005)

likelihood_tsir <- lapply(alphavec, function(x) {
	alpha <- x
	
	tsir_fit <- runtsir(boston, alpha=alpha, sbar=0.051, userYhat=bs_fit$Yhat, regtype="user", nsim=1)
	
	data.frame(
		alpha=x,
		logLik=logLik(tsir_fit$glmfit),
		amount=1,
		se=NA,
		type="conditional"
	)
}) %>%
	bind_rows

likelihood_data <- likelihood_list %>%
	lapply(bind_rows) %>%
	bind_rows(.id="type")

tsir_data <- data.frame(
	alpha=0.975, 
	amount=1,
	type="conditional"
)

global_max <- likelihood_data %>%
	group_by(type) %>%
	filter(logLik==max(logLik))

ggplot(filter(likelihood_data, amount > 0.1)) +
	geom_raster(aes(alpha, amount, fill=logLik)) +
	geom_point(data=tsir_data, aes(alpha, amount), size=4, shape=1) +
	geom_point(data=global_max, aes(alpha, amount), size=4, shape=2) +
	scale_x_continuous(limits=c(0.9, 1), expand=c(0,0)) +
	scale_fill_gradientn(colours=terrain.colors(10)) +
	facet_wrap(~type)

ggplot(filter(likelihood_data, amount==0, alpha < 1)) +
	geom_line(aes(alpha, logLik)) +
	facet_wrap(~type, scale="free")
