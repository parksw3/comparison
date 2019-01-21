library(deSolve)
library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw())

sir <- function(t, state, parameters) {
	with(as.list(c(state, parameters)), {
		I <- exp(logI)
		dS <- - beta * S * I/N
		dlogI <- beta * S/N - gamma
		dR <- gamma * I
		list(c(dS, dlogI, dR))
	})
}

par <- c(beta=2, gamma=1, N=1)
yini <- c(S=1-1e-4, logI=log(1e-4), R=0)

out <- as.data.frame(rk(yini, seq(0, 30, by=0.1), sir, par))

simulate_si <- function(delta=0.1, 
						beta=2, 
						IP=1, ## infectious period
						I0=10,
						N=1e5,
						tmax=30,
						type=c("det", "stoch"),
						nsim=1) {
	
	type <- match.arg(type)
	
	if (type=="det") {
		simfun <- function(x, y) x * y
	} else {
		simfun <- function(x, y) rbinom(1, size=x, prob=y)
	}
	
	if (delta==IP) {
		rec <- 1
	} else if (delta < IP) {
		gamma <- -1/delta * log(1 - delta/IP)
		rec <- 1 - exp(-gamma*delta)
	} else {
		stop("delta must be smaller than IP")
	}
	
	reslist <- vector('list', nsim)
	
	for (i in 1:nsim) {
		Svec <- S <- N-I0
		Ivec <- I <- I0
		tvec <- t <- 0
		
		while (t < tmax) {
			infection <- simfun(S, (1 - exp(-beta * I * delta/N)))
			recovery <- simfun(I, rec)
			S <- S - infection
			I <- infection + I - recovery
			Svec <- c(Svec, S)
			Ivec <- c(Ivec, I)
			t <- t + delta
			tvec <- c(tvec, t)
		}
		
		reslist[[i]] <- data.frame(
			time=tvec,
			S=Svec,
			I=Ivec,
			sim=i,
			type=type,
			delta=delta
		)
	}
	
	do.call("rbind", reslist)
}

delta_vec <- c(0.1, 0.2, 0.5, 1)

det_simlist <- lapply(delta_vec, simulate_si)

det_simdata <- det_simlist %>%
	bind_rows %>%
	gather(key, value, -time, -sim, -type, -delta) %>%
	group_by(key, type, time, delta) %>%
	summarize(mean=mean(value)/1e5) %>%
	mutate(
		delta=factor(delta, labels=paste0("Delta~t==", delta))
	)

ggplot(det_simdata %>% filter(key=="I")) +
	geom_line(data=out, aes(time, exp(logI)), col="gray", lwd=2) +
	geom_line(aes(time, mean, group=delta), lwd=1) +
	ylab("Proportion infected") +
	xlab("Time (generations)") +
	facet_grid(~delta, labeller=label_parsed) +
	theme(
		
	)

stoch_simlist1 <- lapply(delta_vec, simulate_si, nsim=500, type="stoch", N=1000)
stoch_simlist2 <- lapply(delta_vec, simulate_si, nsim=500, type="stoch", N=10000)
stoch_simlist3 <- lapply(delta_vec, simulate_si, nsim=500, type="stoch", N=100000)

stoch_simdata <- stoch_simlist %>% 
	bind_rows %>%
	gather(key, value, -time, -sim, -type, -delta) %>%
	group_by(key, type, time, delta) %>%
	summarize(
		mean=mean(value),
		lwr=quantile(value, 0.025),
		upr=quantile(value, 0.975)
	)

ggplot(stoch_simdata %>% filter(key=="I")) +
	geom_line(aes(time, mean)) +
	geom_ribbon(aes(time, ymin=lwr, ymax=upr), alpha=0.5) +
	facet_wrap(~delta)

stoch_summary <- stoch_simlist %>%
	bind_rows %>%
	group_by(type, sim, delta) %>%
	summarize(
		max=max(I),
		maxtime=head(time[I==max(I)],1),
		zerotime=head(time[I==0], 1),
		final=1-tail(S, 1)/(head(S,1)+10)
	)

ggplot(stoch_summary) + geom_boxplot(aes(as.character(delta), maxtime))
