library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_size = 12,
									 base_family = "Times"))

if (.Platform$OS.type=="windows") {
	windowsFonts(Times=windowsFont("Times"))
} 

## use Dark2 rather than Set1 because colour #6 of Set1 is yellow (ugh/too light)
scale_colour_discrete <- function(...,palette="Dark2") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Dark2") scale_fill_brewer(...,palette=palette)

load("../sim/gillespie_sim.rda")

nsim <- length(reslist)

datalist <- vector('list', nsim)

set.seed(101)
for (i in 1:nsim) {
	rr <- reslist[[i]]
	ii <- c(rr$incidence[1], diff(rr$incidence))
	tt <- rr$time
	
	dd <- as.data.frame(table(ceiling(tt[ii==1])))
	colnames(dd) <- c("time", "incidence")
	dd$time <- as.numeric(as.character(dd$time))
	
	datalist[[i]] <- dd
}

simulate_si <- function(delta=0.1, 
						beta=2, 
						IP=1, ## infectious period
						I0=10,
						N=1e5,
						tmax=20,
						type=c("det", "stoch"),
						nsim=1) {
	
	type <- match.arg(type)
	
	ind <- seq(1, tmax/delta+1, by=1/delta)
	
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
			time=tail(tvec[ind], -1),
			incidence=-diff(Svec[ind]),
			sim=i,
			type=type,
			delta=delta
		)
	}
	
	do.call("rbind", reslist)
}

delta_vec <- c(0.1, 0.2, 0.5, 1)

stoch_simlist <- lapply(delta_vec, simulate_si, nsim=250, type="stoch", N=1e5)

stoch_simdata <- stoch_simlist %>% 
	bind_rows %>%
	gather(key, value, -time, -sim, -type, -delta) %>%
	group_by(key, type, time, delta) %>%
	summarize(
		mean=mean(value),
		lwr=quantile(value, 0.025),
		upr=quantile(value, 0.975)
 	) %>% 
	mutate(
		sim="discrete",
		delta=paste0("Delta~t==", delta, "~generation")
	)

gillespie_simdata <- datalist %>%
	bind_rows(.id="sim") %>%
	filter(time <= 20) %>%
	gather(key, value, -time, -sim) %>%
	group_by(key, time) %>%
	summarize(
		mean=mean(value),
		lwr=quantile(value, 0.025),
		upr=quantile(value, 0.975)
	) %>%
	mutate(sim="gillespie")

gsim <-ggplot(stoch_simdata) +
	geom_line(data=gillespie_simdata, aes(time, mean, col=sim), lwd=1) +
	geom_ribbon(data=gillespie_simdata, aes(time, ymin=lwr, ymax=upr, fill=sim), alpha=0.3) +
	geom_line(aes(time, mean, col=sim), lwd=1) +
	geom_ribbon(aes(time, ymin=lwr, ymax=upr, fill=sim), alpha=0.3) +
	facet_wrap(~delta, labeller = label_parsed, nrow=1) +
	ylab("Incidence") +
	xlab("Generations") +
	theme(
		panel.spacing.x = grid::unit(0, "cm"),
		strip.background = element_blank(),
		legend.position = "top",
		legend.title = element_blank()
	)

ggsave("sir_delta_t_dynamics.pdf", width=6, height=4)
