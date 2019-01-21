library(emdbook)
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

probfun <- function(x, y, tvec=1:30) {
	x <- unname(x); y <- unname(y)
	par <- c(beta=x, gamma=1, N=1)
	yini <- c(S=y-1e-5, logI=log(1e-5), R=1-y)
	
	out <- as.data.frame(rk(yini, tvec, sir, par, hmax=1e-3))
	
	dd <- data.frame(
		time=head(tvec, -1),
		true=1-tail(out$S,-1)/head(out$S,-1),
		exponential=1-exp(-x*head(exp(out$logI), -1)),
		linear=x*head(exp(out$logI), -1),
		R0=as.character(x),
		sus=as.character(y)
	)
	dd
}

R0_vec <- c(1.5, 2, 10, 30)
sus_vec <- c(1, 0.5, 0.1, 0.05)

problist <- apply2d(probfun, x=R0_vec, y=sus_vec, use_plyr = FALSE)

probdata <- problist %>%
	bind_rows %>%
	gather(key, value, -R0, -sus, -time) %>%
	mutate(value=ifelse(value > 1, 1, value),
		   R0=factor(R0, levels=c("1.5", "2", "10", "30"), labels=c("R0=1.5", "R0=2", "R0=10", "R0=30")),
		   sus=factor(sus, levels=c("1", "0.5", "0.1", "0.05"), labels=paste0(c("100%", "50%", "10%", "5%"))))

prob_true <- probdata %>% 
	filter(key=="true")

prob_discrete <- probdata %>%
	filter(key!="true")

gtot <- ggplot(prob_discrete) +
	geom_line(aes(time, value, col=key, lty=key), lwd=1.5) +
	geom_point(data=prob_true, aes(time, value), shape=1) +
	ylab("Probability of infection per generation") +
	xlab("Time (generations)") +
	facet_grid(sus~R0) +
	theme(
		panel.spacing = grid::unit(0, "cm"),
		legend.title = element_blank(),
		legend.position = c(0.085, 0.95)
	)

ggsave("sir_probability_of_infection.pdf", gtot, width=8, height=8)
