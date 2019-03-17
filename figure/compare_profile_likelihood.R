library(tidyr)
library(dplyr)
library(deSolve)
library(ggplot2); theme_set(theme_bw(base_size=12, base_family = "Times"))
library(gridExtra)
library(emdbook)
library(bbmle)

if (.Platform$OS.type=="windows") {
	windowsFonts(Times=windowsFont("Times"))
} 

scale_colour_discrete <- function(...,palette="Dark2") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Dark2") scale_fill_brewer(...,palette=palette)

rprob <- 0.7
theta <- 10

gfun <- function(t, y, pars) {
	with(as.list(c(y, pars)), {
		dS <- - beta * S * I/N
		dI <- beta * S * I/N - gamma * I
		dR <- gamma * I
		
		list(c(dS, dI, dR))
	})
}

gfun_sinusoidal <- function(t, y, pars) {
	with(as.list(c(y, pars)), {
		beta <- b0 * (1 + b1 * cos(2 * pi * t/ 26))
		dS <- mu * (N - S) - beta * S * I/N
		dI <- beta * S * I/N - gamma * I - mu * I
		dR <- gamma * I - mu * R
		dCI <- beta * S * I/N
		dCR <- gamma * I
		list(c(dS, dI, dR, dCI, dCR))
	})
}

objfun <- function(param,
				   data,
				   N=1e5,
				   type=c("incidence", "mortality", "prevalence")) {
	type <- match.arg(type)
	
	init <- ifelse(type=="prevalence", 1, 0)
	
	out <- with(as.list(param), {
		par <- c(beta=exp(log.beta), gamma=exp(log.gamma), N=N)
		I0 <- plogis(logit.I0)
		y <- c(S=N*(1-I0), I=N*I0, R=0)
		
		times <- c(init:20)
		
		out <- as.data.frame(rk4(y, times, gfun, par))
		
		out
	})
	
	rprob <- plogis(param[["logit.rprob"]])
	size <- exp(param[["log.size"]])
	
	yhat <- switch(type,
		incidence=-diff(out$S),
		prevalence=out$I,
		mortality=diff(out$R)
	)
	
	nll <- -sum(dnbinom(data[["prevalence"]], mu=yhat*rprob, size=size, log=TRUE))
	
	nll
}

parnames(objfun) <- c("log.beta", "log.gamma", "logit.rprob", "logit.I0", "log.size")

yini <- c(S=1e5-10, I=10, R=0)
pars <- c(beta=2, gamma=1, N=1e5)
tvec <- 0:20

yini2 <- c(S=0.05*5e6, I=1e-4*5e6, R=0, CI=0, CR=0)
pars2 <- c(mu=1/(50*26), b0=500/26, b1=0.15, gamma=1, N=5e6)
tvec2 <- seq(0, 42, by=1)

out <- as.data.frame(ode(yini, tvec, gfun, pars))
out2 <- as.data.frame(ode(yini2, tvec2, gfun_sinusoidal, pars2))

truedata <- data.frame(
	time=1:20,
	incidence=-diff(out$S),
	prevalence=tail(out$I, -1),
	mortality=diff(out$R)
)

truedata_sinusoidal <- data.frame(
	time=1:42,
	incidence=diff(out2$CI),
	prevalence=tail(out2$I, -1),
	mortality=diff(out2$CR)
)

set.seed(101)
obsdata <- data.frame(
	time=1:20,
	incidence=rbetabinom(20, size=round(-diff(out$S)), prob=rprob, theta=theta),
	prevalence=rbetabinom(20, size=round(tail(out$I, -1)), prob=rprob, theta=theta),
	mortality=rbetabinom(20, size=round(diff(out$R)), prob=rprob, theta=theta)
)

start <- c(log.beta=log(2), log.gamma=log(1), logit.rprob=qlogis(0.7), logit.I0=qlogis(1e-4), log.size=1)

m_inc <- mle2(objfun, start=start, data=list(data=obsdata, type="incidence"))
m_prev <- mle2(objfun, start=start, data=list(data=obsdata, type="prevalence"))
m_mort <- mle2(objfun, start=start, data=list(data=obsdata, type="mortality"))

pp_inc <- profile(m_inc, trace=TRUE, which=1:3, alpha=0.05, continuation="naive")
pp_prev <- profile(m_prev, trace=TRUE, which=1:3, alpha=0.05, continuation="naive")
pp_mort <- profile(m_mort, trace=TRUE, which=1:3, alpha=0.05, continuation="naive")

profile_list <- list(incidence=pp_inc, prevalence=pp_prev, mortality=pp_mort)

profile_data <- profile_list %>%
	lapply(function(pp) {
		prof <- pp@profile
		
		reslist <- list()
		for (i in 1:length(prof)) {
			reslist[[i]] <- data.frame(
				par.vals=prof[[i]]$par.vals[,names(prof)[[i]]],
				z=prof[[i]]$z,
				par.names=names(prof)[[i]]
			)
		}
		
		reslist %>%
			bind_rows
	}) %>%
	bind_rows(.id="key")

truedata2 <- truedata %>%
	gather(key, value, -time) %>%
	mutate(key = sub("(.)", "\\U\\1", key, perl=TRUE))

truedata_sinusoidal2 <- truedata_sinusoidal %>%
	gather(key, value, -time) %>%
	mutate(key = sub("(.)", "\\U\\1", key, perl=TRUE))

g1 <- ggplot(truedata2) +
	geom_line(aes(time, value, col=key, lty=key), lwd=1.1) +
	scale_x_continuous("Time", expand=c(0, 0)) +
	scale_y_continuous("True cases") +
	theme(
		panel.grid = element_blank(),
		legend.title = element_blank(),
		legend.position=c(0.82, 0.78)
	)
	
g2 <- (g1 %+% truedata_sinusoidal2) +
	#scale_y_continuous("Observed cases") +
	theme(
		legend.position="none"
	)

g3 <- ggplot(filter(profile_data, par.names=="log.beta")) +
	geom_vline(xintercept=2, lty=2, lwd=0.5, col='grey') +
	geom_point(aes(exp(par.vals), abs(z)^2, col=key, shape=key)) +
	geom_smooth(aes(exp(par.vals), abs(z)^2, col=key, lty=key), se=FALSE, n=300, method="lm", formula=y~poly(x, 5)) +
	scale_y_continuous("Deviance difference", limits=c(-2.1, 6), expand=c(0,-2)) +
	scale_x_continuous("Transmission rate", limits=c(1.5, 2.2)) +
	theme(
		legend.position="none",
		panel.grid=element_blank()
	)

g4 <- ggplot(filter(profile_data, par.names=="log.gamma")) +
	geom_vline(xintercept=1, lty=2, lwd=0.5, col='grey') +
	geom_point(aes(1/exp(par.vals), abs(z)^2, col=key, shape=key)) +
	geom_smooth(aes(1/exp(par.vals), abs(z)^2, col=key, lty=key), se=FALSE, n=300, method="lm", formula=y~poly(x, 5)) +
	scale_y_continuous("Deviance difference", limits=c(-3.1, 7), expand=c(0,-3)) +
	scale_x_continuous("Mean generation time", limits=c(0.8, 2)) +
	theme(
		legend.position="none",
		panel.grid=element_blank()
	)

g5 <- ggplot(filter(profile_data, par.names=="logit.rprob")) +
	geom_vline(xintercept=0.7, lty=2, lwd=0.5, col='grey') +
	geom_point(aes(plogis(par.vals), abs(z)^2, col=key, shape=key)) +
	geom_smooth(aes(plogis(par.vals), abs(z)^2, col=key, lty=key), se=FALSE, n=300, method="lm", formula=y~poly(x, 5)) +
	scale_y_continuous("Deviance difference", limits=c(-6.1, 10), expand=c(0,-6)) +
	scale_x_continuous("Reporting rate", limits=c(0.4, 1)) +
	theme(
		legend.position="none",
		panel.grid=element_blank()
	)

gfinal <- arrangeGrob(
	arrangeGrob(g1 + ggtitle("A"), g2 + ggtitle("B"), nrow=1),
	arrangeGrob(g3 + ggtitle("C"), g4 + ggtitle("D"), g5 + ggtitle("E"), nrow=1),
	nrow=2
) 

ggsave("compare_profile_likelihood.pdf", gfinal, width=8, height=6)
