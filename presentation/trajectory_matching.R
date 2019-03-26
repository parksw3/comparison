library(deSolve)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)

measles_data <- read.csv("../data/measlesUKUS.csv")

ld <- measles_data %>%
	filter(loc=="LONDON", year > 1955, year < 1960)

betafun_sinusoidal <- function(t, parameters) {
	with(as.list(c(parameters)), {
		b0 * (1 + b1 * cos(2 * pi * t - b2))
	})
}

betafun_term <- function(t, parameters) {
	t <- t %% 1
	
	with(as.list(c(parameters)), {
		if ((t>=7/365&&t<=100/365) || (t>=115/365&&t<=199/365) || (t>=252/365&&t<=300/365) || (t>=308/365&&t<=356/365)) {
			R0 * (1.0+amplitude*0.2411/0.7589)
		} else {
			R0 * (1.0-amplitude)
		}
	})
}

sir <- function(t, state, parameters, betafun) {
	
	beta <- betafun(t, parameters)
	
	with(as.list(c(state, parameters)), {
		I <- exp(logI)
		dS <- mu * (N - S) - beta * S * I/N
		dlogI <- beta * S/N - gamma - mu
		dR <- gamma * I
		dC <- beta * S * I/N
		list(c(dS, dlogI, dR, dC))
	})
}


params1 <- c(mu=1/(50), b0=500, b1=0.15, gamma=26, N=mean(ld$pop), b2=0)
yini1 <- c(S=0.05*mean(ld$pop), logI=log(1e-4*mean(ld$pop)), R=0, C=0)
tvec <- seq(min(ld$decimalYear), 1+max(ld$decimalYear), by=0.001)

out1 <- as.data.frame(ode(yini1, tvec, sir, params1, betafun=betafun_sinusoidal))

params2 <- c(mu=1/(50), R0=500, amplitude=0.5, gamma=26, N=mean(ld$pop))
out2 <- as.data.frame(ode(yini1, tvec, sir, params2, betafun=betafun_term))

params3 <- c(mu=1/50, b0=500, b1=0, gamma=26, N=mean(ld$pop), b2=0)
out3 <- as.data.frame(ode(yini1, tvec, sir, params3, betafun=betafun_sinusoidal))

transdata1 <- data.frame(
	time=seq(0, 1, by=0.001),
	trans=500
)

transdata2 <- data.frame(
	time=seq(0, 1, by=0.001),
	trans=500 * (1 + 0.15 * cos(2 * pi * seq(0, 1, by=0.001)))
)

transdata3 <- data.frame(
	time=seq(0, 1, by=0.001),
	trans=sapply(seq(0, 1, by=0.001), betafun_term, params2)
)

g1 <- ggplot(transdata1) +
	geom_line(aes(time, trans)) +
	scale_x_continuous("Time (years)", expand=c(0, 0), breaks=c(0, 0.2, 0.4, 0.6, 0.8)) +
	scale_y_continuous("Transmission rate") +
	theme(
		panel.grid = element_blank()
	)

g2 <- g1 %+% transdata2

g3 <- g1 %+% transdata3

g4 <- ggplot(out3) +
	geom_point(data=ld, aes(decimalYear, cases), shape=1) +
	geom_line(aes(time, exp(logI) * 0.3)) +
	scale_x_continuous("Time (years)", expand=c(0, 0), limits=c(1956, 1960), breaks=1956:1959) +
	scale_y_continuous("Cases") +
	theme(
		panel.grid = element_blank()
	)

g5 <- g4 %+% out1

g6 <- g5 %+% out2

gtot1 <- arrangeGrob(g1, g4, nrow=1)
gtot2 <- arrangeGrob(g2, g5, nrow=1)
gtot3 <- arrangeGrob(g3, g6, nrow=1)

ggsave("traj1.pdf", gtot1, width=6, height=3)
ggsave("traj2.pdf", gtot2, width=6, height=3)
ggsave("traj3.pdf", gtot3, width=6, height=3)
