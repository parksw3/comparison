library(tsiR)
library(tidyr)
library(dplyr)
library(deSolve)
library(gridExtra)
library(ggplot2); theme_set(theme_bw(base_size = 12,
									 base_family = "Times"))

if (.Platform$OS.type=="windows") {
	windowsFonts(Times=windowsFont("Times"))
} 

## use Dark2 rather than Set1 because colour #6 of Set1 is yellow (ugh/too light)
scale_colour_discrete <- function(...,palette="Dark2") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Dark2") scale_fill_brewer(...,palette=palette)

reconfun <- function(data, epsilon=1e-7) {
	fitdata <- data.frame(
		y=cumsum(data$births),
		x=cumsum(data$cases)
	)
	
	ff <- loess(y~x, data=fitdata)
	
	data.frame(
		Z=ff$residuals,
		rho=1/(coef(ff)[2] + 2 * coef(ff)[3] * 1:260)
	)
}

sir <- function(t, state, parameters) {
	with(as.list(c(state, parameters)), {
		beta <- b0 * (1 + b1 * cos(2 * pi * t/ 26))
		
		mu <- mu0 + mu1 * t
		
		N <- S + I + R
		
		dS <- mu * N - beta * S * I/N - dr * S
		dI <- beta * S * I/N - gamma * I - dr * I
		dR <- gamma * I - dr * R
		dC <- beta * S * I/N
		dB <- mu * N
		list(c(dS, dI, dR, dC, dB), N=N, mu=mu)
	})
}

par1 <- c(mu0=1/(50*26), dr=1/(50*26), mu1=0, b0=500/26, b1=0.15, gamma=1)
par2 <- c(mu0=1/(50*26), mu1=1/260 * 1/(50*26) * 0.5, dr=1/(50*26), b0=500/26, b1=0.15, gamma=1)
yini <- c(S=0.05*5e6, I=1e-4*5e6, R=5e6 - (0.0501) * 5e6, C=0, B=0)
tvec <- seq(0, 260, by=1)

out1 <- as.data.frame(ode(yini, tvec,sir, par1))
out2 <- as.data.frame(ode(yini, tvec,sir, par2))

rho1 <- rep(0.7, 260)
rho2 <- seq(from=0.6, to=0.8, length.out=260)

data1 <- data.frame(
	births=diff(out1$B),
	cases=diff(out1$C) * rho1,
	pop=head(out1$N, -1),
	time=1:260
)

data2 <- data.frame(
	births=diff(out1$B),
	cases=diff(out1$C) * rho2,
	pop=head(out1$N, -1),
	time=1:260
)

data3 <- data.frame(
	births=diff(out2$B),
	cases=diff(out2$C) * rho1,
	pop=head(out2$N, -1),
	time=1:260
)

data4 <- data.frame(
	births=diff(out2$B),
	cases=diff(out2$C) * rho2,
	pop=head(out1$N, -1),
	time=1:260
)

rt1 <- runtsir(data1, regtype="gaussian")
rt2 <- runtsir(data1, regtype="lm")
rt3 <- runtsir(data1, regtype="spline")
rt4 <- runtsir(data1, regtype="lowess")
rt5 <- runtsir(data1, regtype="loess")

tsir_rho_data <- data.frame(
	gaussian=1/rt1$rho,
	lm=1/rt2$rho,
	spline=1/rt3$rho,
	lowess=1/rt4$rho,
	loess=1/rt5$rho,
	time=1:260
) %>%
	gather(key, value, -time)

tsir_Z_data <- data.frame(
	gaussian=rt1$Z[,1],
	lm=rt2$Z,
	spline=rt3$Z,
	lowess=rt4$Z,
	loess=rt5$Z,
	time=1:260
) %>%
	gather(key, value, -time)

g1 <- ggplot(tsir_rho_data) +
	geom_hline(yintercept=0.7, col="grey", lwd=3, alpha=0.7) +
	geom_line(aes(time, value, col=key, lty=key), lwd=1) +
	scale_y_continuous("Reporting rate") +
	scale_x_continuous("Time (generations)") +
	theme(
		legend.position = "none"
	)

g2 <- ggplot(tsir_Z_data) +
	geom_line(data=data.frame(time=1:260, y=head(out1$S - 0.05*5e6, -1)), aes(time, y), col="grey", lwd=3, alpha=0.7) +
	geom_line(aes(time, value, col=key, lty=key), lwd=1) +
	scale_y_continuous("Susceptible dynamics") +
	scale_x_continuous("Time (generations)") +
	theme(
		legend.title = element_blank()
	)

g3 <- arrangeGrob(g1, g2, nrow=1, widths=c(0.45, 0.6))

ggsave("susceptible_reconstruction_tsir.pdf", g3, width=8, height=3)
