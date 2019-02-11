library(splines)
library(tidyr)
library(dplyr)
library(emdbook)
library(deSolve)
library(ggplot2); theme_set(theme_bw(base_size = 12,
									 base_family = "Times"))

if (.Platform$OS.type=="windows") {
	windowsFonts(Times=windowsFont("Times"))
} 

## use Dark2 rather than Set1 because colour #6 of Set1 is yellow (ugh/too light)
scale_colour_discrete <- function(...,palette="Dark2") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Dark2") scale_fill_brewer(...,palette=palette)

sir <- function(t, state, parameters) {
	with(as.list(c(state, parameters)), {
		beta <- b0 * (1 + b1 * cos(2 * pi * t/ 26))
		I <- exp(logI)
		dS <- mu * (N - S) - beta * S * I/N
		dlogI <- beta * S/N - gamma - mu
		dR <- gamma * I
		dC <- beta * S * I/N
		list(c(dS, dlogI, dR, dC))
	})
}

params <- c(mu=1/(50*26), b0=500/26, b1=0.15, gamma=1, N=5e6)
yini <- c(S=0.05*5e6, logI=log(1e-4*5e6), R=0, C=0)
tvec <- seq(0, 260, by=1)
smooth_tvec <- seq(1, 260.01, by=0.1)

out <- as.data.frame(ode(yini, tvec,sir, params))
out_smooth <- as.data.frame(ode(yini, smooth_tvec, sir, params))

true_incidence <- data.frame(
	time=1:260,
	incidence=diff(out$C)
)

rprob <- 0.7
theta <- 10
nsim <- 100

reslist <- vector('list', nsim)

set.seed(101)
for (i in 1:nsim) {
	simulated_incidence <- data.frame(
		time=1:260,
		incidence=rbetabinom(260, size=round(true_incidence$incidence), prob=rprob, theta=theta)
	)
	
	gfit <- MASS::glm.nb(incidence~ns(time, knots=seq(1, 260, by=6)), data=simulated_incidence)
	
	logI <- predict(gfit, newdata=data.frame(time=smooth_tvec))
	
	reslist[[i]] <- data.frame(
		time=head(smooth_tvec, -1),
		grad=diff(logI)/0.1
	)
}

graddata <- reslist %>%
	bind_rows(.id="sim")

true_grad <- data.frame(
	time=out_smooth$time,
	grad=500/26 * (1 + 0.15 * cos(2 * pi * smooth_tvec/26)) * out_smooth$S/(5e6) - 1
)

g <- ggplot(graddata) +
	geom_line(aes(time, grad, group=sim, col="fitted"), alpha=0.1) +
	geom_line(data=true_grad, aes(time, grad, col="true"), lwd=1) +
	scale_y_continuous("Gradient", limits=c(-0.5, 0.6)) +
	scale_x_continuous("Time (generations)", expand=c(0,0)) +
	scale_color_manual(values=c("black", "red")) +
	theme(
		panel.grid = element_blank(),
		legend.position = c(0.9, 0.9),
		legend.title = element_blank(),
		legend.direction="horizontal"
	)

ggsave("gradient_matching_sinusoidal.pdf", g, width=8, height=3)
