library(tidyr)
library(dplyr)
library(emdbook)
library(mgcv)
library(deSolve)
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

par <- c(beta=2, gamma=1, N=1e5)
y <- c(S=1e5-10, logI=log(10), R=0)
time <- seq(0, 20, by=1)
smooth_time <- seq(0, 20, by=0.01)

out <- as.data.frame(ode(y, time, sir, par))
out_smooth <- as.data.frame(ode(y, smooth_time, sir, par))

true_incidence <- data.frame(
	time=1:20,
	incidence=-diff(out$S)
)

rprob <- 0.7
theta <- 10
nsim <- 100

reslist <- vector('list', nsim)

set.seed(101)
for (i in 1:nsim) {
	simulated_incidence <- data.frame(
		time=1:20,
		incidence=rbetabinom(20, size=round(true_incidence$incidence), prob=rprob, theta=theta)
	)
	
	gfit <- gam(incidence~s(time), data=simulated_incidence, family=nb)
	I <- exp(predict(gfit, newdata=data.frame(time=smooth_time)))
	
	reslist[[i]] <- data.frame(
		time=head(smooth_time, -1),
		grad=diff(I)/0.01
	)
}

graddata <- reslist %>%
	bind_rows(.id="sim")

true_grad <- data.frame(
	time=out_smooth$time,
	grad=rprob * (2 * out_smooth$S/1e5 - 1) * exp(out_smooth$logI)
)

g <- ggplot(graddata) +
	geom_line(aes(time, grad, group=sim, col="fitted"), alpha=0.1) +
	geom_line(data=true_grad, aes(time, grad, col="true"), lwd=1) +
	ylab("Gradient") +
	scale_x_continuous("Time (generations)", expand=c(0,0)) +
	scale_color_manual(values=c("black", "red")) +
	theme(
		panel.grid = element_blank(),
		legend.position = c(0.9, 0.9),
		legend.title = element_blank()
	)

ggsave("gradient_matching.pdf", g, width=6, height=4)
