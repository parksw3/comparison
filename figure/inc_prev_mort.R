library(gridExtra)
library(deSolve)
library(ggplot2); theme_set(theme_bw(base_size = 12,
									 base_family = "Times"))

if (.Platform$OS.type=="windows") {
	windowsFonts(Times=windowsFont("Times"))
} 

scale_colour_discrete <- function(...,palette="Dark2") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Dark2") scale_fill_brewer(...,palette=palette)

afun <- function(label,size=8) {
	annotate(geom="text",label=label,x=5,y=Inf,
			 ## http://stackoverflow.com/questions/20083700/how-to-have-annotated-text-style-to-inherit-from-theme-set-options
			 family= theme_get()$text[["family"]],
			 size=size,
			 vjust=1.5,hjust=-0.5)
	## vjust=0.98,hjust=0.02)
}

sir <- function(t, state, parameters) {
	with(as.list(c(state, parameters)), {
		beta <- b0 * (1 + b1 * cos(2 * pi * t/ 26))
		I <- exp(logI)
		dS <- mu * N - beta * S * I/N
		dlogI <- beta * S/N - gamma
		dR <- gamma * I
		dC <- beta * S * I/N
		list(c(dS, dlogI, dR, dC))
	})
}

params <- c(mu=1/(50*26), b0=500/26, b1=0.15, gamma=1, N=5e6)
yini <- c(S=0.05*5e6, logI=log(1e-4*5e6), R=0, C=0)
tvec <- seq(0, 260, by=1)
tvec2 <- seq(0, 260, by=0.5)
smooth_tvec <- seq(1, 260.01, by=0.01)

out <- as.data.frame(ode(yini, tvec, sir, params))
out2 <- as.data.frame(ode(yini, tvec2, sir, params))
out_smooth <- as.data.frame(ode(yini, smooth_tvec, sir, params))

incidence_data <- data.frame(
	time=1:260,
	incidence=diff(out$C)
)

incidence_data2 <- data.frame(
	time=seq(0.5, 260, by=0.5),
	incidence=diff(out2$C)
)

mortality_data <- data.frame(
	time=1:260,
	mortality=diff(out$R)
)

mortality_data2 <- data.frame(
	time=seq(0.5, 260, by=0.5),
	mortality=diff(out2$R)
)

g1 <- ggplot(incidence_data) +
	geom_point(aes(time, incidence, col="incidence")) +
	geom_line(aes(time, incidence, col="incidence"), lwd=1) +
	geom_point(data=mortality_data, aes(time, mortality, col="mortality")) +
	geom_line(data=mortality_data, aes(time, mortality, col="mortality"), lwd=1) +
	geom_line(data=out_smooth, aes(time, exp(logI), col="prevalence"), lwd=1) +
	scale_x_continuous("Time (generations)", expand=c(0,0)) +
	scale_y_continuous("Cases", limits=c(0, 2.7e4), expand=c(0,0)) +
	afun("a") +
	theme(
		legend.title = element_blank(),
		legend.position = c(1-0.18, 0.9),
		legend.direction="horizontal",
		panel.grid=element_blank()
	)

g2 <- ggplot(incidence_data2) +
	geom_point(aes(time, incidence, col="incidence")) +
	geom_line(aes(time, incidence, col="incidence"), lwd=1) +
	geom_point(data=mortality_data2, aes(time, mortality, col="mortality")) +
	geom_line(data=mortality_data2, aes(time, mortality, col="mortality"), lwd=1) +
	geom_line(data=out_smooth, aes(time, exp(logI), col="prevalence"), lwd=1) +
	scale_x_continuous("Time (generations)", expand=c(0,0)) +
	scale_y_continuous("Cases", limits=c(0, 2.7e4), expand=c(0,0)) +
	afun("b") +
	theme(
		legend.title = element_blank(),
		legend.position = "none",
		panel.grid=element_blank()
	)

gtot <- arrangeGrob(g1, g2, nrow=2)

ggsave("inc_prev_mort.pdf", gtot, width=8, height=6)
