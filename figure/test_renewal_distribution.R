library(dplyr)
library(mgcv)
library(pomp)
library(ggplot2); theme_set(theme_bw(base_size = 12,
									 base_family = "Times"))

source("../R/fitfun_pomp_renewal_basis.R")
source("../R/basis.R")


if (.Platform$OS.type=="windows") {
	windowsFonts(Times=windowsFont("Times"))
} 

scale_colour_discrete <- function(...,palette="Dark2") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Dark2") scale_fill_brewer(...,palette=palette)

load("../data/gillespie_sinusoidal_data.rda")

globals <- Csnippet(paste0("double N0=5000000; double Gmean=1; double mu=0.0007692308;"))

k <- seq(0, 26, by=2)

BX <- mgcv::cSplineDes(seq(0, 260, by=0.1)%%26,k)
basis_covar <- as.data.frame(BX)
colnames(basis_covar) <- paste0("B", 1:ncol(basis_covar))
basis_covar$time <- seq(0, 260, by=0.1)

bb <- estimate_basis(500/26 *(1 + 0.15 * cos(1:26 * 2 * pi/26)))

j <- 1

dd <- datalist[[j]][1:(10*26),]

pomp_arg <- make_pomp_renewal_basis()

pomp_arg2 <- append(pomp_arg, list(data=dd[,c("time", "incidence")], globals=globals, t0=0, covar=basis_covar, tcovar="time"))

pomp_model <- do.call(pomp, pomp_arg2)

param1 <- c(bb, rho=0.7, S0=0.05, I0=1e-4, disp=5, Gvar=1)
param2 <- c(bb, rho=0.7, S0=0.05, I0=1e-4, disp=5, Gvar=0.5^2)
param3 <- c(bb, rho=0.7, S0=0.05, I0=1e-4, disp=5, Gvar=0.1^2)

set.seed(101)
ss1 <- replicate(500, {
	ss <- simulate(pomp_model, param=param1, as.data.frame=TRUE)
	data.frame(
		incidence=ss$incidence,
		time=1:(10*26)
	)}, simplify=FALSE) %>%
	bind_rows(.id="sim") %>%
	group_by(time) %>%
	summarize(
		median=median(incidence),
		lwr=quantile(incidence, 0.025),
		upr=quantile(incidence, 0.975)
	) %>%
	mutate(
		Gvar=1
	)

ss2 <- replicate(500, {
	ss <- simulate(pomp_model, param=param2, as.data.frame=TRUE)
	data.frame(
		incidence=ss$incidence,
		time=1:(10*26)
	)}, simplify=FALSE) %>%
	bind_rows(.id="sim") %>%
	group_by(time) %>%
	summarize(
		median=median(incidence),
		lwr=quantile(incidence, 0.025),
		upr=quantile(incidence, 0.975)
	) %>%
	mutate(
		Gvar=0.5
	)

ss3 <- replicate(500, {
	ss <- simulate(pomp_model, param=param3, as.data.frame=TRUE)
	data.frame(
		incidence=ss$incidence,
		time=1:(10*26)
	)}, simplify=FALSE) %>%
	bind_rows(.id="sim") %>%
	group_by(time) %>%
	summarize(
		median=median(incidence),
		lwr=quantile(incidence, 0.025),
		upr=quantile(incidence, 0.975)
	) %>%
	mutate(
		Gvar=0.1
	)

ss_all <- bind_rows(ss1, ss2, ss3) %>%
	mutate(
		Gvar=factor(Gvar,levels=c(0.1, 0.5, 1))
	) %>%
	rename(
		CV=Gvar
	)

g1 <- ggplot(ss_all) +
	geom_ribbon(aes(time, ymin=lwr, ymax=upr, fill=CV), alpha=0.3) +
	geom_line(aes(time, median, col=CV)) +
	scale_x_continuous("Time (generations)", expand=c(0, 0)) +
	scale_y_log10("Cases", limits=c(80, 30000)) +
	theme(
		panel.grid=element_blank(),
		strip.background = element_blank()
	)

ggsave("test_renewal_distribution.pdf", g1, width=8, height=3)
