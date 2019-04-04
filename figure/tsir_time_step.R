library(mgcv)
library(tidyr)
library(dplyr)
library(deSolve)
library(ggplot2); theme_set(theme_bw(base_size=12, base_family = "Times"))
library(gridExtra)
library(emdbook)
source("../R/tsir_util.R")

if (.Platform$OS.type=="windows") {
	windowsFonts(Times=windowsFont("Times"))
} 

scale_colour_discrete <- function(...,palette="Dark2") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Dark2") scale_fill_brewer(...,palette=palette)

rprob <- 0.7
theta <- 10

simulate_data <- function(oo,
						  tstep=5,
						  state="CI",
						  seed=101,
						  rprob=0.7,
						  theta=10) {
	set.seed(seed)
	
	tcut <- seq(0, 365, by=tstep)
	tcut[length(tcut)] <- 365
	
	twhich <- (oo$time %% 365) %in% tcut
	
	oo_subset <- oo[twhich,]
	
	birth <- diff(oo_subset$CB)
	cases <- diff(oo_subset[[state]])
	
	data.frame(
		year=head(floor(oo_subset$time/365),-1),
		birth=rpois(length(birth), birth),
		cases=rbetabinom(length(cases), rprob, round(cases), theta),
		trueS=head(oo_subset$S,-1),
		factor=rep(1:(length(tcut)-1), 100)[1:length(cases)],
		time=head(oo_subset$time/365,-1)
	)
}

simple_tsir <- function(data, N=5e6) {
	regdata <- data.frame(
		y=cumsum(data$birth),
		x=cumsum(data$cases)
	)
	
	lfit <- lm(y~x, data=regdata)
	
	rho <- 1/coef(lfit)[2]
	Z <- residuals(lfit)
	
	Smean <- seq(0.03, 0.2, by=0.001)
	liklist <- rep(NA, length(Smean))
	
	for (i in 1:length(Smean)) {
		fitdata <- data.frame(
			Inew=tail(data$cases,-1)/rho,
			Iprev=head(data$cases,-1)/rho,
			rt=head(data$factor,-1),
			S=head(Smean[i]*N+Z,-1),
			N=N
		)
		
		tempdata <- fitdata[fitdata$rt==max(fitdata$rt),]
		tempdata$rt <- 0
		
		fitdata <- bind_rows(fitdata, tempdata)
		
		fitdata$offterm <- log(fitdata$S) - log(fitdata$N)
		
		gfit <- gam(log(Inew) ~ log(Iprev) + s(rt, bs="cc") + offset(offterm), data=fitdata)
		
		liklist[[i]] <- logLik(gfit)
	}
	
	j <- which.max(liklist)
	
	fitdata <- data.frame(
		Inew=tail(data$cases,-1)/rho,
		Iprev=head(data$cases,-1)/rho,
		rt=head(data$factor,-1),
		S=head(Smean[j]*N+Z,-1),
		N=N
	)
	
	tempdata <- fitdata[fitdata$rt==max(fitdata$rt),]
	tempdata$rt <- 0
	
	fitdata <- bind_rows(fitdata, tempdata)
	
	fitdata$offterm <- log(fitdata$S) - log(fitdata$N)
	
	gfit <- gam(log(Inew) ~ log(Iprev) + s(rt, bs="cc") + offset(offterm), data=fitdata)
	
	tdata <- data.frame(
		time=seq(0, max(fitdata$rt), by=0.1)/max(fitdata$rt),
		raw=exp(predict(gfit, newdata=data.frame(rt=seq(0, max(fitdata$rt), by=0.1), Iprev=1, offterm=0))),
		fix=fixfun(predict(gfit, newdata=data.frame(rt=seq(0, max(fitdata$rt), by=0.1), Iprev=1, offterm=0)), coef(gfit)[2], N, fitdata$Iprev)
	)
	
	list(
		beta=tdata,
		rho=rho,
		Smean=Smean[j],
		gfit=gfit,
		Z=Z,
		alpha=coef(gfit)[2]
	)
}

gfun <- function(t, y, pars) {
	with(as.list(c(y, pars)), {
		beta <- b0 * (1 + b1 * cos(2 * pi * t/365))
		dS <- mu * (N - S) - beta * S * I/N
		dE <- beta * S * I/N - (sigma + mu) * E
		dI <- sigma * E - (gamma + mu) * I
		dR <- gamma * I - mu * R
		
		dCB <- mu * N
		dCE <- beta * S * I/N
		dCI <- sigma * E
		dCR <- gamma * I
		
		list(c(dS, dE, dI, dR, dCB, dCE, dCI, dCR))
	})
}

yini <- c(S=0.05*5e6, E= 1e-4*5e6, I=1e-4*5e6, R=(1-(0.05+2e-4))*5e6, CB=0, CE=0, CI=0, CR=0)
pars <- c(mu=1/(50*365), b0=1250/365, b1=0.15, sigma=1/8, gamma=1/5, N=5e6)
tvec <- seq(0, 60 * 365, by=1)

out <- as.data.frame(ode(yini, tvec, gfun, pars))

out_subset <- out %>% 
	filter(time >= 40*365)

data_5 <- simulate_data(out_subset)
data_7 <- simulate_data(out_subset, tstep=7)
data_13 <- simulate_data(out_subset, tstep=13)
data_14 <- simulate_data(out_subset, tstep=14)

fit_5 <- simple_tsir(data_5)
fit_7 <- simple_tsir(data_7)
fit_13 <- simple_tsir(data_13)
fit_14 <- simple_tsir(data_14)

g_case1 <- ggplot(data_5) +
	geom_line(aes(time, cases)) +
	geom_point(aes(time, cases), size=0.5) +
	scale_x_continuous("Time (years)", expand=c(0, 0)) +
	scale_y_continuous("Cases") +
	theme(
		panel.grid = element_blank(),
		panel.border = element_blank(),
		axis.line = element_line()
	) + 
	ggtitle("5 days")

g_case2 <- (g_case1 %+% data_7) + ggtitle("7 days")
g_case3 <- g_case1 %+% data_13 + ggtitle("13 days")
g_case4 <- g_case1 %+% data_14 + ggtitle("14 days")

trans_default <- ggplot(fit_5$beta) +
	stat_function(fun=function(x) 1250/365 * (1 + 0.15 * cos(2 * pi * x))*5, lwd=2, aes(col="true"), alpha=0.5) +
	geom_line(aes(time, raw, col="raw"), lty=2) +
	scale_x_continuous("Time (years)", expand=c(0,0)) +
	scale_y_continuous("Scaled transmission rate", limits=c(9, 30)) +
	scale_colour_manual(values=c("black", "red", "grey")) +
	theme(
		legend.title = element_blank(),
		panel.grid = element_blank()
	)

g_trans1 <- (trans_default %+% fit_5$beta) +
	geom_line(aes(time, fix * 13/5, col="corrected")) +
	theme(
		legend.position = c(0.5, 0.92),
		legend.key.height = grid::unit(0.1, "cm"),
		legend.direction = "horizontal"
	)

g_trans2 <- (trans_default %+% fit_7$beta) +
	geom_line(aes(time, fix * 13/7, col="corrected"))  +
	theme(
		legend.position = "none"
	)

g_trans3 <- (trans_default %+% fit_13$beta) +
	geom_line(aes(time, fix * 13/13, col="corrected"))  +
	theme(
		legend.position = "none"
	)

g_trans4 <- (trans_default %+% fit_14$beta) +
	geom_line(aes(time, fix * 13/14, col="corrected"))  +
	theme(
		legend.position="none"
	)

gfinal <- arrangeGrob(
	arrangeGrob(g_case1, g_case2, g_case3, g_case4, nrow=1),
	arrangeGrob(g_trans1, g_trans2, g_trans3, g_trans4, nrow=1)
)

ggsave("tsir_time_step.pdf", gfinal, width=10, height=6)
