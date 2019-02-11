library(tidyr)
library(dplyr)
library(gridExtra)
library(splines)
library(scam)
library(mgcv)
library(deSolve)
library(ggplot2); theme_set(theme_bw())

load("../sim/boston_bifurcation_gradient.rda")

measles_data <- read.csv("../data/measlesUKUS.csv")

measles_US <- measles_data %>% 
	filter(country=="US") %>%
	mutate(cases=ifelse(is.na(cases), 0, cases)) %>%
	filter(year >= 1920, year <= 1940)

measles_list <- measles_US %>%
	split(as.character(.$loc))

bs <- measles_list$BOSTON
bs$time <- 1:nrow(bs)

smooth_time <- seq(0, 26, by=0.01)

BX <- mgcv::cSplineDes(head(bs$biweek, -1), seq(0, 26, by=2))
BXsmooth <- mgcv::cSplineDes(smooth_time, seq(0, 26, by=2))
bdata <- as.data.frame(BX)

## susceptible reconstruction
regdata <- data.frame(
	x=cumsum(bs$cases),
	y=cumsum(bs$rec)
)

regfit <- scam(y~s(x, bs="cr"), data=regdata)

Z <- regfit$residuals
rho <- 1/tsiR::derivative(regdata$x, predict(regfit))

gfit <- MASS::glm.nb(cases~ns(bs$time, knots=seq(1, nrow(bs), by=6)), data=bs)

S0vec <- seq(0.03, 0.07, by=0.001)
llvec <- rep(NA, length(S0vec))

for (j in 1:length(S0vec)) {
	fitdata <- data.frame(
		grad=diff(predict(gfit, newdata = data.frame(time=1:(nrow(bs)+1))))+1+1/(50*26),
		logS=head(log(S0vec[j]*bs$pop + Z),-1),
		biweek=head(bs$biweek,-1)
	)
	
	fitdata$offterm <- fitdata$logS - head(log(bs$pop),-1)
	
	fitdata <- append(fitdata, bdata)
	
	lfit <- glm(grad ~ -1 + V1+V2+V3+V4+V5+V6+V7+V8+V9+V10+V11+V12+V13 + offset(offterm), data=fitdata,
				family = gaussian("log"))
	
	llvec[j] <- logLik(lfit)
}

j <- which.max(llvec)

fitdata <- data.frame(
	grad=diff(predict(gfit, newdata = data.frame(time=1:(nrow(bs)+1))))+1+1/(50*26),
	logS=head(log(S0vec[j]*bs$pop + Z),-1),
	biweek=head(bs$biweek,-1)
)	

fitdata$offterm <- fitdata$logS - head(log(bs$pop),-1)

fitdata <- append(fitdata, bdata)

lfit <- glm(grad ~ -1 + V1+V2+V3+V4+V5+V6+V7+V8+V9+V10+V11+V12+V13 + offset(offterm), data=fitdata,
			family = gaussian("log"))

sir <- function(t, state, parameters) {
	with(as.list(c(state, parameters)), {
		logbeta <- mgcv::cSplineDes(t %% 26, seq(0, 26, by=2)) %*% 
			c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13)
		beta <- exp(logbeta) + relR
		I <- exp(logI)
		dS <- mu * N - beta * S * I/N
		dlogI <- beta * S/N - gamma-mu
		dR <- gamma * I
		dC <- beta * S * I/N
		list(c(dS, dlogI, dR, dC))
	})
}

par <- c(coef(lfit), mu=1/(50*26), N=mean(bs$pop), gamma=1, relR=0)
# yini <- c(S=mean(bs$pop)*S0vec[j], logI=log(520), R=0, C=0)

nsim <- 20

set.seed(101)
init <- data.frame(
	S=seq(from=1e-2, to=1e-1, length.out=nsim) * mean(bs$pop),
	logI=log(seq(from=1e-5, to=1e-3, length.out=nsim) * mean(bs$pop)),
	R=0,
	C=0
)

init[] <- apply(init, 2, sample)

rr <- reslist %>%
	bind_rows

rr_clean <- rr %>%
	mutate(prev=round(prev, 11)) %>%
	group_by(R0, sim) %>%
	mutate(
		dup=duplicated(prev)
	) %>%
	filter(!dup) %>%
	mutate(
		period=length(prev),
		oo=order(prev)
	) %>%
	arrange(R0) %>%
	group_by(R0, period) %>%
	mutate(
		dup=duplicated(prev)
	) %>%
	filter(!dup)

rr_2 <- rr_clean %>%
	filter(period == 2) %>%
	filter(R0 > 20, sim != "4", sim != "10") 

rr_3 <- rr_clean %>%
	filter(period == 3) %>%
	filter(R0 > 18)

rr_4 <- rr_clean %>%
	filter(period==4) %>%
	group_by(R0, sim) %>%
	filter(!all(prev > 1e-5), R0 < 20, R0 != 16.1) %>%
	ungroup %>%
	mutate(
		prev=round(prev, -log10(prev)+2)
	) %>%
	group_by(R0) %>%
	mutate(
		oo=as.numeric(cut(prev, c(1e-20, 3e-9, 2e-8, 1e-5, 1)))
	)

g1 <- ggplot(bs) +
	geom_line(aes(decimalYear, cases)) +
	scale_x_continuous("Time (years)", expand=c(0, 0)) +
	scale_y_continuous("Biweekly cases") +
	theme(
		panel.grid = element_blank()
	)

g2 <- ggplot(data.frame(time=smooth_time/26, trans=exp(BXsmooth %*% coef(lfit)))) +
	geom_line(aes(time, trans)) +
	scale_x_continuous("Time (years)", expand=c(0, 0)) +
	scale_y_continuous("Transmission rate/infectious period") +
	theme(
		panel.grid = element_blank()
	)

g3 <- ggplot(rr_2) +
	geom_line(aes(R0, prev, group=interaction(period, oo))) +
	geom_line(data=rr_3, aes(R0, prev, group=interaction(period, oo))) +
	geom_line(data=rr_4, aes(R0, prev, group=interaction(oo))) +
	geom_point(data=rr_clean, aes(R0, prev), shape=".") +
	geom_vline(xintercept=20.6, lty=2) +
	scale_y_log10("Prevalence") +
	scale_x_continuous("Basic reproductive number", expand=c(0, 0)) +
	theme(
		panel.grid = element_blank()
	)

gall <- arrangeGrob(g1, g2, g3, nrow=1)

ggsave("boston_bifurcation_gradient_fig.pdf", gall, width=8, height=3)
