library(dplyr)
library(splines)
library(scam)
library(deSolve)
library(mgcv)
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
		dS <- mu * (N - S) - beta * S * I/N
		dlogI <- beta * S/N - gamma - mu
		dC <- beta * S * I/N
		list(c(dS, dlogI, dC))
	})
}

par <- c(coef(lfit), mu=1/(50*26), N=mean(bs$pop), gamma=1, relR=0)
# yini <- c(S=mean(bs$pop)*S0vec[j], logI=log(520), R=0, C=0)

nsim <- 20

set.seed(101)
init <- data.frame(
	S=seq(from=1e-2, to=1e-1, length.out=nsim) * mean(bs$pop),
	logI=log(seq(from=1e-5, to=1e-3, length.out=nsim) * mean(bs$pop)),
	C=0
)

init[] <- apply(init, 2, sample)

##relative R0
R0vec <- seq(-5, 15, by=0.1)

reslist <- vector('list', length(R0vec))

for (i in 1:length(R0vec)) {
	print(i)
	relR <- R0vec[i]
	
	cc <- par
	cc[["relR"]] <- relR
	
	reslist[[i]] <- lapply(1:nsim, function(x) {
		yini <- unlist(init[x,])
		out <- as.data.frame(rk4(yini, 1:4000, sir, cc))
		oo <- out[seq(3480, by=26, length.out=20),]
		data.frame(
			R0=20.6+relR,
			prev=exp(oo$logI)/mean(bs$pop)
		)
	}) %>%
		bind_rows(.id="sim")
	
	save("reslist", file="boston_bifurcation_gradient.rda")
}

save("reslist", file="boston_bifurcation_gradient.rda")
