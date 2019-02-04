library(splines)
library(mgcv)
load("../sim/gillespie_sim.rda")
load("../data/gillespie_data.rda")

N <- 1e5
nsim <- length(datalist)

fitlist <- vector('list', nsim)

for (i in 1:nsim) {
	print(i)
	dd <- datalist[[i]][1:20,]
	rr <- reslist[[i]]
	
	## assume that we know S exactly
	S <- approx(x=rr$time, y=N-rr$incidence-10, xout=1:20)$y
	
	logI <- log(dd$incidence)
	
	gfit <- lm(log(incidence)~-1 + ns(time, knots=seq(1, 20, by=6)), data=dd)
	
	fitdata <- data.frame(
		grad=diff(predict(gfit, newdata = data.frame(time=1:21)))+1,
		logS=log(S)
	)
	
	fitdata$offterm <- fitdata$logS - log(N)
	
	lfit <- gam(grad ~ offset(offterm), data=fitdata,
				family = gaussian("log"))
	
	logR0 <- coef(lfit)[[1]]
	
	cdata <- data.frame(
		param=c("beta"),
		mean=c(exp(logR0))
	)
	
	fitlist[[i]] <- cdata
}

save("fitlist", file="gradient_fit.rda")
