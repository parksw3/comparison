library(splines)
library(mgcv)
load("../data/tsir_data.rda")

N <- 1e5
nsim <- length(datalist)

fitlist <- vector('list', nsim)

for (i in 1:nsim) {
	print(i)
	dd <- datalist[[i]][1:20,]
	
	I <- dd$incidence/0.7
	Inew <- tail(I, -1)
	Iprev <- head(I, -1)
	
	## assume that we know S exactly
	S <- dd$S
	
	logI <- log(dd$incidence)
	
	gfit <- lm(log(incidence)~-1 + ns(t, knots=seq(1, 20, by=6)), data=dd)
	
	fitdata <- data.frame(
		grad=diff(predict(gfit, newdata = data.frame(t=1:21)))+1,
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

save("fitlist", file="tsir_gradient_fit.rda")
