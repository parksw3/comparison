library(dplyr)
library(tsiR)
measles_data <- read.csv("../data/measlesUKUS.csv")

measles_US <- measles_data %>% 
	filter(country=="US") %>%
	mutate(cases=ifelse(is.na(cases), 0, cases)) %>%
	filter(year >= 1920, year <= 1940)

measles_list <- measles_US %>%
	split(as.character(.$loc))

bs <- measles_list$BOSTON
bs$births <- bs$rec
bs$time <- 1:nrow(bs)

alphavec <- seq(0.9, 1, by=0.005)

rrbase <- runtsir(bs, alpha=0.97, nsim=1, sbar=0.035)

reslist <- vector('list', length(alphavec))

for (i in 1:length(alphavec)) {
	print(i)
	alpha <- alphavec[i]
	
	rr <- runtsir(bs, userYhat=rrbase$Yhat, regtype="user", alpha=alpha, nsim=1, sbar=0.035,
				  inits.fit = TRUE)

	reslist[[i]] <- data.frame(
		alpha=alpha,
		reg=logLik(rr$glmfit)[[1]],
		gof=-sum((rr$res[,1]-bs$cases)^2)
	)
}

resdata <- reslist %>%
	bind_rows

plot(resdata$alpha, resdata$gof, type="l")
