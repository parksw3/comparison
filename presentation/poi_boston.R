library(dplyr)
library(ggplot2); theme_set(theme_bw())

rr <- read.csv("../data/measlesUKUS.csv")

boston <- rr %>%
	filter(loc=="BOSTON") %>%
	filter(year >= 1920, year < 1940) %>%
	rename(
		time=decimalYear,
		cases=cases,
		pop=pop,
		births=rec
	)

regdata <- data.frame(
	Y=cumsum(boston$births),
	X=cumsum(boston$cases)
)

lfit <- lm(Y~X, data=regdata)

fitdata <- data.frame(
	Inew=tail(boston$cases,-1)*coef(lfit)[2]+1,
	Iprev=head(boston$cases,-1)*coef(lfit)[2]+1,
	N=head(boston$pop,-1),
	S=0.07 * mean(boston$pop) + head(lfit$residuals, -1),
	biweek=head(boston$biweek, -1)
)

ffit <- lm(log(Inew) ~ -1 + as.factor(biweek) + log(Iprev) + offset(log(S) - log(N)), data=fitdata)

resdata <- data.frame(
	incidence=fitdata$Iprev,
	poi=exp(predict(ffit))/fitdata$S,
	biweek=fitdata$biweek
)

poi <- lapply(1:26, function(x) {
	data.frame(
		incidence=seq(0, max(fitdata$Iprev[fitdata$biweek==x])+10, by=1),
		poi=exp(coef(ffit)[x]) * seq(0, max(fitdata$Iprev[fitdata$biweek==x])+10, by=1)^tail(coef(ffit), 1)/mean(fitdata$N),
		biweek=x
	)
}) %>%
	bind_rows

ggplot(resdata) +
	geom_point(aes(incidence, poi)) +
	geom_line(data=poi, aes(incidence, poi)) +
	facet_wrap(~biweek, scale="free")

