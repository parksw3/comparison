library(splines)
library(mgcv)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)

measles_data <- read.csv("../data/measlesUKUS.csv")

ld <- measles_data %>%
	filter(loc=="LONDON")

regdata <- data.frame(
	x=cumsum(ld$cases),
	y=cumsum(ld$rec)
)

regfit <- lm(y~x, data=regdata)

Z <- regfit$residuals
rho <- 1/coef(regfit)[[2]]

logI <- log(ld$cases)

gfit <- MASS::glm.nb(cases~ns(decimalYear, df=150), data=ld)

pdata <- data.frame(
	time=ld$decimalYear,
	pp=predict(gfit)
)

g1 <- ggplot(ld) +
	geom_point(aes(decimalYear, log(cases)), shape=1) +
	geom_line(data=pdata, aes(time, pp)) +
	scale_x_continuous("Time (years)", expand=c(0, 0), breaks=seq(1945, 1960, by=5)) +
	scale_y_continuous("log(cases)") +
	theme(
		panel.grid = element_blank()
	)

gdata <- data.frame(
	time=ld$decimalYear,
	gg=(predict(gfit, newdata=data.frame(decimalYear=ld$decimalYear+0.1))-pdata$pp)/0.1
)

g2 <- ggplot(gdata) +
	geom_line(aes(time, gg)) +
	scale_x_continuous("Time (years)", expand=c(0, 0), breaks=seq(1945, 1960, by=5)) +
	scale_y_continuous("Derivative") +
	theme(
		panel.grid = element_blank()
	)

fitdata <- data.frame(
	S=0.035*mean(ld$pop)+Z,
	y=gdata$gg + 26 + 1/50,
	N=ld$pop,
	biweek=ld$biweek
)

tempdata <- fitdata[fitdata$biweek==26,]
tempdata$biweek <- 0

fitdata <- rbind(fitdata, tempdata)

fitdata$offterm <- log(fitdata$S) - log(fitdata$N)

gamfit <- gam(y~s(biweek, bs="cc") + offset(offterm), data=fitdata,
	family = gaussian("log"))

tdata <- data.frame(
	time=seq(0, 26, by=0.01),
	rate=exp(predict(gamfit, newdata=data.frame(biweek=seq(0, 26, by=0.01), offterm=0)))
)

g3 <- ggplot(tdata) +
	geom_line(aes(time/26, rate)) +
	scale_x_continuous("Time (years)", expand=c(0,0), breaks=c(0, 0.2, 0.4, 0.6, 0.8)) +
	scale_y_continuous("Transmission rate") +
	theme(
		panel.grid = element_blank()
	)

ggsave("grad1.pdf", g1, width=6, height=3)
ggsave("grad2.pdf", g2, width=6, height=3)
ggsave("grad3.pdf", g3, width=6, height=3)
