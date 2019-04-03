library(dplyr)
library(ggplot2); theme_set(theme_bw())

likelihood_list <- list()

load("../sim/compare_likelihood_lognormal.rda")

likelihood_list$conditional <- liklist

load("../sim/compare_likelihood_lognormal_process.rda")

likelihood_list$unconditional <- liklist

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

tsir_likelihood <- tsirlist %>%
	lapply(function(x) data.frame(
		amount=1,
		logLik=logLik(x$glmfit),
		type="conditional",
		alpha=x$alpha
	)) %>%
	bind_rows

likelihood_data <- likelihood_list %>%
	lapply(bind_rows) %>%
	bind_rows(.id="type")

global_max <- likelihood_data %>%
	group_by(type) %>%
	filter(logLik==max(logLik))

plot(tsir_likelihood$alpha, tsir_likelihood$logLik)

abline(h=max(tsir_likelihood$logLik)-qchisq(0.95, 1)/2)

g1 <- ggplot(filter(likelihood_data, amount > 0)) +
	geom_raster(aes(alpha, amount, fill=logLik)) +
	stat_contour(data=filter(likelihood_data, type=="conditional"),
				 aes(alpha, amount, z=logLik), 
				 breaks=c(global_max$logLik[1]-qchisq(0.95, 2)/2),
				 col="black", lty=2) +
	stat_contour(data=filter(likelihood_data, type=="unconditional"),
				 aes(alpha, amount, z=logLik), 
				 breaks=c(global_max$logLik[2]-qchisq(0.95, 2)/2),
				 col="black", lty=2) +
	geom_point(data=global_max, aes(alpha, amount), size=2, shape=1) +
	scale_x_continuous(expression(alpha), expand=c(0,0)) +
	scale_y_continuous("Proportion of process variance", expand=c(0, 0), breaks=1:9/10) +
	scale_fill_gradientn(colours=terrain.colors(10)) +
	facet_wrap(~type) +
	theme(
		strip.background = element_blank(),
		panel.spacing = grid::unit(1, "cm")
	)

ggsave("compare_conditional_likelihood.pdf", g1, width=8, height=3)
