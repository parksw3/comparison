library(dplyr)
library(ggplot2); theme_set(theme_bw())

filenames <- c(
	"london_pomp_exp_det_fixed_density",
	"london_pomp_exp_det_fixed_frequency",
	"london_pomp_exp_det_exp_density",
	"london_pomp_exp_det_exp_frequency",
	"london_pomp_tsir_det_fixed_density",
	"london_pomp_tsir_det_fixed_frequency",
	"london_pomp_tsir_det_exp_density",
	"london_pomp_tsir_det_exp_frequency"
)

typedata <- data.frame(
	type=1:8,
	infection=c("exponential", "exponential", "exponential", "exponential", "linear", "linear", "linear", "linear"),
	dependence=c("density", "frequency", "density", "frequency", "density", "frequency", "density", "frequency"),
	generation=c("fixed", "fixed", "exponential", "exponential", "fixed", "fixed", "exponential", "exponential")
)

L <- lapply(paste0("../analysis/", filenames, ".rda"), load, .GlobalEnv)

comblist <- lapply(paste0(filenames, "_list"), get)

combfit <- lapply(comblist, function(x) {
	ll <- sapply(x, "[[", 2)
	x[[which.max(ll)]]$mle
})

BX <- mgcv::cSplineDes(seq(0, 26, by=0.1), seq(0, 26, by=2))

combtrans <- combfit %>% 
	lapply(function(x) {
		data.frame(
			trans=BX %*% coef(x)[1:13],
			time=seq(0, 26, by=0.1)
		)
	}) %>%
	bind_rows(.id="type") %>%
	merge(typedata)

ggplot(combtrans) +
	geom_line(aes(time, trans, col=dependence), lwd=1) +
	facet_grid(generation~infection)

combfoi <- combfit %>% 
	lapply(function(x) {
		BXsim <- mgcv::cSplineDes(ld_pomp_arg$covar$biweek, seq(0, 26, by=2))
		beta <- c(BXsim %*% coef(x)[1:13])
		
		sim <- trajectory(x)
		
		data.frame(
			biweek=ld_pomp_arg$covar$biweek,
			foi=beta*sim["I",,]/x@covar[,"pop"]
		)
	}) %>%
	bind_rows(.id="type") %>%
	merge(typedata)

ggplot(combfoi) +
	geom_boxplot(aes(biweek, foi, fill=dependence, group=interaction(biweek, dependence))) +
	facet_grid(generation~infection)

combsim <- combfit %>% 
	lapply(function(x) {
		BXsim <- mgcv::cSplineDes(ld_pomp_arg$covar$biweek, seq(0, 26, by=2))
		beta <- c(BXsim %*% coef(x)[1:13])
		
		sim <- trajectory(x)
		
		if (dim(sim)[1]==3) {
			dd <- data.frame(
				time=1:length(sim[2,,]),
				cases=sim["C",,]
			)
		} else {
			dd <- data.frame(
				time=1:length(sim[2,,]),
				cases=sim["I",,]
			)
		}
		
		dd$logLik <- logLik(x)
		
		dd$Rsquare <- cor(dd$cases, ld_pomp_arg$data$cases)^2
		
		dd
	}) %>%
	bind_rows(.id="type") %>%
	merge(typedata)

ggplot(combsim) +
	geom_line(aes(time, cases, col=dependence)) + 
	geom_point(data=ld_pomp_arg$data, aes(time, cases)) +
	facet_grid(generation~infection)
