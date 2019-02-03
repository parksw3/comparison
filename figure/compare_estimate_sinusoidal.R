library(dplyr)
library(ggplot2); theme_set(theme_bw())

allcover <- list()
alltrans <- list()

load("../sinusoidal_fit/sinusoidal_tsir_fit_default.rda")

allcover$tsir_default <- bind_rows(fitlist, .id="sim")
alltrans$tsir_default <- bind_rows(translist, .id="sim")

load("../sinusoidal_fit/sinusoidal_tsir_fit_cyclic.rda")

allcover$tsir_cyclic <- bind_rows(fitlist, .id="sim")
alltrans$tsir_cyclic <- bind_rows(translist, .id="sim")

load("../sinusoidal_fit/sinusoidal_gradient_fit.rda")

allcover$gradient <- bind_rows(fitlist, .id="sim")
alltrans$gradient <- bind_rows(translist, .id="sim")

alldata <- bind_rows(allcover, .id="type")

coverdata <- alldata %>%
	filter(param %in% c("R0", "rprob")) %>%
	group_by(param, type) %>%
	summarize(
		cover=mean(coverage)
	)

ggplot(coverdata) +
	geom_rect(xmin=-Inf, xmax=Inf, ymin=0.887, ymax=0.983, alpha=0.1) +
	geom_point(aes(type, cover), size=2) +
	geom_hline(yintercept=0.95, lty=2) +
	scale_y_continuous("Coverage probability", limits=c(0, 1)) +
	facet_wrap(~param, scale="free_x")

estdata <- alldata %>%
	filter(param %in% c("R0", "rprob")) %>%
	group_by(param, type) %>%
	summarize(
		est.mean=mean(mean),
		est.lwr=quantile(mean, 0.025),
		est.upr=quantile(mean, 0.975)
	)

ggplot(estdata) +
	geom_point(aes(type, est.mean), size=2) +
	geom_errorbar(aes(type, ymin=est.lwr, ymax=est.upr), width=0.1, lwd=1) +
	facet_wrap(~param, scale="free")

transdata <- alltrans %>%
	bind_rows(.id="type")

ggplot(transdata) +
	geom_line(aes(time, beta, group=interaction(sim, type)), alpha=0.2) +
	stat_function(fun=function(x) (500/26 * (1 + 0.15 * cos(x * 2 * pi/26))), col=2, lwd=1) +
	facet_wrap(~type)

