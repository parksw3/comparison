library(dplyr)
library(ggplot2); theme_set(theme_bw())

load("../sim/compare_likelihood_lognormal.rda")

tsir_likelihood <- tsirlist %>%
	lapply(function(x) data.frame(
		amount=1,
		logLik=logLik(x$glmfit),
		type="conditional",
		alpha=x$alpha
	)) %>%
	bind_rows

likelihood_data <- liklist %>%
	bind_rows

tsir_max <- tsir_likelihood %>%
	filter(logLik==max(logLik))

likd <- likelihood_data %>%
	group_by(amount) %>%
	filter(logLik==max(logLik)) %>%
	bind_rows(tsir_max) %>%
	filter(amount >= 0.1)

g1 <- ggplot(likd) +
	geom_line(aes(amount, logLik), lwd=5, col="grey", alpha=0.3) +
	geom_point(aes(amount, logLik), shape=1) +
	scale_x_continuous("Proportion of process variance") +
	scale_y_continuous("log-likelihood") +
	theme(
		panel.border = element_blank(),
		panel.grid = element_blank(),
		axis.line = element_line()
	)

ggsave("process.pdf", g1, width=6, height=3)
