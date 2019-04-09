library(dplyr)
library(ggplot2); theme_set(theme_bw())

load("../sinusoidal_fit/sinusoidal_tsir_fit_default.rda")

transdata <- translist %>%
	bind_rows(.id="sim")

gtrans <- ggplot(transdata) +
	geom_line(aes(time, beta, group=sim), alpha=0.2) +
	stat_function(fun=function(x) (500/26 * (1 + 0.15 * cos(x * 2 * pi/26))), col=2, lwd=1) +
	scale_x_continuous("Time (generations)", expand=c(0, 0)) +
	scale_y_continuous("Transmission rates", expand=c(0, 0)) +
	theme(
		strip.background = element_blank(),
		panel.grid = element_blank()
	)

ggsave("tsir_transmission.pdf", gtrans, width=6, height=4)
