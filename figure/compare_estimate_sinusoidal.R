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

templist <- templist2 <- list()
for (i in 0:9) {
	fn <- paste0("../sinusoidal_fit/sinusoidal_pomp_0_1_fit_", i, ".rda")
	
	load(fn)
	
	templist[[i+1]] <- fitlist %>% 
		bind_rows(.id="sim") %>%
		mutate(sim=as.character(as.numeric(sim)+i*10))
	
	templist2[[i+1]] <- translist %>% 
		bind_rows(.id="sim") %>%
		mutate(sim=as.character(as.numeric(sim)+i*10))
}

allcover$pomp <- bind_rows(templist, .id="sim")
alltrans$pomp <- bind_rows(templist2, .id="sim")

transdata <- alltrans %>%
	bind_rows(.id="type")

gtrans <- ggplot(transdata) +
	geom_line(aes(time, beta, group=interaction(sim, type)), alpha=0.2) +
	stat_function(fun=function(x) (500/26 * (1 + 0.15 * cos(x * 2 * pi/26))), col=2, lwd=1) +
	facet_wrap(~type)

ggsave("compare_transmission.pdf", gtrans, width=6, height=4)
