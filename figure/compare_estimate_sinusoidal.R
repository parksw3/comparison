library(dplyr)
library(ggplot2); theme_set(theme_bw())

alltrans <- list()

load("../sinusoidal_fit/sinusoidal_tsir_fit_cyclic.rda")

alltrans$tsir_cyclic <- bind_rows(translist, .id="sim")

load("../sinusoidal_fit/sinusoidal_gradient_fit.rda")

alltrans$gradient <- bind_rows(translist, .id="sim")

templist <- list()
for (i in 0:9) {
	fn <- paste0("../sinusoidal_fit/sinusoidal_pomp_0_1_fit_", i, ".rda")
	
	load(fn)
	
	templist[[i+1]] <- translist %>% 
		bind_rows(.id="sim") %>%
		mutate(sim=as.character(as.numeric(sim)+i*10))
}

alltrans$pomp <- bind_rows(templist, .id="sim")

templist <- list()
for (i in 0:9) {
	fn <- paste0("../sinusoidal_fit/sinusoidal_pomp_renewal_fit_", i, ".rda")
	
	load(fn)
	
	templist[[i+1]] <- translist %>% 
		bind_rows(.id="sim") %>%
		mutate(sim=as.character(as.numeric(sim)+i*10))
}

alltrans$renewal <- bind_rows(templist, .id="sim")

transdata <- alltrans %>%
	bind_rows(.id="type")

gtrans <- ggplot(transdata) +
	geom_line(aes(time, beta, group=interaction(sim, type)), alpha=0.2) +
	stat_function(fun=function(x) (500/26 * (1 + 0.15 * cos(x * 2 * pi/26))), col=2, lwd=1) +
	scale_x_continuous(expand=c(0, 0)) +
	facet_wrap(~type, nrow=1)

ggsave("compare_transmission.pdf", gtrans, width=6, height=4)
