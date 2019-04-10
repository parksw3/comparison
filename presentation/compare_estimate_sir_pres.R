library(tidyr)
library(dplyr)
library(gridExtra)
library(ggplot2); theme_set(theme_bw(base_size = 12,
									 base_family = "Times"))

if (.Platform$OS.type=="windows") {
	windowsFonts(Times=windowsFont("Times"))
} 

scale_colour_discrete <- function(...,palette="Dark2") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Dark2") scale_fill_brewer(...,palette=palette)

classify <- data.frame(
	type=c("trajectory", "gradient", "tsir", "pomp", "renewal_pomp"),
	time=c("Continuous", "Continuous", "Discrete", "Discrete", "Discrete"),
	est=c("Simulation", "Regression", "Regression", "Simulation", "Simulation")
)

allcover_sir <- list()

load("../sir_fit/trajectory_fit.rda")

allcover_sir$trajectory <- bind_rows(fitlist, .id="sim")

templist <- list()
for (i in 0:9) {
	fn <- paste0("../sir_fit/pomp_0_1_fit_", i, ".rda")
	
	load(fn)
	
	templist[[i+1]] <- fitlist %>% 
		bind_rows(.id="sim") %>%
		mutate(sim=as.character(as.numeric(sim)+i*10))
}

allcover_sir$pomp <- bind_rows(templist)

templist <- list()
for (i in 0:9) {
	fn <- paste0("../sir_fit/pomp_renewal_fit_fine_", i, ".rda")
	
	load(fn)
	
	templist[[i+1]] <- fitlist %>%
		bind_rows(.id="sim") %>%
		mutate(sim=as.character(as.numeric(sim)+i*10))
}

allcover_sir$renewal_pomp <- bind_rows(templist)

load("../sir_fit/tsir_fit.rda")

allcover_sir$tsir <- bind_rows(fitlist, .id="sim")

load("../sir_fit/tsir_fit_raw.rda")

tsir_raw <- bind_rows(fitlist, .id="sim")

alldata_sir <- bind_rows(allcover_sir, .id="type") %>%
	mutate(sim="gillespie")

estdata <- alldata_sir %>%
	mutate(param=ifelse(param=="beta", "R0", param)) %>%
	filter(param %in% c("R0")) %>%
	group_by(param, type, sim) %>%
	mutate(error=mean/2) %>%
	summarize(
		bias=mean(error, na.rm=TRUE),
		RMSE=sqrt(mean(error^2, na.rm=TRUE)),
		coverage=mean(coverage, na.rm=TRUE)
	) %>%
	ungroup %>%
	merge(classify) %>%
	mutate(
		type=factor(type, levels=c("trajectory", "gradient", "tsir", "pomp", "renewal_pomp"),
					labels=c("trajectory\nmatching",
							 "gradient\nmatching",
							 "TSIR",
							 "SMC",
							 "Renewal"))
	)

tsir_raw_summary <- tsir_raw %>%
	filter(param %in% c("beta")) %>%
	mutate(error=mean/2) %>%
	summarize(
		bias=mean(error, na.rm=TRUE),
		RMSE=sqrt(mean(error^2, na.rm=TRUE)),
		coverage=mean(coverage, na.rm=TRUE)
	)
	
g1 <- ggplot(estdata) +
	geom_rect(xmin=-Inf, xmax=Inf, ymin=0.887, ymax=0.983, alpha=0.1) +
	geom_hline(yintercept=0.95, lty=1) +
	geom_point(aes(type, coverage), size=3) +
	geom_point(aes("TSIR", tsir_raw_summary$coverage), size=3, shape=2) +
	scale_y_continuous("Coverage probability", limits=c(0, 1)) +
	scale_shape_manual("Time-scale", values=c(1, 16)) +
	scale_colour_discrete("Estimation") +
	# ggtitle("Gillespie simulation") +
	xlab("")

g2 <- ggplot(estdata) +
	geom_hline(yintercept=1, lty=1) +
	geom_point(aes(type, bias), size=3) +
	geom_point(aes("TSIR", tsir_raw_summary$bias), size=3, shape=2) +
	scale_y_continuous("Relative bias") +
	scale_shape_manual("Time-scale", values=c(1, 16)) +
	scale_colour_discrete("Estimation") +
	ggtitle("") +
	xlab("")

ggsave("bias.pdf", g2, width=4, height=2.5)

ggsave("coverage.pdf", g1, width=4, height=2.5)

