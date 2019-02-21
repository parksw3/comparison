library(tidyr)
library(dplyr)
library(gridExtra)
library(ggplot2); theme_set(theme_bw(base_size = 12,
									 base_family = "Times"))

if (.Platform$OS.type=="windows") {
	windowsFonts(Times=windowsFont("Times"))
} 

## use Dark2 rather than Set1 because colour #6 of Set1 is yellow (ugh/too light)
scale_colour_discrete <- function(...,palette="Dark2") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Dark2") scale_fill_brewer(...,palette=palette)

classify <- data.frame(
	type=c("trajectory", "gradient", "tsir", "pomp", "renewal_jags"),
	time=c("Continuous", "Continuous", "Discrete", "Discrete", "Discrete"),
	est=c("Simulation", "Regression", "Regression", "Simulation", "Simulation")
)

allcover_sir <- list()
allcover_tsir <- list()

load("../sir_fit/trajectory_fit.rda")

allcover_sir$trajectory <- bind_rows(fitlist, .id="sim")

load("../sir_fit/gradient_fit.rda")

allcover_sir$gradient <- bind_rows(fitlist, .id="sim")

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
	fn <- paste0("../sir_fit/jags_fit_", i, ".rda")
	
	load(fn)
	
	templist[[i+1]] <- fitlist %>%
		bind_rows(.id="sim") %>%
		mutate(sim=as.character(as.numeric(sim)+i*10))
}

allcover_sir$renewal_jags <- bind_rows(templist)

load("../sir_fit/tsir_fit.rda")

allcover_sir$tsir <- bind_rows(fitlist, .id="sim")

## load TSIR fits

load("../tsir_fit/tsir_trajectory_fit.rda")

allcover_tsir$trajectory <- bind_rows(fitlist, .id="sim")

load("../tsir_fit/tsir_gradient_fit.rda")

allcover_tsir$gradient <- bind_rows(fitlist, .id="sim")

templist <- list()
for (i in 0:9) {
	fn <- paste0("../tsir_fit/tsir_pomp_0_1_fit_", i, ".rda")
	
	load(fn)
	
	templist[[i+1]] <- fitlist %>% 
		bind_rows(.id="sim") %>%
		mutate(sim=as.character(as.numeric(sim)+i*10))
}

allcover_tsir$pomp <- bind_rows(templist)

load("../tsir_fit/tsir_tsir_fit.rda")

allcover_tsir$tsir <- bind_rows(fitlist, .id="sim")

alldata_sir <- bind_rows(allcover_sir, .id="type") %>%
	mutate(sim="gillespie")
alldata_tsir <- bind_rows(allcover_tsir, .id="type") %>%
	mutate(sim="tsir")

alldata <- rbind(alldata_sir, alldata_tsir)

estdata <- alldata %>%
	mutate(param=ifelse(param=="beta", "R0", param)) %>%
	filter(param %in% c("R0")) %>%
	group_by(param, type, sim) %>%
	mutate(error=log(mean/2)) %>%
	summarize(
		bias=mean(error),
		RMSE=sqrt(mean(error^2)),
		coverage=mean(coverage)
	) %>%
	ungroup %>%
	merge(classify) %>%
	mutate(
		type=factor(type, levels=c("trajectory", "gradient", "tsir", "pomp", "renewal_jags"),
					labels=c("trajectory\nmatching",
							 "gradient\nmatching",
							 "tSIR",
							 "POMP",
							 "Renewal\n(JAGS)"))
	)

g1 <- ggplot(estdata) +
	geom_rect(xmin=-Inf, xmax=Inf, ymin=0.887, ymax=0.983, alpha=0.1) +
	geom_hline(yintercept=0.95, lty=1) +
	geom_point(aes(type, coverage, shape=time, col=est), size=5) +
	scale_y_continuous("Coverage probability", limits=c(0, 1)) +
	scale_shape_manual("Time-scale", values=c(1, 16)) +
	scale_colour_discrete("Estimation") +
	ggtitle("Gillespie simulation") +
	xlab("")

g2 <- ggplot(estdata) +
	geom_hline(yintercept=0, lty=1) +
	geom_point(aes(type, bias, shape=time, col=est), size=5) +
	scale_y_continuous("Bias") +
	scale_shape_manual("Time-scale", values=c(1, 16)) +
	scale_colour_discrete("Estimation") +
	ggtitle("") +
	xlab("")
	
g3 <- ggplot(estdata) +
	geom_hline(yintercept=0, lty=1) +
	geom_point(aes(type, RMSE, shape=time, col=est), size=5) +
	scale_y_continuous("RMSE") +
	scale_shape_manual("Time-scale", values=c(1, 16)) +
	scale_colour_discrete("Estimation") +
	ggtitle("") +
	xlab("")

g_legend<-function(a.gplot) {
	tmp <- ggplot_gtable(ggplot_build(a.gplot))
	leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
	legend <- tmp$grobs[[leg]]
	return(legend)
}

mylegend<-g_legend(g1)

nl <- function(x) x + theme(legend.position = "none")

gcomp <- arrangeGrob(
	arrangeGrob(nl(g1) %+% filter(estdata, sim=="gillespie"), 
				nl(g2) %+% filter(estdata, sim=="gillespie"), 
				nl(g3) %+% filter(estdata, sim=="gillespie"),
				nl(g1) %+% filter(estdata, sim=="tsir") + ggtitle("tSIR simulation"), 
				nl(g2) %+% filter(estdata, sim=="tsir"), 
				nl(g3) %+% filter(estdata, sim=="tsir"),
				nrow=2), 
	mylegend,
	nrow=1, widths=c(10, 1)
)

ggsave("compare_estimate_sir.pdf", gcomp, width=12, height=8)
