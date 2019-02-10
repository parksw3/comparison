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
	type=c("trajectory", "gradient", "tsir", "pomp"),
	time=c("Continuous", "Continuous", "Discrete", "Discrete"),
	est=c("Simulation", "Regression", "Regression", "Simulation")
)

allcover <- list()

load("../sir_fit/trajectory_fit.rda")

allcover$trajectory <- bind_rows(fitlist, .id="sim")

load("../sir_fit/gradient_fit.rda")

allcover$gradient <- bind_rows(fitlist, .id="sim")

templist <- list()
for (i in 0:9) {
	fn <- paste0("../sir_fit/pomp_0_1_fit_", i, ".rda")
	
	load(fn)
	
	templist[[i+1]] <- fitlist %>% 
		bind_rows(.id="sim") %>%
		mutate(sim=as.character(as.numeric(sim)+i*10))
}

allcover$pomp <- bind_rows(templist)

load("../sir_fit/tsir_fit.rda")

allcover$tsir <- bind_rows(fitlist, .id="sim")

alldata <- bind_rows(allcover, .id="type")

estdata <- alldata %>%
	filter(param %in% c("beta")) %>%
	group_by(param, type) %>%
	mutate(error=log(mean/2)) %>%
	summarize(
		bias=mean(error),
		RMSE=sqrt(mean(error^2)),
		coverage=mean(coverage)
	) %>%
	ungroup %>%
	merge(classify) %>%
	mutate(
		type=factor(type, levels=c("trajectory", "gradient", "tsir", "pomp"),
					labels=c("trajectory\nmatching",
							 "gradient\nmatching",
							 "tSIR",
							 "POMP"))
	)

g1 <- ggplot(estdata) +
	geom_rect(xmin=-Inf, xmax=Inf, ymin=0.887, ymax=0.983, alpha=0.1) +
	geom_hline(yintercept=0.95, lty=1) +
	geom_point(aes(type, coverage, shape=time, col=est), size=5) +
	scale_y_continuous("Coverage probability", limits=c(0, 1)) +
	scale_shape_manual("Time-scale", values=c(1, 16)) +
	scale_colour_discrete("Estimation") +
	xlab("")

g2 <- ggplot(estdata) +
	geom_hline(yintercept=0, lty=1) +
	geom_point(aes(type, bias, shape=time, col=est), size=5) +
	scale_y_continuous("Bias") +
	scale_shape_manual("Time-scale", values=c(1, 16)) +
	scale_colour_discrete("Estimation") +
	xlab("")
	
g3 <- ggplot(estdata) +
	geom_point(aes(type, RMSE, shape=time, col=est), size=5) +
	scale_y_continuous("RMSE") +
	scale_shape_manual("Time-scale", values=c(1, 16)) +
	scale_colour_discrete("Estimation") +
	xlab("")

g_legend<-function(a.gplot) {
	tmp <- ggplot_gtable(ggplot_build(a.gplot))
	leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
	legend <- tmp$grobs[[leg]]
	return(legend)
}

mylegend<-g_legend(g1)

nl <- function(x) x + theme(legend.position = "none")

gcomp <- arrangeGrob(arrangeGrob(nl(g1), nl(g2), nl(g3), nrow=1), mylegend,
					 nrow=1, widths=c(10, 1))

ggsave("compare_estimate_sir.pdf", gcomp, width=12, height=4)
