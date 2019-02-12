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
	type=c("pomp1", "pomp2", "pomp3", "pomp4"),
	delta=c("0.1", "0.2", "0.5", "1.0")
)

allcover_sir <- list()

templist <- list()
for (i in 0:9) {
	fn <- paste0("../sir_fit/pomp_0_1_fit_", i, ".rda")
	
	load(fn)
	
	templist[[i+1]] <- fitlist %>% 
		bind_rows(.id="sim") %>%
		mutate(sim=as.character(as.numeric(sim)+i*10))
}

allcover_sir$pomp1 <- bind_rows(templist)

templist <- list()
for (i in 0:9) {
	fn <- paste0("../sir_fit/pomp_0_2_fit_", i, ".rda")
	
	load(fn)
	
	templist[[i+1]] <- fitlist %>% 
		bind_rows(.id="sim") %>%
		mutate(sim=as.character(as.numeric(sim)+i*10))
}

allcover_sir$pomp2 <- bind_rows(templist)

templist <- list()
for (i in 0:9) {
	fn <- paste0("../sir_fit/pomp_0_5_fit_", i, ".rda")
	
	load(fn)
	
	templist[[i+1]] <- fitlist %>% 
		bind_rows(.id="sim") %>%
		mutate(sim=as.character(as.numeric(sim)+i*10))
}

allcover_sir$pomp3 <- bind_rows(templist)

templist <- list()
for (i in 0:9) {
	fn <- paste0("../sir_fit/pomp_1_0_fit_", i, ".rda")
	
	load(fn)
	
	templist[[i+1]] <- fitlist %>% 
		bind_rows(.id="sim") %>%
		mutate(sim=as.character(as.numeric(sim)+i*10))
}

allcover_sir$pomp4 <- bind_rows(templist)

alldata <- bind_rows(allcover_sir, .id="type") %>%
	merge(classify)

estdata <- alldata %>%
	filter(param %in% c("beta")) %>%
	group_by(param, type, delta) %>%
	mutate(error=log(mean/2)) %>%
	summarize(
		bias=mean(error),
		RMSE=sqrt(mean(error^2)),
		coverage=mean(coverage, na.rm=TRUE)
	) %>%
	ungroup 

g1 <- ggplot(estdata) +
	geom_rect(xmin=-Inf, xmax=Inf, ymin=0.887, ymax=0.983, alpha=0.1) +
	geom_hline(yintercept=0.95, lty=1) +
	geom_point(aes(delta, coverage), size=5, shape=1) +
	scale_y_continuous("Coverage probability", limits=c(0, 1)) +
	scale_shape_manual("Time-scale", values=c(1, 16)) +
	scale_colour_discrete("Estimation") +
	ggtitle("Gillespie simulation") +
	xlab("")

g2 <- ggplot(estdata) +
	geom_hline(yintercept=0, lty=1) +
	geom_point(aes(delta, bias), size=5, shape=1) +
	scale_y_continuous("Bias") +
	scale_shape_manual("Time-scale", values=c(1, 16)) +
	scale_colour_discrete("Estimation") +
	ggtitle("") +
	xlab("")
	
g3 <- ggplot(estdata) +
	geom_hline(yintercept=0, lty=1) +
	geom_point(aes(delta, RMSE), size=5, shape=1) +
	scale_y_continuous("RMSE") +
	scale_shape_manual("Time-scale", values=c(1, 16)) +
	scale_colour_discrete("Estimation") +
	ggtitle("") +
	xlab("")

gcomp <- arrangeGrob(
	g1, g2, g3,
	nrow=1
)

ggsave("compare_estimate_sir_pomp.pdf", gcomp, width=12, height=4)
