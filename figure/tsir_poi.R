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

load("../sim/gillespie_sim.rda")

nsim <- length(reslist)

N <- 1e5
fitlist <- vector('list', nsim)

for (i in 1:nsim) {
	print(i)
	rr <- reslist[[i]]
	
	ii <- c(rr$incidence[1], diff(rr$incidence))
	tt <- rr$time
	
	dd <- as.data.frame(table(ceiling(tt[ii==1])))
	
	I <- dd$Freq[1:20]
	
	Inew <- tail(I, -1)
	Iprev <- head(I, -1)
	
	S <- head(N - cumsum(I) - 10, -1)
	
	prev <- approx(x=rr$time, y=rr$prevalence+10, xout=1:19)$y
	
	fitlist[[i]] <- data.frame(
		incidence=Iprev,
		Inew=Inew,
		S=S,
		prevalence=prev,
		POI=Inew/S
	)
}

ff <- fitlist %>% 
	bind_rows(.id="sim")

lfit <- glm(log(Inew) ~ log(incidence) + offset(log(S)), data=ff)

plot(lfit, 2)

g1 <- ggplot(ff) +
	geom_point(aes(prevalence, incidence), shape=".") +
	geom_abline(intercept=0, slope=1, lty=2, lwd=1) + 
	scale_x_continuous("Prevalence", expand=c(0, 0)) +
	scale_y_continuous("Incidence", expand=c(0, 0)) +
	theme(
		legend.position = "none",
		panel.grid = element_blank()
	)

g2 <- ggplot(ff) +
	geom_point(aes(prevalence, POI), shape=".")  +
	stat_function(fun=function(x) 1 - exp(- 2 * x/N), aes(col="Hazard"), lwd=1) +
	stat_function(fun=function(x) 3.899021 * x^0.9226/N, aes(col="TSIR"), lwd=1) +
	scale_x_continuous("Prevalence", expand=c(0, 0)) +
	scale_y_continuous("Probability of infection", expand=c(0, 0)) +
	theme(
		legend.position = "none",
		panel.grid = element_blank()
	)

g3 <- ggplot(ff) +
	geom_point(aes(incidence, POI), shape=".") +
	stat_function(fun=function(x) 1 - exp(- 2 * x/N), aes(col="Hazard"), lwd=1) +
	stat_function(fun=function(x) 3.899021 * x^0.9226/N, aes(col="TSIR"), lwd=1) +
	scale_x_continuous("Incidence", expand=c(0, 0)) +
	scale_y_continuous("Probability of infection", expand=c(0, 0)) +
	theme(
		legend.title = element_blank(),
		legend.position = c(0.8, 0.15), 
		panel.grid = element_blank()
	)

gtot <- arrangeGrob(g1, g2, g3, nrow=1)

ggsave("tsir_poi.pdf", gtot, width=8, height=3)
