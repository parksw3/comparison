library(emdbook)
load("../sim/gillespie_sim.rda")

nsim <- length(reslist)
rprob <- 0.7
theta <- 10

datalist <- vector('list', nsim)

set.seed(101)
for (i in 1:nsim) {
	rr <- reslist[[i]]
	ii <- c(rr$incidence[1], diff(rr$incidence))
	tt <- rr$time
	
	dd <- as.data.frame(table(ceiling(tt[ii==1])))
	colnames(dd) <- c("time", "incidence")
	dd$time <- as.numeric(as.character(dd$time))
	
	dd$incidence <- rbetabinom(nrow(dd), size=dd$incidence, prob=rprob, theta=theta)
	
	datalist[[i]] <- dd
}

save("datalist", file="gillespie_data.rda")
