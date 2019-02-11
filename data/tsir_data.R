library(emdbook)
load("../sim/tsir_sim.rda")

nsim <- length(reslist)
rprob <- 0.7
theta <- 10

datalist <- vector('list', nsim)

set.seed(101)
for (i in 1:nsim) {
	rr <- reslist[[i]]
	
	dd <- rr
	
	dd$incidence <- rbetabinom(nrow(dd), size=dd$I, prob=rprob, theta=theta)
	
	datalist[[i]] <- dd
}

save("datalist", file="tsir_data.rda")
