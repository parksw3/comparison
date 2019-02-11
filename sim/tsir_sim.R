source("../R/gillespie.R")

nsim <- 100
reslist <- vector('list', nsim)

set.seed(101 )
for (i in 1:nsim) {
	print(i)
	gg <- tsir.run()
	
	dd <- gg
	
	reslist[[i]] <- dd
}

save("reslist", file="tsir_sim.rda")
