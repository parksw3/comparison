source("../R/gillespie.R")

nsim <- 100
reslist <- vector('list', nsim)

set.seed(101 )
for (i in 1:nsim) {
	print(i)
	gg <- gillespie.run(itmax=2e5, ret="all")
	
	dd <- data.frame(
		time=gg[,1],
		incidence=cumsum(gg[,2]==1),
		prevalence=cumsum((gg[,2]==1) - (gg[,2]==2)),
		recovery=cumsum(gg[,2]==2)
	)
	
	reslist[[i]] <- dd
}

save("reslist", file="gillespie_sim.rda")
