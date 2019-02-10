source("../R/gillespie_sinusoidal.R")

argvals <- commandArgs(trailingOnly=TRUE)
batch_num <- as.numeric(argvals[1])

fn <- paste0("gillespie_sinusoidal_sim_", batch_num, ".rda")

nsim <- 10
reslist <- vector('list', nsim)

set.seed(batch_num)
for (i in 1:nsim) {
	print(i)
	gg <- gillespie.run(itmax=1e8, ret="all")
	
	dd <- data.frame(
		time=gg[,1],
		birth=cumsum(gg[,2]==1),
		incidence=cumsum(gg[,2]==2),
		susceptible=cumsum((gg[,2]==1) - (gg[,2]==2) - (gg[,2]==4)) + 0.05 * 5e6,
		prevalence=cumsum((gg[,2]==2) - (gg[,2]==3) - (gg[,2]==5)) + 1e-4 * 5e6,
		recovery=cumsum(gg[,2]==3)
	)
	
	reslist[[i]] <- dd
	
	save("reslist", file=fn)
}

save("reslist", file=fn)
