library(emdbook)

rprob <- 0.7
theta <- 10

datalist <- vector('list', 100)

set.seed(101)
for (i in 0:9) {
	fn <- paste0("../sim/gillespie_sinusoidal_sim_", i, ".rda")
	
	load(fn)
	
	for (j in 1:10) {
		print(paste(i, j, sep=", "))
		rr <- reslist[[j]]
		rr <- rr[rr$time<260,]
		
		ii <- c(rr$incidence[1], diff(rr$incidence)) ## incidence
		bb <- c(rr$birth[1], diff(rr$birth)) ## birth
		tt <- rr$time
		
		dd <- as.data.frame(table(ceiling(tt[ii==1])))
		colnames(dd) <- c("time", "tincidence")
		dd$time <- as.numeric(as.character(dd$time))
		
		dd$incidence <- rbetabinom(nrow(dd), size=dd$tincidence, prob=rprob, theta=theta)
		dd$birth <- as.data.frame(table(ceiling(tt[bb==1])))[,2]
		dd$biweek <- dd$time %% 26
		dd$biweek[dd$biweek==0] <- 26
		dd$prevalence <- round(approx(x=rr$time, y=rr$prevalence, xout=dd$time)$y) ## true prevalence just because?
		dd$susceptible <- round(approx(x=rr$time, 
									   y=rr$susceptible,
									   xout=dd$time)$y) ## true susceptible just because
		
		datalist[[i*10 + j]] <- dd
		
	}
}

save("datalist", file="gillespie_sinusoidal_data.rda")
