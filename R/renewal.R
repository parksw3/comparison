renewal <- function(R0 = 2, Gmean=1, Gvar=1,
					ell=40, nstep=10,
					reporting=0.7,
					r=10,
					N0=1e5,
					I0=10,
					tmax=20) {
	Gscale <- Gvar/Gmean
	Gshape <- Gmean/Gscale
	
	tstep <- 1:ell
	
	khat <- (tstep/nstep)^(Gshape - 1) * exp(-(tstep/nstep)/Gscale) * 1/nstep
	k <- rep(0, ell)
	
	for (i in 1:ell) {
		k[i] <- R0/N0 * khat[ell-i+1]/sum(khat)
	}
	
	tlength <- tmax * nstep + ell
	
	I <- S <- phi <- pSI <- rep(0, tlength)
	
	I[1:ell] <- 0
	S[1:ell] <- N0
	I[ell] <- I0
	S[ell] <- N0 - I0
	
	phi[1:ell] <- 0
	pSI[1:ell] <- 0
	phi[ell] <- sum(k[1:ell] * I[1:ell])
	pSI[ell] <- 1 - exp(-phi[ell])
	
	for (t in (ell+1):tlength) {
		I[t] <- rbinom(1,S[t-1],pSI[t-1])
		S[t] <- S[t-1] - I[t]
		phi[t] <- sum(k[1:ell] * I[(t-ell+1):t])
		pSI[t] <- 1 - exp(-phi[t])
	}
	
	mu <- p <- incidence <- rep(0, tmax)
	
	for (i in 1:tmax) {
		mu[i] <- sum(I[((ell+1):(ell+nstep)) + (i-1) * nstep]) + 1e-10
		p[i] <- r/(r + mu[i] * reporting)
		incidence[i] <- rnbinom(1, r, p[i])
	}
	
	list(
		I=I,
		mu=mu,
		incidence=incidence
	)
}
