make_pomp_tsir <- function(process=c("deterministic", "poisson", "nbinom"),
						   density=FALSE) {
	process <- match.arg(process)
	
	foi <- ifelse(density, "foi = Beta * S * pow(I, alpha)/N;\n", "foi = Beta * S * pow(I, alpha);\n")
	
	pp <- switch (process,
		deterministic = "infection = foi;\n",
		poisson = "infection = rpois(foi);\n",
		nbinom = "infection = rnbinom_mu(I+1e-10, foi);\n"
	)
	
	rprocess <- paste0(
		"double foi, infection;\n",
		foi,
		pp,
		"if (infection > S) infection = S;\n",
		"S = B + S - infection;\n",
		"I = infection;"
	)
	
	initlz <- paste0(
		"S=nearbyint(S0);\n",
		"I=nearbyint(I0);\n"
	) 
	
	dmeas <- paste0(
		"lik = dnbinom_mu(cases, disp, I*rho, give_log);"
	)
	
	rmeas <- paste0(
		"cases = rnbinom_mu(disp, I*rho);"
	)
	
	toEst <- paste0(
		"Tdisp=log(disp);\n",
		"Talpha=log(alpha);\n"
	)
	
	fromEst <- paste0(
		"Tdisp=exp(disp);\n",
		"Talpha=exp(alpha);\n"
	)
	
	pomp_arg <- list(
		times="time",
		rprocess=pomp::discrete.time.sim(pomp::Csnippet(rprocess), delta.t=1),
		dmeasure=pomp::Csnippet(dmeas),
		rmeasure=pomp::Csnippet(rmeas),
		initializer = pomp::Csnippet(initlz),
		toEstimationScale=pomp::Csnippet(toEst),
		fromEstimationScale=pomp::Csnippet(fromEst),
		statenames=c("S", "I"),
		paramnames=c("disp", "alpha")
	)

	pomp_arg
}
