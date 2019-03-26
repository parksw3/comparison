make_pomp_tsir_lognormal <- function(density=FALSE) {
	foi <- ifelse(density, "foi = Beta * S * pow(I, alpha)/N;\n", "foi = Beta * S * pow(I, alpha);\n")
	
	rprocess <- paste0(
		"double foi, infection;\n",
		foi,
		"infection = rnorm(log(foi), sigma_p);\n",
		"logI=infection;\n",
		"I = exp(logI);\n"
	)
	
	initlz <- paste0(
		"I=I0;\n",
		"logI=log(I0);\n"
	) 
	
	dmeas <- paste0(
		"lik = dnorm(log(cases), logI + log(rho), sigma_o, give_log);"
	)
	
	rmeas <- paste0(
		"cases = exp(rnorm(logI + log(rho), sigma_o));"
	)
	
	pomp_arg <- list(
		times="time",
		rprocess=pomp::discrete.time.sim(pomp::Csnippet(rprocess), delta.t=1),
		dmeasure=pomp::Csnippet(dmeas),
		rmeasure=pomp::Csnippet(rmeas),
		initializer = pomp::Csnippet(initlz),
		statenames=c("I", "logI"),
		paramnames=c("alpha", "sigma_p", "sigma_o")
	)

	pomp_arg
}
