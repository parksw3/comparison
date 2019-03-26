make_pomp_tsir_lognormal_process <- function(density=FALSE) {
	foi <- ifelse(density, "foi = Beta * S * pow(I, alpha)/N;\n", "foi = Beta * S * pow(I, alpha);\n")
	
	rprocess <- paste0(
		"double foi, infection;\n",
		foi,
		"infection = rnorm(log(foi), sigma_p);\n",
		"if (exp(infection) > S) infection = log(S - 1);\n",
		"logI=infection;\n",
		"I = exp(logI);\n",
		"S = B + S - I;"
	)
	
	initlz <- paste0(
		"S=S0;\n",
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
		statenames=c("I", "logI", "S"),
		paramnames=c("alpha", "sigma_p", "sigma_o")
	)

	pomp_arg
}
