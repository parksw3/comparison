make_pomp_renewal <- function(ell=60, nstep=10) {
	
	states <- paste0("I", 1:ell)
	
	mu <- paste0("(", paste(paste0("I", (ell-nstep+1):ell),collapse = "+"), ")")
	
	rprocess <- paste0(
		"double phi, pSI, infection, khatsum, Gshape, Gscale;\n",
		"double khat[", ell, "], k[", ell ,"];\n",
		"Gscale = Gvar/Gmean;\n",
		"Gshape = Gmean/Gscale;\n",
		paste(paste0("khat[", 1:ell, "]=pow(",(1:ell/nstep), ",(Gshape - 1)) * exp(-(", 1:ell/nstep, ")/Gscale);"), collapse="\n"),"\n",
		paste0("khatsum=",paste(paste0("khat[", 1:ell, "]"), collapse="+"),
			  ";\n"),
		paste(paste0("k[", 1:ell, "] = R0/N0 * khat[", ell:1, "]/khatsum;"), collapse = "\n"), "\n",
		paste0("phi=", paste(paste0("k[", 1:ell, "]*I", 1:ell), collapse="+"), ";\n"),
		"pSI=1-exp(-phi);\n",
		"infection=rbinom(S, pSI);\n",
		ifelse(nstep > 1, paste(paste0(head(states, -1), "=", tail(states, -1),";"), collapse="\n"), ""), 
		"\n",
		paste0("I", ell, "=infection;\n"),
		paste0("S=S-infection;\n")
	)
	
	initlz <- paste0(
		"S=nearbyint(N0*(1-I0));\n",
		paste(paste0("I", 1:(ell-1), "=0;"), collapse="\n"), "\n",
		"I", ell,"=nearbyint(N0*I0);\n"
	) 
	
	dmeas <- paste0(
		"lik = dnbinom_mu(incidence, disp, ", mu, "*rho, give_log);"
	)
	
	rmeas <- paste0(
		"incidence = rnbinom_mu(disp, ", mu, "*rho);"
	)
	
	toEst <- paste0(
		"TR0=log(R0);\n",
		"TGvar=log(Gvar);\n",
		"TI0=logit(I0);\n",
		"Trho=logit(rho);\n",
		"Tdisp=log(disp);"
	)
	
	fromEst <- paste0(
		"TR0=exp(R0);\n",
		"TGvar=exp(Gvar);\n",
		"TI0=expit(I0);\n",
		"Trho=expit(rho);\n",
		"Tdisp=exp(disp);"
	)
	
	pomp_arg <- list(
		times="time",
		rprocess=pomp::euler.sim(pomp::Csnippet(rprocess), delta.t=1/nstep),
		dmeasure=pomp::Csnippet(dmeas),
		rmeasure=pomp::Csnippet(rmeas),
		initializer = pomp::Csnippet(initlz),
		toEstimationScale=pomp::Csnippet(toEst),
		fromEstimationScale=pomp::Csnippet(fromEst),
		statenames=c("S", states),
		paramnames=c("R0", "Gvar", "I0", "rho", "disp")
	)

	pomp_arg
}

