make_pomp_renewal_basis <- function(ell=60, nstep=10) {
	
	states <- paste0("I", 1:ell)
	
	mu <- paste0("(", paste(paste0("I", (ell-nstep+1):ell),collapse = "+"), ")")
	
	rprocess <- paste0(
		"double phi, pSI, moveS, infection, khatsum, Gshape, Gscale, R0, birth;\n",
		"double khat[", ell, "], k[", ell ,"];\n",
		"Gscale = Gvar/Gmean;\n",
		"Gshape = Gmean/Gscale;\n",
		"R0 = b1*B1+b2*B2+b3*B3+b4*B4+b5*B5+b6*B6+b7*B7+b8*B8+b9*B9+b10*B10+b11*B11+b12*B12+b13*B13;\n",
		"birth = rpois(mu * N0 * dt);\n",
		paste(paste0("khat[", 1:ell, "]=pow(",(1:ell/nstep), ",(Gshape - 1)) * exp(-(", 1:ell/nstep, ")/Gscale);"), collapse="\n"),"\n",
		paste0("khatsum=",paste(paste0("khat[", 1:ell, "]"), collapse="+"),
			  ";\n"),
		paste(paste0("k[", 1:ell, "] = R0/N0 * khat[", ell:1, "]/khatsum;"), collapse = "\n"), "\n",
		paste0("phi=", paste(paste0("k[", 1:ell, "]*I", 1:ell), collapse="+"), ";\n"),
		"pSI=1-exp(-(phi+mu*dt));\n",
		"moveS=rbinom(S, pSI);\n",
		"infection=rbinom(moveS, phi/(phi + mu*dt));\n",
		ifelse(nstep > 1, paste(paste0(head(states, -1), "=", tail(states, -1),";"), collapse="\n"), ""), 
		"\n",
		paste0("I", ell, "=infection;\n"),
		paste0("S=birth+S-moveS;\n")
	)
	
	initlz <- paste0(
		"S=nearbyint(N0*S0);\n",
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
		"Tb1=log(b1);\n",
		"Tb2=log(b2);\n",
		"Tb3=log(b3);\n",
		"Tb4=log(b4);\n",
		"Tb5=log(b5);\n",
		"Tb6=log(b6);\n",
		"Tb7=log(b7);\n",
		"Tb8=log(b8);\n",
		"Tb9=log(b9);\n",
		"Tb10=log(b10);\n",
		"Tb11=log(b11);\n",
		"Tb12=log(b12);\n",
		"Tb13=log(b13);\n",
		"TGvar=log(Gvar);\n",
		"TS0=logit(S0);\n",
		"TI0=logit(I0);\n",
		"Trho=logit(rho);\n",
		"Tdisp=log(disp);"
	)
	
	fromEst <- paste0(
		"Tb1=exp(b1);\n",
		"Tb2=exp(b2);\n",
		"Tb3=exp(b3);\n",
		"Tb4=exp(b4);\n",
		"Tb5=exp(b5);\n",
		"Tb6=exp(b6);\n",
		"Tb7=exp(b7);\n",
		"Tb8=exp(b8);\n",
		"Tb9=exp(b9);\n",
		"Tb10=exp(b10);\n",
		"Tb11=exp(b11);\n",
		"Tb12=exp(b12);\n",
		"Tb13=exp(b13);\n",
		"TGvar=exp(Gvar);\n",
		"TS0=expit(S0);\n",
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
		paramnames=c(
			"b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8", "b9", "b10", "b11", "b12", "b13",
			"Gvar", "S0", "I0", "rho", "disp"
		)
	)

	pomp_arg
}
