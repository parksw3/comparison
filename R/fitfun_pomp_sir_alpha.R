rprocess <- pomp::Csnippet("
	double infection, recovery;
	
	infection = rbinom(S, 1 - exp(-beta*pow(I,alpha)/N*dt));	

	recovery = rbinom(I, (1-exp(-Gamma*dt)));

	S += - infection;
	I += infection - recovery;
	C += infection;
")

initlz <- pomp::Csnippet("
	S = nearbyint(N*(1-I0));
	I = nearbyint(N*I0);
	C = 0;
")

dmeas <- pomp::Csnippet("
	lik = dnbinom_mu(incidence, disp, C*rho, give_log);
")

rmeas <- pomp::Csnippet("
	incidence = rnbinom_mu(disp, C*rho);
")

toEst <- pomp::Csnippet("
	Tbeta=log(beta);
	TI0=logit(I0);
	Trho=logit(rho);
	Tdisp=log(disp);
	Talpha=log(alpha);
")

fromEst <- pomp::Csnippet("
	Tbeta=exp(beta);
	TI0=expit(I0);
	Trho=expit(rho);
	Tdisp=exp(disp);
	Talpha=exp(alpha);
")

pomp_arg <- list(
	times="time",
	rprocess=pomp::euler.sim(rprocess, delta.t=1/10),
	dmeasure=dmeas,
	rmeasure=rmeas,
	initializer = initlz,
	toEstimationScale=toEst,
	fromEstimationScale=fromEst,
	statenames=c("S", "I", "C"),
	zeronames=c("C"),
	paramnames=c("beta", "I0", "rho", "disp", "alpha")
)
