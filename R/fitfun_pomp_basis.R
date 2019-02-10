rprocess <- pomp::Csnippet("
	double beta, birth;
	double rate[4], trans[4];

	beta = b1*B1+b2*B2+b3*B3+b4*B4+b5*B5+b6*B6+b7*B7+b8*B8+b9*B9+b10*B10+b11*B11+b12*B12+b13*B13;
	
	birth = rpois(mu * N * dt);

	rate[0] = beta * I/N;
	rate[1] = mu;
	rate[2] = Gamma;
	rate[3] = mu;

	// transitions between classes
  	reulermultinom(2,S,&rate[0],dt,&trans[0]);
	reulermultinom(2,I,&rate[2],dt,&trans[2]);

	S += birth - trans[0] - trans[1];
	I += trans[0] - trans[2] - trans[3];
	C += trans[0];
")

initlz <- pomp::Csnippet("
	S = nearbyint(N*S0);
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
	Tb1=log(b1);
	Tb2=log(b2);
	Tb3=log(b3);
	Tb4=log(b4);
	Tb5=log(b5);
	Tb6=log(b6);
	Tb7=log(b7);
	Tb8=log(b8);
	Tb9=log(b9);
	Tb10=log(b10);
	Tb11=log(b11);
	Tb12=log(b12);
	Tb13=log(b13);
	TS0=logit(S0);
	TI0=logit(I0);
	Trho=logit(rho);
	Tdisp=log(disp);
")

fromEst <- pomp::Csnippet("
	Tb1=exp(b1);
	Tb2=exp(b2);
	Tb3=exp(b3);
	Tb4=exp(b4);
	Tb5=exp(b5);
	Tb6=exp(b6);
	Tb7=exp(b7);
	Tb8=exp(b8);
	Tb9=exp(b9);
	Tb10=exp(b10);
	Tb11=exp(b11);
	Tb12=exp(b12);
	Tb13=exp(b13);
	TS0=expit(S0);
	TI0=expit(I0);
	Trho=expit(rho);
	Tdisp=exp(disp);
")

pomp_arg <- list(
	times="time",
	rprocess=pomp::euler.sim(rprocess, delta.t=1/10),
	dmeasure=dmeas,
	rmeasure=rmeas,
	initializer = initlz,
	tcovar="time",
	toEstimationScale=toEst,
	fromEstimationScale=fromEst,
	statenames=c("S", "I", "C"),
	zeronames=c("C"),
	paramnames=c("b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8", "b9", "b10", "b11", "b12", "b13",
				 "S0", "I0", "rho", "disp")
)
