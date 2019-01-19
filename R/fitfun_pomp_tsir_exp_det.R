rprocess <- pomp::Csnippet("
	double beta, lambda, infection, theta, births, recovery;
	
	beta = b1*B1+b2*B2+b3*B3+b4*B4+b5*B5+b6*B6+b7*B7+b8*B8+b9*B9+b10*B10+b11*B11+b12*B12+b13*B13;

	births = rec * dt;

	if (I == 0) {
		lambda=0;
		infection=0;
	} else {
		lambda = beta * S * pow(I, alpha)/pop * dt;
		infection = lambda;
	}

	recovery = I * (1-exp(-gamma*dt));

	if (infection > S) infection = S;

	S += births - infection;
	I += infection - recovery;
	C += infection;
")

initlz <- pomp::Csnippet("
	S = nearbyint(S0);
	I = nearbyint(I0);
	C = 0;
")

dmeas <- pomp::Csnippet("
	lik = dnbinom_mu(cases, disp, C*rho, give_log);
")

rmeas <- pomp::Csnippet("
	cases = rnbinom_mu(disp, C*rho);
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
						Talpha=log(alpha);
						Tgamma=log(gamma);
						Tm=log(m);
						TS0=log(S0);
						TI0=log(I0);
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
						  Talpha=exp(alpha);
						  Tgamma=exp(gamma);
						  Tm=exp(m);
						  TS0=exp(S0);
						  TI0=exp(I0);
						  Trho=expit(rho);
						  Tdisp=exp(disp);
						  ")

pomp_arg <- list(
	times="time",
	rprocess=pomp::euler.sim(rprocess, delta.t=1/2),
	dmeasure=dmeas,
	rmeasure=rmeas,
	initializer = initlz,
	tcovar="index",
	toEstimationScale=toEst,
	fromEstimationScale=fromEst,
	statenames=c("S", "I", "C"),
	zeronames=c("C"),
	paramnames=c("b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8", "b9", "b10", "b11",
		     	 "b12", "b13", "alpha", "gamma", "m", "S0", "I0", "rho", "disp")
)
