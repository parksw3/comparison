rprocess <- pomp::Csnippet("
	double beta, lambda, infection, theta, births;
	
	if (biweek == 1) 
		beta = b1;
	else if (biweek == 2)
		beta = b2;
	else if (biweek == 3)
		beta = b3;
	else if (biweek == 4)
		beta = b4;
	else if (biweek == 5)
		beta = b5;
	else if (biweek == 6)
		beta = b6;
	else if (biweek == 7)
		beta = b7;
	else if (biweek == 8)
		beta = b8;
	else if (biweek == 9)
		beta = b9;
	else if (biweek == 10)
		beta = b10;
	else if (biweek == 11)
		beta = b11;
	else if (biweek == 12)
		beta = b12;
	else if (biweek == 13)
		beta = b13;
	else if (biweek == 14)
		beta = b14;
	else if (biweek == 15)
		beta = b15;
	else if (biweek == 16)
		beta = b16;
	else if (biweek == 17)
		beta = b17;
	else if (biweek == 18)
		beta = b18;
	else if (biweek == 19)
		beta = b19;
	else if (biweek == 20)
		beta = b20;
	else if (biweek == 21)
		beta = b21;
	else if (biweek == 22)
		beta = b22;
	else if (biweek == 23)
		beta = b23;
	else if (biweek == 24)
		beta = b24;
	else if (biweek == 25)
		beta = b25;
	else if (biweek == 26)
		beta = b26;
	
	theta = rpois(m);
	births = rpois(rec);

	if (I == 0) {
		lambda=0;
		infection=0;
	} else {
		lambda = beta * S * pow((I + theta), alpha)/pop;
		infection = rnbinom_mu(I, lambda);
	}

	if (infection > S) infection = S;

	S = S + births - infection;
	I = infection;
")

initlz <- pomp::Csnippet("
	S = nearbyint(S0);
	I = nearbyint(I0);
")

dmeas <- pomp::Csnippet("
	lik = dnbinom_mu(cases, disp, I*rho, give_log);
")

rmeas <- pomp::Csnippet("
	cases = rnbinom_mu(disp, I*rho);
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
	Tb14=log(b14);
	Tb15=log(b15);
	Tb16=log(b16);
	Tb17=log(b17);
	Tb18=log(b18);
	Tb19=log(b19);
	Tb20=log(b20);
	Tb21=log(b21);
	Tb22=log(b22);
	Tb22=log(b22);
	Tb23=log(b23);
	Tb24=log(b24);
	Tb25=log(b25);
	Tb26=log(b26);
	Talpha=log(alpha);
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
	Tb14=exp(b14);
	Tb15=exp(b15);
	Tb16=exp(b16);
	Tb17=exp(b17);
	Tb18=exp(b18);
	Tb19=exp(b19);
	Tb20=exp(b20);
	Tb21=exp(b21);
	Tb22=exp(b22);
	Tb22=exp(b22);
	Tb23=exp(b23);
	Tb24=exp(b24);
	Tb25=exp(b25);
	Tb26=exp(b26);
	Talpha=exp(alpha);
	Tm=exp(m);
	TS0=exp(S0);
	TI0=exp(I0);
	Trho=expit(rho);
	Tdisp=exp(disp);
")

pomp_arg <- list(
	times="time",
	rprocess=discrete.time.sim(rprocess),
	dmeasure=dmeas,
	rmeasure=rmeas,
	initializer = initlz,
	tcovar="index",
	toEstimationScale=toEst,
	fromEstimationScale=fromEst,
	statenames=c("S", "I"),
	paramnames=c("b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8", "b9", "b10", "b11",
		     	 "b12", "b13", "b14", "b15", "b16", "b17", "b18", "b19", "b20", "b21",
			     "b22", "b23", "b24", "b25", "b26", "alpha", "m", "S0", "I0", "rho", "disp")
)
