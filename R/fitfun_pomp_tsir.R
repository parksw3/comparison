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
	lik = dbetabinom(cases, I, rho, disp, 0);
")

rmeas <- pomp::Csnippet("
	cases = rbetabinom(I, rho, disp);
")

pomp_arg <- list(
	times="time",
	rprocess=discrete.time.sim(rprocess),
	dmeasure=dmeas,
	rmeasure=rmeas,
	initializer = initlz,
	tcovar="index",
	statenames=c("S", "I"),
	paramnames=c("b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8", "b9", "b10", "b11",
		     	 "b12", "b13", "b14", "b15", "b16", "b17", "b18", "b19", "b20", "b21",
			     "b22", "b23", "b24", "b25", "b26", "alpha", "m", "S0", "I0", "rho", "disp")
)
