## adapted from https://rpubs.com/bbolker/SIRgillespie
ratefun <- function(x,param,t) {
	with(as.list(c(x,param)),{
		beta <- b0 * (1 + b1 * cos(t/26 * 2 * pi))
		c(
			birth=mu*N,
			inf=beta*S*I/N,  ## scale infection by pop size
			recover=gamma*I,
			deathS=mu*S,
			deathI=mu*I)
	})
}

transfun <- function(x,w) {
	switch(w,
		   x + c(1, 0),
		   x + c(-1,1),
		   x + c(0,-1),
		   x + c(-1,0),
		   x + c(0,-1)
	)   
}

gillespie.run <- function(param=c(mu=1/(50*26), b0=500/26, b1=0.15, gamma=1, N=5e6),
						  S0=0.05 * 5e6,
						  I0=1e-4 * 5e6,
						  itmax=1e7,
						  tmax=26*10,
						  ret=c("final","all"),
						  seed) {
	ret <- match.arg(ret)
	if (ret=="all") {
		rmat <- matrix(NA,nrow=itmax,ncol=2,
					   dimnames=list(NULL,c("t","trans")))
	}
	
	if (!missing(seed)) set.seed(seed)
	
	x <- c(S=S0,I=I0)
	it <- 1
	t <- 0
	trans <- c(0,0)
	while (x["I"]>0 & it<=itmax & t<tmax) {
		r <- ratefun(x,param,t)
		t <- t+rexp(1,rate=sum(r))
		w <- sample(length(r),size=1,prob=r)
		x <- transfun(x,w)
		if (ret=="all") rmat[it,] <- c(t,w)
		it <- it+1
	}
	
	if (ret=="all") return(rmat[!is.na(rmat[,1]),])
	return(c(x,t=t,it=it))
}
