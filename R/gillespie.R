## adapted from https://rpubs.com/bbolker/SIRgillespie
ratefun <- function(x,param,t,
					betafun,
					rhofun) {
	beta <- betafun(t, param)
	
	birth <- rhofun(x, t, param)
	
	with(as.list(c(x,param)),{
		c(
			birth=birth,
			inf=beta*S*(I+nu)/N,  ## scale infection by pop size
			prog=sigma*E,
			recover=gamma*I)
	})
}

transfun <- function(x,w) {
	switch(w,
		   x + c(1,0,0,0),
		   x + c(-1,1,0,0),
		   x + c(0,-1,1,0), 
		   x + c(0,0,-1,1)
	)   
}

gillespie.run <- function(trans.param,
				epi.param=c(sigma=365/8, gamma=365/5, nu=1e-6,N=1e7),
				S0=0.035 * 1e7,
				E0=0,
				I0=100,
				R0=0,
				itmax=1e3,
				ret=c("final","all"),
				betafun, rhofun,
				seed) {
	ret <- match.arg(ret)
	if (ret=="all") {
		rmat <- matrix(NA,nrow=itmax,ncol=2,
					   dimnames=list(NULL,c("t","trans")))
	}
	
	if (!missing(seed)) set.seed(seed)
	
	x <- c(S=S0,E=E0,I=I0,R=R0)
	it <- 1
	t <- 0
	trans <- c(0,0,0,0)
	param <- c(trans.param, epi.param)
	while (x["I"]>0 & it<=itmax) {
		r <- ratefun(x,param,t,betafun,rhofun)
		t <- t+rexp(1,rate=sum(r))
		w <- sample(length(r),size=1,prob=r)
		x <- transfun(x,w)
		if (ret=="all") rmat[it,] <- c(t,w)
		it <- it+1
	}
	
	if (ret=="all") return(rmat[!is.na(rmat[,1]),])
	return(c(x,t=t,it=it))
}
