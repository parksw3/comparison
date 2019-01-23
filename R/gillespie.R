## adapted from https://rpubs.com/bbolker/SIRgillespie
ratefun <- function(x,param) {
	with(as.list(c(x,param)),{
		c(
			inf=beta*S*I/N,  ## scale infection by pop size
			recover=gamma*I)
	})
}

transfun <- function(x,w) {
	switch(w,
		   x + c(-1,1),
		   x + c(0,-1)
	)   
}

gillespie.run <- function(param=c(beta=2, gamma=1, N=1e5),
						  S0=1e5-10,
						  I0=10,
						  itmax=1e5,
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
	while (x["I"]>0 & it<=itmax) {
		r <- ratefun(x,param)
		t <- t+rexp(1,rate=sum(r))
		w <- sample(length(r),size=1,prob=r)
		x <- transfun(x,w)
		if (ret=="all") rmat[it,] <- c(t,w)
		it <- it+1
	}
	
	if (ret=="all") return(rmat[!is.na(rmat[,1]),])
	return(c(x,t=t,it=it))
}
