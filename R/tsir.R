tsir.run <- function(param=c(beta=2, N=1e5, alpha=0.97),
					 S0=1e5-10,
					 I0=10,
					 tmax=20) {
	with(as.list(param), {
		rmat <- matrix(0, ncol=3, nrow=tmax)
		x <- c(S=S0,I=I0)
		
		for (t in 1:tmax) {
			inf <- rnbinom(n=1, mu=beta * x[["S"]] * x[["I"]]^alpha/N, size=x[["I"]] + 1e-10) 
			
			inf <- min(inf, x[["S"]])
			
			x <- c(S=x[["S"]]-inf, I=inf)
			rmat[t,] <- c(t=t, x)	
		}
		
		
		return(setNames(as.data.frame(rmat), c("t", "S", "I")))
	})
}
