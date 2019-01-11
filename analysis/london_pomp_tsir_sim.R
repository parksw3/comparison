data <- ld
nsim <- 10
res <- matrix(0,length(data$cases),nsim)
Sres <- matrix(0,length(data$cases),nsim)
period <- ld$biweek
alpha <- coef(m)[["alpha"]]
beta <- coef(m)[1:26]
pop <- ld$pop

S_start <- coef(m)[["S0"]]
I_start <- coef(m)[["I0"]]

method <- "negbin"
epidemics <- "cont"
add.noise.sd <- 0
adj.rho <- 1/coef(m)[["rho"]]

for(ct in 1:nsim){
	
	S <- rep(0,length(data$cases))
	I <- rep(0,length(data$cases))
	
	S[1] <- S_start
	I[1] <- I_start
	
	for (t in 2:(nrow(data))){
		
		lambda <- min(S[t-1],unname(beta[period[t-1]] * S[t-1] * (I[t-1])^alpha/pop[t-1]))
		
		#if(lambda < 1 || is.nan(lambda) == T){lambda <- 0}
		if(is.nan(lambda) == T){lambda <- 0}
		
		if(method == 'deterministic'){
			I[t] <- lambda * rnorm( n = 1, mean = 1, sd=mul.noise.sd)
			if(I[t] < 0 && lambda >= 0 ){
				warning('infected overflow  -- reduce multiplicative noise sd')
			}
		}
		if(method == 'negbin'){
			I[t] <- rnbinom(n=1,mu=lambda,size=I[t-1]+1e-10)
		}
		if(method == 'pois'){
			I[t] <- rpois(n=1,lambda=lambda)
		}
		if(epidemics == 'cont'){
			I[t] <- I[t]
		}
		if(epidemics == 'break'){
			
			t0s <- epitimes(data,threshold)$start
			if(t %in% t0s){
				I[t] <- adj.rho[t]*data$cases[t]
			}
		}
		S[t] <- max(S[t-1] + data$rec[t-1] - I[t] + rnorm(n=1,mean=0,sd=add.noise.sd),0)
		
		if(S[t] < 0 && (S[t-1] + data$births[t-1] - I[t]) >0 ){
			warning('susceptible overflow  -- reduce additive noise sd')
		}
	}
	res[,ct] <- I / adj.rho
	Sres[,ct] <- S
}

plot(ld$cases)
lines(res[,3], col=2)
