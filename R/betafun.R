# trans.param <- c(c1=26, c2=30, c3=27, c4=26, c5=30, c6=28, c7=25, c8=20, c9=20, c10=23, c11=36, c12=35, c13=28) * gamma

# rhofun <- function(x, t, param) 1/50 * (1e7 - x[["S"]])


betafun_cs <- function(t, param) {
	with(as.list(param), {
		knots <- seq(0, 1, by=1/13)
		coef <- c(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13)
		t <- t%%1
		
		BX <- cSplineDes(t,knots)
		
		c(BX %*% coef)
	})
}

betafun_cos <- function(t, param) {
	with(as.list(param), {
		b0 * (1 + b1 * cos(t*2*pi))
	})
}

