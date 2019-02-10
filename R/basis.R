estimate_basis <- function(x,
						   k=seq(0, 26, by=2)) {
	BX <- mgcv::cSplineDes(0:26,k)
	basis_covar <- as.data.frame(BX)
	colnames(basis_covar) <- paste0("B", 1:ncol(basis_covar))
	
	basis_covar$y <- c(tail(x, 1), x)
	
	cc <- coef(lm(y~-1 + ., data=basis_covar))
	names(cc) <- paste0("b", 1:length(cc))
	
	cc
}
