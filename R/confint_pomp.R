confint_pomp <- function(pomp_object,
						 prob=0.95,
						 par=c(1,2),
						 delta=0.01,
						 rwsd_arg,
						 trace=FALSE,
						 seed=101) {
	set.seed(seed)
	
	cc <- coef(pomp_object)
	
	ll_max <- logmeanexp(replicate(10,logLik(pfilter(pomp_object,Np=1000))),se=TRUE)
	
	maxdata <- data.frame(as.list(cc), ll=ll_max[1], ll.se=ll_max[2], ldiff=0)
	
	maxdiff <- qchisq(prob, 1)/2
	cutoff <- c(-maxdiff, maxdiff)
	
	reslist <- vector('list', length(par))
	
	for (p in par) {
		ldiff <- 0
		
		uplist <- list()
		i <- 1
		
		while (abs(ldiff) < maxdiff) {
			## up
			cc2 <- cc
			cc2[p] <- cc[p] + delta * i
			
			mprof <- mif2(
				pomp_object,
				Nmif=50,
				start=cc2,
				Np=1000,
				cooling.fraction.50=0.95,
				rw.sd=do.call(rw.sd, rwsd_arg[-p]),
				transform=TRUE) %>%
				continue(Nmif=50, cooling.fraction=0.8) %>%
				continue(Nmif=50, cooling.fraction=0.6) %>%
				continue(Nmif=50, cooling.fraction=0.2) %>%
				continue(Nmif=50, cooling.fraction=0.1)
			
			ll_prof <- logmeanexp(replicate(10,logLik(pfilter(mprof,Np=1000))),se=TRUE)
			
			ldiff <- ll_max[1] - ll_prof[1]
			
			uplist[[i]] <- data.frame(
				as.list(coef(mprof)),
				ll=ll_prof[1],
				ll.se=ll_prof[2],
				ldiff=ldiff
			)
			
			if(trace) print(uplist[[i]])
			
			i <- i + 1
		}
		
		ldiff <- 0
		
		downlist <- list()
		i <- 1
		
		while (abs(ldiff) < maxdiff) {
			## down
			cc2 <- cc
			cc2[p] <- cc[p] - delta * i
			
			mprof <- mif2(
				pomp_object,
				Nmif=50,
				start=cc2,
				Np=1000,
				cooling.fraction.50=0.95,
				rw.sd=do.call(rw.sd, rwsd_arg[-p]),
				transform=TRUE) %>%
				continue(Nmif=50, cooling.fraction=0.8) %>%
				continue(Nmif=50, cooling.fraction=0.6) %>%
				continue(Nmif=50, cooling.fraction=0.2) %>%
				continue(Nmif=50, cooling.fraction=0.1)
			
			ll_prof <- logmeanexp(replicate(10,logLik(pfilter(mprof,Np=1000))),se=TRUE)
			
			ldiff <- ll_prof[1] - ll_max[1]
			
			downlist[[i]] <- data.frame(
				as.list(coef(mprof)),
				ll=ll_prof[1],
				ll.se=ll_prof[2],
				ldiff=ldiff
			)
			
			if(trace) print(downlist[[i]])
			
			i <- i + 1
		}
		
		pdata <- c(uplist, downlist, list(maxdata)) %>% 
			bind_rows
		
		pdata <- pdata[order(pdata[,p]),]
		
		sfit <- spline(x=pdata[,p], y=pdata$ldiff)
		
		tt <- approx(sfit$y, sfit$x, xout=cutoff)$y
		
		reslist[[p]] <- data.frame(
			par=names(cc)[p],
			lwr=sort(tt)[1],
			upr=sort(tt)[2]
		)
	}
	
	resdata <- reslist %>% bind_rows
	
	resdata
}

## higher cooling fraction
confint_pomp2 <- function(pomp_object,
						 prob=0.95,
						 par=c(1,2),
						 delta=0.01,
						 rwsd_arg,
						 trace=FALSE,
						 seed=101) {
	set.seed(seed)
	
	cc <- coef(pomp_object)
	
	ll_max <- logmeanexp(replicate(10,logLik(pfilter(pomp_object,Np=1000))),se=TRUE)
	
	maxdata <- data.frame(as.list(cc), ll=ll_max[1], ll.se=ll_max[2], ldiff=0)
	
	maxdiff <- qchisq(prob, 1)/2
	cutoff <- c(-maxdiff, maxdiff)
	
	reslist <- vector('list', length(par))
	
	for (p in par) {
		ldiff <- 0
		
		uplist <- list()
		i <- 1
		
		while (abs(ldiff) < maxdiff) {
			## up
			cc2 <- cc
			cc2[p] <- cc[p] + delta * i
			
			mprof <- mif2(
				pomp_object,
				Nmif=100,
				start=cc2,
				Np=1000,
				cooling.fraction.50=0.95,
				rw.sd=do.call(rw.sd, rwsd_arg[-p]),
				transform=TRUE) %>%
				continue(Nmif=100, cooling.fraction=0.8) %>%
				continue(Nmif=100, cooling.fraction=0.6)
			
			ll_prof <- logmeanexp(replicate(10,logLik(pfilter(mprof,Np=1000))),se=TRUE)
			
			ldiff <- ll_max[1] - ll_prof[1]
			
			uplist[[i]] <- data.frame(
				as.list(coef(mprof)),
				ll=ll_prof[1],
				ll.se=ll_prof[2],
				ldiff=ldiff
			)
			
			if(trace) print(uplist[[i]])
			
			i <- i + 1
		}
		
		ldiff <- 0
		
		downlist <- list()
		i <- 1
		
		while (abs(ldiff) < maxdiff) {
			## down
			cc2 <- cc
			cc2[p] <- cc[p] - delta * i
			
			mprof <- mif2(
				pomp_object,
				Nmif=100,
				start=cc2,
				Np=1000,
				cooling.fraction.50=0.95,
				rw.sd=do.call(rw.sd, rwsd_arg[-p]),
				transform=TRUE) %>%
				continue(Nmif=100, cooling.fraction=0.8) %>%
				continue(Nmif=100, cooling.fraction=0.6)
			
			ll_prof <- logmeanexp(replicate(10,logLik(pfilter(mprof,Np=1000))),se=TRUE)
			
			ldiff <- ll_prof[1] - ll_max[1]
			
			downlist[[i]] <- data.frame(
				as.list(coef(mprof)),
				ll=ll_prof[1],
				ll.se=ll_prof[2],
				ldiff=ldiff
			)
			
			if(trace) print(downlist[[i]])
			
			i <- i + 1
		}
		
		pdata <- c(uplist, downlist, list(maxdata)) %>% 
			bind_rows
		
		pdata <- pdata[order(pdata[,p]),]
		
		sfit <- spline(x=pdata[,p], y=pdata$ldiff)
		
		tt <- approx(sfit$y, sfit$x, xout=cutoff)$y
		
		reslist[[p]] <- data.frame(
			par=names(cc)[p],
			lwr=sort(tt)[1],
			upr=sort(tt)[2]
		)
	}
	
	resdata <- reslist %>% bind_rows
	
	resdata
}

