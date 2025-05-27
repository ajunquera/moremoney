#setwd("C:/Users/alvar/UAB/OneDrive - Universitat Autònoma de Barcelona/Methodology/Causal inference/Regression discontinuity/Distributional effects/2018 Qu Yoon supplementary/code")

previouswd <- getwd()

setwd("C:/Users/afernan5/Nextcloud/AFJ/AXL/inputcode/2018quyoon_code")

## prepared in March, 2017
## A collection of functions to
## (i) estimate conditional quantile process
## (ii) QTE and its uniform confidence band
## (iii) bandwidth selection
## (iv) Score and Wald tests

source("qte_rdd_funcs.R")

## 1. Functions to estimate the QTE and its uniform band

rdd.rqpro <- function (x,y,tt=c(0.2,0.8),m=9,x.eval,bandw,method=1,kr=3){
	# main functions to estimate conditonal quantile process
	# method=1 implements Procedure 1
	# method=2 implements Procedure 2
	# method=3 estimate conditional quantiles without monotonicity enforcement
	# kr=1 uses Gaussian kernel, kr=2 uses Cauchy kernel
	# kr=3 uses Epanechnikov kernel (default option), kr=4 uses truncated Gaussian kernel
	n  <- length(y)
	taus <- seq(tt[1],tt[2],length.out=m)
	taus <- sort(unique(c(taus,0.5)))
	h.tau <- band.covt(bandw,taus)
	if(method==1){est <- lprq.pro.rdd(x,y,taus,x.eval,h.tau,mon=1,kr)$b0}
	if(method %in% c(2,3)){est <- lprq.pro.rdd(x,y,taus,x.eval,h.tau,mon=0,kr)$b0}
	if(method==2){est <- rearrange2(est,taus,n)}
	return(list(taus = taus, Q0 = est))
}

rdd.qte <- function(x,y,d,x.eval,alpha=0.9,tt,m,bandw,kr=3,bias=0,eql=0){
	# main function to estimate QTE and its confidence bands
	# provide uniform and pointwise bands ('uci' and 'pci' respectively)
	taus <- seq(tt[1],tt[2],length.out=m)
	taus <- sort(unique(c(taus,0.5)))
	taul <- ext.tau(taus)
	ind  <- taul%in%taus
	h.taul<- band.covt(bandw,taul)
	Q.0.p <- lprq.pro.rdd(x[d==1],y[d==1],taul,x.eval,h.tau=h.taul,mon=0,kr)
	Q.0.m <- lprq.pro.rdd(x[d==0],y[d==0],taul,x.eval,h.tau=h.taul,mon=0,kr)
	fqx.p <- condf(Q.0.p$b0,taul,5000,ad=2,kr)
	fqx.m <- condf(Q.0.m$b0,taul,5000,ad=2,kr)
	fx <- marf(x,x.eval)
	ba <- kxb(x,y,x.eval,taus,length(y),d.x=1,h.taul[ind],kr)
	Q.2.p <- local.poly(y[d==1],x[d==1],x.eval,taus,h.taul[ind],order=2)$Q.2
	Q.2.m <- local.poly(y[d==0],x[d==0],x.eval,taus,h.taul[ind],order=2)$Q.2
	if(bias==0){
	out <- rdd.ci(d,alpha,Q.0.p$b0[ind],Q.0.m$b0[ind],fqx.p[ind],fqx.m[ind],fx,Q.2.p,Q.2.m,h.taul[ind],taus,ba,length(y),n.sim=3000,bias=0,eql=0)
	}
	if(bias==1 & eql==0){
	out <- rdd.ci(d,alpha,Q.0.p$b0[ind],Q.0.m$b0[ind],fqx.p[ind],fqx.m[ind],fx,Q.2.p,Q.2.m,h.taul[ind],taus,ba,length(y),n.sim=3000,bias=1,eql=0)
	}
	if(bias==1 & eql==1){
	Q.2.p.eql <- rep(mean(Q.2.p),length(Q.2.p))
	Q.2.m.eql <- rep(mean(Q.2.m),length(Q.2.m))
	out <- rdd.ci(d,alpha,Q.0.p$b0[ind],Q.0.m$b0[ind],fqx.p[ind],fqx.m[ind],fx,Q.2.p.eql,Q.2.m.eql,h.taul[ind],taus,ba,length(y),n.sim=3000,bias=1,eql=1)
	}
	return(list(taus = taus, qte = out$qte, uci = out$uci, pci = out$pci))
}

## 2. Bandwidth selectors

rdd.bandwidth <- function(x,y,d,x.eval,tt,m,method,val=NULL,band=NULL,kr=3,br=0){
	# method=1 produces the leave-one-out CV bandwidth
	# method=2 produces the interior point leave-one-out CV bandwidth
	# method=3 produce the interior point optimal bandwidth
	# method=4 produce the boundary point optimal bandwidth
	# method=5 produce the IK bandwidth
	# val specifies the set of potential values for the CV bandwidths, used in grid search
	# band specifies the bandwidths to estimate the second order derivatives
	# if band is not set by a user, it automatically calculates the bandwidths
	# band[1] is for h^{int}, band[2:3] is for h^{bdy}, and band[4:5] is for h^{ik}
	# kr=1 uses Gaussian kernel, kr=2 uses Cauchy kernel
	# kr=3 uses Epanechnikov kernel (default option), kr=4 uses truncated Gaussian kernel
	# if br=1, interior point optimal bandwidth is robust to the potential break in regression function
	if(is.null(val))
		stop("Cross validation methods require the list of possible bandwidth values.")
	if(max(method %in% 4) & !is.null(band) & length(band)<3)
		stop("the length of band must be 3 or larger.")
	if(max(method %in% 5) & !is.null(band) & length(band)<5)
		stop("the length of band must be 5.")
	if(max(method %in% c(1,3,4,5))){
		h.tau1 <- cv.med.rd(x,y,x.eval,val,xl=0.5,kr)$h.cv
	}
	if(max(method %in% c(1,3,4,5))){
		taus <- seq(tt[1],tt[2],length.out=m)
		taus <- sort(unique(c(taus,0.5)))
		taul <- ext.tau(taus)
		ind  <- taul%in%taus
		where <- taul%in%0.5
		fx    <- marf(x,x.eval)
		h.tau <- band.covt(h.tau1,taul)
	}
	if(max(method %in% 2)){
		h.tau2 <- cv.med.rd.int(x,y,x.eval,val,xl=0.5,kr)$h.cv
	}
	if(max(method %in% 3)){
		Q.0 <- lprq.pro.rdd(x,y,taul,x.eval,h.tau=h.tau,mon=1,kr)
		fqx <- condf(Q.0$b0,taul,5000,ad=2,kr)
		if(max(method %in% 3) & is.null(band)){h2 <- min(abs(max(x)-x.eval),abs(x.eval-min(x)))}
		else{h2 <- band[1]}
		h2  <- min(abs(max(x)-x.eval),abs(x.eval-min(x)))
		if(br==0){Q.2 <- local.poly(y,x,x.eval,taus=0.5,bandw3=h2,order=3)$Q.2}
		if(br==1){Q.2 <- local.poly.2(y,x,x.eval,taus=0.5,bandw3=h2,order=3)$Q.2}
		h.tau3 <- band.opt(length(y),fx,fqx[where],Q.2,kr)$h.opt
	}
	if(max(method %in% c(4,5))){
		Q.0.p <- lprq.pro.rdd(x[d==1],y[d==1],taul,x.eval,h.tau=h.tau,mon=1,kr)
		Q.0.m <- lprq.pro.rdd(x[d==0],y[d==0],taul,x.eval,h.tau=h.tau,mon=1,kr)
		fqx.p <- condf(Q.0.p$b0,taul,5000,ad=2,kr)
		fqx.m <- condf(Q.0.m$b0,taul,5000,ad=2,kr)
	}
	if(max(method %in% 4)){
		if(max(method %in% 4) & is.null(band)){h2.p <- abs(max(x)-x.eval); h2.m <- abs(x.eval-min(x))}
		else{h2.p <- band[2]; h2.m <- band[3]}
		Q.2.p <- local.poly(y[d==1],x[d==1],x.eval,taus=0.5,bandw3=h2.p,order=3)$Q.2
		Q.2.m <- local.poly(y[d==0],x[d==0],x.eval,taus=0.5,bandw3=h2.m,order=3)$Q.2
		h.tau.p <- band.opt.rd(length(y[d==1]),fx,fqx.p[where],Q.2.p)$h.opt
		h.tau.m <- band.opt.rd(length(y[d==0]),fx,fqx.m[where],Q.2.m)$h.opt
		h.tau4  <- min(h.tau.p,h.tau.m)
	}
	if(max(method %in% 5)){
		if(max(method %in% 5) & is.null(band)){h2.p <- 0.5*abs(max(x)-x.eval); h2.m <- 0.5*abs(x.eval-min(x))}
		else{h2.p <- band[4]; h2.m <- band[5]}
		Q.2.p <- local.poly(y[d==1],x[d==1],x.eval,taus=0.5,bandw3=h2.p,order=2)$Q.2
		Q.2.m <- local.poly(y[d==0],x[d==0],x.eval,taus=0.5,bandw3=h2.m,order=2)$Q.2
		h.tau5 <- band.IK.3(x,d,x.eval,length(y),fx,fqx.p[where],fqx.m[where],Q.2.p,Q.2.m,h2.p,h2.m,order=2)$h.ik
	}
	bandwidths <- vector("list",5)
	names(bandwidths) <- c("hcv","hcvi","hint","hbdy","hik")
	loc <- (1:5) %in% method
	for(i in 1:5){
		if(loc[i]==1) {bandwidths[[i]] <- get(paste("h.tau",i,sep=""))}
	}
	return(bandwidths)
}

# 3. Score tests

Score <- function(x,y,x.eval,alpha=c(0.9,0.95),tt,m,bandw,kr=3,left=0){
	# main function for the score test
	# if left=0 (or 1), use observations from the right (or left) side of the cutoff point
	taus <- seq(tt[1],tt[2],length.out=m)
	taus <- sort(unique(c(taus,0.5)))
	h.tau<- band.covt(bandw,taus)
	Q.0  <- lprq.pro.rdd(x,y,taus,x.eval,h.tau,mon=1,kr)
	ba   <- kxb(x,y,x.eval,taus,length(y),d.x=1,h.tau,kr)
	stat <- score.stat(x,y,xt=x.eval,Q.0,h.tau,taus,ba,length(y),left)
	crit <- uci.rd(x,xt=x.eval,alpha=alpha,ba,h.tau,taus,length(y),n.sim=3000,left)
	return(list(test = stat, crit = crit))
}

# 4. Wald tests

Wald <- function(x,y,d,x.eval,alpha=c(0.9,0.95),tt,m,bandw,bandw2,kr,test.type=c(1,2,3),sign.opt=1,eql){
	# main function for Wald tests
	# test.type=1 for significance hypothesis, 2 for homogeneiry, and 3 for unambiguity hypothesis
	# if sign.opt=1, the unambiguity hypothesis tests if effects are uniformly positive
	# if sign.opt=2, the unambiguity hypothesis tests if effects are uniformly negative
	taus  <- seq(tt[1],tt[2],length.out=m)
	taus  <- sort(unique(c(taus,0.5)))
	h.tau <- band.covt(bandw,taus)
	taul  <- ext.tau(taus)
	ind   <- taul%in%taus
	h.taul<- band.covt(bandw,taul)
	Q.0.p <- lprq.pro.rdd(x[d==1],y[d==1],taul,x.eval,h.tau=h.taul,mon=1,kr)
	Q.0.m <- lprq.pro.rdd(x[d==0],y[d==0],taul,x.eval,h.tau=h.taul,mon=1,kr)
	fqx.p <- condf(Q.0.p$b0,taul,5000,ad=2,kr)
	fqx.m <- condf(Q.0.m$b0,taul,5000,ad=2,kr)
	fx  <- marf(x,x.eval)
	ba  <- kxb(x,y,x.eval,taus,length(y),d.x=1,h.tau,kr)
	h2  <- band.covt(bandw2,taus)
	ba2 <- kxb(x,y,x.eval,taus,length(y),d.x=1,h2,kr)

	# bias estimation without equality constraints
	Q.2.p <- local.poly(y[d==1],x[d==1],x.eval,taus,h2,order=2)$Q.2
	Q.2.m <- local.poly(y[d==0],x[d==0],x.eval,taus,h2,order=2)$Q.2

	# bias estimation with equality constraints
	Q.2.e.p <- rep(mean(Q.2.p),length(Q.2.p))
	Q.2.e.m <- rep(mean(Q.2.m),length(Q.2.m))

	# Wald test without bias correction
	wald.s <- wald.sim.15(d,Q.0.p$b0[ind],Q.0.m$b0[ind],fqx.p[ind],fqx.m[ind],fx,Q.2.p,Q.2.m,h.tau,taus,ba,length(y),n.sim=3000)
	test <- wald.test.15(alpha=alpha,wald.s$Q.diff,wald.s$G.diff,wald.s$bias.cor,fqx.p[ind],fqx.m[ind],bis=0,h.tau,taus,length(y),n.sim=3000,test.type,sign.opt)
	stat <- test$test.stat
	crit <- test$cr.value

	# Robust Wald test, bias corrected
	if(eql==0){wald.r <- wald.sim.rbs.15(d,Q.0.p$b0[ind],Q.0.m$b0[ind],fqx.p[ind],fqx.m[ind],fx,Q.2.p,Q.2.m,h.tau,h2,taus,ba,ba2,length(y),n.sim=3000)}
	if(eql==1){wald.r <- wald.sim.rbs.eql.15(d,Q.0.p$b0[ind],Q.0.m$b0[ind],fqx.p[ind],fqx.m[ind],fx,Q.2.e.p,Q.2.e.m,h.tau,h2,taus,ba,ba2,length(y),n.sim=3000)}
	test.r  <- wald.test.15(alpha=alpha,(wald.r$Q.diff-wald.r$bias.cor),wald.r$G.diff,rep(0,length(taus)),fqx.p[ind],fqx.m[ind],bis=0,h.tau,taus,length(y),n.sim=3000,test.type,sign.opt)
	stat.r  <- test.r$test.stat
	crit.r  <- test.r$cr.value

	return(list(wald.test = stat, wald.crit = crit, wald.robust.test = stat.r, wald.robust.crit = crit.r))
}

## 6. Functions used in the empirical application.
## Use these functions when data is big.
## But be aware that these functions can be slower.

# 6.1. Bandwidth estimation

rdd.bandwidth.app <- function(x,y,d,x.eval,tt,m,method,val=NULL,band=NULL,kr,br=0){
	# method=1 produces the leave-one-out CV bandwidth
	# method=2 produces the interior point leave-one-out CV bandwidth
	# method=3 produce the interior point optimal bandwidth
	# method=4 produce the boundary point optimal bandwidth
	# method=5 produce the IK bandwidth
	if(max(method %in% c(1,2)) & is.null(val))
		stop("Cross validation methods require the list of possible bandwidth values.")

	if(max(method %in% 4) & !is.null(band) & length(band)<3)
		stop("the length of band must be 3 or larger.")

	if(max(method %in% 5) & !is.null(band) & length(band)<5)
		stop("the length of band must be 5.")

	if(max(method %in% c(1,3,4,5))){
		h.tau1 <- cv.rd.emp(x,y,x.eval,val,xl=0.5,kr)$h.cv
	}
	if(max(method %in% c(1,3,4,5))){
		taus <- seq(tt[1],tt[2],length.out=m)
		taus <- sort(unique(c(taus,0.5)))
		taul <- ext.tau(taus)
		ind  <- taul%in%taus
		where <- taul%in%0.5
		fx    <- marf(x,x.eval)
		h.tau <- band.covt(h.tau1,taul)
	}
	if(max(method %in% 2)){
		h.tau2 <- cv.rd.emp.int(x,y,x.eval,val,xl=0.5,kr)$h.cv
	}
	if(max(method %in% 3)){
		Q.0 <- lprq.pro.rdd(x,y,taul,x.eval,h.tau=h.tau,mon=0,kr)
		fqx <- condf(sort(Q.0$b0),taul,5000,ad=2,kr)
		if(max(method %in% 3) & is.null(band)){h2 <- min(abs(max(x)-x.eval),abs(x.eval-min(x)))}
		else{h2 <- band[1]}
		h2  <- min(abs(max(x)-x.eval),abs(x.eval-min(x)))
		if(br==0){Q.2 <- local.poly(y,x,x.eval,taus=0.5,bandw3=h2,order=3)$Q.2}
		if(br==1){Q.2 <- local.poly.2(y,x,x.eval,taus=0.5,bandw3=h2,order=3)$Q.2}
		h.tau3 <- band.opt(length(y),fx,fqx[where],Q.2,kr)$h.opt
	}
	if(max(method %in% c(4,5))){
		Q.0.p <- lprq.pro.rdd(x[d==1],y[d==1],taul,x.eval,h.tau=h.tau,mon=0,kr)
		Q.0.m <- lprq.pro.rdd(x[d==0],y[d==0],taul,x.eval,h.tau=h.tau,mon=0,kr)
		fqx.p <- condf(sort(Q.0.p$b0),taul,5000,ad=2,kr)
		fqx.m <- condf(sort(Q.0.m$b0),taul,5000,ad=2,kr)
	}
	if(max(method %in% 4)){
		if(max(method %in% 4) & is.null(band)){h2.p <- abs(max(x)-x.eval); h2.m <- abs(x.eval-min(x))}
		else{h2.p <- band[2]; h2.m <- band[3]}
		Q.2.p <- local.poly(y[d==1],x[d==1],x.eval,taus=0.5,bandw3=h2.p,order=3)$Q.2
		Q.2.m <- local.poly(y[d==0],x[d==0],x.eval,taus=0.5,bandw3=h2.m,order=3)$Q.2
		h.tau.p <- band.opt.rd(length(y[d==1]),fx,fqx.p[where],Q.2.p)$h.opt
		h.tau.m <- band.opt.rd(length(y[d==0]),fx,fqx.m[where],Q.2.m)$h.opt
		h.tau4  <- min(h.tau.p,h.tau.m)
	}
	if(max(method %in% 5)){
		if(max(method %in% 5) & is.null(band)){h2.p <- 0.5*abs(max(x)-x.eval); h2.m <- 0.5*abs(x.eval-min(x))}
		else{h2.p <- band[4]; h2.m <- band[5]}
		Q.2.p <- local.poly(y[d==1],x[d==1],x.eval,taus=0.5,bandw3=h2.p,order=2)$Q.2
		Q.2.m <- local.poly(y[d==0],x[d==0],x.eval,taus=0.5,bandw3=h2.m,order=2)$Q.2
		h.tau5 <- band.IK.3.big(x,d,x.eval,length(y),fx,fqx.p[where],fqx.m[where],Q.2.p,Q.2.m,h2.p,h2.m)$h.ik
	}
	bandwidths <- vector("list",5)
	names(bandwidths) <- c("hcv","hcvi","hint","hbdy","hik")
	loc <- (1:5) %in% method
	for(i in 1:5){
		if(loc[i]==1) {bandwidths[[i]] <- get(paste("h.tau",i,sep=""))}
	}
	return(bandwidths)
}

# 6.2. Uniform (and pointwise) confidence intervals (used in application)

rdd.qte.app <- function(x,y,d,x.eval,alpha=0.9,tt,m,bandw,kr,bias=0,eql=0){
	# main function for the confidence bands
	# provide QTE and its uniform and pointwise bands
	taus <- seq(tt[1],tt[2],length.out=m)
	taus <- sort(unique(c(taus,0.5)))
	taul <- ext.tau(taus)
	ind  <- taul%in%taus
	h.taul<- band.covt(bandw,taul)
	Q.0.p <- lprq.pro.rdd(x[d==1],y[d==1],taul,x.eval,h.tau=h.taul,mon=0,kr)
	Q.0.m <- lprq.pro.rdd(x[d==0],y[d==0],taul,x.eval,h.tau=h.taul,mon=0,kr)
	fqx.p <- condf(Q.0.p$b0,taul,5000,ad=2,kr)
	fqx.m <- condf(Q.0.m$b0,taul,5000,ad=2,kr)
	fx <- marf(x,x.eval)
	ba <- kxb(x,y,x.eval,taus,length(y),d.x=1,h.taul[ind],kr)
	Q.2.p <- local.poly(y[d==1],x[d==1],x.eval,taus,h.taul[ind],order=2)$Q.2
	Q.2.m <- local.poly(y[d==0],x[d==0],x.eval,taus,h.taul[ind],order=2)$Q.2
	if(bias==0){
	out <- rdd.ci.app(d,alpha,Q.0.p$b0[ind],Q.0.m$b0[ind],fqx.p[ind],fqx.m[ind],fx,Q.2.p,Q.2.m,h.taul[ind],taus,ba,length(y),n.sim=3000,bias=0,eql=0)
	}
	if(bias==1 & eql==0){
	out <- rdd.ci.app(d,alpha,Q.0.p$b0[ind],Q.0.m$b0[ind],fqx.p[ind],fqx.m[ind],fx,Q.2.p,Q.2.m,h.taul[ind],taus,ba,length(y),n.sim=3000,bias=1,eql=0)
	}
	if(bias==1 & eql==1){
	Q.2.p.eql <- rep(mean(Q.2.p),length(Q.2.p))
	Q.2.m.eql <- rep(mean(Q.2.m),length(Q.2.m))
	out <- rdd.ci.app(d,alpha,Q.0.p$b0[ind],Q.0.m$b0[ind],fqx.p[ind],fqx.m[ind],fx,Q.2.p.eql,Q.2.m.eql,h.taul[ind],taus,ba,length(y),n.sim=3000,bias=1,eql=1)
	}
	return(list(taus = taus, qte = out$qte, uci = out$uci, pci = out$pci))
}

# 6.3. Score and Wald tests (used in application with a large size data)

Score.app <- function(x,y,x.eval,alpha=c(0.9,0.95),tt,m,bandw,kr,left=0){
	# main function for the score test
	taus <- seq(tt[1],tt[2],length.out=m)
	taus <- sort(unique(c(taus,0.5)))
	h.tau<- band.covt(bandw,taus)
	Q.0  <- lprq.pro.rdd(x,y,taus,x.eval,h.tau,mon=1,kr)
	ba   <- kxb(x,y,x.eval,taus,length(y),d.x=1,h.tau,kr)
	stat <- score.stat(x,y,xt=x.eval,Q.0,h.tau,taus,ba,length(y),left)
	crit <- uci.rd.big(x,xt=x.eval,alpha,ba,h.tau,taus,length(y),n.sim=5000,left)
	return(list(test = stat, crit = crit))
}

Wald.app <- function(x,y,d,x.eval,alpha=c(0.9,0.95),tt,m,bandw,bandw2,kr,test.type=c(1,2,3),sign.opt=1,eql){
	# main function for Wald tests
	taus  <- seq(tt[1],tt[2],length.out=m)
	taus  <- sort(unique(c(taus,0.5)))
	h.tau <- band.covt(bandw,taus)
	taul  <- ext.tau(taus)
	ind   <- taul%in%taus
	h.taul<- band.covt(bandw,taul)
	Q.0.p <- lprq.pro.rdd(x[d==1],y[d==1],taul,x.eval,h.tau=h.taul,mon=1,kr)
	Q.0.m <- lprq.pro.rdd(x[d==0],y[d==0],taul,x.eval,h.tau=h.taul,mon=1,kr)
	fqx.p <- condf(Q.0.p$b0,taul,5000,ad=2,kr)
	fqx.m <- condf(Q.0.m$b0,taul,5000,ad=2,kr)
	fx  <- marf(x,x.eval)
	ba  <- kxb(x,y,x.eval,taus,length(y),d.x=1,h.tau,kr)
	h2  <- band.covt(bandw2,taus)
	ba2 <- kxb(x,y,x.eval,taus,length(y),d.x=1,h2,kr)

	# bias estimation without equality constraints
	Q.2.p <- local.poly(y[d==1],x[d==1],x.eval,taus,h2,order=2)$Q.2
	Q.2.m <- local.poly(y[d==0],x[d==0],x.eval,taus,h2,order=2)$Q.2

	# bias estimation with equality constraints
	Q.2.e.p <- rep(mean(Q.2.p),length(Q.2.p))
	Q.2.e.m <- rep(mean(Q.2.m),length(Q.2.m))

	# Wald test without bias correction
	wald.s <- wald.sim.big(d,Q.0.p$b0[ind],Q.0.m$b0[ind],fqx.p[ind],fqx.m[ind],fx,Q.2.p,Q.2.m,h.tau,taus,ba,length(y),n.sim=5000,testing=1)
	test  <- wald.test.15(alpha=alpha,wald.s$Q.diff,wald.s$G.diff,wald.s$bias.cor,fqx.p[ind],fqx.m[ind],bis=0,h.tau,taus,length(y),n.sim=5000,test.type,sign.opt)
	stat <- test$test.stat
	crit <- test$cr.value

	# Robust Wald test, bias corrected
	if(eql==0){
		wald.r <- wald.sim.big.rbs(d,Q.0.p$b0[ind],Q.0.m$b0[ind],fqx.p[ind],fqx.m[ind],fx,Q.2.p,Q.2.m,h.tau,taus,ba,length(y),n.sim=5000,testing=1)
		}
	if(eql==1){
		wald.r <- wald.sim.big.rbs.eql(d,Q.0.p$b0[ind],Q.0.m$b0[ind],fqx.p[ind],fqx.m[ind],fx,Q.2.e.p,Q.2.e.m,h.tau,h2,taus,ba,ba2,length(y),n.sim=5000,testing=1)
		}
	test.r  <- wald.test.15(alpha=alpha,(wald.r$Q.diff-wald.r$bias.cor),wald.r$G.diff,rep(0,length(taus)),fqx.p[ind],fqx.m[ind],bis=0,h.tau,taus,length(y),n.sim=5000,test.type,sign.opt)
	stat.r  <- test.r$test.stat
	crit.r  <- test.r$cr.value

	return(list(wald.test = stat, wald.crit = crit, wald.robust.test = stat.r, wald.robust.crit = crit.r))
}

setwd(previouswd)
