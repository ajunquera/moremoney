
## Prepared in March, 2017
## A collection of supporting functions
## They will be called to implement the main functions in 'qte_rdd.R'

require(quantreg)

## 1. Functions to estimate the QTE and its uniform band

# 1.1. Function to Estimate conditional quantiles

lprq.pro.rdd <- function (x,y,tau,x.eval,h.tau,mon=TRUE,kr=3) {
	# rdd.rqpro() calls this function
	n  <- length(y)
	fv <- NULL
	z  <- x - x.eval
	if(kr==1){wx <- dnorm((1/h.tau)%x%z)}
	if(kr==2){wx <- dcauchy((1/h.tau)%x%z)}
	if(kr==3){wx <- depa((1/h.tau)%x%z)}
	if(kr==4){wx <- dtnorm((1/h.tau)%x%z)}
	if(kr==3){
		hk <- which.max(h.tau)		
		ind.0 <- (wx[(1+(hk-1)*n):(hk*n)]==0)
		y  <- y[ind.0==0]
		z  <- z[(ind.0==0)]
		wx <- wx[rep(ind.0,length(h.tau))==0]
		}
	if(mon==1){r <- rq.pro(y, z, taus = tau, w = wx)
			fv <- r$coef}
	if(mon==0){r <- NULL
		for(k in 1:length(tau)){
			rt <- rq(y ~ z, tau = tau[k], w = as.vector(wx[(1+(k-1)*length(y)):(k*length(y))]))$coef
			r  <- cbind(r,as.numeric(rt))
		}
		fv <- cbind(fv,r)}
		b0 <- fv[1,]
		b1 <- fv[2,]
	return(list(b0 = b0, b1 = b1))
}

rq.pro <- function(y, x, taus, w=NULL){
	# lprq.pro.rdd() calls this function	
    if(is.matrix(x)==1){
    n = nrow(x)
    p = ncol(x)}
    if(is.matrix(x)==0){
    n = length(x)
    p = 1}
    m = length(taus)
    pp = p+1
    
    if(length(w)==0){
    Y = rep(y, m)
    x = cbind(rep(1,n),x)
    D = diag(m)
    D[lower.tri(D)] = 1 
    X = D %x% x
    }
    
    if(length(w)>0){
    Y = as.vector(y*w)
    x = cbind(rep(1,n),x)
    D = diag(m)
    D[lower.tri(D)] = 1 
    X = (D %x% x)*as.vector(w)
    }
    
    sX = as.matrix.csr(X)
    K = m*pp
   	R <- array(0,c((m-1),(pp*m)))
	for(i in 1:(m-1)){
	  for(j in 1:(pp*m)){
		if(j==(i*pp+1)) R[i,j] <- 1
	  }
	}

    sR = as.matrix.csr(R)
    r2 = rep(0, nrow(R))
       
    tau.v = rep(taus, each=n)  
    rhs = t(X)%*%(1-tau.v)     
    coeff1 =  myrq.fit.sfnc(sX, Y, sR, r2, tau=tau.v, rhs=rhs, nsubmax = 1e+6, tmpmax=1e+6)$coef
    #coeff <- matrix(coeff1,ncol=pp)
    coeff = coeff1[1:K]
    gamma.m = matrix(coeff, ncol=m)
    D = diag(m)
    D[upper.tri(D)]=1
    bhat = gamma.m %*% D
  return(list(coef = bhat))   # col for different tau's, and row for x
}

# 1.2. Function to estimate a uniform band

rdd.ci <- function(d,alpha,Qp,Qm,fxp,fxm,f0,Q2.p,Q2.m,bandw,taus,ba,n,n.sim=3000,bias=0,eql=0){
	# uniform and pointwise confidence bands
	# rdd.qte() calls this function
	fxa  <- 0.5*(fxp+fxm)
	if(bias==0){ob <- wald.sim.15(d,Qp,Qm,fxp,fxm,f0,Q2.p,Q2.m,bandw,taus,ba,n,n.sim,testing=0)}
	if(bias==1 & eql==0){ob <- wald.sim.rbs.15(d,Qp,Qm,fxp,fxm,f0,Q2.p,Q2.m,bandw,bandw,taus,ba,ba,n,n.sim,testing=0)}
	if(bias==1 & eql==1){ob <- wald.sim.rbs.eql.15(d,Qp,Qm,fxp,fxm,f0,Q2.p,Q2.m,bandw,bandw,taus,ba,ba,n,n.sim,testing=0)}
	
	Q.d <- ob$Q.diff
	G.d <- ob$G.diff
	if(bias==0) {bias.adj <- rep(0,length(taus))}
	if(bias==1) {bias.adj <- ob$bias.cor}

	zz <- apply(abs(G.d),1,max)
	za <- quantile(zz,prob=alpha)
	sig <- za/(sqrt(n*bandw)*fxa)
	uci <- cbind((Q.d-bias.adj) - sig, (Q.d-bias.adj) + sig)  # uniform bands
	
	zss <- apply(G.d,2,quantile,c((1-alpha)/2,1-(1-alpha)/2))
	cc  <- t(zss)/(sqrt(n*bandw)*fxa)	
	pci <- cbind((Q.d-bias.adj) + cc[,1], (Q.d-bias.adj) + cc[,2])	# pointwise bounds

	return(list(qte = (Q.d-bias.adj), uci = uci, pci = pci))
}

rdd.ci.app <- function(d,alpha,Qp,Qm,fxp,fxm,f0,Q2.p,Q2.m,bandw,taus,ba,n,n.sim=3000,bias=0,eql=0){
	# uniform and pointwise confidence bands
	# rdd.qte.app() calls this function
	fxa  <- 0.5*(fxp+fxm)
	if(bias==0){ob <- wald.sim.big(d,Qp,Qm,fxp,fxm,f0,Q2.p,Q2.m,bandw,taus,ba,n,n.sim,testing=0)}
	if(bias==1 & eql==0){ob <- wald.sim.big.rbs(d,Qp,Qm,fxp,fxm,f0,Q2.p,Q2.m,bandw,taus,ba,n,n.sim,testing=0)}
	if(bias==1 & eql==1){ob <- wald.sim.big.rbs.eql(d,Qp,Qm,fxp,fxm,f0,Q2.p,Q2.m,bandw,bandw,taus,ba,ba,n,n.sim,testing=0)}
	
	Q.d <- ob$Q.diff
	G.d <- ob$G.diff
	if(bias==0) {bias.adj <- rep(0,length(taus))}
	if(bias==1) {bias.adj <- ob$bias.cor}

	zz <- apply(abs(G.d),1,max)
	za <- quantile(zz,prob=alpha)
	sig <- za/(sqrt(n*bandw)*fxa)
	uci <- cbind((Q.d-bias.adj) - sig, (Q.d-bias.adj) + sig)  # uniform bands
	
	zss <- apply(G.d,2,quantile,c((1-alpha)/2,1-(1-alpha)/2))
	cc  <- t(zss)/(sqrt(n*bandw)*fxa)	
	pci <- cbind((Q.d-bias.adj) + cc[,1], (Q.d-bias.adj) + cc[,2])	# pointwise bounds

	return(list(qte = (Q.d-bias.adj), uci = uci, pci = pci))
}

## 2. Bandwidth selectors

# 2.1. CV bandwidths - h^{cv} and h^{cvi}

cv.med.rd <- function(x,y,x.eval,val,xl,kr){
	# cv bandwidth treating the cutoof point as a boundary point
	dd  <- (x>x0)
	cut <- c(quantile((x[dd==1]-x.eval),prob=xl),quantile((x[dd==0]-x.eval),prob=(1-xl)))
	index <- which((x-x.eval <= cut[1]) & (x-x.eval >= cut[2]))
	cri <- array(0,c(length(val),1))
	for (i in 1:length(val)){
		h.tau <- val[i]
		cri.h <- NULL
		for (j in index){
			xx=x[-j]
			yy=y[-j]
			if(x[j]>=x.eval){sgn <- ((xx >= x.eval) & (xx>x[j]))}
			if(x[j]< x.eval){sgn <- ((xx <  x.eval) & (xx<x[j]))}
			Qy <- lprq.pro.rdd(xx[sgn],yy[sgn],tau=0.5,x[j],h.tau,mon=0,kr)$b0
			cri.h <- c(cri.h,abs(y[j]-Qy))
        	}
		cri[i] <- mean(cri.h)
	}
	return(list(h.cv = val[which.min(cri)], cand = cbind(val,cri)))
}

cv.med.rd.int <- function(x,y,x.eval,val,xl,kr){
	# cv bandwidth treating the cutoof point as an interior point
	cut   <- quantile(abs(x-x.eval),probs=xl)
	index <- which(abs(x-x.eval) <= cut)
	cri   <- array(0,c(length(val),1))
	for (i in 1:length(val)){
		h.tau <- val[i]
		cri.h <- NULL
		for (j in index){
			xx=x[-j]
			yy=y[-j]
			Qy  <- lprq.pro.rdd(xx,yy,tau=0.5,x[j],h.tau,mon=0,kr)$b0
			cri.h <- c(cri.h,abs(y[j]-Qy))
        }
		cri[i] <- mean(cri.h)
	}
	return(list(h.cv = val[which.min(cri)], cand = cbind(val,cri)))
}

# 2.2. Functions for the optimal bandwidths

band.opt.rd <- function(n,fx,fqx,Q.2){
	# boundary point bandwidth for RDD
	h.med <- (((0.5*0.5*4.497)/(fx*(fqx^2)*0.0134*(Q.2^2)))^{1/5})*n^{-1/5}		
	return(list(h.opt = h.med))	
}

band.opt <- function(n,fx,fqx,trq,kr=3){
# interior point optimal bandwidth 
	if(kr==1){h.med <- (((1*0.25*(0.282^{1}))/(1*fx*(fqx^2)*(trq^2)))^{1/(4+1)})*(n^{-1/(4+1)})}
	if(kr==3){h.med <- (((1*0.25*(0.6^{1}))/((0.2^2)*fx*(fqx^2)*(trq^2)))^{1/(4+1)})*(n^{-1/(4+1)})}
	if(kr==4){h.med <- (((1*0.25*(0.282095^{1}))/((0.9999813^2)*fx*(fqx^2)*(trq^2)))^{1/(4+1)})*(n^{-1/(4+1)})}
	return(list(h.opt = as.numeric(h.med)))	
}

# 2.3. Function for IK bandwidth

band.IK.3 <- function(x,d,x0,n,fx,fqx.p,fqx.m,Q.2.p,Q.2.m,hp,hm,order=2){
	# IK bandwidth selection
	z <- x - x0
	r <- array(0,c(2,1))
	for(k in 1:2){
		if(k==1){hh = hp; fqx = fqx.p; dd = d}
		if(k==2){hh = hm; fqx = fqx.m; dd = 1-d}
		wx <- depa(z/hh)
			if(order==2) {zt <- cbind(1,(z/hh),((z/hh)^2))}
			if(order==3) {zt <- cbind(1,(z/hh),((z/hh)^2),((z/hh)^3))}
			Nf <- (t(zt) %*% diag(dd*wx) %*% zt)/(n*hh)
			Mf <- (t(zt) %*% diag(dd*wx*wx) %*% zt)/(n*hh)
			r[k]  <- 3*((n*(hh^{5}))^{-1})*((solve(Nf) %*% Mf %*% solve(Nf))[3,3])/(fqx^2)
		}
	rp <- r[1]
	rm <- r[2]
	con1 <- 0.5*0.5*4.497*((fqx.p^{-2})+(fqx.m^{-2}))
	con2 <- 0.0134*fx*(((Q.2.p-Q.2.m)^{2})+rp+rm)
	h.med <- (con1/con2)^{1/5}*(n^{-1/5})
	return(list(h.ik = as.numeric(h.med)))	
}

## 3. Function to estimate the bias term

local.poly <- function(y,x,x.eval,taus,bandw3,order=3){
	# estimate local polynomial quantile regression
	# order=2 for quadratic, 3 for cubic, 4 for quartic models
	z <- x - x.eval
	n <- length(z)	
	wx<- depa((1/bandw3)%x%z)		
	if(order==2){lc.model <- y~z+I(z^2)}
	if(order==3){lc.model <- y~z+I(z^2)+I(z^3)}
	if(order==4){lc.model <- y~z+I(z^2)+I(z^3)+I(z^4)}
	coe <- NULL
	for(qq in 1:length(taus)){
		coe   <- cbind(coe,rq(lc.model, tau = taus[qq], w = as.vector(wx[(1+(qq-1)*n):(qq*n)]))$coef)
	}
	b.0   <- coe[1,]
	b.1st <- coe[2,]
	b.2nd <- coe[3,]
	if(order==3) {b.3rd <- coe[4,]}
	if(order==4) {b.3rd <- coe[4,]; b.4th <- coe[5,]}
	if(order==2) {return(list(Q = b.0, Q.1 = b.1st, Q.2 = (2*b.2nd)))}
	if(order==3) {return(list(Q = b.0, Q.1 = b.1st, Q.2 = (2*b.2nd), Q.3 = (6*b.3rd)))}
	if(order==4) {return(list(Q = b.0, Q.1 = b.1st, Q.2 = (2*b.2nd), Q.3 = (6*b.3rd), Q.4 = (24*b.4th)))}
}

local.poly.2 <- function(y,x,x.eval,taus,bandw3,order=3){
	# a variation of local.poly()
	# use observations from both sides of the cutoff but allow a break in level
	n  <- length(y)
	zs <- x - x.eval
	wx <- depa((1/bandw3)%x%zs)
	ds <- (zs >= 0)
	if(order==1){lc.model <- y~ds+zs}
	if(order==2){lc.model <- y~ds+zs+I(zs^2)}
	if(order==3){lc.model <- y~ds+zs+I(zs^2)+I(zs^3)}
	if(order==4){lc.model <- y~ds+zs+I(zs^2)+I(zs^3)+I(zs^4)}
	coe <- NULL
	for(qq in 1:length(taus)){
		coe   <- cbind(coe,rq(lc.model, tau = taus[qq], w = as.vector(wx[(1+(qq-1)*n):(qq*n)]))$coef)
	}
	b.0.m <- coe[1,]
	b.0.p <- coe[1,]+coe[2,]
	b.1st <- coe[3,]
	if(order%in%c(2,3,4)) {b.2nd <- coe[4,]}
	if(order%in%c(3,4)) {b.3rd <- coe[5,]}
	if(order%in%c(4)) {b.4th <- coe[6,]}
	if(order==1){return(list(Qm = b.0.m, Qp = b.0.p, Q.1 = b.1st))}
	if(order==2){return(list(Qm = b.0.m, Qp = b.0.p, Q.1 = b.1st, Q.2 = (2*b.2nd)))}
	if(order==3){return(list(Qm = b.0.m, Qp = b.0.p, Q.1 = b.1st, Q.2 = (2*b.2nd), Q.3 = (6*b.3rd)))}
	if(order==4){return(list(Qm = b.0.m, Qp = b.0.p, Q.1 = b.1st, Q.2 = (2*b.2nd), Q.3 = (6*b.3rd), Q.4 = (24*b.4th)))}
}

## 4. Furcntions for the Score test

score.stat <- function(x,y,xt,Q.est,bandw,taus,ba,n,left=0){
	# test statistic for the significance hypothesis
	# use Epanechnikov kernel
	# if left=0 (or 1), use observations from the right (or left) side of the cutoff point
	ss <- array(0,c(length(taus),1))
	for(i in 1:length(taus)){
		t <- taus[i]
		uhat <- y - Q.est$b0[i] - (ba$zi)*(Q.est$b1[i])
		if(left==1) {pa2 <- sum((t-(uhat<0))*ba$K.x[,i]*(x < xt))/sqrt(n*bandw[i])}
		if(left==0) {pa2 <- sum((t-(uhat<0))*ba$K.x[,i]*(x >= xt))/sqrt(n*bandw[i])}
		ss[i] <- abs(pa2)
	}
	return(max(ss))
}

uci.rd <- function(x,xt,alpha,ba,bandw,taus,n,n.sim,left=0){
	# critical value for the significance hypothesis
	# use Epanechnikov kernel
	# if left=0 (or 1), use observations from the right (or left) side of the cutoff point
	zs  <- array(0,c(n.sim,length(taus)))
	umat <- matrix(runif(n*n.sim),ncol=n.sim)
	for(i in 1:length(taus)){
		t  <- taus[i]
		for(j in 1:n.sim){
			uu <- umat[,j]
			if(left==1) {pa2 <- sum((t-(uu-t<0))*ba$K.x[,i]*((x < xt)-0.5+(15/16)*(ba$zi/bandw[i])))/sqrt(n*bandw[i])}
			if(left==0) {pa2 <- sum((t-(uu-t<0))*ba$K.x[,i]*((x >=xt)-0.5-(15/16)*(ba$zi/bandw[i])))/sqrt(n*bandw[i])}			
			zs[j,i] <- abs(pa2)
		}
	}
	zz <- apply(zs,1,max)
	z.a <- quantile(zz,probs=alpha)
	return(z.a)
}

## 5. Functions for Wald tests

wald.test.15 <- function(alpha,Q.d,g.in,bias.adj,fxp,fxm,bis,bandw,taus,n,n.sim,test.type,sign.opt){
	# performs Wald tests
	# if sign.opt=1, the unambiguity hypothesis tests if effects are uniformly positive
	# if sign.opt=2, the unambiguity hypothesis tests if effects are uniformly negative
	fxa <- 0.5*(fxp+fxm)
	if(bis==0) {bias.adj <- rep(0,length(taus))}
	tsall <- NULL
	zaall <- NULL
	Q.in <- sqrt(n*bandw)*fxa*(Q.d - bias.adj)
	if(1%in%test.type){  # test for the Treatment Significance
		ts <- max(abs(Q.in))  		# test statistics
		zz <- apply(abs(g.in),1,max)
		za <- quantile(zz,probs=alpha)		# critical values
			tsall <- c(tsall,ts)
			zaall <- cbind(zaall,za)
		}
	if(2%in%test.type){  # test for the Treatment Homegeneity
		m2 <- sqrt(n*bandw)*fxa
		ts <- max(abs(Q.in - (m2/mean(m2))*rep(mean(Q.in),length(taus))))	# test statistic
		m2.mat <- matrix(rep((m2/mean(m2)),n.sim),byrow=TRUE,nrow=n.sim)
		zmean <- apply((g.in),1,mean)
		zz <- apply(abs(g.in - m2.mat*matrix(rep(zmean,length(taus)),nrow=n.sim)),1,max)
		za <- quantile(zz,probs=alpha)				# critical values
			tsall <- c(tsall,ts)
			zaall <- cbind(zaall,za)
		}
	if(3%in%test.type){  # test for the Treatment Unambiguity
		if(sign.opt==1){		
			ts <- max(abs((Q.in <= 0)*Q.in))  		# test statistics
			zz <- apply(abs((g.in <= 0)*g.in),1,max)
		}
		if(sign.opt==2){
			ts <- max(abs((Q.in >= 0)*Q.in))  		# test statistics
			zz <- apply(abs((g.in >= 0)*g.in),1,max)			
		}
		za <- quantile(zz,probs=alpha)				# critical values
			tsall <- c(tsall,ts)
			zaall <- cbind(zaall,za)
		}
	return(list(test.stat = tsall, cr.value = zaall))
}

wald.sim.15 <- function(d,Qp,Qm,fxp,fxm,f0,Q2.p,Q2.m,bandw,taus,ba,n,n.sim,testing=1){
	# calculates critical values by simulation
	# estimate f0 and N_{b}(x) simultaneously
	# Wald tests with no bias correction
	fxa <- 0.5*(fxp+fxm) 
	if (testing==1){
		rtp  <- rep(1,length(taus))
		rtm  <- rep(1,length(taus))
    }
	if (testing==0){
		rtp<-fxa/fxp
      	rtm<-fxa/fxm
    }
	zs  <- array(0,c(n.sim,length(taus),2))
	umat <- matrix(runif(n*n.sim),ncol=n.sim)
	for(i in 1:length(taus)){
		t  <- taus[i]
		x1 <- cbind(1,(ba$zi/bandw[i]))
		l1 <- cbind(rep(1,n),rep(1,n))
		p1.p <- solve((t(x1) %*% diag(d*(ba$K.x[,i])) %*% x1)/(n*bandw[i]))
		p1.m <- solve((t(x1) %*% diag((1-d)*(ba$K.x[,i])) %*% x1)/(n*bandw[i]))
		for(j in 1:n.sim){
			uu <- umat[,j]
			pp <- (t-(uu-t<0))*(ba$K.x[,i])
			p2.p <- c(crossprod((pp*d),x1[,1]),crossprod((pp*d),x1[,2]))/sqrt(n*bandw[i])
			p2.m <- c(crossprod((pp*(1-d)),x1[,1]),crossprod((pp*(1-d)),x1[,2]))/sqrt(n*bandw[i])
			zs[j,i,1] <- rtp[i]*(p1.p %*% p2.p)[1]
			zs[j,i,2] <- rtm[i]*(p1.m %*% p2.m)[1]
		}
	}
	g1 <- zs[,,1]
	g2 <- zs[,,2]
	Q.d <- Qp - Qm
	bias.adj <- (1/2)*(-0.1157)*(Q2.p - Q2.m)*(bandw^2)
	g.in <- (g1 - g2)

	return(list(Q.diff = Q.d, G.diff = g.in, bias.cor = bias.adj))
}

wald.sim.rbs.15 <- function(d,Qp,Qm,fxp,fxm,f0,Q2.p,Q2.m,bandw,bandw2,taus,ba,ba2,n,n.sim,testing=1){
	# calculates critical values by simulation
	# Robust Wald test
    fxa <-0.5*(fxp+fxm)
	if (testing==1){
		rtp <- rep(1,length(taus))
		rtm <- rep(1,length(taus))
    }
	if (testing==0){
		rtp<-fxa/fxp
      	rtm<-fxa/fxm
    }
	zs  <- array(0,c(n.sim,length(taus),2))
	umat <- matrix(runif(n*n.sim),ncol=n.sim)
	for(i in 1:length(taus)){
		t  <- taus[i]
		x1 <- cbind(1,(ba$zi/bandw[i]))
		x2 <- cbind(1,(ba2$zi/bandw2[i]),((ba2$zi/bandw2[i])^2))
		p1.p <- solve((t(x1) %*% diag(d*(ba$K.x[,i])) %*% x1)/(n*bandw[i]))
		p1.m <- solve((t(x1) %*% diag((1-d)*(ba$K.x[,i])) %*% x1)/(n*bandw[i]))
		q1.p <- solve((t(x2) %*% diag(d*(ba2$K.x[,i])) %*% x2)/(n*bandw2[i]))
		q1.m <- solve((t(x2) %*% diag((1-d)*(ba2$K.x[,i])) %*% x2)/(n*bandw2[i]))
		for(j in 1:n.sim){
			uu <- umat[,j]
			pp <- (t-(uu-t<0))*(ba$K.x[,i])
			p2.p <- c(crossprod((pp*d),x1[,1]),crossprod((pp*d),x1[,2]))/sqrt(n*bandw[i])
			p2.m <- c(crossprod((pp*(1-d)),x1[,1]),crossprod((pp*(1-d)),x1[,2]))/sqrt(n*bandw[i])
			qq   <- (t-(uu-t<0))*(ba2$K.x[,i])
			q2.p <- c(crossprod((qq*d),x2[,1]),crossprod((qq*d),x2[,2]),crossprod((qq*d),x2[,3]))/sqrt(n*bandw2[i])
			q2.m <- c(crossprod((qq*(1-d)),x2[,1]),crossprod((qq*(1-d)),x2[,2]),crossprod((qq*(1-d)),x2[,3]))/sqrt(n*bandw2[i])
			zs[j,i,1] <- rtp[i]*(p1.p %*% p2.p)[1] - rtp[i]*((bandw[i]/bandw2[i])^{5/2})*(-0.1157)*(q1.p %*% q2.p)[3]
			zs[j,i,2] <- rtm[i]*(p1.m %*% p2.m)[1] - rtm[i]*((bandw[i]/bandw2[i])^{5/2})*(-0.1157)*(q1.m %*% q2.m)[3]
		}
		}
	g1 <- zs[,,1]
	g2 <- zs[,,2]
	Q.d <- Qp - Qm
	bias.adj <- (1/2)*(-0.1157)*(Q2.p - Q2.m)*(bandw^2)
	g.in <- (g1 - g2)

	return(list(Q.diff = Q.d, G.diff = g.in, bias.cor = bias.adj))
}

wald.sim.rbs.eql.15 <- function(d,Qp,Qm,fxp,fxm,f0,Q2.p,Q2.m,bandw,bandw2,taus,ba,ba2,n,n.sim,testing=1){
	# calculates critical values by simulation
	# equality constrainted Robust Wald tests 
    fxa <- 0.5*(fxp+fxm) 
	if (testing==1){
		rtp <- rep(1,length(taus))
		rtm <- rep(1,length(taus))
    }
	if (testing==0){
		rtp <- fxa/fxp 
      	rtm <- fxa/fxm 
    }
	zs  <- array(0,c(n.sim,length(taus),2))
	p1.p <- array(0,c(2,2,length(taus))); p1.m <- p1.p
	q1.p <- array(0,c(3,3,length(taus))); q1.m <- q1.p
	for(i in 1:length(taus)){
		t  <- taus[i]
		x1 <- cbind(1,(ba$zi/bandw[i]))
		x2 <- cbind(1,(ba2$zi/bandw2[i]),((ba2$zi/bandw2[i])^2))
		p1.p[,,i] <- solve((t(x1) %*% diag(d*(ba$K.x[,i])) %*% x1)/(n*bandw[i]))
		p1.m[,,i] <- solve((t(x1) %*% diag((1-d)*(ba$K.x[,i])) %*% x1)/(n*bandw[i]))
		q1.p[,,i] <- solve((t(x2) %*% diag(d*(ba2$K.x[,i])) %*% x2)/(n*bandw2[i]))
		q1.m[,,i] <- solve((t(x2) %*% diag((1-d)*(ba2$K.x[,i])) %*% x2)/(n*bandw2[i]))		
		}
	for(j in 1:n.sim){
		uu  <- runif(n)
		D1 <- array(0,c(length(taus),2))
		D2 <- array(0,c(length(taus),2))
		for(i in 1:length(taus)){
			t <- taus[i]
			x1 <- cbind(1,(ba$zi/bandw[i]))
			x2 <- cbind(1,(ba2$zi/bandw2[i]),((ba2$zi/bandw2[i])^2))
			pp <- (t-(uu-t<0))*(ba$K.x[,i])
			p2.p <- c(crossprod((pp*d),x1[,1]),crossprod((pp*d),x1[,2]))/sqrt(n*bandw[i])
			p2.m <- c(crossprod((pp*(1-d)),x1[,1]),crossprod((pp*(1-d)),x1[,2]))/sqrt(n*bandw[i])
			qq   <- (t-(uu-t<0))*(ba2$K.x[,i])
			#q2.p <- c(crossprod((qq*d),x2[,1]),crossprod((qq*d),x2[,2]),crossprod((qq*d),x2[,3]))/(sqrt(n*bandw2[i])*fxp[i]*(bandw2[i]^{5/2}))
			q2.p <- c(crossprod((qq*d),x2[,1]),crossprod((qq*d),x2[,2]),crossprod((qq*d),x2[,3]))/(sqrt(n*bandw2[i])*(bandw2[i]^{5/2}))
			 q2.m <- c(crossprod((qq*(1-d)),x2[,1]),crossprod((qq*(1-d)),x2[,2]),crossprod((qq*(1-d)),x2[,3]))/(sqrt(n*bandw2[i])*(bandw2[i]^{5/2}))
			D1[i,1] <- rtp[i]*(p1.p[,,i] %*% p2.p)[1]
			D1[i,2] <- rtm[i]*(p1.m[,,i] %*% p2.m)[1]			
			D2[i,1] <- rtp[i]*(q1.p[,,i] %*% q2.p)[3]
			D2[i,2] <- rtm[i]*(q1.m[,,i] %*% q2.m)[3]
		}
		zs[j,,1] <- D1[,1] - mean(D2[,1])*(bandw^{5/2})*(-0.1157)
		zs[j,,2] <- D1[,2] - mean(D2[,2])*(bandw^{5/2})*(-0.1157)	
	}
	g1 <- zs[,,1]
	g2 <- zs[,,2]
	Q.d <- Qp - Qm
	bias.adj <- (1/2)*(-0.1157)*(Q2.p - Q2.m)*(bandw^2)
	g.in <- (g1 - g2)

	return(list(Q.diff = Q.d, G.diff = g.in, bias.cor = bias.adj))
}

# 6. Functions for the Score and Wald tests when data is big (used in application)

uci.rd.big <- function(x,xt,alpha,ba,bandw,taus,n,n.sim,left=1){
	zs  <- array(0,c(n.sim,length(taus)))
	for(j in 1:n.sim){
	uu <- runif(n)
	for(i in 1:length(taus)){
		t  <- taus[i]
		if(left==1) {pa2 <- sum((t-(uu-t<0))*ba$K.x[,i]*((x < xt)-0.5+(15/16)*(ba$zi/bandw[i])))/sqrt(n*bandw[i])}
		if(left==0) {pa2 <- sum((t-(uu-t<0))*ba$K.x[,i]*((x >=xt)-0.5-(15/16)*(ba$zi/bandw[i])))/sqrt(n*bandw[i])}			
		zs[j,i] <- abs(pa2)
		}
	}
	zz <- apply(zs,1,max)
	z.a <- quantile(zz,probs=alpha)
	return(z.a)
}

wald.sim.big <- function(d,Qp,Qm,fxp,fxm,f0,Q2.p,Q2.m,bandw,taus,ba,n,n.sim,testing){
	fxa <- 0.5*(fxp+fxm) 
	if(testing==1){
		rtp <- rep(1,length(taus))
		rtm <- rep(1,length(taus))		
	}
	if(testing==0){
		rtp <- fxa/fxp
		rtm <- fxa/fxm		
	}	
	zs  <- array(0,c(n.sim,length(taus),2))
	for(k in 1:2){
	for(j in 1:n.sim){
		uu <- runif(n)
		for(i in 1:length(taus)){
			t  <- taus[i]
			if(k==1) {zs[j,i,k] <- rtp[i]*sum((t-(uu-t<0))*(6.736842-12.63158*(ba$zi/bandw[i]))*d*ba$K.x[,i])/(sqrt(n*bandw[i])*f0)}
			if(k==2) {zs[j,i,k] <- rtm[i]*sum((t-(uu-t<0))*(6.736842+12.63158*(ba$zi/bandw[i]))*(1-d)*ba$K.x[,i])/(sqrt(n*bandw[i])*f0)}
		}
	}
	}
	g1 <- zs[,,1]
	g2 <- zs[,,2]
	Q.d <- Qp - Qm
	bias.adj <- (1/2)*(-0.1157)*(Q2.p - Q2.m)*(bandw^2)
	g.in <- (g1 - g2)

	return(list(Q.diff = Q.d, G.diff = g.in, bias.cor = bias.adj))
}

wald.sim.big.rbs <- function(d,Qp,Qm,fxp,fxm,f0,Q2.p,Q2.m,bandw,taus,ba,n,n.sim,testing){
	# robust Wald procedure
	fxa <- 0.5*(fxp+fxm) 
	if(testing==1){
		rtp <- rep(1,length(taus))
		rtm <- rep(1,length(taus))		
	}
	if(testing==0){
		rtp <- fxa/fxp
		rtm <- fxa/fxm		
	}	
	zs  <- array(0,c(n.sim,length(taus),2))
	for(k in 1:2){
	for(j in 1:n.sim){
		uu <- runif(n)
		for(i in 1:length(taus)){
			t  <- taus[i]
			if(k==1) {zs[j,i,k] <- rtp[i]*sum((t-(uu-t<0))*(14.16061-66.62360*(ba$zi/bandw[i])+64.11594*((ba$zi/bandw[i])^{2}))*d*ba$K.x[,i])/(sqrt(n*bandw[i])*f0)}
			if(k==2) {zs[j,i,k] <- rtm[i]*sum((t-(uu-t<0))*(14.16606+66.66379*(ba$zi/bandw[i])+64.16390*((ba$zi/bandw[i])^{2}))*(1-d)*ba$K.x[,i])/(sqrt(n*bandw[i])*f0)}
		}
	}
	}
	g1 <- zs[,,1]
	g2 <- zs[,,2]
	Q.d <- Qp - Qm
	bias.adj <- (1/2)*(-0.1157)*(Q2.p - Q2.m)*(bandw^2)
	g.in <- (g1 - g2)

	return(list(Q.diff = Q.d, G.diff = g.in, bias.cor = bias.adj))
}

wald.sim.big.rbs.eql <- function(d,Qp,Qm,fxp,fxm,f0,Q2.p,Q2.m,bandw,bandw2,taus,ba,ba2,n,n.sim,testing){
	# Equality constraint is imposed on bias estimation. 
	fxa <- 0.5*(fxp+fxm) 
	if(testing==1){
		rtp <- rep(1,length(taus))
		rtm <- rep(1,length(taus))		
	}
	if(testing==0){
		rtp <- fxa/fxp
		rtm <- fxa/fxm		
	}	
	zs  <- array(0,c(n.sim,length(taus),2))
	for(j in 1:n.sim){
		uu  <- runif(n)
		zs1 <- array(0,c(length(taus),2))
		zs2 <- array(0,c(length(taus),2))
		for(i in 1:length(taus)){
			t <- taus[i]
			zs1[i,1] <- rtp[i]*sum((t-(uu-t<0))*(6.736842-12.63158*(ba$zi/bandw[i]))*d*ba$K.x[,i])/(sqrt(n*bandw[i])*f0) 
			zs2[i,1] <- rtp[i]*sum((t-(uu-t<0))*(64.16-466.66*(ba2$zi/bandw2[i])+554.16*((ba2$zi/bandw2[i])^{2}))*d*ba2$K.x[,i])/(sqrt(n*bandw2[i])*f0*(bandw2[i]^{5/2}))
			zs1[i,2] <- rtm[i]*sum((t-(uu-t<0))*(6.736842+12.63158*(ba$zi/bandw[i]))*(1-d)*ba$K.x[,i])/(sqrt(n*bandw[i])*f0)
			zs2[i,2] <- rtm[i]*sum((t-(uu-t<0))*(64.16+466.66*(ba2$zi/bandw2[i])+554.16*((ba2$zi/bandw2[i])^{2}))*(1-d)*ba2$K.x[,i])/(sqrt(n*bandw2[i])*f0*(bandw2[i]^{5/2}))
		}
		zs[j,,1] <- zs1[,1] - mean(zs2[,1])*(bandw^{5/2})*(-0.1157)
		zs[j,,2] <- zs1[,2] - mean(zs2[,2])*(bandw^{5/2})*(-0.1157)
	}
	g1 <- zs[,,1]
	g2 <- zs[,,2]
	Q.d <- Qp - Qm
	bias.adj <- (1/2)*(-0.1157)*(Q2.p - Q2.m)*(bandw^2)
	g.in <- (g1 - g2)

	return(list(Q.diff = Q.d, G.diff = g.in, bias.cor = bias.adj))
}

# 7. Functions for bandwidth (used in application)

rd.cv.emp <- function(x,y,x.eval,xl,tau,kr){
		h.med <- cv.rd.emp(x,y,x.eval,xl,kr)$h.cv
		h.tau <- h.med*((2*tau*(1-tau)/(pi*dnorm(qnorm(tau))^{2}))^{1/(4+1)})
	return(h.tau)	
}

rd.cv.emp.int <- function(x,y,x.eval,xl,tau,kr){ # interior point cv
		h.med <- cv.rd.emp.int(x,y,x.eval,xl,kr)$h.cv
		h.tau <- h.med*((2*tau*(1-tau)/(pi*dnorm(qnorm(tau))^{2}))^{1/(4+1)})
	return(h.tau)	
}

# Cross-validation for the pilot bandwidth at the median

cv.rd.emp <- function(x,y,x.eval,val,xl,kr){
	# cv bandwidth treating the cutoff as a boundary point
	x.u <- sort(unique(x))
	wdt <- quantile(abs(x-x.eval),probs=xl)
	x.m <- x.u[(abs(x.u-x.eval)<wdt)]	
	cri <- array(0,c(length(val),1))
	for (i in 1:length(val)){
		h.tau <- val[i]
		cri.h <- NULL
		for (j in 1:length(x.m)){
			xe <- x.m[j]
			if(xe >= x.eval){sgn <- ((x.u > x.eval) & (x.u > xe) & (x.u <= xe+h.tau))}
			if(xe <  x.eval){sgn <- ((x.u < x.eval) & (x.u < xe) & (x.u >= xe-h.tau))}
			xx <- x[(x%in%x.u[sgn])]
			yy <- y[(x%in%x.u[sgn])]
			Qy <- lprq.pro.rdd(xx,yy,tau=0.5,xe,h.tau,mon=0,kr)$b0
			cri.h <- c(cri.h,sum(abs(y[(x==xe)]-Qy)))
        }
	cri[i] <- mean(cri.h)
	}
	return(list(h.cv = val[which.min(cri)], cand = cbind(val,cri)))
}

cv.rd.emp.int <- function(x,y,x.eval,val,xl,kr){
	# cv bandwidth treating the cutoff as an interior point
	x.u <- sort(unique(x))
	wdt <- quantile(abs(x-x.eval),probs=xl)
	x.m <- x.u[(abs(x.u-x.eval)<wdt)]
	cri <- array(0,c(length(val),1))
	for (i in 1:length(val)){
		h.tau <- val[i]
		cri.h <- NULL
		for (j in 1:length(x.m)){
			xe <- x.m[j]
			xx <- x[-j]
			yy <- y[-j]
			Qy <- lprq.pro.rdd(xx,yy,tau=0.5,xe,h.tau,mon=0,kr)$b0
			cri.h <- c(cri.h,abs(y[j]-Qy))
        }
	cri[i] <- mean(cri.h)
	}
	return(list(h.cv = val[which.min(cri)], cand = cbind(val,cri)))
}

band.IK.3.big <- function(x,d,x0,n,fx,fqx.p,fqx.m,Q.2.p,Q.2.m,hp,hm){
	# adapted from band.IK.3()
	rp <- (3*266.6277)/(n*(hp^5)*(fqx.p^2)*fx)
	rm <- (3*266.6277)/(n*(hm^5)*(fqx.m^2)*fx)
	con1 <- 0.5*0.5*4.497*((fqx.p^{-2})+(fqx.m^{-2}))
	con2 <- 0.0134*fx*(((Q.2.p-Q.2.m)^{2})+rp+rm)
	h.med <- (con1/con2)^{1/5}*(n^{-1/5})
	return(list(h.ik = as.numeric(h.med)))	
}

# 8. Miscellaneous Functions

kxb <- function(x,y,x.eval,tau,n,d.x,h.tau,kr=1){
	# produce weights and centred running variable
	if(is.matrix(x)==0){
		z  <- x - x.eval
	    	if(kr==1){wx <- dnorm((1/h.tau)%x%z)}
    		if(kr==2){wx <- dcauchy((1/h.tau)%x%z)}
	    if(kr==3){wx <- depa((1/h.tau)%x%z)}
	    if(kr==4){wx <- dtnorm((1/h.tau)%x%z)}
		K.x <- matrix(wx,nrow=n,ncol=length(tau),byrow=FALSE)
		}
	if(is.matrix(x)==1){
		z    <- x - matrix(rep(x.eval,each=n),ncol=length(x.eval))
		if(kr==1){wx <- dnorm((1/h.tau)%x%z[,1])*dnorm((1/h.tau)%x%z[,2])}
		if(kr==2){wx <- dcauchy((1/h.tau)%x%z[,1])*dcauchy((1/h.tau)%x%z[,2])}
		if(kr==3){wx <- depa((1/h.tau)%x%z[,1])*depa((1/h.tau)%x%z[,2])}
		if(kr==4){wx <- dtnorm((1/h.tau)%x%z[,1])*dtnorm((1/h.tau)%x%z[,2])}
		K.x  <- matrix(wx,nrow=n,ncol=length(tau),byrow=FALSE)
		}
	return(list(K.x = K.x, zi = z))
}

condf <- function(condQ,tt,nm,ad=2,kr=1){
	# conditional density estimator	
	ts <- c(0.00001,tt,0.99999)
	Qs <- sort(c(qnorm(0.00001),condQ,qnorm(0.99999)))
	if(kr==1){
		con.Y<- approx(x=ts,y=Qs,xout=runif(nm),method="linear",rule=2)$y
		cond <- density(con.Y, kernel = "gaussian", adjust=ad)
		den  <- approx(x=cond$x,y=cond$y,xout=condQ,method="linear",rule=2)$y
		}
	if(kr==2){
		a1 <- approx(x=ts,y=Qs,xout=runif(nm),method="linear",rule=2)$y
		h1 <- bw.nrd(a1)
		den<- sapply(condQ,est.f.c,a1,h1)
		}	
	if(kr==3){
		con.Y<- approx(x=ts,y=Qs,xout=runif(nm),method="linear",rule=2)$y
		cond <- density(con.Y, kernel = "epanechnikov", adjust=ad)
		den  <- approx(x=cond$x,y=cond$y,xout=condQ,method="linear",rule=2)$y
		}
	if(kr==4){
		a1 <- approx(x=ts,y=Qs,xout=runif(nm),method="linear",rule=2)$y
		h1 <- bw.nrd(a1)*ad
		den<- sapply(condQ,est.f.tn,a1,h1)
		}
	return(den)
}

depa <- function(xx,loc=0,scale=1){
	# Epanechnikov kernel
      nx=(xx-loc)/scale
	return((scale^{-1})*(3/4)*(1-nx^2)*(abs(nx)<1))
}

marf <- function(xx,xeval){
	# estimating marginal density	
	if(is.matrix(xx)==0){
		mar.x <- density(xx)
		den   <- mar.x$y[which.min(abs(mar.x$x - xeval))]
		}
	if(is.matrix(xx)==1){
		bw    <- h.select(x,y,method="cv")
		mar.x <- bkde2D(x, range.x = list(c(0,1),c(0,1)), bandwidth=bw)
		den   <- mar.x$fhat[which.min(abs(mar.x$x1-xeval[1])),which.min(abs(mar.x$x2-xeval[2]))]
		}
	return(den)
}

band.covt <- function(h.med,ttau){
	# use Yu and Jones (1998) rule
	# calculate bandwidths at the specified quantiles from the median bandwidth
	h.tau <- h.med*((2*ttau*(1-ttau)/(pi*dnorm(qnorm(ttau))^{2}))^{1/(4+1)})
	return(h.tau)
}

ext.tau <- function(tau){
	# extend the quantile range to estimate conditional densities
	te <- tau[1]
	tau.add <- seq(te*0.1,te*0.6,length.out=3)
	tau.ext <- unique(sort(c(tau,seq(0.1,0.9,by=0.05),tau.add,1-tau.add)))
	return(tau.ext)
}

rearrange2 <- function(Q,tt,n){
	# rearrangement
	o   <- sort(unique(c(seq(tt[1],tt[length(tt)],length.out=n),tt)))
	Qpe <- approx(x=tt,y=Q,xout=o,method="linear",rule=2)$y
	Qre <- approx(x=o,y=sort(Qpe),xout=tt,method="linear",rule=2)$y
	return(Qre)
	}

dtgen.new <- function(x,d,model,type,ns){
	# generate response variables, used in simulation
	# type=1 then under the null, type=2 then under the alternative
	n.sam<- length(x)	
	rank <- runif(n.sam)
	rb1  <- qnorm(rank)
	if(model==1){
		y <- 1 + x + (0.5+0.3*x)*rb1
		if(type==2) {
		lp <- ns*0.5		
		yp <- 1 + x + (0.5 + 0.3*x)*(rb1+3.3*lp*pmax(rank-0.2, 0))
		y[d==1] <- yp[d==1]
		}
	}
	if(model==2){
		y   <- (0.5 + x + (x^2) + sin(pi*x-1.0)) + (x+1.25)*rb1
		if(type==2) {
		lp <- ns*1.25
		yp <- (0.5+x+(x^2)+sin(pi*x-1.0))+(x+1.25)*(rb1+1.3*lp*pmax(rank-0.2, 0))
		y[d==1] <- yp[d==1]
		}
	}
	if(model==3){
		y <- (1-d)*(0.48+1.27*x+7.18*(x^2)+20.21*(x^3)+21.54*(x^4)+7.33*(x^5)) + d*(0.48+0.84*x-3*(x^2)+7.99*(x^3)-9.01*(x^4)+3.16*(x^5)) + 0.1295*rb1
	if(type==2) {
		lp <- ns*0.1295
		yp <- 0.48+0.84*x-3*(x^2)+7.99*(x^3)-9.01*(x^4)+3.16*(x^5) + 0.1295*(rb1+12.9*lp*pmax(rank-0.2, 0))
		y[d==1] <- yp[d==1]
		}
	}
	if(model==4){
		y   <- (1-d)*3*(x^2)+d*4*(x^2) + 0.1295*rb1
		if(type==2) {
		lp <- ns*0.1295
		yp <- 4*(x^2)+0.1295*(rb1+12.9*lp*pmax(rank-0.2, 0))
		y[d==1] <- yp[d==1]
		}
	}
	return(y)
}

## 9. Functions from Bondell, Reich, and Wang (2010, Biometrika) 
## available at http://www4.stat.ncsu.edu/~bondell/Software/NoCross/NoCrossQuant.R

rq.no.cross = function(y, x, taus)
{
    require(quantreg)
     
    if (length(taus)<2)
	  stop("At least 2 quantile levels should be specified. If only using a single quantile, you should use rq.")

    taus = sort(taus)
    if (max(taus)>=1)
	  stop("All quantile levels should be between 0 and 1, not including the boundary.")
    if (min(taus)<=0)
	  stop("All quantile levels should be between 0 and 1, not including the boundary.")

    n = nrow(x)
    p = ncol(x)
    m = length(taus)
    pp = p+1
    Y = rep(y, m)
    xtemp = x

    x = matrix(0,nrow=n,ncol=p)
    shifts = apply(xtemp,2,min)
    scalings = rep(0,p)
    for (i in 1:p)
    {
    x[,i] = xtemp[,i] - shifts[i]
	scalings[i] = max(x[,i])
	x[,i] = x[,i]/scalings[i]
    }

    x = cbind(rep(1,n),x)

    D = diag(m)
    D[lower.tri(D)] = 1 
    X = D %x% x
    X2 = cbind(X, -X)
    sX = as.matrix.csr(X2)
    K = m*pp
    R1 = (diag(m) %x% rep(1,1) %x% t(c(1, rep(0, p))))[-(1:1),]
    R2 = rep(1,1)
    for (j in 1:p)
    {
    	R2 = cbind(R2, rep(1,1) %x% (diag(1) %x% rep(1,1)))
    }

    R2 = (diag(m) %x% R2)[-(1:1),]
    
    sR = as.matrix.csr(rbind(diag(2*K), cbind(R1, -R2)))
    r2 = rep(0, 2*K + (m-1)*(1))
       
    tau.v = rep(taus, each=n)
    rhs = t(X2)%*%(1-tau.v)
    coeff1 =  myrq.fit.sfnc(sX, Y, sR, r2, tau=tau.v, rhs=rhs, tmpmax=100000)$coef
    coeff = coeff1[1:K]-coeff1[-(1:K)]
    gamma.m = matrix(coeff, ncol=m)
    D = diag(m)
    D[upper.tri(D)]=1
    bhat.temp = gamma.m %*% D
    cov.bhat = array(0,c(pp,pp,m))
    se.bhat = matrix(0,ncol=m,nrow=pp)
    transform.mat = rbind(c(1,-shifts/scalings),cbind(rep(0,p),diag(as.vector(1/scalings))))
    bhat = transform.mat%*%bhat.temp
    for (j in 1:m)
    {
	cov.bhat[,,j] = se.constr(x,y,bhat.temp[,j],taus[j])
   	cov.bhat[,,j] = transform.mat%*%cov.bhat[,,j]%*%t(transform.mat)
	se.bhat[,j] = sqrt(diag(cov.bhat[,,j])) 
    }
	
    vars<-c("intercept",paste("x",1:p,sep=""))
    rownames(bhat)<-rownames(se.bhat)<-vars
    colnames(bhat)<-colnames(se.bhat)<-taus
    dimnames(cov.bhat)<-list(vars,vars,taus) 
	
    constr.qr.fit = NULL
    constr.qr.fit$bhat = bhat
    constr.qr.fit$se.bhat = se.bhat
    constr.qr.fit$cov.bhat = cov.bhat
    constr.qr.fit$taus = taus
    
    return(constr.qr.fit)
}


myrq.fit.sfnc <- function (x, y, R, r, tau, rhs, nsubmax, tmpmax, nnzlmax, cachsz = 64, 
    small = 1e-08, maxiter = 100, warn.mesg = TRUE) 
{
    require(quantreg)
    y <- -y
    r <- -r
    n1 <- length(y)
    m <- x@dimension[2]
    if (n1 != x@dimension[1]) 
        stop("The design matrix A1' and response vector y are not compatible")
    n2 <- length(r)
    if (n2 != R@dimension[1]) 
        stop("The constraint matrix A2' and the constraint right-hand-side are not compatible")
    maxn1n2 <- max(n1, n2)
    u <- rep(1, length = n1)
    if(length(tau)==1)     x1 <- rep(1 - tau, length = n1)
    else x1=1-tau
    x2 <- rep(1, length = n2)
    wwm <- vector("numeric", 6 * m)
    wwm[1:m] <- rhs
    nnzx <- x@ia[x@dimension[1] + 1] - 1
    nnzR <- R@ia[R@dimension[1] + 1] - 1
    nnzdmax <- max(nnzx, nnzR)
    iwmax <- 7 * m + 3
    ao1 <- t(x)
    ao2 <- t(R)
    e <- ao1 %*% x
    g <- ao2 %*% R
    h <- e + g
    nnzemax <- e@ia[e@dimension[1] + 1] - 1
    nnzgmax <- g@ia[g@dimension[1] + 1] - 1
    nnzhmax <- h@ia[h@dimension[1] + 1] - 1
    if (missing(nnzlmax)) 
        nnzlmax <- 4 * nnzdmax
    if (missing(nsubmax)) 
        nsubmax <- nnzhmax
    if (missing(tmpmax)) 
        tmpmax <- 6 * m
    s <- u - x1
    chol.o <- chol(e, tmpmax = tmpmax, nsubmax = nsubmax, nnzlmax = nnzlmax)
    b <- backsolve(chol.o, ao1 %*% y)
    r1 <- y - x %*% b
    z1 <- ifelse(abs(r1) < small, (r1 * (r1 > 0) + small), r1 * 
        (r1 > 0))
    w <- z1 - r1
    z2 <- rep(1, n2)
    wwn1 <- matrix(0, n1, 10)
    wwn1[, 1] <- z1
    wwn1[, 2] <- w
    wwn2 <- matrix(0, n2, 7)
    wwn2[, 2] <- z2
    srqfnc.o <- .Fortran("srqfnc", n1 = as.integer(n1), m = as.integer(m), 
        nnzx = as.integer(nnzx), x = as.double(x@ra), jx = as.integer(x@ja), 
        ix = as.integer(x@ia), ao1 = as.double(ao1@ra), jao1 = as.integer(ao1@ja), 
        iao1 = as.integer(ao1@ia), n2 = as.integer(n2), nnzR = as.integer(nnzR), 
        R = as.double(R@ra), jR = as.integer(R@ja), iR = as.integer(R@ia), 
        ao2 = as.double(ao2@ra), jao2 = as.integer(ao2@ja), iao2 = as.integer(ao2@ia), 
        nnzdmax = as.integer(nnzdmax), d = double(nnzdmax), jd = integer(nnzdmax), 
        id = integer(m + 1), dsub = double(nnzhmax + 1), jdsub = integer(nnzhmax + 
            1), nnzemax = as.integer(nnzemax), e = as.double(e@ra), 
        je = as.integer(e@ja), ie = as.integer(e@ia), nnzgmax = as.integer(nnzgmax), 
        g = double(nnzgmax), jg = integer(nnzgmax), ig = integer(m + 
            1), nnzhmax = as.integer(nnzhmax), h = double(nnzhmax), 
        jh = integer(nnzhmax), ih = integer(m + 1), nsubmax = as.integer(nsubmax), 
        lindx = integer(nsubmax), xlindx = integer(m + 1), nnzlmax = as.integer(nnzlmax), 
        lnz = double(nnzlmax), xlnz = integer(m + 1), iw = integer(m * 
            5), iwmax = as.integer(iwmax), iwork = integer(iwmax), 
        xsuper = integer(m + 1), tmpmax = as.integer(tmpmax), 
        tmpvec = double(tmpmax), maxn1n2 = as.integer(maxn1n2), 
        ww1 = double(maxn1n2), wwm = as.double(wwm), wwn1 = as.double(wwn1), 
        wwn2 = as.double(wwn2), cachsz = as.integer(cachsz), 
        level = as.integer(8), x1 = as.double(x1), x2 = as.double(x2), 
        s = as.double(s), u = as.double(u), c1 = as.double(y), 
        c2 = as.double(r), sol = as.double(b), small = as.double(small), 
        ierr = integer(1), maxiter = as.integer(maxiter), time = double(7), 
        PACKAGE = "quantreg")[c("sol", "ierr", "maxiter", "time")]
    ierr <- srqfnc.o$ierr
    if (ierr == 13) 
        stop("Increase nnzh.factor")
    if (!(ierr == 0) && warn.mesg) 
        warning(sfnMessage(ierr))
    list(coef = -srqfnc.o$sol, ierr = ierr, it = srqfnc.o$maxiter, 
        time = sum(srqfnc.o$time))
}

se.constr = function(x,y,coef,tau)
{
        require(quantreg)
        n = nrow(x)
        p = ncol(x)
        h = bandwidth.rq(tau, n)
        if (tau + h > 1) 
            stop("tau + h > 1:  error in summary.rq")
        if (tau - h < 0) 
            stop("tau - h < 0:  error in summary.rq")
        uhat = c(y - x %*% coef)
        h = (qnorm(tau + h) - qnorm(tau - h)) * min(sqrt(var(uhat)), 
            (quantile(uhat, 0.75) - quantile(uhat, 0.25))/1.34)
        f = dnorm(uhat/h)/h
        fxxinv = diag(p)
        fxxinv = backsolve(qr(sqrt(f) * x)$qr[1:p, 1:p, drop = FALSE], 
            fxxinv)
        fxxinv = fxxinv %*% t(fxxinv)
        cov = tau * (1 - tau) * fxxinv %*% crossprod(x) %*% fxxinv
        cov
}
