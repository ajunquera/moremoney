bwd1o1_probando <- rdd.bandwidth(running1, outcome1, d1,
                        x.eval = 0, tt = c(0.25,0.75), m=3, method = c(3,4,5), 
                        val=c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08), kr = 3)
# 9:34, 9:54


bwd2preo1ik <- rdd.bandwidth(running2rmt1, outcomepre1_d2, d2rmt1,
                             x.eval = 0, tt = c(0.25,0.75), m=3, method = c(3,4,5), 
                             val=c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08), kr = 3)

# rdd.bandwidth ------------
h.tau1_afj <- cv.med.rd(running2rmt1,outcomepre1_d2,0,c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08),xl=0.5,3)$h.cv # AQUÍ ESTÁ EL PROBLEMA !!!!!!!!!
taus_afj <- seq(0.25,0.75,length.out=3)
taus_afj <- sort(unique(c(taus_afj,0.5)))
taul_afj <- ext.tau(taus_afj)
#ind  <- taul%in%taus
where_afj <- taul_afj%in%0.5
fx_afj <- marf(running2rmt1,0)
h.tau <- band.covt(h.tau1,taul)

Q.0.p_afj <- lprq.pro.rdd(running2rmt1[d2rmt1==T],outcomepre1_d2[d2rmt1==T],taul_afj,0,h.tau=h.tau,mon=1,kr)
Q.0.m_afj <- lprq.pro.rdd(running2rmt1[d2rmt1==F],outcomepre1_d2[d2rmt1==F],taul_afj,0,h.tau=h.tau,mon=1,kr)
fqx.p_afj <- condf(Q.0.p_afj$b0,taul,5000,ad=2,kr)
fqx.m_afj <- condf(Q.0.m_afj$b0,taul,5000,ad=2,kr)

h2.p_afj <- abs(max(running2rmt1)-0)
h2.m_afj <- 0.5*abs(0-min(running2rmt1))

Q.2.p_afj <- local.poly(outcomepre1_d2[d2rmt1==T],running2rmt1[d2rmt1==T],0,taus=0.5,bandw3=h2.p_afj,order=2)$Q.2
Q.2.m_afj <- local.poly(outcomepre1_d2[d2rmt1==F],running2rmt1[d2rmt1==F],0,taus=0.5,bandw3=h2.m_afj,order=2)$Q.2		

h.tau5 <- band.IK.3(running2rmt1,d2rmt1,0,length(outcomepre1_d2),
                    fx_afj,fqx.p[where_afj],fqx.m[where_afj], # faltan fqx.p y fqx.m
                    Q.2.p_afj,Q.2.m_afj,h2.p_afj,h2.m_afj,order=2)$h.ik

# cv.med.rd -----------
dd_afj  <- (running2rmt1>0)
cut_afj <- c(quantile((running2rmt1[dd_afj==1]-0),prob=0.5),quantile((running2rmt1[dd_afj==0]-0),prob=(1-0.5)))
index_afj <- which((running2rmt1-0 <= cut_afj[1]) & (running2rmt1-0 >= cut_afj[2]))
cri_afj <- array(0,c(length(c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08)),1))

val_afj <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08)

for (i in 1:length(val_afj)){
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

## Bucle manual para i
h.tau_1 <- val_afj[1]
cri.h_1 <- NULL

for (j in index_afj){
  xx=running2rmt1[-j]
  yy=outcomepre1_d2[-j]
  if(running2rmt1[j]>=0){sgn <- ((xx >= 0) & (xx>running2rmt1[j]))}
  if(running2rmt1[j]< 0){sgn <- ((xx <  0) & (xx<running2rmt1[j]))}
  Qy <- lprq.pro.rdd(xx[sgn],yy[sgn],tau=0.5,running2rmt1[j],h.tau_1,mon=0,3)$b0
  cri.h_1 <- c(cri.h_1,abs(outcomepre1_d2[j]-Qy))
  
  mensaje <- paste("Iteración", j, "completada con éxito.")
  cat(mensaje, "\n")
}

### Bucle manual para j
xx_1=running2rmt1[-1]
yy_1=outcomepre1_d2[-1]
if(running2rmt1[1]>=0){sgn <- ((xx_1 >= 0) & (xx_1>running2rmt1[1]))}
if(running2rmt1[1]< 0){sgn <- ((xx_1 <  0) & (xx_1<running2rmt1[1]))}
Qy_1 <- lprq.pro.rdd(xx_1[sgn],yy_1[sgn],tau=0.5,running2rmt1[1],h.tau_1,mon=0,3)$b0
cri.h_1 <- c(cri.h_1,abs(outcomepre1_d2[1]-Qy_1)) # works for j=1

### Observando por qué hay un problema en la iteración 3074
xx_prob=running2rmt1[-3074]
yy_prob=outcomepre1_d2[-3074]
if(running2rmt1[3074]>=0){sgn <- ((xx_prob >= 0) & (xx_prob>running2rmt1[3074]))}
if(running2rmt1[3074]< 0){sgn <- ((xx_prob <  0) & (xx_prob<running2rmt1[3074]))}
Qy_prob <- lprq.pro.rdd(xx_prob[sgn],yy_prob[sgn],tau=0.5,running2rmt1[3074],h.tau_1,mon=0,3)$b0 # aquí está el problema
cri.h_1 <- c(cri.h_1,abs(outcomepre1_d2[1]-Qy_1)) # works for j=1

## lprq.pro.rdd --------------
x = xx_prob[sgn]
y = yy_prob[sgn]

n_pro  <- length(y)
fv_pro <- NULL
z_pro  <- x - 0

wx <- depa((1/h.tau_1)%x%z_pro)

if(kr==3){
  hk <- which.max(h.tau_1)		
  ind.0 <- (wx[(1+(hk-1)*n_pro):(hk*n_pro)]==0)
  y  <- yy_prob[ind.0==0]
  z  <- z_pro[(ind.0==0)]
  wx <- wx[rep(ind.0,length(h.tau_1))==0]
}


if(mon==0){r <- NULL
for(k in 1:length(tau)){
  rt <- rq(y ~ z, tau = 0.5, w = as.vector(wx[(1+(1-1)*length(y)):(1*length(y))]))$coef
  r  <- cbind(r,as.numeric(rt))
}
fv <- cbind(fv,r)}
b0 <- fv[1,]
b1 <- fv[2,]
return(list(b0 = b0, b1 = b1))


#### Observando por qué NO hay problema en 3073
### Observando por qué hay un problema en la iteración 3074
xx_noprob=running2rmt1[-3073]
yy_noprob=outcomepre1_d2[-3073]
if(running2rmt1[3073]>=0){sgn_np <- ((xx_noprob >= 0) & (xx_noprob>running2rmt1[3073]))}
if(running2rmt1[3073]< 0){sgn_np <- ((xx_noprob <  0) & (xx_noprob<running2rmt1[3073]))}

Qy_prob <- lprq.pro.rdd(xx_prob[sgn],yy_prob[sgn],tau=0.5,running2rmt1[3074],h.tau_1,mon=0,3)$b0 # aquí está el problema
cri.h_1 <- c(cri.h_1,abs(outcomepre1_d2[1]-Qy_1)) # works for j=1

# lprq.pro.rdd
x_np = xx_noprob[sgn_np]
y_np = yy_noprob[sgn_np]

n_nopro  <- length(y_np)
fv_nopro <- NULL
z_nopro  <- x_np - 0

wx_np <- depa((1/h.tau_1)%x%z_nopro)

hk <- which.max(h.tau_1)		
ind.0_np <- (wx_np[(1+(hk-1)*n_nopro):(hk*n_nopro)]==0)
y_np  <- yy_noprob[ind.0_np==0]
z_np  <- z_pro[(ind.0==0)]
wx <- wx[rep(ind.0,length(h.tau_1))==0]


# función bw con método que asegura monotonicidad -----------
rdd.bandwidth.int <- function(x,y,d,x.eval,tt,m,method,val=NULL,band=NULL,kr=3,br=0){
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
    h.tau1 <- cv.med.rd.int(x,y,x.eval,val,xl=0.5,kr)$h.cv
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

rdd.bandwidth.AFJ <- function(x,y,d,x.eval,tt,m,method,val=NULL,band=NULL,kr=3,br=0){
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
    h.tau1 <- cv.med.rd.AFJ(x,y,x.eval,val,xl=0.5,kr)$h.cv # Aquí hay cambio en la función
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

cv.med.rd.AFJ <- function(x,y,x.eval,val,xl,kr){
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
      Qy <- lprq.pro.rdd(xx[sgn],yy[sgn],tau=0.5,x[j],h.tau,mon=T,kr)$b0 # HE CAMBIADO mon=0 por mon=1
      cri.h <- c(cri.h,abs(y[j]-Qy))
      
      mensa <- paste("Iteración", j, "para valor", i, "completada con éxito.")
      cat(mensa, "\n")
    }
    cri[i] <- mean(cri.h)
  }
  return(list(h.cv = val[which.min(cri)], cand = cbind(val,cri)))
}


rdd.bandwidth.prev <- function(x,y,d,x.eval,tt,m,method,val=NULL,band=NULL,kr=3,br=0, prev){
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
  if(max(method %in% 4) & !is.null(band) & length(band)<3)
    stop("the length of band must be 3 or larger.")
  if(max(method %in% 5) & !is.null(band) & length(band)<5)
    stop("the length of band must be 5.")
  if(max(method %in% c(1,3,4,5))){
    h.tau1 <- prev # Aquí hay cambio en la función. Se incorpora prev para indicar la hcv que resultaría si la hubieramos calculado antes. También he eliminado las dos primeras líneas. Si se añade prev no necesitamos val.
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
