## Prepared in March, 2017
## this file replicates the simulation exercise
## running it over different models, sample sizes, bandwidth choices replicates results in the simulation section

rm(list=ls())
source("qte_rdd.R")

s.sam <- c(500,1000,2000)	# sample sizes
supp  <- c(-1,1)				# support of x
x0    <- 0		# cut-off point
krn   <- 3		# choice of the kernel function, 3 for Epanechnikov kernel
type  <- 1		# 1 for the null, 2 for the alternative
ns    <- c(0.3,0.6,1.0,2.0)	# for the power calculations, varies under the alternative. sets the size of deviation from the null

L <- 2000		# number of simulation repetitions

model <- 1	# sets which model to use, Models 1 to 4

if(type==1){ num <- length(s.sam) } 	# under H0
if(type==2){ num <- length(ns) }		# under H1

for(k in 1:1){	
# loop over sample sizes (under the null), or over deviations from the null (under the alternative)

R <- array(0,c(4,2))		# rejection ratios for Score and Wald tests (for 2 nominal levels)
R2 <- array(0,c(4,2))		# rejection ratios for robust Wald tests
R3 <- array(0,c(4,2))		# rejection ratios for the robust Wald tests with constraint

if(type==1){n.sam <- s.sam[k]}
if(type==2){n.sam <- s.sam[2]}

A <- array(0,c(L,3));  A.w <- array(0,c(L,3,3)); A.w2 <- A.w;	A.w3 <- A.w
# collects test statistics and critical values. 
# A: score test, A.w: Wald test, A.w2: Robust Wald, A.w3: Robust Wald with equality constraint
B <- array(0,c(L,5))		# saves bandwidths

for(j in 1:L){	# beginning of the inner-loop

# 1. Data generation

if(model%in%c(1,2)){x <- runif(n.sam,min=supp[1],max=supp[2])}
if(model%in%c(3,4)){x <- 2*rbeta(n.sam,shape1=2,shape2=4,ncp=0)-1}
d <- (x>x0)
y <- dtgen.new(x,d,model,type,ns[k])

# 2. Bandwidth selection

H <- rdd.bandwidth(x,y,d,x0,tt=c(0.2,0.8),m=9,method=c(1,2,3,4,5),seq(0.1,0.5,length.out=15),kr=krn)
h.med <- H$hcv
h.med <- min(max(h.med,0.1),0.5)

# 3. QTE and confidence bands

Qy <- rdd.rqpro(x,y,tt=c(0.2,0.8),m=9,x.eval=x0,bandw=h.med,method=1,kr=krn)$Q0		# quantile process

cband <- rdd.qte(x,y,d,x.eval=x0,alpha=0.9,tt=c(0.2,0.8),m=9,bandw=h.med,kr=krn,bias=0,eql=0)	# QTE and bands

# 4. Score test

test <- Score(x,y,x.eval=x0,alpha=c(0.9,0.95),tt=c(0.2,0.8),m=9,bandw=H$hcvi,kr=krn,left=0)
stat.s <- test$test
crit.s <- test$crit

# 5. Wald and Wald Robust tests

test1 <- Wald(x,y,d,x.eval=x0,alpha=c(0.9,0.95),tt=c(0.2,0.8),m=9,bandw=h.med,bandw2=h.med,kr=krn,test.type=c(1,2,3),sign.opt=1,eql=0)

test2 <- Wald(x,y,d,x.eval=x0,alpha=c(0.9,0.95),tt=c(0.2,0.8),m=9,bandw=h.med,bandw2=h.med,kr=krn,test.type=c(1,2,3),sign.opt=1,eql=1)

stat.w <- test1$wald.test
crit.w <- test1$wald.crit
stat.r <- test1$wald.robust.test
crit.r <- test1$wald.robust.crit
stat.e <- test2$wald.robust.test
crit.e <- test2$wald.robust.crit

# 6. Save outcomes and calculates rejection ratios

A[j,1] <- stat.s			# Score test
A[j,2:3] <- crit.s
A.w[j,1,] <- stat.w		# Wald test
A.w[j,2:3,] <- crit.w
A.w2[j,1,] <- stat.r		# Robust Wald tesrt
A.w2[j,2:3,] <- crit.r
A.w3[j,1,] <- stat.e		# Robust Wald with equality constraint
A.w3[j,2:3,] <- crit.e
B[j,] <- pmin(pmax(c(H$hcv,H$hint,H$hbdy,H$hik,H$hcvi),0.1),0.5)		# bandwidths

print("j = "); print(j)
}	# end of the inner-loop

# Save rejection ratios for Score test
	R[1,1] <- sum(A[,1] > A[,2])/L
	R[1,2] <- sum(A[,1] > A[,3])/L
# Save rejection ratios for Wald tests
for(ii in 1:3){
	R[(1+ii),1] <- sum(A.w[,1,ii] > A.w[,2,ii])/L
	R[(1+ii),2] <- sum(A.w[,1,ii] > A.w[,3,ii])/L
	R2[(1+ii),1] <- sum(A.w2[,1,ii] > A.w2[,2,ii])/L
	R2[(1+ii),2] <- sum(A.w2[,1,ii] > A.w2[,3,ii])/L
	R3[(1+ii),1] <- sum(A.w3[,1,ii] > A.w3[,2,ii])/L
	R3[(1+ii),2] <- sum(A.w3[,1,ii] > A.w3[,3,ii])/L	
}

# save outcomes
if(type==1){save(R,R2,R3,A,A.w,A.w2,A.w3,B,krn,L,type,file=paste("outcome","_type",type,"_model",model,"_n",n.sam,".rda",sep="")) }
if(type==2){save(R,R2,R3,A,A.w,A.w2,A.w3,B,krn,L,type,ns,file=paste("outcome","_type",type,"_model",model,"_n",n.sam,"_ns",ns[k],".rda",sep="")) }
}