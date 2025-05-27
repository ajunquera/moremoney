## Prepared in March, 2017
## Replicate the analysis using Card, Chetty, and Weber (2007) data
## Tested on a MAC desktop with 8GB ram and 3.2GHz processor; each result (say obtaining a uniform confidence band or a test and its critical value) takes about 2-12 hours to obtain using a single core.
## Different segments of the code can be run independently.

rm(list=ls())
require("foreign")
source("qte_rdd.R")

run.bandw <- 1		# if 1, estimate bandwidths, and if 0, do not estimate the bandwidth
run.qte   <- 0		# if 1, estimate QTE and confidence bands
run.test  <- 0		# if 1, do tests

## 1. Bandwidth estimation

if(run.bandw==1){

# 1.1. Severance pay (sp)

data <- read.dta("mar13_2015.dta")	# read data
drop <- is.na(data$noneduration) + is.na(data$tenure_cat_sp)		# strictly positive if missing occurs
data <- data[(drop==0),]	# keep complete observations
dur  <- data$noneduration	# unemployment duration (in days)
ten  <- data$tenure_cat_sp	# job tenure at the firm from which the worker got laid off (in months)

# sample restriction, includes only individuals who worked at least one month in the past five years at a firm different from the one from which they were just laid off
y  <- dur[(data$emp_prior==1)]	
x  <- ten[(data$emp_prior==1)]	
x0 <- 35.5		# cutoff
d  <- (x>x0)		# treatment assignment

h5.bdy <- rdd.bandwidth.app(x,y,d,x.eval=x0,tt=c(0.2,0.8),m=25,method=c(1,4,5),val=(3:10),band=c((max(x)-min(x)),rep(0.5*(max(x)-min(x)),2),(0.5*(max(x)-x0)),0.5*(x0-min(x))),kr=3)	# estimate boundary point bandwidth

h5.int <- rdd.bandwidth.app(x,y,d,x.eval=x0,tt=c(0.2,0.8),m=25,method=c(2,3),val=(3:10),band=(max(x)-min(x)),kr=3,br=1)	# estimate interior point bandwidths

save(h5.bdy,h5.int,file="Bandwidth_sp.rda")	# save the outcome for the severance pay

# 1.2. Extended benefits (eb)

drop <- is.na(data$tenure_cat_eb)	# strictly positive if missing occurs
data <- data[(drop==0),]
dur  <- data$noneduration		# unemployment duration (in days)
ten  <- data$tenure_cat_eb		# number of months worked at any firms in the past five years (in months)

y <- dur[(data$emp_prior==1)]	# sample restriction
x <- ten[(data$emp_prior==1)]
d <- (x>x0)

h5.bdy <- rdd.bandwidth.app(x,y,d,x.eval=x0,tt=c(0.2,0.8),m=25,method=c(1,4,5),val=(3:10),band=c((max(x)-min(x)),rep(0.5*(max(x)-min(x)),2),(0.5*(max(x)-x0)),0.5*(x0-min(x))),kr=3)

h5.int <- rdd.bandwidth.app(x,y,d,x.eval=x0,tt=c(0.2,0.8),m=25,method=c(2,3),val=(3:10),band=(max(x)-min(x)),kr=3,br=1)

save(h5.bdy,h5.int,file="Bandwidth_eb.rda")	# save outcomes for the extended benefits
}

## 2. QTE and Confidence bands

if(run.qte==1){

h.med <- 4.5		# bandwidth at the median. change it in order to try other bandwidths

# 2.1. Severance pay 

data <- read.dta("mar13_2015.dta")
drop <- is.na(data$noneduration)+is.na(data$tenure_cat_sp)
data <- data[(drop==0),]
dur  <- data$noneduration
ten  <- data$tenure_cat_sp

y  <- dur[(data$emp_prior==1)]
x  <- ten[(data$emp_prior==1)]
# further restriction in sample, includes only individuals who worked more than four months in the past five years at a firm different from the one from which they were just laid off
#y <- dur[(data$emp_prior4==1)]
#x <- ten[(data$emp_prior4==1)]
x0 <- 35.5
d  <- (x>x0)

qte.sp <- rdd.qte.app(x,y,d,x.eval=x0,alpha=0.9,tt=c(0.2,0.8),m=25,bandw=h.med,kr=3,bias=0,eql=0)	# QTE and its uniform bands without bias correction

qte.sp2 <- rdd.qte.app(x,y,d,x.eval=x0,alpha=0.9,tt=c(0.2,0.8),m=25,bandw=h.med,kr=3,bias=1,eql=0)	# QTE and its uniform bands with bias correction

qte.sp3 <- rdd.qte.app(x,y,d,x.eval=x0,alpha=0.9,tt=c(0.2,0.8),m=25,bandw=h.med,kr=3,bias=1,eql=1)	# QTE and its uniform bands with equality constrained bias correction

save(qte.sp,qte.sp2,qte.sp3,h.med,file=paste("QTE_sp","_h",h.med,".rda",sep=""))

# 2.2. Extended benefits

drop <- is.na(data$tenure_cat_eb)
data <- data[(drop==0),]
dur  <- data$noneduration
ten  <- data$tenure_cat_eb

y  <- dur[(data$emp_prior==1)]
x  <- ten[(data$emp_prior==1)]
#y <- dur[(data$emp_prior4==1)]
#x <- ten[(data$emp_prior4==1)]
d  <- (x>x0)

qte.eb <- rdd.qte.app(x,y,d,x.eval=x0,alpha=0.9,tt=c(0.2,0.8),m=25,bandw=h.med,kr=3,bias=0,eql=0)	# QTE and uniform bands without bias correction

qte.eb2 <- rdd.qte.app(x,y,d,x.eval=x0,alpha=0.9,tt=c(0.2,0.8),m=25,bandw=h.med,kr=3,bias=1,eql=0)	# QTE and uniform bands with bias correction

qte.eb3 <- rdd.qte.app(x,y,d,x.eval=x0,alpha=0.9,tt=c(0.2,0.8),m=25,bandw=h.med,kr=3,bias=1,eql=1)	# QTE and uniform bands with equality constrained bias correction

save(qte.eb,qte.eb2,qte.eb3,h.med,file=paste("QTE_eb","_h",h.med,".rda",sep=""))
}

## 3. Score and Wald tests

if(run.test==1){

h.med <- 4.5		# bandwidth at the median
alevel = c((1:99/100),(991:999/1000))

# 3.1. Severance pay

data <- read.dta("mar13_2015.dta")
drop <- is.na(data$noneduration)+is.na(data$tenure_cat_sp)
data <- data[(drop==0),]
dur  <- data$noneduration
ten  <- data$tenure_cat_sp

y  <- dur[(data$emp_prior==1)]
x  <- ten[(data$emp_prior==1)]
x0 <- 35.5
d  <- (x>x0)

score.S <- Score.app(x,y,x.eval=x0,alpha=alevel,tt=c(0.2,0.8),m=25,bandw=h.med,kr=3,left=1)		# Score test

test.S <- Wald.app(x,y,d,x.eval=x0,alpha=alevel,tt=c(0.2,0.8),m=25,h.med,h.med,kr=3,test.type=c(1,2,3),sign.opt=1,eql=0)		# Wald test and Robust Wald test with quantile-by-quantile bias estimation

test.S2 <- Wald.app(x,y,d,x.eval=x0,alpha=alevel,tt=c(0.2,0.8),m=25,h.med,h.med,kr=3,test.type=c(1,2,3),sign.opt=1,eql=1)		# Wald test and Robust Wald test with equality constrained bias estimation

save(score.S,test.S,test.S2,file=paste("test_sp","_h",h.med,".rda",sep=""))

# 3.2. Extended benefits

drop <- is.na(data$tenure_cat_eb)
data <- data[(drop==0),]
dur  <- data$noneduration
ten  <- data$tenure_cat_eb

y <- dur[(data$emp_prior==1)]
x <- ten[(data$emp_prior==1)]
d <- (x>x0)

score.E <- Score.app(x,y,x.eval=x0,alpha=alevel,tt=c(0.2,0.8),m=25,bandw=h.med,kr=3,left=1)

test.E <- Wald.app(x,y,d,x.eval=x0,alpha=alevel,tt=c(0.2,0.8),m=25,h.med,h.med,kr=3,test.type=c(1,2,3),sign.opt=1,eql=0)

test.E2 <- Wald.app(x,y,d,x.eval=x0,alpha=alevel,tt=c(0.2,0.8),m=25,h.med,h.med,kr=3,test.type=c(1,2,3),sign.opt=1,eql=1)

save(score.E,test.E,test.E2,file=paste("test_eb","_h",h.med,".rda",sep=""))
}