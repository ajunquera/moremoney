#' Regression Discontinuity (RD) Design for Categorical Outcomes
#'
#' This function gives the optimal bandwidth, point estimates and robust confidence 
#' intervals for the RD estimand, when the outcome variable is categorical. The methods use
#' a nonparametric multinomial logit (MNL) transformation.
#' @param DAT The data, n by (J+2). n is the sample size and J+1 is the number of categories for the 
#' outcome. The first column of DAT is the running variable, the rest J+1 columns are binary dummies 
#' for each of the J+1 categories (each row sums to one). [If you are familiar with MATLAB, 
#' the last J+1 columns follows the format in the command mnrfit] The matrix DAT should have 
#' been sorted (ascendingly) according to the first column.
#' @param c The cutoff
#' @param H0_t J by 1. The null hypothesis when you perform a t-test.
#' @param H0_R q by J, where q is number of restrictions. H0_R and H0_r are restriction
#' constants as in the null hypothesis when you perform a Wald-test. E.g. H0_R=I_J 
#' if we want a joint significance test, H0_R=[1,-1] if we want a test of equivalence for the
#' treatment effects across two categories.
#' @param H0_r q by 1. See above.
#' @param level Confidence interval level. e.g. level=0.95.
#' @return The result is a list, which contains following results.
#' 
#' h_opt: Optimal bandwidth (a scalar). 
#' 
#' ATE: Point estimates (J by 1). 
#' 
#' ttest: (Standard/Non-robust) t-test statistics and p-values (J by 2). 
#' 
#' ttest_rob: Robust t-test statistics and p-values. 
#' 
#' 
#' ci: (Standard/Non-robust) Confidence interval. 
#' 
#' ci_rob: Robust confidence interval. 
#' 
#' Waldtest: (Standard/Non-robust) Wald test statistic and the p-value (1 by 2). 
#' 
#' Waldtest_rob: Robust Wald statistic.
#' @keywords RD
#' @export
#' @details The details can be found in Xu (2017+). Introductions to regression discontinuity
#' designs can be found in Imbens and Lemieux (2008) and Cattaneo and Escanciano (2017).
#' 
#' [1] Cattaneo, M.D., and J.C. Escanciano (2017): "Introduction: Regression Discontinuity Designs,"
#' in Regression Discontinuity Designs: Theory and Applications, \emph{Advances in Econometrics}, volume 38.
#' 
#' [2] Imbens, G.W., and T. Lemieux (2008): "Regression Discontinuity Designs: a Guide to 
#' Practice," \emph{Journal of Econometrics}, 142, 615-635.
#' 
#' [3] Xu, K.-L. (2017+): "Regression Discontinuity with Categorical Outcomes," \emph{Journal of Econometrics}, in press.
#' @examples
#' #Read your data into R, naming it as DAT (see above for the format of DAT)
#' #To illustrate, in the following I use the data in birth.weight.data (40750 obs)
#' 
#' >library("rd.categorical")
#' >data("birth.weight.data")
#' >DAT=birth.weight.data
#' >J=length(DAT[1,])-2
#' >c=1.5             #your cutoff value
#' >g=rd.mnl(DAT, c, matrix(0,J,1),diag(J),matrix(0,J,1),0.95)
#' 
#' >g$h_opt           #Chosen bandwidth
#' [1] 0.1788468
#' 
#' >g$ATE             #RD estimates
#'            [,1]
#'[1,] -0.01844537
#'[2,]  0.02043593
#'
#' >g$ttest_rob       #95% Robust t-tests
#'          [,1]         [,2]
#'[1,] -4.279627 1.872071e-05
#'[2,]  3.879861 1.045162e-04
#'
#' >g$Waldtest_rob    #95% Robust Wald joint-significance test
#'          [,1]        [,2]
#'[1,]  32.87078 7.28116e-08
#'
#'@author Ke-Li Xu, Indiana University, \email{keli.xu15@gmail.com}

#Author's note: from rd_mnl_4.m (in MATLAB, 2/14/2017)
#     **Revised August 24, 2017. (Now works for all J. A typo is also corrected)**

rd.mnl<-function(DAT,c,H0_t,H0_R,H0_r,level)
{
  n=length(DAT[,1])
  J=length(DAT[1,])-2
  
  dat=DAT[,1:2] #dat is tailored for the command "multinom"
  
  for (j in 1:J) {
     for (i in 1:n) {if (DAT[i,j+2]==1) dat[i,2]=j+1}
  }
  
  for (i in 1:n) {if (dat[i,1]>=c) break}
  
  L=dat[1:(i-1),]  # the left sample
  R=dat[i:n,]      # the right sample
  
  n_f=length(L[,1])
  n_r=length(R[,1])
  
  g=matrix(0,J,1)
  g2=matrix(0,J,1)
  g3=matrix(0,J,1)
  mu=matrix(0,J,1)
  Gamma=matrix(0,J,J)
  
  library(nnet)
  ss=relevel(factor(L[,2]),ref=J+1)
  #ss=relevel(factor(L[,2]),ref=3)
  res=multinom(ss ~ L[,1]+I(L[,1]^2)+I(L[,1]^3))
  res=summary(res)$coefficients
  res=matrix(res,nrow=J,ncol=4)
  
  sum=1
  
  sum
  for (j in 1:J) {
    g[j]=res[j,1]+res[j,2]*c+res[j,3]*c^2+res[j,4]*c^3
    g2[j]=3*2*res[j,4]*c+2*res[j,3] # 2nd derivative
    g3[j]=3*2*res[j,4] # 3rd derivative
    sum=sum+exp(g[j])
  }
  
  sum
  
  for (j in 1:J) {
    mu[j]=exp(g[j])/sum
  }
  
  for (j in 1:J){
    for (k in 1:J){
      if (j==k) Gamma[j,k]=mu[j]*(1-mu[j])
      else Gamma[j,k]=-mu[j]*mu[k]
    }
    
  }
  
  den=density(DAT[,1],from=c,to=c)
  den_c=den$y[1]
  
  K0=2.702^5 # the kernel constant (flat kernel is used) (for bandwidth to estimate mu)
  K2=3.557^7 # the kernel constant (flat kernel is used) (for bandwidth to estimate g2, the 2nd derivative of g)
  h_mu_f=K0*sum(diag(Gamma))/(n*den_c*(norm(Gamma%*%g2,type="F"))^2)
  h_mu_f=h_mu_f^(1/5)
  h_g2_f=K2*sum(diag(solve(Gamma)))/(n*den_c*(norm(g3,type="F"))^2)
  h_g2_f=h_g2_f^(1/7)
  
  ss=relevel(factor(R[,2]),ref=J+1)
  res=multinom(ss ~ R[,1]+I(R[,1]^2)+I(R[,1]^3))
  res=summary(res)$coefficients
  res=matrix(res,nrow=J,ncol=4)
  
  sum=1
  for (j in 1:J) {
    g[j]=res[j,1]+res[j,2]*c+res[j,3]*c^2+res[j,4]*c^3
    g2[j]=3*2*res[j,4]*c+2*res[j,3] # 2nd derivative
    g3[j]=3*2*res[j,4] # 3rd derivative
    sum=sum+exp(g[j])
  }
  
  for (j in 1:J) mu[j]=exp(g[j])/sum
  
  
  Gamma=matrix(0,J,J)
  for (j in 1:J) {
    for (k in 1:J) {
      if (j==k) Gamma[j,k]=mu[j]*(1-mu[j])
      else Gamma[j,k]=-mu[j]*mu[k]
    }
  }
  
  h_mu_r=K0*sum(diag(Gamma))/(n*den_c*(norm(Gamma%*%g2,type="F"))^2)
  h_mu_r=h_mu_r^(1/5)
  h_g2_r=K2*sum(diag(Gamma^(-1)))/(n*den_c*(norm(g3,type="F"))^2)
  h_g2_r=h_g2_r^(1/7)
  
  for (i in 1:n_f) {if (L[i,1]>=c-h_mu_f) break}
  L1=L[i:n_f,] # the left sample used to estimate mu
  for (i in 1:n_f) {if (L[i,1]>=c-h_g2_f) break}
  L2=L[i:n_f,] # the left sample used to estimate the 2nd derivative of g
  
  for (i in 1:n_r) {if (R[i,1]>=c+h_mu_r) break}
  R1=R[1:i,] # the right sample used to estimate mu
  for (i in 1:n_r) {if (R[i,1]>=c+h_g2_r) break}
  R2=R[1:i,] # the right sample used to estimate the 2nd derivative of g
  
  # now we use these 4 subsamples to estimate Gamma_f (using L1), Gamma_r (using R1), g2_f (using L2) and g2_r (using R2).
  
  ss=relevel(factor(L1[,2]),ref=J+1)
  res=multinom(ss ~ L1[,1])
  res=summary(res)$coefficients
  res=matrix(res,nrow=J,ncol=2)
  
  sum=1
  for (j in (1:J)) {
    g[j]=res[j,1]+res[j,2]*c
    sum=sum+exp(g[j])
  }
  
  for (j in 1:J) mu[j]=exp(g[j])/sum
  
  Gamma_f=matrix(0,J,J)
  for (j in 1:J) {
    for (k in 1:J) {
      if (j==k) Gamma_f[j,k]=mu[j]*(1-mu[j])
      else Gamma_f[j,k]=-mu[j]*mu[k]
    }
  }
  
  ss=relevel(factor(R1[,2]),ref=J+1)
  res=multinom(ss ~ R1[,1])
  res=summary(res)$coefficients
  res=matrix(res,nrow=J,ncol=2)
  
  sum=1
  for (j in 1:J)
  {
    g[j]=res[j,1]+res[j,2]*c;
    sum=sum+exp(g[j]);
    mu[j]=exp(g[j]);
  }
  
  for (j in 1:J) {mu[j]=exp(g[j])/sum}
  
  Gamma_r=matrix(0,J,J)
  for (j in 1:J) {
    for (k in 1:J) {
      if (j==k) Gamma_r[j,k]=mu[j]*(1-mu[j])
      else Gamma_r[j,k]=-mu[j]*mu[k]
    }
  }
  
  g2_f=matrix(0,J,1)
  g2_r=matrix(0,J,1)
  ss=relevel(factor(L2[,2]),ref=J+1)
  res=multinom(ss ~ L2[,1]+I(L2[,1]^2))
  res=summary(res)$coefficients
  res=matrix(res,nrow=J,ncol=3)
  
  for (j in 1:J) {g2_f[j]=2*res[j,3]}
  
  ss=relevel(factor(R2[,2]),ref=J+1)
  res=multinom(ss ~ R2[,1]+I(R2[,1]^2))
  res=summary(res)$coefficients
  res=matrix(res,nrow=J,ncol=3)
  
  for (j in 1:J) {g2_r[j]=2*res[j,3]}
  
  # now to calculate one final product: optimal bandwidth
  
  h_opt=2.702^5*sum(diag(Gamma_f+Gamma_r))/(n*den_c*(norm(Gamma_r%*%g2_r-Gamma_f%*%g2_f,type="F"))^2)
  h_opt=h_opt^(0.2)
  
  #pick the subsample within the h-neighborhood
  for (i in 1:n_f) {
    if (L[i,1]>=c-h_opt) break;
  }
  
  L3=L[i:n_f,]  # the left sample used to estimate mu
  
  for (i in 1:n_r) {
    if (R[i,1]>=c+h_opt) break
  }
  
  R3=R[1:i,] # the right sample used to estimate mu
  
  library(nnet)
  ss=relevel(factor(L3[,2]),ref=J+1)
  te_f=multinom(ss ~ L3[,1])
  res_f=summary(te_f)$coefficients
  res_f=matrix(res_f,nrow=J,ncol=2)
  
  ss=relevel(factor(R3[,2]),ref=J+1)
  te_r=multinom(ss ~ R3[,1])
  res_r=summary(te_r)$coefficients
  res_r=matrix(res_r,nrow=J,ncol=2)
  #te=multinom(ss ~ L[,1]+I(L[,1]^2)+I(L[,1]^3)) #results are somewhat different from MATLAB
  
  g_f=numeric()
  g_r=numeric()
  mu_f=numeric()
  mu_r=numeric()
  ATE=matrix(0,J,1)
  sum_f=1;
  sum_r=1;
  for (j in 1:J) {
    g_f[j]=res_f[j,1]+res_f[j,2]*c
    g_r[j]=res_r[j,1]+res_r[j,2]*c
    sum_f=sum_f+exp(g_f[j]);
    sum_r=sum_r+exp(g_r[j]);
  }
  
  for (j in 1:J) {
    mu_f[j]=exp(g_f[j])/sum_f;
    mu_r[j]=exp(g_r[j])/sum_r;
    ATE[j]=mu_r[j]-mu_f[j];
  }
  
  # Now do inference
  
  q=length(H0_R[,1]) #number of restrictions
  
  c0=-1/6
  c1=4
  
  Gamma_f=matrix(0,J,J)
  for (j in 1:J) {
    for (k in 1:J) {
      if (j==k) Gamma_f[j,k]=mu_f[j]*(1-mu_f[j])
      else Gamma_f[j,k]=-mu_f[j]*mu_f[k]
    }
  }
  
  Gamma_r=matrix(0,J,J)
  for (j in 1:J) {
    for (k in 1:J) {
      if (j==k) Gamma_r[j,k]=mu_r[j]*(1-mu_r[j])
      else Gamma_r[j,k]=-mu_r[j]*mu_r[k] 
    }
  }
  
  B=c0*h_opt^2*(Gamma_r%*%g2_r-Gamma_f%*%g2_f)/2
  V=c1*(Gamma_r+Gamma_f)/den_c
  
  ttest=matrix(0,J,2)
  ttest_rob=matrix(0,J,2)
  ci=matrix(0,J,2)
  ci_rob=matrix(0,J,2)
  Waldtest=matrix(0,1,2)
  Waldtest_rob=matrix(0,1,2)
  
  Waldtest[1]=n*h_opt*t(H0_R%*%(ATE-B)-H0_r)%*%solve(H0_R%*%V%*%t(H0_R))%*%(H0_R%*%(ATE-B)-H0_r)
  Waldtest[2]=1-pchisq(Waldtest[1],q) #p-value

  lv=qnorm(1-(1-level)/2) #lv=1.645, if level=0.9
  
  for (j in 1:J) {
    ttest[j,1]=sqrt(n*h_opt)*(ATE[j]-H0_t[j]-B[j])/sqrt(V[j,j])
    ttest[j,2]=2*pnorm(-abs(ttest[j,1])) #p-value
    ci[j,1]=ATE[j]-B[j]-lv*sqrt(V[j,j])/sqrt(n*h_opt)
    ci[j,2]=ATE[j]-B[j]+lv*sqrt(V[j,j])/sqrt(n*h_opt)
  }

 #robust tests (higher order correction terms in the variance estimates)

  VR=(Gamma_r*(c1+20*(h_opt/h_g2_r)^5+10*(h_opt/h_g2_r)^3)+Gamma_f*(c1+20*(h_opt/h_g2_f)^5+10*(h_opt/h_g2_f)^3))/den_c
  
  Waldtest_rob[1]=n*h_opt*t(H0_R%*%(ATE-B)-H0_r)%*%solve(H0_R%*%VR%*%t(H0_R))%*%(H0_R%*%(ATE-B)-H0_r)
  Waldtest_rob[2]=1-pchisq(Waldtest_rob[1],q) #p-value

  for (j in 1:J) {
    ttest_rob[j,1]=sqrt(n*h_opt)*(ATE[j]-H0_t[j]-B[j])/sqrt(VR[j,j])
    ttest_rob[j,2]=2*pnorm(-abs(ttest_rob[j,1])) #p-value
    ci_rob[j,1]=ATE[j]-B[j]-lv*sqrt(VR[j,j])/sqrt(n*h_opt)
    ci_rob[j,2]=ATE[j]-B[j]+lv*sqrt(VR[j,j])/sqrt(n*h_opt)
  }
  
  #Reporting all results
  list(h_opt=h_opt,ATE=ATE,ttest=ttest,ttest_rob=ttest_rob,ci=ci,ci_rob=ci_rob,Waldtest=Waldtest,Waldtest_rob=Waldtest_rob)
  }

#Old notes for people who did not install the rd.categorical package

#For empirical researchers

#Step 1: Input the data, naming it as "DAT". (e.g. DAT has dimension 40750 by 4, as in the CA infant birth weight example)
#Step 2: Define the cutoff c (e.g. c=1.5, as in the infant birth weight example)
#Step 3: Run all the lines above (i.e. define a function "rd_mnl_final" in R).
#Step 4: Run the following:

#J=length(DAT[1,])-2
#rd.mnl(DAT,c,matrix(0,J,1),diag(J),matrix(0,J,1),0.90) #90% t-tests and joint significance Wald test

# Run the following, if a test for equivalent ATE across categories is desired, only for J>=2:
#EQ=matrix(0,J-1,J)
#for (i in 1:(J-1)) {EQ[i,i]=1; EQ[i,i+1]=-1}

#rd.mnl(DAT,c,matrix(0,J,1),EQ,0,0.90) #90% t-tests and equivalence Wald test

#Step 5: Read the results
