# Simulating the data
set.seed(1234)
library(mvtnorm)
n=800
nbre1=400
nbre2=600
rZIGP<-function(n,rho=stop("no rho arg"),kesai=stop("no kesai arg"),lambda=stop("no 
         lambda arg"))
{
  #check if parameters are valid
  if (rho<0) { return ("rho has to be in [0,1]!")}
  if (rho>1) { return ("rho has to be in [0,1]!")}
  #inversion method
  x<- double(n)
  u<-runif(n,0,1)
  upper<-max(u)
  s<-double(1000)
  #P(X=0)
  p<-rho+(1-rho)*exp(-lambda)
  s[1]<-p
  if (upper>0){
    rekursive<-FALSE
    i<-1
    while (s[i]<upper){
      #P(X=i)
      if (rekursive==FALSE){
        p<-(1-rho)*lambda*(lambda+kesai*i)^(i-1)/exp(lgamma(i+1))*exp(-(lambda+kesai*i))}
      if(p=="NaN") {s[i]=2
      p=0}
      if (p==Inf){
        rekursive<-TRUE
        log.p.alt<-log((1-rho)*lambda*(lambda+kesai*(i-1))^((i-1)-1)/exp(lgamma(i-1+1))*exp(-(lambda+kesai*(i-1))))
      }
      if (rekursive==TRUE){
        log.p<-log((lambda+kesai*(i-1))*exp(-kesai)/i*(1+kesai/(lambda+kesai*(i-1)))^(i-1))+log.p.alt
        log.p.alt<-log.p
        p<-exp(log.p)
      }
      if (ceiling(i/1000)==floor(i/1000)){
        temp<-double(1000)
        s<-c(s,temp)
      }
      s[i+1]<-s[i]+p 
      i<-i+1
    }
  }
  for (j in 1:length(u)){
    i<-1
    while (u[j]>s[i]){
      i<-i+1
    }
    x[j]<-i-1
  }
  return(x)
}
nosim=500

nob	=nbre1-1
ini	= 200
nmax  	= nob + ini
y	= matrix(0,ncol=nmax,nrow=nosim)
x=      matrix(0,ncol=nmax,nrow=nosim)
xt1    = matrix(0,ncol=nob,nrow=nosim)
yt1	= matrix(0,ncol=nob,nrow=nosim)
om_s	= 0.6
alpha_s	= 0.9	
beta_s	= -0.2
rho_s       =0
kesai_s     =0.1
r_s         =0.5
lam_s1=matrix(0,ncol=nob,nrow=nosim)
true=c(om_s,alpha_s,beta_s,rho_s,kesai_s,r_s)
lam=NULL
lamr=NULL
lam[1]	= 5
lamr[1]     =(1-kesai_s)/(1-rho_s)*lam[1]
y[,1]	= 2.0
for (isi in 1:nosim)
{
  x[isi,]<-rnorm(nmax,0,1)
}

for(isi in 1:nosim)
{
  #print(isi)
  for (t in 2:nmax)
  { 
    lam[t] = exp(om_s+ alpha_s*log(y[isi,t-1]+1) + beta_s*log(lam[t-1])+x[isi,t]*r_s)
    lamr[t]= (1-kesai_s)/(1-rho_s)*lam[t]
    y[isi,t]	  = rZIGP(1,rho_s,kesai_s, lamr[t])
  }
  lam_s1[isi,]=lam[(ini+1):nmax]
}
for(isi in 1:nosim){
  yt1[isi,]=y[isi,-1:-ini]
}
for(isi in 1:nosim){
  xt1[isi,]=x[isi,-1:-ini]
}

nob	= nbre2-nbre1
nmax  	= nob + ini
y	= matrix(0,ncol=nmax,nrow=nosim)
x=      matrix(0,ncol=nmax,nrow=nosim)
xt2    = matrix(0,ncol=nob,nrow=nosim)
yt2	= matrix(0,ncol=nob,nrow=nosim)
om_s	= 0.2
alpha_s	= 0.9	
beta_s	= -0.2
rho_s       =0.6
kesai_s     =0.1
r_s         =0.1
lam_s2=matrix(0,ncol=nob,nrow=nosim)
true=c(om_s,alpha_s,beta_s,rho_s,kesai_s,r_s)
lam=NULL
lamr=NULL
lam[1]	= 2
lamr[1]     =(1-kesai_s)/(1-rho_s)*lam[1]
y[,1]	= 2.0
for (isi in 1:nosim)
{
  x[isi,]<-rnorm(nmax,0,1)
}

for(isi in 1:nosim)
{
  #print(isi)
  for (t in 2:nmax)
  { 
    lam[t] = exp(om_s+ alpha_s*log(y[isi,t-1]+1) + beta_s*log(lam[t-1])+x[isi,t]*r_s)
    lamr[t]= (1-kesai_s)/(1-rho_s)*lam[t]
    y[isi,t]	  = rZIGP(1,rho_s,kesai_s, lamr[t])
  }
  lam_s2[isi,]=lam[(ini+1):nmax]
}
for(isi in 1:nosim){
  yt2[isi,]=y[isi,-1:-ini]
}
for(isi in 1:nosim){
  xt2[isi,]=x[isi,-1:-ini]
}


nob	= n-nbre2+1
nmax  	= nob + ini
y	= matrix(0,ncol=nmax,nrow=nosim)
x=      matrix(0,ncol=nmax,nrow=nosim)
xt3    = matrix(0,ncol=nob,nrow=nosim)
yt3	= matrix(0,ncol=nob,nrow=nosim)
om_s	= 0.6
alpha_s	= 0.9	
beta_s	= -0.2
rho_s       =0
kesai_s     =0.1
r_s         =0.5
lam_s3=matrix(0,ncol=nob,nrow=nosim)
true=c(om_s,alpha_s,beta_s,rho_s,kesai_s,r_s)
lam=NULL
lamr=NULL
lam[1]	= 5
lamr[1]     =(1-kesai_s)/(1-rho_s)*lam[1]
y[,1]	= 2.0
for (isi in 1:nosim)
{
  x[isi,]<-rnorm(nmax,0,1)
}

for(isi in 1:nosim)
{
  #print(isi)
  for (t in 2:nmax)
  { 
    lam[t] = exp(om_s+ alpha_s*log(y[isi,t-1]+1) + beta_s*log(lam[t-1])+x[isi,t]*r_s)
    lamr[t]= (1-kesai_s)/(1-rho_s)*lam[t]
    y[isi,t]	  = rZIGP(1,rho_s,kesai_s, lamr[t])
  }
  lam_s3[isi,]=lam[(ini+1):nmax]
}
for(isi in 1:nosim){
  yt3[isi,]=y[isi,-1:-ini]
}
for(isi in 1:nosim){
  xt3[isi,]=x[isi,-1:-ini]
}

yt<-cbind(yt1,yt2,yt3)
xt<-cbind(xt1,xt2,xt3)
lam_s<-cbind(lam_s1,lam_s2,lam_s3)


saveRDS(yt, file="yt_RSm1.rds",version=2)
saveRDS(xt, file="xt_RSm1.rds",version=2)
saveRDS(lam_s, file="lam_RSm1.rds",version=2)

