# Simulation code
#install.packages("mvtnorm")
#install.packages("foreach")
#install.packages("doParallel")
library(mvtnorm)
library(foreach)
library(doParallel)

## the  R function source required to compute the parameter estimation for full ALG model
mc.theta<-function(yt,xt,M,nob,nosim){
  result = NULL
  S1	 = matrix(0,nosim,6)
  S2	 = matrix(0,nosim,6)
  S3	 = matrix(0,nosim,6)
  S4	 = matrix(0,nosim,6)
  S5	 = matrix(0,nosim,6)
  S6	 = matrix(0,nosim,6)
  S7	 = matrix(0,nosim,6)
  S8     =matrix(0,nosim,6)
  r1	 = NULL
  r2	 = NULL
  r3     	= NULL
  r4     =NULL
  r5    =NULL
  r6    =NULL
  r7    =NULL
  r8    =NULL
  
  for (isi in 1:nosim)
  {
    # Setting up starting values
    om	=	0.5
    alpha	=	0.1
    beta	=	0.1
    rho   =   0.6
    kesai =   0.1
    r     =   0.3
    lam   =   yt[1]
    lam_o	=	NULL
    lamr_o=      NULL
    lam_n	=	NULL
    lamr_n=     NULL
    star=NULL
    star_lam=NULL
    lam_o[1]	=	 2
    lamr_o[1]=(1-kesai)/(1-rho)*lam_n[1]
    lam_n[1] =	2
    lamr_n[1]=(1-kesai)/(1-rho)*lam_n[1]
    count	=	0
    #****************************************************************New Part
    
    # *******************Random walk Metropolis
    
    ### set stepsize 
    
    step_om	=	0.3
    step_al	= 	0.3	
    step_beta  =	0.3
    step_rho   =      0.1
    step_kesai =      0.1
    step_r     =      0.3
    step_lam   =      1
    ind	=   	5000   ## burin-in
    draws = matrix(0,M,7)
    lam_t <-NULL
    lam_e<-NULL
    #result  = matrix(0,M-ind,3)
    lam0=mean(yt)
    for (i in 1:M){
      lik1=0
      lik2=0
      lik_lam1=0
      lik_lam2=0
      lognor1=0
      lognor2=0
      old=c(om,alpha,beta,rho,kesai,r,lam)
      lam_o[1]=lam
      
      for (t in 2:nob)
      { 
        lam_o[t] 	= exp(om  + alpha*log(yt[t-1]+1) + beta*log(lam_o[t-1])+r*xt[t])
        lamr_o[t]   = (1-kesai)/(1-rho)*lam_o[t]
        if (yt[t]==0) {lik1= lik1+log(rho+(1-rho)*exp(-lamr_o[t]))}
        else         {lik1= lik1+log(1-rho)+log(lamr_o[t])+(yt[t]-1)*log(lamr_o[t]+kesai*yt[t])-(lamr_o[t]+kesai*yt[t])-sum(log(1:yt[t]))}
        
      }
      
      if (i <= ind)	
      {
        repeat 
        {	 star[1]	= om +  rnorm(1,0,0.2)*step_om
        {break}				
        }
        repeat 
        {
          star[2]	=  alpha +  rnorm(1,0,0.2)*step_al
          star[3]        =  beta  +   rnorm(1,0,0.2)*step_beta
          if ( (star[2]>0 )&(abs(star[3]) <1) & (abs(star[2]+ star[3])<1) ) { break}
          
        }
        repeat 
        {	 star[4]	= rho +  rnorm(1,0,0.2)*step_rho
        if ( star[4] < 1 & star[4]>0 ){break}	#om_sta = om_n ;				
        }
        repeat 
        {	 star[5]	= kesai +  rnorm(1,0,0.2)*step_kesai
        if ( star[5] < 1 & star[5]>0 ){break}	#om_sta = om_n ;				
        }
        repeat 
        {	 star[6]	= r +  rnorm(1,0,0.2)*step_r
        if ( star[6] < 1 & star[6]>0 ){break}	#om_sta = om_n ;				
        }
        repeat 
        {	 star[7]	= lam0 +  rnorm(1,0,0.2)*step_lam
        if ( star[7]>0 ){break}	#om_sta = om_n ;				
        }
      }
      else
      {
        repeat
        {
          star=rmvnorm(1, IK_mean, IK_cov)
          if ((star[2] >0 )&( abs(star[3]) <1) & (abs(star[2]+ star[3])<1)&star[4]>0&star[4]<1&star[5]<1&star[5]>0&star[6]<1&star[6]>0&star[7]>0)
            
          { break}
        }
        lognor2=dmvnorm(star,IK_mean,IK_cov,log=T)
        lognor1=dmvnorm(old,IK_mean,IK_cov,log=T)
        
      }
      lam_n[1]=star[7]
      for (t in 2:nob)
      { 
        lam_n[t] 	= exp(star[1]  + star[2]*log(yt[t-1]+1) + star[3]*log(lam_n[t-1])+star[6]*xt[t])
        lamr_n[t]   =(1-star[5])/(1-star[4])*lam_n[t]
        if (yt[t]==0) {lik2= lik2+log(star[4]+(1-star[4])*exp(-lamr_n[t]))}
        else         {lik2= lik2+log(1-star[4])+log(lamr_n[t])+(yt[t]-1)*log(lamr_n[t]+star[5]*yt[t])-(lamr_n[t]+star[5]*yt[t])-sum(log(1:yt[t]))}
        
        
      }
      
      lik_lam1=dnorm(lam,lam0,0.2,log=T)
      lik_lam2=dnorm(star[7],lam0,0.2,log=T)
      lik = lognor1-lognor2 + lik2- lik1+lik_lam2-lik_lam1 ## the first two terms are from the proposals
      u 	=	 runif(1)
      if (log(u)<lik) 
      {
        om      	=  star[1]
        alpha 	=  star[2]
        beta    	=  star[3]
        rho         =  star[4]
        kesai       =  star[5]
        r           =  star[6]
        lam         =  star[7]
        count	= count+1
        
      }
      draws[i,] = c(om ,alpha,beta,rho,kesai,r,lam)
      lam_t[1]=draws[i,7]
      for (t in 2:nob)
      { 
        lam_t[t] 	= exp(draws[i,1]  + draws[i,2]*log(yt[t-1]+1) + draws[i,3]*log(lam_t[t-1])+draws[i,6]*xt[t])
      }
      lam_e[i]=lam_t[nob]
      
      if(i == ind ) {IK_mean=c(mean(draws[501:ind,1]),mean(draws[501:ind,2]),mean(draws[501:ind,3]),mean(draws[501:ind,4]),mean(draws[501:ind,5]),mean(draws[501:ind,6]),mean(draws[501:ind,7]))
      IK_cov=cov(draws[501:ind,])
      
      }
      
      
    }
    ############################# PLOT
    
    names		 = c("om","alpha","beta","rho","kesai","r","lam","lam_nob")
    
    
    ######################################################################
    #############################
    
    MCMC=(ind+1):M
    S1[isi,1]=0
    S1[isi,2]=mean(draws[MCMC,1])
    S1[isi,3]=median(draws[MCMC,1])
    S1[isi,4]=sd(draws[MCMC,1])
    S1[isi,5]=quantile(draws[MCMC,1],0.025)
    S1[isi,6]=quantile(draws[MCMC,1],0.975)
    
    S2[isi,1]=0
    S2[isi,2]=mean(draws[MCMC,2])
    S2[isi,3]=median(draws[MCMC,2])
    S2[isi,4]=sd(draws[MCMC,2])
    S2[isi,5]=quantile(draws[MCMC,2],0.025)
    S2[isi,6]=quantile(draws[MCMC,2],0.975)
    
    S3[isi,1]=0
    S3[isi,2]=mean(draws[MCMC,3])
    S3[isi,3]=median(draws[MCMC,3])
    S3[isi,4]=sd(draws[MCMC,3])
    S3[isi,5]=quantile(draws[MCMC,3],0.025)
    S3[isi,6]=quantile(draws[MCMC,3],0.975)
    
    S4[isi,1]=0
    S4[isi,2]=mean(draws[MCMC,4])
    S4[isi,3]=median(draws[MCMC,4])
    S4[isi,4]=sd(draws[MCMC,4])
    S4[isi,5]=quantile(draws[MCMC,4],0.025)
    S4[isi,6]=quantile(draws[MCMC,4],0.975)
    
    S5[isi,1]=0
    S5[isi,2]=mean(draws[MCMC,5])
    S5[isi,3]=median(draws[MCMC,5])
    S5[isi,4]=sd(draws[MCMC,5])
    S5[isi,5]=quantile(draws[MCMC,5],0.025)
    S5[isi,6]=quantile(draws[MCMC,5],0.975)
    
    S6[isi,1]=0
    S6[isi,2]=mean(draws[MCMC,6])
    S6[isi,3]=median(draws[MCMC,6])
    S6[isi,4]=sd(draws[MCMC,6])
    S6[isi,5]=quantile(draws[MCMC,6],0.025)
    S6[isi,6]=quantile(draws[MCMC,6],0.975)
    
    S7[isi,1]=0
    S7[isi,2]=mean(draws[MCMC,7])
    S7[isi,3]=median(draws[MCMC,7])
    S7[isi,4]=sd(draws[MCMC,7])
    S7[isi,5]=quantile(draws[MCMC,7],0.025)
    S7[isi,6]=quantile(draws[MCMC,7],0.975)
    
    S8[isi,1]=0
    S8[isi,2]=mean(lam_e[MCMC])
    S8[isi,3]=median(lam_e[MCMC])
    S8[isi,4]=sd(lam_e[MCMC])
    S8[isi,5]=quantile(lam_e[MCMC],0.025)
    S8[isi,6]=quantile(lam_e[MCMC],0.975)
  }
  #r0=true
  r1=apply(S1,2,mean)
  r2=apply(S2,2,mean)
  r3=apply(S3,2,mean)
  r4=apply(S4,2,mean)
  r5=apply(S5,2,mean)
  r6=apply(S6,2,mean)
  r7=apply(S7,2,mean)
  r8=apply(S8,2,mean)
  result=round(rbind(r1,r2,r3,r4,r5,r6,r7,r8),4)
  colnames(result) <- c("true value","mean","median","std","P025","P975")
  rownames(result) <- names
  return(result)
}

no_cores<-detectCores()
cl<-makeCluster(no_cores)
registerDoParallel(cl)
cv=c(10,2.8625,1.7789,1.5683,1.3370,0.6539)

## Data input
yt<-readRDS("yt_RSm1.rds")
xt<-readRDS("xt_RSm1.rds")
lam_s<-readRDS("lam_RSm1.rds")
nosim=nrow(yt)
noleng=ncol(yt)
nstart=300 # initial time  period
nob=nstart+1 #  start forecast point

lam_t<-NULL
lam_r<-NULL
lamr_t<-NULL
lam_f<-matrix(0,ncol=noleng-nob,nrow=nosim)
lam_e<-matrix(0,ncol=noleng-nob+1,nrow=nosim)
lam025_e<-matrix(0,ncol=noleng-nob+1,nrow=nosim)
lam975_e<-matrix(0,ncol=noleng-nob+1,nrow=nosim)
choice<-matrix(0,ncol=noleng-nob+1,nrow=nosim)
est<-array(0,dim=c(8,noleng-nob+1,nosim))
est025<-array(0,dim=c(8,noleng-nob+1,nosim))
est975<-array(0,dim=c(8,noleng-nob+1,nosim))

r<-1/2
s1=30
s2=50
s3=70
s4=90
s5=120
s6=140
s=c(s1,s2,s3,s4,s5,s6)
M=10000 ## MCMC iteration

# main code
save<-foreach(isi=1:nosim,.packages='mvtnorm') %dopar%
  {
    nob=nstart+1
    result=list(matrix(0,ncol=8,nrow=(noleng-nob+1)),matrix(0,ncol=8,nrow=(noleng-nob+1)),matrix(0,ncol=8,nrow=(noleng-nob+1)),0,0)
    for(k in 1:(noleng-nob+1))
    {
      i=1
      while(i<=5)
      {
        temp=mc.theta(t(as.matrix(yt[isi,(nob-s[i]+1):nob])),t(as.matrix(xt[isi,(nob-s[i]+1):nob])),M,s[i],1)
        theta_o1=temp[,2]
        theta025_o1=temp[,5]
        theta975_o1=temp[,6]  
        #############calculate L(I,thetatuta)  
        
        temp=mc.theta(t(as.matrix(yt[isi,(nob-s[i+1]+1):nob])),t(as.matrix(xt[isi,(nob-s[i+1]+1):nob])),M,s[i+1],1)
        theta_o2=temp[,2]
        theta025_o2=temp[,5]
        theta975_o2=temp[,6]
        
        lik2=0 
        for( t in 1:s[i+1]){
          lam_t[t] 	= ifelse(t==1,theta_o2[7],exp(theta_o2[1]  + theta_o2[2]*log(yt[isi,nob-s[i+1]+t-1]+1) + theta_o2[3]*log(lam_t[t-1])+theta_o2[6]*xt[isi,nob-s[i+1]+t]))
          lam_r[t]   =(1-theta_o2[5])/(1-theta_o2[4])*lam_t[t]
          if (yt[isi,nob-s[i+1]+t]==0) {lik2= lik2+log(theta_o2[4]+(1-theta_o2[4])*exp(-lam_r[t]))}
          else         {lik2= lik2+log(1-theta_o2[4])+log(lam_r[t])+(yt[isi,nob-s[i+1]+t]-1)*log(lam_r[t]+theta_o2[5]*yt[isi,nob-s[i+1]+t])-(lam_r[t]+theta_o2[5]*yt[isi,nob-s[i+1]+t])-sum(log(1:yt[isi,nob-s[i+1]+t]))}
        }
        
        ##################calculate  Lthetahat
        theta_h1=theta_o1
        theta_h2=theta_o2
        lik_h1=0
        for( t in 1:s[i+1]){
          lam_t[t]=ifelse(t==1,theta_o2[7],exp(theta_h1[1]+ theta_h1[2]*log(yt[isi,nob-s[i+1]+t-1]+1) + theta_h1[3]*log(lam_t[t-1])+theta_h1[6]*xt[isi,nob-s[i+1]+t]))
          lam_r[t]   =(1-theta_h1[5])/(1-theta_h1[4])*lam_t[t]
          if (yt[isi,nob-s[i+1]+t]==0) {lik_h1= lik_h1+log(theta_h1[4]+(1-theta_h1[4])*exp(-lam_r[t]))}
          else         {lik_h1= lik_h1+log(1-theta_h1[4])+log(lam_r[t])+(yt[isi,nob-s[i+1]+t]-1)*log(lam_r[t]+theta_h1[5]*yt[isi,nob-s[i+1]+t])-(lam_r[t]+theta_h1[5]*yt[isi,nob-s[i+1]+t])-sum(log(1:yt[isi,nob-s[i+1]+t]))}
          
        }
        
        ###############calculete diff
        diff<-NULL
        diff[i]=(abs(lik2-lik_h1))^r
        if(diff[i]<cv[i+1])  {i=i+1}
        else {break}
      }
      lag=i
      result[[1]][k,]=theta_o1
      result[[2]][k,]=theta025_o1
      result[[3]][k,]=theta975_o1
      result[[4]][k]=lag
      result[[5]][k]=lam_s[isi,nob]
      nob=nob+1
    }
    result
  }
saveRDS(save, "result_RSm1.rds")
stopCluster(cl)
stopImplicitCluster()

