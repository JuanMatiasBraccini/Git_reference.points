library(mvtnorm)

#note: K and LINF must be in same order or magnitude

#growth pars
SCALER.mag=1000
SCALER.mag1=125
SCALER.mag2=10000
SCALER.mag3=5000

#whiskery
K.mean.W=0.369
Linf.mean.W=120.7/SCALER.mag1


#Dusky
K.mean.D=.0367
Linf.mean.D=374.4/SCALER.mag2


#Sandbar
K.mean.S=.04
Linf.mean.S=244.2/SCALER.mag3

#Gummy              #Walker 2010
K.mean.G=.123
Linf.mean.G=201.9/SCALER.mag


#create dummy age-length data
fn.get.sigma=function(Age,NN,LO,LINF,KK,Jit,SCALE,TITULO)
{
  LO=LO/SCALE
  age=sort(rep(0:Age,NN))
  mid.FL.fem=LO+(LINF-LO)*(1-exp(-KK*age))
  mid.FL.fem.scat=jitter(mid.FL.fem,Jit)
  plot(age,mid.FL.fem*SCALE,type='l',main=TITULO,ylab="length (cm)")
  points(age,mid.FL.fem.scat*SCALE)
  Dummy=data.frame(Age=age,FL=mid.FL.fem.scat)
  
  #fit traditional vonB growth model
  Start.pars=list(Linf=jitter(LINF,10),K=jitter(KK,10))
  fit=nls(FL~LO+(Linf-LO)*(1-exp(-K*Age)),data=Dummy,start=Start.pars)
  sigma=vcov(fit)
  return(sigma=sigma)
}

par(mfcol=c(2,2))
Sigma.G=fn.get.sigma(20,5,33.5,Linf.mean.G,K.mean.G,30,SCALER.mag,"gummy")

Sigma.W=fn.get.sigma(19,5,25,Linf.mean.W,K.mean.W,1000,SCALER.mag1,"whiskery") 
Sigma.D=fn.get.sigma(55,5,75.3,Linf.mean.D,K.mean.D,100,SCALER.mag2,"dusky")
Sigma.S=fn.get.sigma(29,5,42.5,Linf.mean.S,K.mean.S,15,SCALER.mag3,"sandbar")

# Sigma.G=fn.get.sigma(20,5,33.5,Linf.mean.G,K.mean.G,3,SCALER.mag,"gummy")
# Sigma.W=fn.get.sigma(19,5,25/SCALER.mag,Linf.mean.W,K.mean.W,200,"whiskery") 
# Sigma.D=fn.get.sigma(55,5,75.3/SCALER.mag,Linf.mean.D,K.mean.D,20,"dusky")
# Sigma.S=fn.get.sigma(29,5,42.5/SCALER.mag,Linf.mean.S,K.mean.S,2.5,"sandbar")


#Store all var-covar matrices in a list for using elsewhere
SIGMA=list(gummy=Sigma.G,whiskery=Sigma.W,dusky=Sigma.D,sandbar=Sigma.S)



#Plot
Plot.fn=function(k.mean,Linf.mean,SIG,SCALE)
{
  aa=rmvnorm(1000,mean=c(k.mean,Linf.mean),sigma=SIG)
  par(mfcol=c(2,1))
  hist(aa[,1],xlab="k",main=paste("CV=",round(sd(aa[,1])/mean(aa[,1]),3)))
  hist(aa[,2]*SCALE,xlab="Linf (m)",main=paste("CV=",round(sd(aa[,2])/mean(aa[,2]),3)))
  
}

Plot.fn(K.mean.G,Linf.mean.G,Sigma.G,SCALER.mag)
Plot.fn(K.mean.W,Linf.mean.W,Sigma.W,SCALER.mag1)
Plot.fn(K.mean.D,Linf.mean.D,Sigma.D,SCALER.mag2)
Plot.fn(K.mean.S,Linf.mean.S,Sigma.S,SCALER.mag3)




# #Plot
# Plot.fn=function(k.mean,Linf.mean,SIG)
# {
#   store=data.frame(k=NA,Lin=NA)
#   for (i in 1:1000){
#     growth.pars.sim=mvrnorm(1,mu=c(Linf.mean,k.mean),Sigma=SIG,tol=1e-10)
#     if(growth.pars.sim[1]<=0.01)    #repeat until sensible pars obtained
#     { repeat 
#     {
#       growth.pars.sim=mvrnorm(1,mu=c(k.mean,Linf.mean),Sigma=sigma,tol=1e-10)
#       if(growth.pars.sim[1]>0)break
#     }
#     }
#     store[i,]=growth.pars.sim
#   }
#   k.sim=growth.pars.sim[1]
#   Linf.sim=growth.pars.sim[2]
#   par(mfcol=c(2,1))
#   hist(store[,1],main="linf")
#   hist(store[,2],main="k")
# }
# 
# Plot.fn(K.mean.G,Linf.mean.G,Sigma.G)
# Plot.fn(K.mean.W,Linf.mean.W,Sigma.W)
# Plot.fn(K.mean.D,Linf.mean.D,Sigma.D)
# Plot.fn(K.mean.S,Linf.mean.S,Sigma.S)




