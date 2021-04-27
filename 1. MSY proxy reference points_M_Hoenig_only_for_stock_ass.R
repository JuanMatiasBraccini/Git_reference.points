#notes:re run H calculation for gummy and whiskery sharks for stock assessments



#VERY IMPORTANT: check if simulations make sense (growth in particular) by tracking
                # selectivity and growth curve. Resample if too much variability

rm(list=ls(all=TRUE))
library(triangle)
library(TeachingDemos)      #for grey scale
library(mvtnorm)      #for multivariate normal pdf
library(PBSmapping)
handl_OneDrive=function(x)paste('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias',x,sep='/')

#source sigmas for var-covar matrix
source(handl_OneDrive("Analyses/Reference Points/5.Derive var_covar Matrix.r"))

#source indirect estimation of M
#M.gummy=0.186  #Table 20 Simpfendorfer et al 1996
M.gummy=0.283      #Walker et al 2000
M.whiskery=0.27  #Simpfendorfer et al 2000



#---DATA SECTION----
species.names <- c("gummy","whiskery")
N.sp=length(species.names)


#---PARAMETERS SECTION----

#Management decisions
Prop.FMsy=0.75  #management parameter. Proportion of Fmsy (as used by NMFS)


#Life history par vectors (Gummy,Whiskery)
  #Min values
Min.max.age=c(16,15)
Min.breed.freq=c(1,0.5)
Min.age.mat=c(4,6)
Min.fec=c(1,4)

  #Max values
Max.max.age=c(21,20) #put +1 because using floor in random sample
Max.breed.freq=c(1,0.5)
Max.age.mat=c(6,8) #put +1 because using floor in random sample
Max.fec=c(31,28)


species.list <-vector("list", length(species.names))
names(species.list) <- species.names

pars.names <- c("max.age", "M","fec","breed.freq","age.mat")
pars.list <- vector("list", length(pars.names))
names(pars.list) <- pars.names

#Fill in species list of pars
for (i in 1:N.sp)
{
  species.list[[i]]=list(max.age=c(Min.max.age[i],Max.max.age[i]),fec=c(Min.fec[i],Max.fec[i]),
                    breed.freq=c(Min.breed.freq[i],Max.breed.freq[i]),
                    age.mat=c(Min.age.mat[i],Max.age.mat[i]))
}


#Average water temperature
Temperature=c(18,18)

#Growth pars (female)
Linf.g=201.9
Linf.w=120.7 #Simfendorfer et al 2000, age and growth
K.g=0.123
K.w=0.369
to.g=-1.55
to.w=-0.6


#Size at birth
Size.birth=c(33.5,25)

#FL to TL pars             
b.g=NA
b.w=8.891
a.g=NA
a.w=1.046

#TL to TW pars 
bwt.g=0.927e-6  #(reported as 0.927e-9 but I changed scale to match others)
bwt.w=0.0000163
awt.g=3.206
awt.w=2.733


#Selectivity
  #gummy                #Walker 2010
theta1=186
theta2=36695
mesh= 6.5 #(in inches, = 16.5 cm)
alphabeta.g=theta1*mesh
beta.g=-0.5*(theta1*mesh-((theta1^2)*(mesh^2)+4*theta2)^0.5)
alpha.g=alphabeta.g/beta.g

   #whiskery                #Simpfendorfer & Unsworth 1998
 alpha.w=64.01339
 beta.w=18.53164
 alphabeta.w=alpha.w*beta.w
 


#Put pars in data frame for loop
Growth.pars=data.frame(Linf=c(Linf.g,Linf.w),K=c(K.g,K.w),to=c(to.g,to.w))
Length.conv=data.frame(b=c(b.g,b.w),a=c(a.g,a.w))
Weight.conv=data.frame(bwt=c(bwt.g,bwt.w),awt=c(awt.g,awt.w))
Sel.pars=data.frame(alphabeta=c(alphabeta.g,alphabeta.w),
                    alpha=c(alpha.g,alpha.w),
                    beta=c(beta.g,beta.w))


#---PROCEDURE SECTION----

#1. Priors
  #1.1. Natural mortality

  #1.2 Max age
A.max=function(MIN) MIN

  #1.3 Fecundity
FEC=function(spec.)
{
  if (spec.=="gummy")Fec=1.12*exp(-3.033+0.00396*total.length*10)  #add 10 scaler as it was calculated in mm
  if (spec.=="whiskery")
  {
    #Fec=0.314*mid.FL.fem-17.8    #underestimates Fec because mid.FL only gets to Linf (i.e. 120)!!
    Fec=rep(round(rtriangle(1, a=Min.fec[2], b=Max.fec[2], c=mean(c(Min.fec[2],Max.fec[2])))),
          length(mid.FL.fem))
    
    if(USE.MEAN=="YES")Fec=rep(mean(c(Min.fec[2],Max.fec[2])),length(mid.FL.fem))
    if(USE.MEAN=="Linear")Fec=round(seq(Min.fec[2],Max.fec[2],length.out=length(mid.FL.fem)))
  }
  
  if (spec.=="dusky")
  {
    Fec=rep(round(rtriangle(1, a=Min.fec[3], b=Max.fec[3], c=mean(c(Min.fec[3],Max.fec[3])))),
            length(mid.FL.fem))
    if(USE.MEAN=="YES")Fec=rep(mean(c(Min.fec[3],Max.fec[3])),length(mid.FL.fem))
    if(USE.MEAN=="Linear")Fec=round(seq(Min.fec[3],Max.fec[3],length.out=length(mid.FL.fem)))
  }
  if (spec.=="sandbar")
  {
    Fec=rep(round(rtriangle(1, a=Min.fec[4], b=Max.fec[4], c=mean(c(Min.fec[4],Max.fec[4])))),
            length(mid.FL.fem))
    if(USE.MEAN=="YES")Fec=rep(mean(c(Min.fec[4],Max.fec[4])),length(mid.FL.fem))
    if(USE.MEAN=="Linear")Fec=round(seq(Min.fec[4],Max.fec[4],length.out=length(mid.FL.fem)))
  }
  
  return(Fec)
}

  #1.4 Reproductive periodicity
REP.PER=function(MIN,MAX)runif(1,MIN,MAX)

  #1.5 age at maturity
AGE.MAT=function(MIN,MAX)floor(runif(1,MIN,MAX))  


  #1.6 growth pars (multivariate normal distribution)
GROWTH=function(Linf.mean,k.mean,sigma)
{
  growth.pars.sim=rmvnorm(1,mean=c(Linf.mean,k.mean),sigma=sigma)
  if(growth.pars.sim[2]<=0.001)    #repeat until sensible pars obtained
  { repeat 
    {
      growth.pars.sim=rmvnorm(1,mean=c(Linf.mean,k.mean),sigma=sigma)
      if(growth.pars.sim[2]>0.001)break
    }
  }
  return(growth.pars.sim)
}


#2. Spawning potential ratio method (Brooks et al 2010)
  #2.1. Analytically-derived SPR.mer
SPR=function(max.age,M,fec,breed.freq,sex.ratio,age.mat,sel.age,F.mult=0,Plus.gp="Y",MAT)
{
  #age
  age=0:max.age
  
  #survivorship
  surv=exp(-M)
  
  #fecundity  
  #fecundity=rep(fec*breed.freq*sex.ratio,length(age))
  fecundity=fec*breed.freq*sex.ratio
  
  #maturity
  if(MAT=="knife")maturity=ifelse(age>=age.mat,1,0)   
  if(MAT=="ogive")maturity=plogis(age,age.mat,1)
    
  #if maximum age is a terminal group
  if(!(Plus.gp=="Y"))
  {
    phi.o=0.0         #unexploited spawners per recruit
    cum.survive=1.0   #cumulative survival
    z=0.0
    for (i in 2:(max.age)  )
    {
      z=M[i] + F.mult*sel.age[i]
      z.ts=(M[i]+F.mult*sel.age[i])*spawn.time
      phi.o=phi.o+cum.survive*fecundity[i]*maturity[i]*exp(-z.ts)
      cum.survive=cum.survive*exp(-z )
    }
  }
  
  #if maximum age is plus group
  if(Plus.gp=="Y")
  {
    phi.o=0.0
    cum.survive=1.0
    z=0.0
    for (i in 2:(max.age)  )
    {
      z=M[i] + F.mult*sel.age[i]
      z.ts=(M[i]+F.mult*sel.age[i])*spawn.time
      phi.o=phi.o+cum.survive*fecundity[i]*maturity[i]*exp(-z.ts)
      cum.survive=cum.survive*exp(-z )
    }
    #plus group  
    z= M[max.age+1] + F.mult*sel.age[max.age+1]
    z.ts=(M[max.age+1]+F.mult*sel.age[max.age+1])*spawn.time
    phi.o=phi.o + fecundity[max.age+1]*maturity[max.age+1]*cum.survive*exp(-z.ts)/( 1- exp(-z ) )
  }
  
  #maximum lifetime reproductive rate at low density
  alpha=phi.o*surv[1]
  
  #steepness
  h=alpha/(4+alpha)
  
  #spawning potential ratio at maximum excess recruitment (MER) (Beverton-Holt relationship)
  SPR.mer=1/alpha^0.5

  #optimal depletionlevel (i.e.depletion at MER, the proportional reduction from unexploited level)
  Dep.MER=((alpha^0.5)-1)/(alpha-1) 
  
  #overfished threshold
  p=max(c(0.5,(1-mean(M))))

  if(h<0.2)return(Warning="Bad joint priors")
  if(h>=0.2)return(list(phi.o=phi.o,alpha=alpha,h=h,SPR.mer=SPR.mer,Dep.MER=Dep.MER,p=p,
                        fecundity=fecundity,maturity=maturity))
}

Stock.depletion=function(ab.index,init.depl,dep.mer,p)
{
  #scaled current index
  Scaled.relative.index=init.depl*(ab.index[length(ab.index)]/ab.index[1])
  
  #Current depletion level
  Depletion=Scaled.relative.index/dep.mer
  
  #Stock status
  if(Scaled.relative.index<(p*dep.mer)){Stock.status="Overfished"}else
  {Stock.status="Not overfished"}      

  return(Stock.status)
}


  #2.2 Find Fmer thru numerical approximation
#note: Equilibrium Per recruit function
EqPerRec=function(Fe=0)
{
  #survivorship
  surv=exp(-mm)    
  
  #note: age, fecundity and maturity are set outside the function (in loop)
  
  #eggs
  phieggs=0.0
  
  #survival
  cum.survive=1.0
  z=0.0
  for (i in 2:(A.MAX)  )
  {
    z=mm[i] + Fe*Sel.A[i]
    z.ts=(mm[i]+Fe*Sel.A[i])*spawn.time
    phieggs=phieggs+cum.survive*fecundity[i]*maturity[i]*exp(-z.ts)
    cum.survive=cum.survive*exp(-z )
  }
  #plus group  
  z= mm[A.MAX+1] + Fe*Sel.A[A.MAX+1]
  z.ts=(mm[A.MAX+1]+Fe*Sel.A[A.MAX+1])*spawn.time
  phieggs=phieggs + fecundity[A.MAX+1]*maturity[A.MAX+1]*cum.survive*exp(-z.ts)/( 1- exp(-z ) )
  
  
  #Spawning potential ratio
  SPR=phieggs/phieggs0
  
  #Objective function
  epsilon=(Observed-SPR)^2   				
    
  return(list(phieggs0=phieggs0,phieggs=phieggs,epsilon=epsilon,SPR=SPR))
}

  #2.3 Find SPR for given depletion level thru numerical approximation
Find.SPR=function(Spr)
{
  DEPLET=((((1/Spr)^2)^0.5)-1)/(((1/Spr)^2)-1)
  
  #Objective function
  epsilon=(Observed-DEPLET)^2   				#objective function
    
  return(list(epsilon=epsilon,Depletion=DEPLET))
}

  #2.4 Find B from a given SPR level (B-H function)
Find.B<-function(recr.par, R0, spr, spr0)
{
  sprF=spr/spr0
  ssb=spr0*R0*(sprF*recr.par - 1.0)/(recr.par-1.0)
  return(ssb)
}



#3. Zhou et al (2012) FMSY proxy
Zhou=function(M) FMSY=0.41*M


#4. Standard per recruit analysis 
Per.recruit=function(Linf,K,to,b,a,M,bwt,awt,age.mat,fec,breed.freq,fishing)
{
  N=length(age)
  
    #length at mid age
  FL=Linf*(1-exp(-K*(age+0.5-to)))
  TL=b+a*FL
  
    #mortality at age
  Fat.age=fishing*Sel.A
  Z.at.age=Fat.age+M

    #numbers at age
  Rel.numbers=vector(length=N)
  Rel.numbers[1]=1
  for(i in 2:N)Rel.numbers[i]=Rel.numbers[i-1]*exp(-Z.at.age[i])
  
    #harvest fraction
  harvest=(Fat.age/Z.at.age)*Rel.numbers*(1-exp(-Z.at.age))
 
    #weigth at age
  Wt=bwt*TL^awt
  
    #fecundity at age (considering Maturity and breeding freq)
  fecundity=ifelse(age>=age.mat,fec*breed.freq,0)
  
    #YPR
  YPR=Wt*harvest
  
    #Spawning stock biomass per recruit                           
  SSB=ifelse(age>=age.mat,Wt*Rel.numbers,0)
  
    #Eggs per recruit    (female)                                         
  Eggs=fecundity*sex.ratio*Rel.numbers
  
    #Objective functions
  Epsilon=-1*sum(YPR)
  Epsilon.mer=-1*sum(Eggs)
  
    #Totals per recruit
    return(list(N.per.rec=sum(Rel.numbers),Yield.per.rec=sum(YPR),SSB.per.rec=sum(SSB),
            Eggs.per.rec=sum(Eggs),Epsilon=Epsilon,Epsilon.mer=Epsilon.mer))
  
}


#5. Selectivity function (note: Fork length in mm)
Select=function(alphabeta,alpha,beta)
{
  sel.fem=((mid.FL.fem*10/alphabeta)^alpha)*(exp(alpha-(mid.FL.fem*10/beta)))
}


#---MAIN SECTION----
N.sim=10000

#1. Spawning Potential Ratio and derived methods
Species.SPR.MER=Species.B.MER=Species.F.MER=Species.F.est.Conv=Species.M=Species.Sel=Species.Fec=
  Species.len=Species.F.40=Species.F.30=Species.F.zhou=Species.Trigged=vector("list", length=N.sp)
names(Species.SPR.MER)=names(Species.B.MER)=names(Species.F.MER)=names(Species.F.est.Conv)=
  names(Species.M)=names(Species.Sel)=names(Species.Fec)=names(Species.F.40)=names(Species.F.30)=
  names(Species.F.zhou)=names(Species.Trigged)=names(Species.len)=species.names

#Monte Carlo simulations
#note: obtain stock rec pars and biological reference points

#add or remove variability in growth pars
Growth.var="YES"
#Growth.var="NO"

#use Mean Fecundity for whiskery, dusky and sandbar
#USE.MEAN="YES"   #for using mean
#USE.MEAN="Linear" #for using linear relation
USE.MEAN="NO"     #for using random sample

b.f <- function(Fstart)
{
  temp.spr = EqPerRec(Fstart)$SPR    
  temp.SSB =Find.B(spr.temp$alpha,R0=1,spr=temp.spr,spr0=1)  #R0 and spr0 =1 as it's on relative scale
  yy<-abs(temp.SSB - obs.biom )
  return(yy)
}
Fstart=0.1

system.time(for (a in 1:N.sp)
{
  spec.=species.names[a]
  WT=Weight.conv[a,]  
  GROW=Growth.pars[a,]
  TO=Growth.pars[a,3]
  SIG=SIGMA[[a]]
  Lo=Size.birth[a]
  AA=species.list[[a]]$max.age
  FF=species.list[[a]]$fec
  BF=species.list[[a]]$breed.freq
  b.fem=Length.conv[a,1]
  a.fem=Length.conv[a,2]
  alphabeta.S=Sel.pars[a,1]
  alpha.S=Sel.pars[a,2]
  beta.S=Sel.pars[a,3]
  sex.ratio=0.5
  AMat=species.list[[a]]$age.mat
  Temper=Temperature[a]
  r=1
  spawn.time = 0 # specify time of the year when spawning (or pupping) occurs as a fraction beteween 0 and 1
  
  w.g=Wght.G[a]
  w.ng=Wght.noG[a] 
  
  SPR.out=f.out=b.out=Store.M=Store.Fec=Convergence=Store.Sel=Store.len=f.out.30=f.out.40=
    f.zhou=A.Max=vector("list", length=N.sim)
  
  Trigged=rep(NA,N.sim)
  
  #Threshold
  for (j in 1:N.sim)
  {
    #1. draw random samples of input parameters
    A.MAX=A.max(AA[1])
    age=0:A.MAX
    
    if(Growth.var=="YES") GR=GROWTH(GROW[[1]],GROW[[2]],SIG)
    if(Growth.var=="NO") GR=matrix(c(GROW[[1]],GROW[[2]]),ncol=2)
    
    Rep=REP.PER(BF[2],BF[1])
    A.MAT=AGE.MAT(AMat[1],AMat[2])    
    
    if(!spec.=="gummy")
    {
      mid.FL.fem=Lo+(GR[1]-Lo)*(1-exp(-GR[2]*age))  #modified version
      total.length=b.fem+a.fem*mid.FL.fem
    }
    
    if(spec.=="gummy")
    {
      total.length=Lo+(GR[1]-Lo)*(1-exp(-GR[2]*age)) #gummy pars are fn(total length)
      mid.FL.fem=total.length
    }
    
    Fec=FEC(spec.)
    
    #put a cap on gummy Fecundity to avoid predicting beyond data range
    if(spec.=="gummy")
    {
      Fec=ifelse(Fec>Max.fec[1],Max.fec[1],Fec)
      Fec=ifelse(Fec<0,0,Fec)
    }
    
    if(spec.=="gummy") mm=rep(M.gummy,length(age))
    if(spec.=="whiskery") mm=rep(M.whiskery,length(age))
    Sel.A=Select(alphabeta.S,alpha.S,beta.S)
    
    
    #2. calculate SPR.mer quantities 
    spr.temp=SPR(A.MAX,mm,Fec,Rep,sex.ratio,A.MAT,Sel.A,0,"Y","ogive")
    
    if(length(spr.temp)>1)
    {
      #3. Calculate Threshold (i.e. MER)reference points
      #3.1.1 Biomass (Bmer)
      BTHR=spr.temp$Dep.MER
    }
    
    
    #3. Repeat if nonsense outputs (e.g. h<0.2, negative Biomass RP), i.e. obtain sensible joint priors 
    Trigger=0
    if (length(spr.temp)==1|(BTHR*(1/Prop.FMsy))<0|(BTHR*Prop.FMsy)<0|BTHR<0)Trigger=1
    Trigged[j]=Trigger
    if(Trigger==1)    
    { repeat 
    {
      #3.1 draw random samples
      A.MAX=A.max(AA[1])
      age=0:A.MAX
      
      if(Growth.var=="YES") GR=GROWTH(GROW[[1]],GROW[[2]],SIG)
      if(Growth.var=="NO") GR=matrix(c(GROW[[1]],GROW[[2]]),ncol=2)
      
      Rep=REP.PER(BF[2],BF[1])
      A.MAT=AGE.MAT(AMat[1],AMat[2])    
      
      if(!spec.=="gummy")
      {
        mid.FL.fem=Lo+(GR[1]-Lo)*(1-exp(-GR[2]*age))  #modified version
        total.length=b.fem+a.fem*mid.FL.fem
      }
      
      if(spec.=="gummy")
      {
        total.length=Lo+(GR[1]-Lo)*(1-exp(-GR[2]*age)) #gummy pars are fn(total length)
        mid.FL.fem=total.length
      }
      
      Fec=FEC(spec.)
      
      #3.2 put a cap on gummy Fecundity to avoid predicting beyond data range
      if(spec.=="gummy")
      {
        Fec=ifelse(Fec>Max.fec[1],Max.fec[1],Fec)
        Fec=ifelse(Fec<0,0,Fec)
      }
      
      
      if(spec.=="gummy") mm=rep(M.gummy,length(age))
      if(spec.=="whiskery") mm=rep(M.whiskery,length(age))
      Sel.A=Select(alphabeta.S,alpha.S,beta.S)
      
      
      #3.3 calculate SPR.mer quantities
      spr.temp=SPR(A.MAX,mm,Fec,Rep,sex.ratio,A.MAT,Sel.A,0,"Y","ogive")
      
      if(length(spr.temp)>1)
      {
        #3.3.1 Calculate reference points
        
        #3.3.1.1 Threshold (i.e. MER)    
        
        #Biomass (Bmer)
        BTHR=spr.temp$Dep.MER       
      }
      
      #break if sensible joint prior
      Trigger=0
      
      if (length(spr.temp)==1|(BTHR*(1/Prop.FMsy))<0|(BTHR*Prop.FMsy)<0|BTHR<0)Trigger=1
      
      if(Trigger==0)break
    }
    }
    
    #store quantities of interest
    SPR.out[[j]]=spr.temp
    Store.M[[j]]=mm
    Store.Sel[[j]]=Sel.A
    Store.len[[j]]=mid.FL.fem
    Store.Fec[[j]]=Fec
    b.out[[j]]=BTHR
    A.Max[[j]]=A.MAX
  }
  
  
  BTHR=unlist(b.out)
  
  #Biomass limit (Blim = p x BTHR)
  b.limt=Prop.FMsy*BTHR      
  b.limt=BTHR-(mean(BTHR)-mean(b.limt))
  
  #Biomass target (Blim = p x BTHR)
  b.targt=(1/Prop.FMsy)*BTHR 
  b.targt=BTHR+(mean(b.targt)-mean(BTHR))
  
  
  for (j in 1:N.sim)
  {
    
    #3.1.2. F (Fmer)
    #Set up useful vectors
    spr.temp=SPR.out[[j]]
    A.MAX=A.Max[[j]]
    age=0:A.MAX
    mm=Store.M[[j]]
    Sel.A=Store.Sel[[j]]
    fecundity=spr.temp$fecundity
    maturity=spr.temp$maturity
    lx=vector(length=length(age))    #no fishing
    lz=vector(length=length(age))    #with fishing
    Observed=spr.temp$SPR.mer
    phieggs0=spr.temp$phi.o
    fn.yield=function(Fe) EqPerRec(Fe=Fe)$epsilon
    
    #numerical approximation
    Fmer.fit <- nlminb( start=Fstart, objective=fn.yield, lower=0.0001, upper=5.0,
                        control=list(eval.max=500, iter.max=500  )   )
    FTHR=Fmer.fit$par
    FTHR.Conv=Fmer.fit$convergence
    
    #Store convergence code
    Convergence[[j]]=c(FTHR.Conv=FTHR.Conv)
    
    #3.2.2 F limit
    obs.biom <- b.limt[j]   
    Flim.fit = nlminb(start=Fstart, objective=b.f, lower=0, upper=3)
    FLim = Flim.fit$par
    Flim.Conv = Flim.fit$convergence
    
    
    #3.3. F Target
    obs.biom <- b.targt[j]   
    FTar.fit = nlminb(start=Fstart, objective=b.f, lower=0, upper=3)
    FTar = FTar.fit$par
    FTar.Conv = FTar.fit$convergence
    
    
    #     #5. Calculate F reference points for arbitrary SPR
    #     Observed=0.40 #target (Tsai et al 2011 page 1388)
    #     F40.fit=optimize(fn.yield,lower=0.0001,upper=5)   #numerical approximation
    #     f.out.40[[j]]=F40.fit$minimum
    #     
    #     Observed=0.30 #limit
    #     F30.fit=optimize(fn.yield,lower=0.0001,upper=5)   #numerical approximation
    #     f.out.30[[j]]=F30.fit$minimum
    #     
    #     #6. Calculate Fmsy proxy from Zhou et al 2012
    #     f.zhou[[j]]=Zhou(mean(mm))  
    
    
    #store
    b.out[[j]]=c(BTHR[j],b.limt[j],b.targt[j])
    f.out[[j]]=c(FTHR,FLim,FTar)
    
  }
  
  
  #STORE  
  Species.SPR.MER[[a]]=SPR.out
  Species.M[[a]]=Store.M
  Species.Fec[[a]]=Store.Fec
  Species.Sel[[a]]=Store.Sel
  Species.len[[a]]=Store.len
  Species.B.MER[[a]]=b.out
  Species.F.MER[[a]]=f.out
  Species.F.est.Conv[[a]]=do.call(rbind,Convergence)
  Species.F.40[[a]]=f.out.40
  Species.F.30[[a]]=f.out.30
  Species.F.zhou[[a]]=f.zhou
  Species.Trigged[[a]]=Trigged
})


#Extract reference points (MSY proxies)
B.lim=FMER=BMER=B.thre=B.tar=F40=FZHOU=F30=Steep=SprMer=vector("list", length=N.sp)
names(B.lim) =names(BMER) =names(FMER)= names(B.tar)=names(B.thre)=names(Steep)=names(SprMer)<-species.names
for (i in 1:N.sp)
{
  Blim=Bthre=Btar=f40=fZhou=f30=H=SPRMER=data.frame(title=rep(NA,N.sim))
  Fmer=data.frame(Fmer=rep(NA,N.sim),Flim=rep(NA,N.sim),Ftar=rep(NA,N.sim))
  Bmer=data.frame(Bmer=rep(NA,N.sim),Blim=rep(NA,N.sim),Btar=rep(NA,N.sim))
  for (j in 1:N.sim)
  {
    #Biomass    

    Bmer[j,]=Species.B.MER[[i]][[j]]

    #steepness
    H[j,1]=Species.SPR.MER[[i]][[j]]$h
    
    #SPR.mer
    SPRMER[j,1]=Species.SPR.MER[[i]][[j]]$SPR.mer
    
    #F
  #  f30[j,1]=Species.F.30[[i]][[j]]
    Fmer[j,]=Species.F.MER[[i]][[j]]
  #  f40[j,1]=Species.F.40[[i]][[j]]
  #  fZhou[j,1]=Species.F.zhou[[i]][[j]]
  }
  
  #B.thre[[i]]=Bthre
  #B.lim[[i]]=Blim
  BMER[[i]]=Bmer
  FMER[[i]]=Fmer
  B.tar[[i]]=Btar
 # FZHOU[[i]]=fZhou
#  F40[[i]]=f40
#  F30[[i]]=f30
  Steep[[i]]=H
  SprMer[[i]]=SPRMER
}


#ACA
Steep.list=Steep
for (i in 1:N.sp)
{
  this=Steep[[i]][,1]
  Steep.list[[i]]=round(c(median=median(this),sixty.per=quantile(this,probs=0.6),
                          CI95=sort(this)[c(floor(0.025*length(this)),ceiling(0.975*length(this)))]),3)
}
Steep.list
M.gummy
M.whiskery
