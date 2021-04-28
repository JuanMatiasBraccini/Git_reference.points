###########################MSY PROXY REFERENCE POINTS###########################

#notes: This script implements different MSY proxy reference points

# FIRST: the Spawning Potential Ratio reference points derived in Brooks
#       et al (2010) for assessment of data poor fisheries using a biological reference point 
#       which is compared to abundance data ("overfished state").
#       To determine if "Overfishing" is occurring (page 169), need current estimate of F (tagging, etc)

#       Also, due to the relationship between steepness (lower bound=0.2) and alpha, the Maximum
#       Reproductive Rate at low density (lower bound=1), this approach can also be used to evaluate
#       if life history parameter combinations satisfy the lower bounds

#       Finally, by reparameterising the S-R relationship in terms of alpha, 
#       we avoid the need for estimating steepness!!

#       Gummy shark not included as whiskery is used as representative of this life history group

#       assumptions:   for steepness:  if using maturity schedule, assumes that reproduction occurs immediately once
#                               reaching maturing for first age at maturity
#                       for reference point: in addition, these assumes constant catchability and selectivity thru time



#       required:  life history parameters and an index of abundance for stock status assessment
#                  life history data only for steepness

#      F.mer is estimated by finding F that yields SPR.mer
#      F.lim is estimated by finding F that yields SPR.lim, which is estimated by finding the SPR.lim
#           that yields p*optimal depletion
#      F.targ is estimated by finding F that yields SPR.tar, which is estimated by finding the SPR.tar
#           that yields 1.2*optimal depletion

# SECOND: Zhou et al (2012) FMSY proxy


# THIRD:  Per recruit analysis



#VERY IMPORTANT: check if simulations make sense (growth in particular) by tracking
                # selectivity and growth curve. Resample if too much variability

rm(list=ls(all=TRUE))
library(triangle)
library(TeachingDemos)      #for grey scale
library(mvtnorm)      #for multivariate normal pdf
library(PBSmapping)
if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')

#source sigmas for var-covar matrix
source(handl_OneDrive("Analyses/Reference Points/5.Derive var_covar Matrix.r"))

#source indirect estimation of M
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/M_indirect.R"))



#---DATA SECTION----

#Abundance index
species.names <- c("gummy","whiskery", "dusky","sandbar")
abundance.species.list <-depletion.species.list<- vector("list", length(species.names))
names(abundance.species.list) <-names(depletion.species.list) <- species.names
N.sp=length(species.names)

  #abundance series
abundance.species.list$whiskery=c(3.76,0.81)    #first and last year of cpue data.



  #initial depletion level of the abundance series
depletion.species.list$Whiskery=0.95       #dummy. Use real value. Note: if catchability/selectivity change, depletion 
                                                                    # will change if inferred from cpue


#Empirical F at age estimates (McAuley et al 2007)
Dusky.F.emp=c(.21,.145,.05,.03,.03,.025,.022,.02,.03,.01,0)
names(Dusky.F.emp)=0:10 #age classes studied
Sandb.F.emp=c(0,0.06,.28,.08,.04,.03,.04,.03,.025,.022,.02,.01,0)
names(Sandb.F.emp)=c("0-3","3-6","6-9","9-12","12-15","15-18","18-24",25:30)
Empirc.F=list(Dusky=Dusky.F.emp,Sandbar=Sandb.F.emp)


#Maturity ogive
  #Whiskery
# FL.w=c(92,96,100,104,108,112,116,120,124,128,132)
# Prop.M.w=c(0,0,0,0.35,0.4,0.25,0.58,0.7,0.8,1,1)
# fit.Mat=function(DAT,theta)fit=nls(Prop.M~1/(1+exp(-log(19)*((FL-L50)/(L95-L50)))),data=DAT,start=theta)
# Whis.Fit=fit.Mat(DAT=data.frame(Prop.M=Prop.M.w,FL=FL.w),theta=c(L50=115,L95=125))



#---PARAMETERS SECTION----

#Management decisions

Prop.FMsy=0.75  #management parameter. Proportion of Fmsy (as used by NMFS)
#Prop.FMsy=0.65  #sensitivity
#Prop.FMsy=0.85


#Life history par vectors (Gummy,Whiskery,Dusky,Sandbar)
  #Min values
Min.max.age=c(16,15,40,30)
Min.breed.freq=c(1,0.5,0.5,0.5)
Min.age.mat=c(4,6,26,13)
Min.fec=c(1,4,2,4)

  #Max values
Max.max.age=c(21,20,56,40) #put +1 because using floor in random sample
Max.breed.freq=c(1,0.5,0.333,0.5)
Max.age.mat=c(6,8,36,20) #put +1 because using floor in random sample
Max.fec=c(31,28,18,10)

  #Mean values
#Whisk.Fec.relation=0.314*(100:140)-17.8 #Simfendorfer & Unsworth 1998
#Mean.fec=c(mean(Whisk.Fec.relation),9.9,7.7) #mean Sandbar fec from (Table 4 Baremore & Hale 2010) 
#Mean.fec=c(mean(Whisk.Fec.relation),9.9,6.5)
#Sd.fec=c(sd(Whisk.Fec.relation),2.7,2.3)  #sd=2 fec is CV=30% as found for dusky



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
Temperature=c(18,18,18,24)

#Growth pars (female)
Linf.g=201.9
#Linf.w=128.2  #Simfendorfer et al 2000, tag-recapture
Linf.w=120.7 #Simfendorfer et al 2000, age and growth
Linf.d=374.4
Linf.s=244.2

K.g=0.123
#K.w=0.288
K.w=0.369
K.d=.0367
K.s=.040

to.g=-1.55
to.w=-0.6
to.d=-3.3
to.s=-4.8

#Size at birth
Size.birth=c(33.5,25,75.3,42.5)

#FL to TL pars             
b.g=NA
b.w=8.891
b.d=4.000
b.s=5.8188
a.g=NA
a.w=1.046
a.d=1.100
a.s=1.11262


#TL to TW pars 
bwt.g=0.927e-6  #(reported as 0.927e-9 but I changed scale to match others)
bwt.w=0.0000163
bwt.d=1.2334e-5
bwt.s=0.000002
awt.g=3.206
awt.w=2.733
awt.d=2.855
awt.s=3.2069

#Selectivity
  #gummy                #Walker 2010
theta1=186
theta2=36695
mesh= 6.5 #(in inches, = 165 mm)
alphabeta.g=theta1*mesh
beta.g=-0.5*(theta1*mesh-((theta1^2)*(mesh^2)+4*theta2)^0.5)
alpha.g=alphabeta.g/beta.g

mesh= 7 #(in inches, = 178 mm)
alphabeta.g7=theta1*mesh
beta.g7=-0.5*(theta1*mesh-((theta1^2)*(mesh^2)+4*theta2)^0.5)
alpha.g7=alphabeta.g7/beta.g7


   #whiskery                #Simpfendorfer & Unsworth 1998
 # alpha.w=64.01339
 # beta.w=18.53164
 # alphabeta.w=alpha.w*beta.w   #old, used in paper as it was the values used in Simpfendorfer for stock assessment
 
 
 theta1=173.70
 theta2=26415
 mesh= 6.5 #(in inches, = 16.5 cm)
 alphabeta.w=theta1*mesh
 beta.w=-0.5*(theta1*mesh-((theta1^2)*(mesh^2)+4*theta2)^0.5)
 alpha.w=alphabeta.w/beta.w
 
 mesh= 7 #(in inches)
 alphabeta.w7=theta1*mesh
 beta.w7=-0.5*(theta1*mesh-((theta1^2)*(mesh^2)+4*theta2)^0.5)
 alpha.w7=alphabeta.w7/beta.w7
 
 
   #dusky                #Simpfendorfer & Unsworth 1998
 theta1=130.13
 theta2=29237
 mesh= 6.5 #(in inches, = 16.5 cm)
 alphabeta.d=theta1*mesh
 beta.d=-0.5*(theta1*mesh-((theta1^2)*(mesh^2)+4*theta2)^0.5)
 alpha.d=alphabeta.d/beta.d

 mesh= 7 #(in inches)
 alphabeta.d7=theta1*mesh
 beta.d7=-0.5*(theta1*mesh-((theta1^2)*(mesh^2)+4*theta2)^0.5)
 alpha.d7=alphabeta.d7/beta.d7
 
 
   #sandbar                #McAuley et al 2005
 theta1=135.58
 theta2=117001
 mesh= 6.5
 alphabeta.s=theta1*mesh
 beta.s=-0.5*(theta1*mesh-((theta1^2)*(mesh^2)+4*theta2)^0.5)
 alpha.s=alphabeta.s/beta.s

 mesh= 7
 alphabeta.s7=theta1*mesh
 beta.s7=-0.5*(theta1*mesh-((theta1^2)*(mesh^2)+4*theta2)^0.5)
 alpha.s7=alphabeta.s7/beta.s7
 

#Put pars in data frame for loop
Growth.pars=data.frame(Linf=c(Linf.g,Linf.w,Linf.d,Linf.s),K=c(K.g,K.w,K.d,K.s),to=c(to.g,to.w,to.d,to.s))
Length.conv=data.frame(b=c(b.g,b.w,b.d,b.s),a=c(a.g,a.w,a.d,a.s))
Weight.conv=data.frame(bwt=c(bwt.g,bwt.w,bwt.d,bwt.s),awt=c(awt.g,awt.w,awt.d,awt.s))
Sel.pars=data.frame(alphabeta=c(alphabeta.g,alphabeta.w,alphabeta.d,alphabeta.s),
                    alpha=c(alpha.g,alpha.w,alpha.d,alpha.s),
                    beta=c(beta.g,beta.w,beta.d,beta.s))


#---PROCEDURE SECTION----

#1. Priors
  #1.1. Natural mortality
#source external function 
Wght.G=c(1,1,1,1)    # weight given to methods using growth parameters (high uncertainty for whiskery so give lower weight)
Wght.noG=c(1,10,1,1)  # weight given to methods not using growth parameters

  #1.2 Max age
A.max=function(MIN,MAX,MODE) floor(rtriangle(1, a=MIN, b=MAX, c=MODE))

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
    A.MAX=A.max(AA[1],AA[2],AA[1])
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
    
    mm=M.fun(A.MAX,GR[2],GR[1],Temper,A.MAT,TO,WT[1,1],WT[1,2],w.g,w.ng)
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
      A.MAX=A.max(AA[1],AA[2],AA[1])
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
      
      
      mm=M.fun(A.MAX,GR[2],GR[1],Temper,A.MAT,TO,WT[1,1],WT[1,2],w.g,w.ng)
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


#Model checks
  #1 .extract times priors had to be resampled (0= no resampling, 1= resampling)
TAB.Trigger=lapply(Species.Trigged,table)
names(TAB.Trigger)=species.names


  #2. check convergence (0= converged)
TAB.converge=lapply(Species.F.est.Conv,table)
names(TAB.converge)=species.names
print(TAB.converge)

  #3. check steepness and SPRmer
check.cols=1:4
par(mfcol=c(1,1),las=1,mai=c(.8,1,.075,.15),omi=c(.1,.1,.1,.05),mgp=c(1.5,0.95,0))
one=density(Steep[[1]][,1],adjust=2,na.rm =T,from=0.2,to=1)
two=density(Steep[[2]][,1],adjust=2,na.rm =T,from=0.2,to=1)
three=density(Steep[[3]][,1], adjust=2,na.rm =T,from=0.2,to=1)
four=density(Steep[[4]][,1], adjust=2,na.rm =T,from=0.2,to=1)
plot(one, type="l",lty=1,xlab="",lwd=2.5, yaxs="i",xaxs="i",ylab="",main="",
       xlim=c(0.2,1),ylim=c(0, max(one$y,two$y,three$y,four$y)),cex.axis=1.5,cex.lab=1.65,col=check.cols[1])
lines(two, lty=1,col=check.cols[2],lwd=2.5)
lines(three, lty=1,col=check.cols[3],lwd=2.5)
lines(four, lty=1,col=check.cols[4],lwd=2.5)
legend("topright",species.names,lty=1,cex=1.5,bty="n",col=check.cols,lwd=2)
mtext("Density",side=2,line=-1.5,font=1,las=0,cex=1.3,outer =T)
mtext("Steepness",side=1,line=2.1,font=1,las=0,cex=1.3)


one=density(SprMer[[1]][,1],adjust=2,na.rm =T,from=0.2,to=1)
two=density(SprMer[[2]][,1],adjust=2,na.rm =T,from=0.2,to=1)
three=density(SprMer[[3]][,1], adjust=2,na.rm =T,from=0.2,to=1)
four=density(SprMer[[4]][,1], adjust=2,na.rm =T,from=0.2,to=1)
plot(one, type="l",lty=1,xlab="",lwd=2.5, yaxs="i",xaxs="i",ylab="",main="",
     xlim=c(0.2,1),ylim=c(0, max(one$y,two$y,three$y)),cex.axis=1.5,cex.lab=1.65,col=check.cols[1])
lines(two, lty=1,col=check.cols[2],lwd=2.5)
lines(three, lty=1,col=check.cols[3],lwd=2.5)
lines(four, lty=1,col=check.cols[4],lwd=2.5)
legend("topright",species.names,lty=1,cex=1.5,bty="n",col=check.cols,lwd=2)
mtext("Density",side=2,line=-1.5,font=1,las=0,cex=1.3,outer =T)
mtext("SprMer",side=1,line=2.1,font=1,las=0,cex=1.3)



  #4.check variability in growth, selectivity and M
par(mfcol=c(4,3),las=1,mai=c(.6,.5,.075,.15),omi=c(.1,.1,.1,.05),mgp=c(1.5,0.6,0))
fn.plot.1=function(DATA,LEG)
{  plot(DATA[[1]],col="transparent",ylab="",xlab="")
   for (j in 1:N.sim)  lines(  DATA[[j]])
   legend("topleft",LEG,bty='n',cex=.9)
}
for(p in 1:N.sp)fn.plot.1(Species.Sel[[p]],"Select")
for(p in 1:N.sp)fn.plot.1(Species.len[[p]],"Growth")
for(p in 1:N.sp)
  {
  fn.plot.1(Species.M[[p]],"M")
  legend("topright",species.names[p],cex=1.5,bty="n")
  }


#Stock satus assessment
# Whiskery.stock.status=Stock.depletion(abundance.species.list$Whiskery,depletion.species.list$Whiskery,
#                                       Species.SPR.MER$Whiskery$Dep.MER,Species.SPR.MER$Whiskery$p)





#---REPORT SECTION----

SP.names=c("gummy shark","whiskery shark","dusky shark","sandbar shark")

setwd(handl_OneDrive("Analyses/Reference Points/Outputs"))
#setwd("C:/Matias/Analyses/Reference Points/Outputs/p_0.65")
#setwd("C:/Matias/Analyses/Reference Points/Outputs/p_0.85")


write.table(do.call(rbind,TAB.Trigger),"Trigged.times.csv",sep=",")


#1. Natural mortality
# Function for filling in ages to max age for all iterations
add.missing.age=function(Data,A)
{
  test=vector('list',length=length(Data))
  for(i in 1:length(Data))
  {
    datos=Data[[i]]
    if(length(datos)<A)
    {
      extra=A-length(datos)
      datos=c(datos,rep(NA,extra))
    }
    test[[i]]=datos
  }
  dummy=do.call(rbind,test)
  return(dummy)
}

tiff(file="Figure1.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
par(mfcol=c(4,2),las=1,mai=c(.6,.6,.075,.1),omi=c(.1,.1,.1,.05),mgp=c(1.5,0.95,0))
for (j in 1:N.sp)
{
  mSim=add.missing.age(Species.M[[j]],(Max.max.age[j]+1))
  boxplot(mSim,col="gray60",ylab="",xaxt="n",cex.axis=1.5,cex.lab=1.5,outline=F,
          border="gray10",varwidth=T,xlim=c(0,Max.max.age[j]))
          
  
  axis(1,at=seq(1,Max.max.age[j]+1,1),labels=F,tck=-0.02)
  axis(1,at=seq(1,Max.max.age[j],2),labels=seq(0,Max.max.age[j]-1,2),tck=-0.04,cex.axis=1.5)
#  legend("topright",species.names[j],cex=1.5,bty="n")

  ##Calculate median and export value for population dynamics
  colMedian=function(x) median(x,na.rm=T)
  Median.M=apply(mSim,2,colMedian)
  Median.M=Median.M[!is.na(Median.M)]
  write.csv(cbind(Median.M,Age=seq(0,Max.max.age[j]-1,1)),
      paste(handl_OneDrive("Data/Population dynamics/Parameter inputs for models/"),
        names(Species.M)[j],".M_at_age.csv",sep=""),row.names=F)
}
mtext(expression("     Natural mortality  " (year^-1)),side=2,line=-1.25,font=1,las=0,cex=1.2,outer=T)
mtext("Age",side=1,line=-2,font=1,las=0,cex=1.2,outer=T)
#dev.off()


#2. Selectivity

#2.1 Box plots

#tiff(file="Figure2.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
#par(mfcol=c(3,1),las=1,mai=c(.6,.6,.075,.15),omi=c(.1,.1,.1,.05),mgp=c(1.5,0.95,0))
for (j in 1:N.sp)
{
  SelSim=add.missing.age(Species.Sel[[j]],(Max.max.age[j]+1))
  boxplot(SelSim,col="gray60",ylab="",xaxt="n",cex.axis=1.5,cex.lab=1.5,
          varwidth=T,outline=F,border="gray10",xlim=c(0,Max.max.age[j]))
  if(j==3)mtext("                          Relative selectivity",side=2,line=3,font=1,las=0,cex=1.2,outer=F)
  axis(1,at=seq(1,Max.max.age[j]+1,1),labels=F,tck=-0.02)
  axis(1,at=seq(1,Max.max.age[j],2),labels=seq(0,Max.max.age[j]-1,2),tck=-0.04,cex.axis=1.5)
  if(!j==2)legend("topright",SP.names[j],cex=1.75,bty="n")
  if(j==2)legend("bottomright",SP.names[j],cex=1.75,bty="n")
  
  
  ##Calculate median and export value for population dynamics
  Median.Sel=apply(SelSim,2,colMedian)
  Median.Sel=Median.Sel[!is.na(Median.Sel)]
  write.csv(cbind(Median.Sel,Age=seq(0,Max.max.age[j]-1,1)),
            paste(handl_OneDrive("Data/Population dynamics/Parameter inputs for models/"),
                  names(Species.M)[j],".Selectivity_at_age.csv",sep=""),row.names=F)
  
  
}

mtext("Age",side=1,line=-2,font=1,las=0,cex=1.2,outer=T)
dev.off()




  #2.2 Mean behaviour
ext=seq(130,200,10)
fn.plot.sel=function(species,age,Lo,Linf,K,alphabeta,alpha,beta)
{
  
  mid.FL.fem=Lo+(Linf-Lo)*(1-exp(-K*age))  #use modified version                
  sel.fem=((mid.FL.fem*10/alphabeta)^alpha)*(exp(alpha-(mid.FL.fem*10/beta)))
  # plot(age,mid.FL.fem,pch=19,cex=1.5)
  if(!(species =="sandbar shark"))
    {
      plot(mid.FL.fem,sel.fem,xlim=c(Lo,200), yaxs="i",xaxs="i",type='l',lwd=2,xlab="",
           ylab="",cex.axis=1.5,ylim=c(0,1.1))
    }
  if(species =="sandbar shark")
    {
      plot(mid.FL.fem,sel.fem,xlim=c(Lo,200), yaxs="i",xaxs="i",type='l',lwd=2,xlab="Fork length (cm)",
                       ylab="",cex.lab=2,cex.axis=1.5,ylim=c(0,1.1))
    }
  if (species =="whiskery shark")
  {
    sel.fem.ext=((c(mid.FL.fem,ext)*10/alphabeta)^alpha)*(exp(alpha-(c(mid.FL.fem,ext)*10/beta)))
    lines(c(mid.FL.fem,ext),sel.fem.ext,col=1,lwd=2)
  }
  if(!(species =="sandbar shark"))
  {    
      plot(age,sel.fem,type='l', yaxs="i",xaxs="i",lwd=2,xlab="",ylab="",cex.axis=1.5,ylim=c(0,1.1))
  }

  if(species =="sandbar shark")
  {
    plot(age,sel.fem,type='l', yaxs="i",xaxs="i",lwd=2,xlab="Age",ylab="",cex.lab=2,cex.axis=1.5,ylim=c(0,1.1))
  }
  legend("topright",species,cex=1.5,bty="n")
}
par(mfrow=c(4,2),las=1,mai=c(.7,.6,.075,.15),omi=c(.1,.3,.1,.05))
for (a in 1:N.sp)
{  
  age=seq(0,max(species.list[[a]]$max.age),by=.1)
  GR=Growth.pars[a,];ab=Sel.pars[a,1];al=Sel.pars[a,2];be=Sel.pars[a,3]
  fn.plot.sel(species.names[a],age,Size.birth[a],GR[1,1],GR[1,2],ab,al,be)
}
mtext("Relative selectivity",side=2,line=0,font=1,las=0,cex=1.3,outer =T)  


#Check fecundity
par(mfcol=c(4,1),las=1,mai=c(.6,.6,.075,.1),omi=c(.1,.1,.1,.05),mgp=c(1.5,0.95,0))
for (j in 1:N.sp)
{
  FecSim=add.missing.age(Species.Fec[[j]],(Max.max.age[j]+1))
  boxplot(FecSim,col="gray60",ylab="",xaxt="n",cex.axis=1.5,cex.lab=1.5,outline=F,
          border="gray10",varwidth=T,xlim=c(0,Max.max.age[j]))
  axis(1,at=seq(1,Max.max.age[j]+1,1),labels=F,tck=-0.02)
  axis(1,at=seq(1,Max.max.age[j],2),labels=seq(0,Max.max.age[j]-1,2),tck=-0.04,cex.axis=1.5)
  legend("topleft",SP.names[j],cex=1.75,bty="n")
  
}
mtext("Fecundity",side=2,line=-1.25,font=1,las=0,cex=1.2,outer=T)
mtext("Age",side=1,line=-2,font=1,las=0,cex=1.2,outer=T)



#3. SPR reference points

 #3.1. B reference points  
tiff(file="Figure2.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
par(mfcol=c(4,1),las=1,mai=c(.6,.5,.075,.15),omi=c(.1,.1,.1,.05),mgp=c(1.5,0.6,0))
for (i in 1:N.sp)
{
  one=density(BMER[[i]][,2],adjust=2,na.rm =T,from=0,to=1)
  two=density(BMER[[i]][,1],adjust=2,na.rm =T,from=0,to=1)
  three=density(BMER[[i]][,3], adjust=2,na.rm =T,from=0,to=1)
  plot(one, type="l",lty=1,col=1,xlab="",lwd=2.5, yaxs="i",xaxs="i",yaxt='n',
  ylab="",main="",xlim=c(0,.7),ylim=c(0, max(one$y,two$y,three$y)*1.05),cex.axis=1.75,cex.lab=1.65)
  lines(two, lty=1,col="gray60",lwd=2.5)
  lines(three, lty=5,col="gray40",lwd=2.5)
  legend("topleft",SP.names[i],cex=1.75,bty="n")
  if(i==1)legend("topright",c("Limit","Threshold","Target"),lty=c(1,1,5),
                 cex=1.75,bty="n",col=c(1,"gray60","gray40"),lwd=2)
}
mtext("Density",side=2,line=-1.5,font=1,las=0,cex=1.5,outer =T)
mtext("Relative biomass",side=1,line=2.1,font=1,las=0,cex=1.5)
dev.off()




 #3.2. F reference points
tiff(file="Figure3.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
par(mfcol=c(4,1),las=1,mai=c(.6,.5,.075,.15),omi=c(.1,.1,.1,.05),mgp=c(1.5,0.6,0))
XLim=c(.4,.4,.3,.1)
for (i in 1:N.sp)
{
  one=density(FMER[[i]][,2],adjust=2,na.rm =T,from=0,to=1)
  two=density(FMER[[i]][,1],adjust=2,na.rm =T,from=0,to=1)
  three=density(FMER[[i]][,3], adjust=2,na.rm =T,from=0,to=1)
  plot(one, type="l",lty=1,col=1,xlab="",lwd=2.5, yaxs="i",xaxs="i",yaxt='n',
       ylab="",main="",xlim=c(0,XLim[i]),ylim=c(0, max(one$y,two$y,three$y)*1.05),cex.axis=1.75,cex.lab=1.65)
  lines(two, lty=1,col="gray60",lwd=2.5)
  lines(three, lty=5,col="gray40",lwd=2.5)
  legend("topright",SP.names[i],cex=1.75,bty="n")
  if(i==1)legend("topleft",c("Limit","Threshold","Target"),lty=c(1,1,5),
                 cex=1.75,bty="n",col=c(1,"gray60","gray40"),lwd=2)
}
mtext("Density",side=2,line=-1.5,font=1,las=0,cex=1.5,outer =T)
mtext(expression("Fishing mortality  " (year^-1)),side=1,line=2.6,font=1,las=0,cex=1.5)
dev.off()



#Extract median and 95% CI
BLIM.list=BTHRE.list=BTARG.list=FMER.list=FLIM.list=FTAR.fin.list=FZHOU.list=
  F40.list=F30.list=SprMer.list=Steep.list=Steep.list.percentile=vector('list',length=N.sp)
for (i in 1:N.sp)
{

    this=BMER[[i]][,2]
    BLIM.list[[i]]=round(c(median=median(this),sixty.per=quantile(this,probs=0.6),
                           CI95=sort(this)[c(floor(0.025*length(this)),ceiling(0.975*length(this)))]),3)
    this=BMER[[i]][,1]  
    BTHRE.list[[i]]=round(c(median=median(this),sixty.per=quantile(this,probs=0.6),
                            CI95=sort(this)[c(floor(0.025*length(this)),ceiling(0.975*length(this)))]),3)
    this=BMER[[i]][,3] 
    BTARG.list[[i]]=round(c(median=median(this),sixty.per=quantile(this,probs=0.6),
                            CI95=sort(this)[c(floor(0.025*length(this)),ceiling(0.975*length(this)))]),3)
   
  this=FMER[[i]][,1]
  FMER.list[[i]]=round(c(median=median(this),sixty.per=quantile(this,probs=0.6),
                         CI95=sort(this)[c(floor(0.025*length(this)),ceiling(0.975*length(this)))]),3)
  
  this=FMER[[i]][,2]
  FLIM.list[[i]]=round(c(median=median(this),sixty.per=quantile(this,probs=0.6),
                         CI95=sort(this)[c(floor(0.025*length(this)),ceiling(0.975*length(this)))]),3)
  
  this=FMER[[i]][,3]
  FTAR.fin.list[[i]]=round(c(median=median(this),sixty.per=quantile(this,probs=0.6),
                             CI95=sort(this)[c(floor(0.025*length(this)),ceiling(0.975*length(this)))]),3)
  
#   this=FZHOU[[i]][,1]
#   FZHOU.list[[i]]=round(c(median=median(this),sixty.per=quantile(this,probs=0.6),
#                           CI95=sort(this)[c(floor(0.025*length(this)),ceiling(0.975*length(this)))]),3)
  
#   this=F40[[i]][,1]
#   F40.list[[i]]=round(c(median=median(this),sixty.per=quantile(this,probs=0.6),
#                         CI95=sort(this)[c(floor(0.025*length(this)),ceiling(0.975*length(this)))]),3)
  
#   this=F30[[i]][,1]
#   F30.list[[i]]=round(c(median=median(this),sixty.per=quantile(this,probs=0.6),
#                         CI95=sort(this)[c(floor(0.025*length(this)),ceiling(0.975*length(this)))]),3)


    this=Steep[[i]][,1]
    Steep.list[[i]]=round(c(median=median(this),sixty.per=quantile(this,probs=0.6),
                            CI95=sort(this)[c(floor(0.025*length(this)),ceiling(0.975*length(this)))]),3)
    Steep.list.percentile[[i]]=quantile(this,probs=seq(0,1,0.05))

    #export steepness for populatin dynamics
    Steep.pop.dyn=round(c(mean=mean(this),sd=sd(this)),3)
    write.csv(data.frame(t(Steep.pop.dyn)),paste(handl_OneDrive("Data/Population dynamics/Parameter inputs for models/"),
                names(Species.M)[i],".Steepness.csv",sep=""),row.names=F)


    this=SprMer[[i]][,1]
    SprMer.list[[i]]=round(c(median=median(this),sixty.per=quantile(this,probs=0.6),
                            CI95=sort(this)[c(floor(0.025*length(this)),ceiling(0.975*length(this)))]),3)
    


}
write.csv(do.call(rbind,BLIM.list),"SPR.B.lim.csv",row.names=species.names)
write.csv(do.call(rbind,BTHRE.list),"SPR.B.thre.csv",row.names=species.names)
write.csv(do.call(rbind,BTARG.list),"SPR.B.tar.csv",row.names=species.names)
write.csv(do.call(rbind,FMER.list),"SPR.F.mer.csv",row.names=species.names)
write.csv(do.call(rbind,FLIM.list),"SPR.F.lim.csv",row.names=species.names)
write.csv(do.call(rbind,FTAR.fin.list),"SPR.F.tar.fin.csv",row.names=species.names)
#write.csv(do.call(rbind,F40.list),"SPR.F.40.csv",row.names=species.names)
#write.csv(do.call(rbind,F30.list),"SPR.F.30.csv",row.names=species.names)
#write.csv(do.call(rbind,FZHOU.list),"SPR.F.zhou.csv",row.names=species.names)

write.csv(do.call(rbind,SprMer.list),"SPR.mer.csv",row.names=species.names)
write.csv(do.call(rbind,Steep.list),"h.csv",row.names=species.names)
write.csv(do.call(rbind,Steep.list.percentile),"h_percentiles.csv",row.names=species.names)





#4. Distributions overlap

  #4.1 Calculate densities overlaps and probabilities of being over BRP
Dens.over=function(LIM,THR,TAR,Estim,WHAT,Ley)
{
  #Calculate original densities
  Dens.lim=density(LIM,adjust=3,na.rm =T)
  Dens.th=density(THR,adjust=3,na.rm =T)
  Dens.tar=density(TAR,adjust=3,na.rm =T)
  Estimate=density(Estim,adjust=3,na.rm =T)
  
  #Calculate prob larger than median
  Lim.med=median(LIM)
  Thr.med=median(THR)
  Tar.med=quantile(TAR,prob=0.2)
  fn.Prob.larger<<-function(DAT,Point.Est.BRP)
  {
    Fn <- ecdf(DAT)
    Prob=1-Fn(Point.Est.BRP)
    return(Prob)
  }
  #B.Lim.prob=round(fn.Prob.larger(Estim,Lim.med),2) #median
  #B.Thr.prob=round(fn.Prob.larger(Estim,Thr.med),2)
  #B.Tar.prob=round(fn.Prob.larger(Estim,Tar.med),2)#20% percentile
  
  #Calculate prob smaller than median
  fn.Prob.smaller<<-function(DAT,Point.Est.BRP)
  {
    Fn <- ecdf(DAT)
    Prob=Fn(Point.Est.BRP)
    return(Prob)
  }
  B.Lim.prob=round(fn.Prob.smaller(Estim,Lim.med),2) #median
  B.Thr.prob=round(fn.Prob.smaller(Estim,Thr.med),2)
  B.Tar.prob=round(fn.Prob.smaller(Estim,Tar.med),2)#20% percentile
  
  
  
  #Redo pairwise densities within same range
  RANGO=range(Dens.lim$y, Dens.th$y,Dens.tar$y,Estimate$y)
  fn.redo=function(A,B,C,D,YRANGE,LEG,LEG1,relative.to,whatprob,prob.larger,MED)
  {
    #plot
    LWD=1.75
    lim <- range(A$x, B$x)
    BRP=density(C,adjust=3,na.rm =T,from=lim[1], to=lim[2])
    estimate=density(D,adjust=3,na.rm =T,from=lim[1], to=lim[2])  
    plot(BRP, col=1, ylim=YRANGE,main="",lwd=LWD,ylab="",xlab="", yaxs="i",xaxs="i",cex.axis=1.5,xlim=c(0,1))
    COL=rgb(70, 70, 70, alpha=95, maxColorValue=255)
    #polygon(BRP$x, pmin(BRP$y, estimate$y), col= COL) #overlap
    #lines(estimate, col="gray50",lwd=2.5)
    
    #plot reference point
    polygon(x=c(BRP$x, rev(BRP$x)), y=c(rep(0,length(BRP$x)), rev(BRP$y)), border=1, density=10,
            col=1,lwd=LWD)
    
    #plot biomass
    polygon(x=c(BRP$x, rev(BRP$x)), y=c(rep(0,length(estimate$x)), rev(estimate$y)), border="gray40", density=10,
            col="gray40", angle=135, lty=6,lwd=LWD)
    
    #plot alpha overlap
    IDs=match(round(MED,3),round(estimate$x,3))   #match MED in density
    lines(c(BRP$x[IDs],BRP$x[IDs]),c(0,BRP$y[IDs]), col=1,lwd=LWD)  #plot median
    polygon(x=c(estimate$x[1:IDs], rev(estimate$x[1:IDs])), y=c(rep(0,IDs), rev(estimate$y[1:IDs])), 
            border=COL, col=COL)
    
    #calculate overlapping area
    NN=length(BRP$x)
    BRP.pol=data.frame(PID=rep(1,(NN)),POS=1:(NN),X=c(BRP$x),Y=c(BRP$y))
    Estim.pol=data.frame(PID=rep(1,(NN)),POS=1:(NN),X=c(estimate$x),Y=c(estimate$y))
    Over.pol=data.frame(PID=rep(1,(NN)),POS=1:(NN),X=BRP$x, Y=pmin(BRP$y, estimate$y))
    area.Overlap=calcArea (Over.pol, rollup = 3)
    if(relative.to=="Estim")AREA=calcArea (Estim.pol, rollup = 3)
    if(relative.to=="Tar")AREA=calcArea (BRP.pol, rollup = 3)
    Per.over=area.Overlap*100/AREA
    Per.over=abs(round(Per.over$area))
    
    #  polygon(Over.pol$X, Over.pol$Y, col=rgb(1,0,1,0.2))
    
    greek <- c(LEG,LEG1,"tau","alpha")
    other_stuff <- c('','',paste('=',Per.over,'%'),paste('=',prob.larger))
    .expressions <- mapply(sprintf, greek, other_stuff,MoreArgs = list(fmt = '%s~"%s"'))        
    legend_expressions <-parse(text = .expressions)
    
    
    legend("topright",legend=legend_expressions,
           lty=c(1,6,0,0),cex=1.25,bty="n",col=c(1,"gray40"),
           fill=c("transparent","transparent", "transparent","transparent"),lwd=LWD,merge = TRUE,
           border=c("transparent","transparent", "transparent","transparent"))
    #           fill=c("transparent","transparent", COL,"transparent"),lwd=2,merge = TRUE,
    #           border=c("transparent","transparent", 1,"transparent"))
    
    #    abline(0,0,col="white",lwd=2)
    
    #    return(Overlap=Per.over)
  }
  
  Lim.Est=fn.redo(Dens.lim,Estimate,LIM,Estim,RANGO,"Limit",WHAT,"Estim","Limit",B.Lim.prob,Lim.med)     #Limit Vs Estimate
  legend("topleft",Ley,bty="n",cex=1.5)
  Thr.Est=fn.redo(A=Dens.th,B=Estimate,C=THR,D=Estim,YRANGE=RANGO,LEG="Threshold",LEG1=WHAT,
                  relative.to="Estim",whatprob="Threshold",prob.larger=B.Thr.prob,MED=Thr.med)     #Threshold Vs Estimate
  Tar.Est=fn.redo(Dens.tar,Estimate,TAR,Estim,RANGO,"Target",WHAT,"Tar","Target",B.Tar.prob,Tar.med)     #Target Vs Estimate
  
  return(list(Lim.Est=Lim.Est,Thr.Est=Thr.Est,Tar.Est=Tar.Est))
}


tiff(file="Figure4.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff

par(mfcol=c(3,2),las=1,mai=c(.6,.5,.075,.15),omi=c(.1,.1,.1,.05),mgp=c(1.5,0.6,0))
#happy place
Dummy=rnorm(10000, mean=0.55, sd=0.03)
B.overlap.nice=Dens.over(LIM=BMER[[3]][,2],THR=BMER[[3]][,1],TAR=BMER[[3]][,3],Estim=Dummy,WHAT="Bc",Ley="A")

#nasty place
Dummy2=rnorm(10000, mean=0.30, sd=0.03)
B.overlap.nasty=Dens.over(LIM=BMER[[3]][,2],THR=BMER[[3]][,1],TAR=BMER[[3]][,3],Estim=Dummy2,WHAT="Bc",Ley="B")

mtext("Relative biomass",1,line=-1,outer=T,cex=1.5)
mtext("Density",2,line=-1,outer=T,las=3,cex=1.5)
dev.off()
  
  #4.2 Calculate probabilty of being above or below point estimate
B.Lim.prob.above=round(fn.Prob.larger(Dummy,median(BMER[[3]][,2])),2)
B.Thr.prob.above=round(fn.Prob.larger(Dummy,median(BMER[[3]][,1])),2)
B.Tar.prob.above=round(fn.Prob.larger(Dummy,quantile(BMER[[3]][,3],prob=0.2)),2)

B.Lim.prob.below=round(fn.Prob.smaller(Dummy,median(BMER[[3]][,2])),2)
B.Thr.prob.below=round(fn.Prob.smaller(Dummy,median(BMER[[3]][,1])),2)
B.Tar.prob.below=round(fn.Prob.smaller(Dummy,quantile(BMER[[3]][,3],prob=0.2)),2)

Fig4.prob.above=cbind(B.Lim.prob.above,B.Thr.prob.above,B.Tar.prob.above)
Fig4.prob.below=cbind(B.Lim.prob.below,B.Thr.prob.below,B.Tar.prob.below)
  
write.csv(Fig4.prob.above,"Fig4.prob.above.csv")
write.csv(Fig4.prob.below,"Fig4.prob.below.csv")



#5. SENSITIVITY TEST
#Manually change one parameter at a time to see effect. Don't change gummy, they are ok

#Life history par vectors (Gummy,Whiskery,Dusky,Sandbar)
#Min values
Min.max.age=c(16,15,40,30)  #NEW (less variable for dusky)
#Min.max.age=c(16,15,40,30)
Min.breed.freq=c(1,0.5,0.5,0.5)
Min.age.mat=c(4,6,26,13)  #NEW (less variable for dusky, sandbar)
#Min.age.mat=c(4,6,26,13)
Min.fec=c(1,14,2,10)      #NEW (less variable for whiskery, dusky, sandbar)
#Min.fec=c(1,4,2,4)

#Max values
Max.max.age=c(21,20,56,40) #put +1 because using floor in random sample
Max.breed.freq=c(1,0.5,0.333,0.5)  #NEW (less variable for dusky)
#Max.breed.freq=c(1,0.5,0.333,0.5)
Max.age.mat=c(6,8,36,13) #put +1 because using floor in random sample
Max.fec=c(31,28,18,10)   


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
Temperature=c(18,18,18,24)

#Growth pars (female)
Linf.g=201.9
#Linf.w=120.7 #Simfendorfer et al 2000, age and growth
Linf.w=160   #Last & Stevens max size                 NEW
Linf.d=374.4
Linf.s=244.2

K.g=0.123
#K.w=0.288
K.w=0.369
K.d=.0367
K.s=.040

to.g=-1.55
to.w=-0.6
to.d=-3.3
to.s=-4.8

#Size at birth
Size.birth=c(33.5,25,75.3,42.5)

#FL to TL pars             
b.g=NA
b.w=8.891
b.d=4.000
b.s=5.8188
a.g=NA
a.w=1.046
a.d=1.100
a.s=1.11262


#TL to TW pars 
bwt.g=0.927e-6  #(reported as 0.927e-9 but I changed scale to match others)
bwt.w=0.0000163
bwt.d=1.2334e-5
bwt.s=0.000002
awt.g=3.206
awt.w=2.733
awt.d=2.855
awt.s=3.2069

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

#dusky                #Simpfendorfer & Unsworth 1998
theta1=130.13
theta2=29237
mesh= 6.5 #(in inches, = 16.5 cm)
alphabeta.d=theta1*mesh
beta.d=-0.5*(theta1*mesh-((theta1^2)*(mesh^2)+4*theta2)^0.5)
alpha.d=alphabeta.d/beta.d

#sandbar                #McAuley et al 2005
theta1=135.58
theta2=117001
mesh= 6.5
alphabeta.s=theta1*mesh
beta.s=-0.5*(theta1*mesh-((theta1^2)*(mesh^2)+4*theta2)^0.5)
alpha.s=alphabeta.s/beta.s


#Put pars in data frame for loop
Growth.pars=data.frame(Linf=c(Linf.g,Linf.w,Linf.d,Linf.s),K=c(K.g,K.w,K.d,K.s),to=c(to.g,to.w,to.d,to.s))
Length.conv=data.frame(b=c(b.g,b.w,b.d,b.s),a=c(a.g,a.w,a.d,a.s))
Weight.conv=data.frame(bwt=c(bwt.g,bwt.w,bwt.d,bwt.s),awt=c(awt.g,awt.w,awt.d,awt.s))
Sel.pars=data.frame(alphabeta=c(alphabeta.g,alphabeta.w,alphabeta.d,alphabeta.s),
                    alpha=c(alpha.g,alpha.w,alpha.d,alpha.s),
                    beta=c(beta.g,beta.w,beta.d,beta.s))


#---PROCEDURE SECTION----

#1. Priors
#1.1. Natural mortality
#source external function 
Wght.G=c(1,1,1,1)    # weight given to methods using growth parameters (high uncertainty for whiskery so give lower weight)
Wght.noG=c(1,10,1,1)  # weight given to methods not using growth parameters




#---MAIN SECTION----
N.sim=100

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
    f.zhou=vector("list", length=N.sim)
  #names(f.out)=c("threshold","limit","target")
  Trigged=rep(NA,N.sim)
  for (j in 1:N.sim)
  {
    #1. draw random samples of input parameters
    A.MAX=A.max(AA[1],AA[2],AA[1])
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
    
    
    mm=M.fun(A.MAX,GR[2],GR[1],Temper,A.MAT,TO,WT[1,1],WT[1,2],w.g,w.ng)
    Sel.A=Select(alphabeta.S,alpha.S,beta.S)
    
    
    #2. calculate SPR.mer quantities 
    spr.temp=SPR(A.MAX,mm,Fec,Rep,sex.ratio,A.MAT,Sel.A,0,"Y","ogive")
    
    if(length(spr.temp)>1)
    {
      #3. Calculate reference points
      
      #3.1. Threshold (i.e. MER)    
      
      #3.1.1 Biomass (Bmer)
      BTHR=spr.temp$Dep.MER
      
      #3.1.2. F (Fmer)
      #Set up useful vectors
      fecundity=spr.temp$fecundity
      maturity=spr.temp$maturity
      lx=vector(length=length(age))    #no fishing
      lz=vector(length=length(age))    #with fishing
      Observed=spr.temp$SPR.mer
      phieggs0=spr.temp$phi.o
      fn.yield=function(Fe) EqPerRec(Fe=Fe)$epsilon
      
      #numerical approximation
      Fstart=0.1
      Fmer.fit <- nlminb( start=Fstart, objective=fn.yield, lower=0.0001, upper=5.0,
                          control=list(eval.max=500, iter.max=500  )   )
      FTHR=Fmer.fit$par
      FTHR.Conv=Fmer.fit$convergence
      
      
      #3.2. Limits 
      
      #3.2.1. Biomass (Blim = p x BTHR)
      B.lim=Prop.FMsy*BTHR
      
      #3.2.2 F
      obs.biom <- B.lim   
      b.f <- function(Fstart)
      {
        temp.spr = EqPerRec(Fstart)$SPR    
        temp.SSB =Find.B(spr.temp$alpha,R0=1,spr=temp.spr,spr0=1)  #R0 and spr0 =1 as it's on relative scale
        yy<-abs(temp.SSB - obs.biom )
        return(yy)
      }
      Flim.fit = nlminb(start=Fstart, objective=b.f, lower=0, upper=3)
      FLim = Flim.fit$par
      Flim.Conv = Flim.fit$convergence
      
      
      #3.3. Targets
      
      #3.3.1. Biomass (Blim = p x BTHR)
      B.tar=(1/Prop.FMsy)*BTHR
      
      #3.3.2 F
      obs.biom <- B.tar   
      FTar.fit = nlminb(start=Fstart, objective=b.f, lower=0, upper=3)
      FTar = FTar.fit$par
      FTar.Conv = FTar.fit$convergence
      
      
    }
    
    
    #3. Repeat if nonsense outputs (e.g. h<0.2, negative Biomass RP), i.e. obtain sensible joint priors 
    Trigger=0
    if (length(spr.temp)==1|B.tar<0|B.lim<0|BTHR<0)Trigger=1
    Trigged[j]=Trigger
    if(Trigger==1)    
    { repeat 
    {
      #3.1 draw random samples
      A.MAX=A.max(AA[1],AA[2],AA[1])
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
      
      
      mm=M.fun(A.MAX,GR[2],GR[1],Temper,A.MAT,TO,WT[1,1],WT[1,2],w.g,w.ng)
      Sel.A=Select(alphabeta.S,alpha.S,beta.S)
      
      
      #3.3 calculate SPR.mer quantities
      spr.temp=SPR(A.MAX,mm,Fec,Rep,sex.ratio,A.MAT,Sel.A,0,"Y","ogive")
      
      if(length(spr.temp)>1)
      {
        #3.3.1 Calculate reference points
        
        #3.3.1.1 Threshold (i.e. MER)    
        
        #Biomass (Bmer)
        BTHR=spr.temp$Dep.MER
        
        #F (Fmer)
        #Set up useful vectors
        fecundity=spr.temp$fecundity
        maturity=spr.temp$maturity
        lx=vector(length=length(age))    #no fishing
        lz=vector(length=length(age))    #with fishing
        Observed=spr.temp$SPR.mer
        phieggs0=spr.temp$phi.o
        fn.yield=function(Fe) EqPerRec(Fe=Fe)$epsilon
        
        #numerical approximation
        Fmer.fit <- nlminb( start=0.1 , objective=fn.yield, lower=0.0001, upper=5.0,
                            control=list(eval.max=500, iter.max=500  )   )
        FTHR=Fmer.fit$par
        FTHR.Conv=Fmer.fit$convergence
        
        
        #3.2. Limits 
        
        #3.2.1. Biomass (Blim = p x BTHR)
        B.lim=Prop.FMsy*BTHR
        
        #3.2.2 F
        obs.biom <- B.lim   
        b.f <- function(Fstart)
        {
          temp.spr = EqPerRec(Fstart)$SPR    
          temp.SSB =Find.B(spr.temp$alpha,R0=1,spr=temp.spr,spr0=1)  #R0 and spr0 =1 as it's on relative scale
          yy<-abs(temp.SSB - obs.biom )
          return(yy)
        }
        Flim.fit = nlminb(start=Fstart, objective=b.f, lower=0, upper=3)
        FLim = Flim.fit$par
        Flim.Conv = Flim.fit$convergence
        
        
        #3.3. Targets
        
        #3.3.1. Biomass (Blim = p x BTHR)
        B.tar=(1/Prop.FMsy)*BTHR
        
        #3.3.2 F
        obs.biom <- B.tar   
        FTar.fit = nlminb(start=Fstart, objective=b.f, lower=0, upper=3)
        FTar = FTar.fit$par
        FTar.Conv = FTar.fit$convergence         
        
        
      }
      
      #break if sensible joint prior
      Trigger=0
      if (length(spr.temp)==1|B.tar<0|B.lim<0|BTHR<0)Trigger=1
      
      
      if(Trigger==0)break
    }
    }
    
    
    #store quantities of interest
    SPR.out[[j]]=spr.temp
    Store.M[[j]]=mm
    Store.Sel[[j]]=Sel.A
    Store.len[[j]]=mid.FL.fem
    Store.Fec[[j]]=Fec
    
    
    
    #4.5 Store all Bs and Fs
    b.out[[j]]=c(BTHR,B.lim,B.tar)
    f.out[[j]]=c(FTHR,FLim,FTar)
    
    #4.6 Store convergence code
    Convergence[[j]]=c(FTHR.Conv=FTHR.Conv)
    
    
    #5. Calculate F reference points for arbitrary SPR
    Observed=0.40 #target (Tsai et al 2011 page 1388)
    F40.fit=optimize(fn.yield,lower=0.0001,upper=5)   #numerical approximation
    f.out.40[[j]]=F40.fit$minimum
    
    Observed=0.30 #limit
    F30.fit=optimize(fn.yield,lower=0.0001,upper=5)   #numerical approximation
    f.out.30[[j]]=F30.fit$minimum
    
    #6. Calculate Fmsy proxy from Zhou et al 2012
    f.zhou[[j]]=Zhou(mean(mm))  
    
    
  }
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
FMER=BMER=vector("list", length=N.sp)
names(BMER) =names(FMER)= species.names
for (i in 1:N.sp)
{
  Fmer=data.frame(Fmer=rep(NA,N.sim),Flim=rep(NA,N.sim),Ftar=rep(NA,N.sim))
  Bmer=data.frame(Bmer=rep(NA,N.sim),Blim=rep(NA,N.sim),Btar=rep(NA,N.sim))
  for (j in 1:N.sim)
  {
    Bmer[j,]=Species.B.MER[[i]][[j]]
    Fmer[j,]=Species.F.MER[[i]][[j]]
  }
  
  BMER[[i]]=Bmer
  FMER[[i]]=Fmer
}


#5.1. F reference points
#tiff(file="Figure.F.sensit.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
par(mfcol=c(4,1),las=1,mai=c(.6,.5,.075,.15),omi=c(.1,.1,.1,.05),mgp=c(1.5,0.6,0))
XLim=c(.4,.4,.3,.06)
for (i in 1:N.sp)
{
  one=density(FMER[[i]][,2],adjust=2,na.rm =T,from=0,to=1)
  two=density(FMER[[i]][,1],adjust=2,na.rm =T,from=0,to=1)
  three=density(FMER[[i]][,3], adjust=2,na.rm =T,from=0,to=1)
  plot(one, type="l",lty=1,col=1,xlab="",lwd=2.5, yaxs="i",xaxs="i",yaxt='n',
       ylab="",main="",xlim=c(0,XLim[i]),ylim=c(0, max(one$y,two$y,three$y)*1.05),cex.axis=1.75,cex.lab=1.65)
  lines(two, lty=1,col="gray60",lwd=2.5)
  lines(three, lty=5,col="gray40",lwd=2.5)
  legend("topright",SP.names[i],cex=1.75,bty="n")
  if(i==1)legend("topleft",c("Limit","Threshold","Target"),lty=c(1,1,5),
                 cex=1.75,bty="n",col=c(1,"gray60","gray40"),lwd=2)
}
mtext("Density",side=2,line=-1.5,font=1,las=0,cex=1.5,outer =T)
mtext(expression("Fishing mortality  " (year^-1)),side=1,line=2.6,font=1,las=0,cex=1.5)
#dev.off()





#################################################LIZ BROOK'S###################
#2.2 Spawners per recruit (Liz Brooks)
s.per.recr<-function(nages,fec.age,mat.age,M.age, F.mult, sel.age, spawn.time )
{
  
  spr=0.0
  cum.survive=1.0
  z=0.0
  for (i in 1:(nages-1)  ) {
    z=M.age[i] + F.mult*sel.age[i]
    z.ts=(M.age[i]+F.mult*sel.age[i])*spawn.time
    spr=spr+cum.survive*fec.age[i]*mat.age[i]*exp(-z.ts)
    cum.survive=cum.survive*exp(-z )
    
  }
  
  z= M.age[nages] + F.mult*sel.age[nages]
  z.ts=(M.age[nages]+F.mult*sel.age[nages])*spawn.time
  spr=spr + fec.age[nages]*mat.age[nages]*cum.survive*exp(-z.ts)/( 1- exp(-z ) )
  
  return(spr)
  
}

#2.3 Yield per recruit  (Liz Brooks)
ypr<-function(nages, wgt.age, M.age, F.mult, sel.age )
{
  
  yield=0.0
  cum.survive=1.0
  z=0.0
  
  for (i in 1:(nages-1)  ) {
    z=M.age[i] + F.mult*sel.age[i]
    yield=yield + wgt.age[i]*F.mult*sel.age[i]*(1-exp(-z) )*cum.survive/z
    cum.survive=cum.survive*exp(-z)
  }
  
  z= M.age[nages] + F.mult*sel.age[nages]
  yield=yield + wgt.age[nages]*F.mult*sel.age[nages]*cum.survive/z
  
  return(yield)
  
}

#2.4 Equilibrium SSB (includes Ricker option) (Liz Brooks)
ssb.eq<-function(recr.par, R0, spr, spr0, is.steepness=T, SRtype )
{
  #SRtype=1 is BH
  #SRtype=2 is Ricker
  
  sprF=spr/spr0
  
  if (SRtype==1) {
    if (is.steepness==T)  alpha.hat <- 4*recr.par/(1-recr.par)
    if (is.steepness==F)  alpha.hat <- recr.par
    
    ssb=spr0*R0*(sprF*alpha.hat - 1.0)/(alpha.hat-1.0)
  }
  
  return(ssb)
  
}



############################################





########################################## EXAMPLE USED TO SET UP SCRIPT ############################
#1.1 life history
#Dusky (Brooks et al 2010)
#ages
age.min=0
age.max=40
age=age.min:age.max
r=1     #age at recruitment
first=match(r,age)
last=match(age.max,age)

Numb.pups=7.1     #mean number of pups
cycle.length=3    #reproductive cycle length
sex.ratio=0.5     #embryo sex ratio
FEc= Numb.pups*sex.ratio/cycle.length   #annual female fecundity


#fecundity at age (number of female pups per year)
fec=c(rep(0,21),rep(FEc,length(age)-21))                                                        #NEW

#proportion mature at age
mu=c(rep(0,11),.001,.002,.006,.013,.03,.062,.124,.226,.37,.535,.687,.803,.881,.929,.958,.975,.985,
     .991,.994,.996,.998,.998,.999,.999,rep(1,6))

#survivorship
surv=c(.78,.865,.875,.883,.889,.894,.899,.903,.907,.91,.913,.915,.918,.919,.922,.923,.925,.926,.927,.928,
       .93,.93,.932,.933,.934,.934,.935,.935,.936,.937,.937,.938,.939,.939,.939,.94,.94,.94,.94,.941,
       .942)

m=-log(surv)

#compounded survival
surv.compounded=cumprod(c(1,surv[first:(last)]))                                              #NEW

#unexploited spawners per recruit
phi.o=sum(fec*mu*surv.compounded)                                                #NEW 

#maximum lifetime reproductive rate at low density
alpha=phi.o*surv[1]

#steepness
h=alpha/(4+alpha)

#Spawning potential ratio at maximum excess recruitment (MER) (Beverton-Holt relationship)
SPR.mer=1/alpha^0.5


#Reference point assessment

#1. Stock status
#Depletion at MER (Proportional reduction from unexploited level)
Reference.depletion.MER=((alpha^0.5)-1)/(alpha-1) #this is the optimal level of depletion (life history dependent!)

#Abundance index
Index.1974=2.197; Index.2003=0.24; Depletion.in.1974=0.8
Scaled.relative.index.2003=Depletion.in.1974*(Index.2003/Index.1974)


#Overfished threshold
#this is some proportion of Reference.depletion.MER
p=max(c(0.5,(1-max(m))))

#Stock status
if(Scaled.relative.index.2003<(p*Reference.depletion.MER)){Stock.status="Overfished"}else
{Stock.status="Not overfished"}      



#}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}

###############
### added may 2014 by liz
### comparing predicted analytic MER (which assumes maturity=1 from the first age and 
###     selectivity =1 for all ages) to numerically estimated MER
###  
wd<-("C:/liz_temp/papers/matias_sharks/liz_R/")
setwd(wd)
### dummy parameter values (recycled from above)
#Dusky (Brooks et al 2010)
#ages
age.min=0
age.max=40
age=age.min:age.max
nages <- length(age) 
r=1     #age at recruitment
first=match(r,age)
last=match(age.max,age)

Numb.pups=7.1     #mean number of pups
cycle.length=3    #reproductive cycle length
sex.ratio=0.5     #embryo sex ratio
FEc= Numb.pups*sex.ratio/cycle.length   #annual female fecundity
spawn.time <- 0 # value between 0 and 1 to decrement survival until timing  of spawning

#some fake selectivity parameters to make a dome
a50.up <- 1
slope.up <- 2
a50.down <- 4
slope.down <- 0.4
sel1 <- ( 1/(1+exp(-(age-a50.up)/slope.up) )  )*(1- 1/(1++exp(-(age-a50.down)/slope.down) ) )
sel <- sel1/max(sel1)

#weight at age (dummy values for a and b in l-w eqn)
Linf <- 374
K <- 0.0367
t0 <- -3.3
wgt.age <- 1.5e-6*(Linf*(1-exp(-K*(age[2:nages]-t0)) )  )^2.8

#fecundity at age (number of female pups per year)
fec=c(rep(0,21),rep(FEc,length(age)-21))                                                        #NEW

#proportion mature at age
mu=c(rep(0,11),.001,.002,.006,.013,.03,.062,.124,.226,.37,.535,.687,.803,.881,.929,.958,.975,.985,
     .991,.994,.996,.998,.998,.999,.999,rep(1,6))

#survivorship
surv=c(.78,.865,.875,.883,.889,.894,.899,.903,.907,.91,.913,.915,.918,.919,.922,.923,.925,.926,.927,.928,
       .93,.93,.932,.933,.934,.934,.935,.935,.936,.937,.937,.938,.939,.939,.939,.94,.94,.94,.94,.941,
       .942)

m=-log(surv)

#compounded survival
surv.compounded=cumprod(c(1,surv[first:(last)]))                                              #NEW

#unexploited spawners per recruit
phi.o=sum(fec*mu*surv.compounded)                                                #NEW 

#maximum lifetime reproductive rate at low density
alpha=phi.o*surv[1]

#steepness
steepness=alpha/(4+alpha)



#-------Spawners per recruit -----------------------------
s.per.recr<-function(nages,fec.age,mat.age,M.age, F.mult, sel.age, spawn.time ) {
  
  spr=0.0
  cum.survive=1.0
  z=0.0
  for (i in 1:(nages-1)  ) {
    z=M.age[i] + F.mult*sel.age[i]
    z.ts=(M.age[i]+F.mult*sel.age[i])*spawn.time
    spr=spr+cum.survive*fec.age[i]*mat.age[i]*exp(-z.ts)
    cum.survive=cum.survive*exp(-z )
    
  }
  
  z= M.age[nages] + F.mult*sel.age[nages]
  z.ts=(M.age[nages]+F.mult*sel.age[nages])*spawn.time
  spr=spr + fec.age[nages]*mat.age[nages]*cum.survive*exp(-z.ts)/( 1- exp(-z ) )
  
  return(spr)
  
}

#-------Yield per recruit -----------------------------
ypr<-function(nages, wgt.age, M.age, F.mult, sel.age ) {
  
  yield=0.0
  cum.survive=1.0
  z=0.0
  
  for (i in 1:(nages-1)  ) {
    z=M.age[i] + F.mult*sel.age[i]
    yield=yield + wgt.age[i]*F.mult*sel.age[i]*(1-exp(-z) )*cum.survive/z
    cum.survive=cum.survive*exp(-z)
  }
  
  z= M.age[nages] + F.mult*sel.age[nages]
  yield=yield + wgt.age[nages]*F.mult*sel.age[nages]*cum.survive/z
  
  return(yield)
  
}
#------------------------------------
#-------Equilibrium SSB ------------------------



ssb.eq<-function(recr.par, R0.BH, spr, spr0, is.steepness=T ) {
  
  if (is.steepness==T)  alpha.BH <- 4*recr.par/(1-recr.par)
  if (is.steepness==F)  alpha.BH <- recr.par
  sprF=spr/spr0
  
  ssb=spr0*R0.BH*(sprF*alpha.BH - 1.0)/(alpha.BH-1.0)
  
  
  return(ssb)
  
}
#------------------------------------  
bev.holt.alpha<-function(S,R0,recr.par,spr0, is.steepness=T){
  
  if (is.steepness==T)  alpha.BH <- 4*recr.par/(1-recr.par)
  if (is.steepness==F)  alpha.BH <- recr.par
  
  
  y=rep(0,length(S))
  y=R0*S*alpha.BH/(R0*spr0 +(alpha.BH-1.0)*S)
  return(y)
  
} 
##============================================================================================ 
##----  Calculate Analytical SPR.MER, B.MER, and numerically solve for F.MER ----  #

R0 <- 1   #put unfished recruitment at 1 (everything on relative scale)
sr0 <- 1  #unfished spawners per recruit at 1 (relative scale)
sr0.calc <- s.per.recr(nages=(nages-1), fec.age=fec[first:(last)], mat.age=rep(1, (nages-1)), M.age= m[first:(last)], F.mult=0, sel.age=rep(1, (nages-1)), spawn.time=spawn.time)  


#1. Spawning potential ratio at maximum excess recruitment (MER) (Beverton-Holt relationship)
SPR.MER=1/alpha^0.5


#2. Depletion at MER (Proportional reduction from unexploited level)
B.MER=((alpha^0.5)-1)/(alpha-1) #this is the optimal level of depletion (life history dependent!)

#3. Solve numerically for F.MER  
#(technically, should fix maturity=1 and selectivity=1 at all ages)
# also, yield is maximized in NUMBERS, not WEIGHT, so calculating YPR with a vector of 1 for the weight at age 
F.start=0.12
t.spr <- SPR.MER   # the SPR that we want to match by finding an F that produces it

spr.f <- function(F.start) {
  abs(s.per.recr(nages=(nages-1), fec.age=fec[first:(last)], mat.age=rep(1, (nages-1)), M.age= m[first:(last)], F.mult=F.start, sel.age=rep(1, (nages-1)), spawn.time=spawn.time)/sr0.calc - t.spr )
}
yyy <- nlminb(start=F.start, objective=spr.f, lower=0, upper=3)

F.MER <- yyy$par  #find F that matches SPR.MER



#4. Yield per recruit at MER (in numbers)
#    calculated with F that matches SPR
t.ypr<- ypr((nages-1), wgt.age=rep(1, (nages-1)), M.age=m[first:(last)],  F.mult=F.MER, sel.age=rep(1, (nages-1)) )   
t.ssb <-ssb.eq( recr.par=alpha, R0=R0, spr=t.spr, spr0=sr0, is.steepness=F) 
Y.MER <- t.ypr*t.ssb/t.spr



#5. Let B.Lim = p*B.MER
p <- 0.75
B.Lim <- p*B.MER


#6. calculate corresponding F.MER that produces B.Lim

t.b.ref <- B.Lim   # define a level of B/B0 that we want to match with an F
b.f <- function(F.start) {
  temp.spr = s.per.recr(nages=(nages-1), fec.age=fec[first:(last)], mat.age=rep(1, (nages-1)), M.age= m[first:(last)], F.mult=F.start, sel.age=rep(1, (nages-1)), spawn.time=spawn.time)/sr0.calc    
  temp.SSB = ssb.eq( recr.par=alpha, R0=R0, spr=temp.spr, spr0=sr0, is.steepness=F)
  yy<-abs(temp.SSB - t.b.ref )
  return(yy)
}
yyy <- nlminb(start=F.start, objective=b.f, lower=0, upper=3)

F.Lim <-  yyy$par

#7. Yield per recruit at B.Lim (in numbers)
t.ypr<- ypr((nages-1), wgt.age=rep(1, (nages-1)), M.age=m[first:(last)],  F.mult=F.Lim, sel.age=rep(1, (nages-1)) )   
SPR.Lim <- s.per.recr(nages=(nages-1), fec.age=fec[first:(last)], mat.age=rep(1, (nages-1)), M.age= m[first:(last)], F.mult=F.Lim, sel.age=rep(1, (nages-1)), spawn.time=spawn.time)/sr0.calc    

Y.Lim <- t.ypr*B.Lim/SPR.Lim


#8.  Let B.Tar = 1/p*B.MER
B.Tar <- 1/p*B.MER

#9. calculate F.Tar that produces B.Tar
t.b.ref <- B.Tar   # define a level of B/B0 that we want to match with an F
b.f <- function(F.start) {
  temp.spr = s.per.recr(nages=(nages-1), fec.age=fec[first:(last)], mat.age=rep(1, (nages-1)), M.age= m[first:(last)], F.mult=F.start, sel.age=rep(1, (nages-1)), spawn.time=spawn.time)/sr0.calc    
  temp.SSB = ssb.eq( recr.par=alpha, R0=R0, spr=temp.spr, spr0=sr0, is.steepness=F)
  yy<-abs(temp.SSB - t.b.ref )
  return(yy)
}
yyy <- nlminb(start=F.start, objective=b.f, lower=0, upper=3)

F.Tar <-  yyy$par


#10. Calculate YPR at B.Tar
t.ypr<- ypr((nages-1), wgt.age=rep(1, (nages-1)), M.age=m[first:(last)],  F.mult=F.Tar, sel.age=rep(1, (nages-1)) )   
SPR.Tar <- s.per.recr(nages=(nages-1), fec.age=fec[first:(last)], mat.age=rep(1, (nages-1)), M.age= m[first:(last)], F.mult=F.Tar, sel.age=rep(1, (nages-1)), spawn.time=spawn.time)/sr0.calc    

Y.Tar <- t.ypr*B.Tar/t.spr

MER.analytic <- matrix(NA, nrow=12, ncol=1)
rownames(MER.analytic) <- c("SPR.MER", "B.MER", "F.MER", "Y.MER",
                            "SPR.Lim", "B.Lim", "F.Lim", "Y.Lim",
                            "SPR.Tar", "B.Tar", "F.Tar", "Y.Tar")
MER.analytic[1] <-  SPR.MER
MER.analytic[2] <-  B.MER
MER.analytic[3] <-  F.MER
MER.analytic[4] <-  Y.MER
MER.analytic[5] <-  SPR.Lim
MER.analytic[6] <-  B.Lim
MER.analytic[7] <-  F.Lim
MER.analytic[8] <-  Y.Lim
MER.analytic[9] <-  SPR.Tar
MER.analytic[10] <-  B.Tar
MER.analytic[11] <-  F.Tar
MER.analytic[12] <-  Y.Tar



##============================================================================================ 
##----  Calculate all MER reference points NUMERICALLY ---------------------------------------


R0 <- 1   #put unfished recruitment at 1 (everything on relative scale)
sr0 <- 1  #unfished spawners per recruit at 1 (relative scale)
sr0.calc <- s.per.recr(nages=(nages-1), fec.age=fec[first:(last)], mat.age=mu[first:(last)], M.age= m[first:(last)], F.mult=0, sel.age=sel, spawn.time=spawn.time)  

# 1.Numeric F.MER
# Given observed maturity, selectivity, find F that maximizes yield
#  (still use weight at age = vector of 1s since maximizing yield in numbers)

F.start=0.12

get.yield.f.min <- function(F.start) {
  
  temp.spr = s.per.recr(nages=(nages-1), fec.age=fec[first:(last)], mat.age=mu[first:(last)],  M.age= m[first:(last)], F.mult=F.start, sel.age=sel[first:(last)], spawn.time=spawn.time)/sr0.calc   
  temp.ypr = ypr(nages=(nages-1), wgt.age=rep(1, (nages-1)), M.age=m[first:(last)],  F.mult=F.start, sel.age=sel[first:(last)] )
  temp.SSB = ssb.eq( recr.par=alpha, R0=R0, spr=temp.spr, spr0=sr0, is.steepness=F)
  yield = temp.ypr*temp.SSB/temp.spr           #harvest in weight
  yy=-1*yield
  return(yy)  
}   # end get.yield.f.min function
F.nlmin <- nlminb( start=F.start , objective=get.yield.f.min, lower=0.001, upper=3.0,
                   control=list(eval.max=500, iter.max=500  )   )


F.MER.num <- F.nlmin$par  #optimize yield in numbers

#2. Numeric SPR.MER
SPR.MER.num <- s.per.recr(nages=(nages-1), fec.age=fec[first:(last)], mat.age=mu[first:(last)],  M.age= m[first:(last)], F.mult=F.MER.num, sel.age=sel[first:(last)], spawn.time=spawn.time)/sr0.calc   

#3. Numeric B.MER
B.MER.num <- ssb.eq( recr.par=alpha, R0=R0, spr=SPR.MER.num, spr0=sr0, is.steepness=F)

#4. Numeric Y.MER
Y.MER.num <- (B.MER.num/SPR.MER.num)*ypr(nages=(nages-1), wgt.age=rep(1, (nages-1)), M.age=m[first:(last)],  F.mult=F.MER.num, sel.age=sel[first:(last)] )



#5. Next, define B.Lim as p*B.MER
B.Lim.num <- p*B.MER.num


#6. Find F.Lim.num that produces B.Lim.num
t.b.ref <- B.Lim.num   # define a level of B/B0 that we want to match with an F
b.f <- function(F.start) {
  temp.spr = s.per.recr(nages=(nages-1), fec.age=fec[first:(last)], mat.age=mu[first:(last)], M.age= m[first:(last)], F.mult=F.start, sel.age=sel[first:(last)], spawn.time=spawn.time)/sr0.calc    
  temp.SSB = ssb.eq( recr.par=alpha, R0=R0, spr=temp.spr, spr0=sr0, is.steepness=F)
  yy<-abs(temp.SSB - t.b.ref )
  return(yy)
}
yyy <- nlminb(start=F.start, objective=b.f, lower=0, upper=3)

F.Lim.num <-  yyy$par

#7. Yield per recruit at B.Lim (in numbers)
t.ypr<- ypr((nages-1), wgt.age=rep(1, (nages-1)), M.age=m[first:(last)],  F.mult=F.Lim.num, sel.age=sel[first:(last)] )   
SPR.Lim.num <- s.per.recr(nages=(nages-1), fec.age=fec[first:(last)], mat.age=mu[first:(last)], M.age= m[first:(last)], F.mult=F.Lim.num, sel.age=sel[first:(last)], spawn.time=spawn.time)/sr0.calc    

Y.Lim.num <- t.ypr*B.Lim.num/SPR.Lim.num


#8.  Let B.Tar = 1/p*B.MER
B.Tar.num <- 1/p*B.MER

#9. calculate F.Tar that produces B.Tar
t.b.ref <- B.Tar.num   # define a level of B/B0 that we want to match with an F
b.f <- function(F.start) {
  temp.spr = s.per.recr(nages=(nages-1), fec.age=fec[first:(last)], mat.age=mu[first:(last)], M.age= m[first:(last)], F.mult=F.start, sel.age=sel[first:(last)], spawn.time=spawn.time)/sr0.calc    
  temp.SSB = ssb.eq( recr.par=alpha, R0=R0, spr=temp.spr, spr0=sr0, is.steepness=F)
  yy<-abs(temp.SSB - t.b.ref )
  return(yy)
}
yyy <- nlminb(start=F.start, objective=b.f, lower=0, upper=3)

F.Tar.num <-  yyy$par


#10. Calculate YPR at B.Tar
t.ypr<- ypr((nages-1), wgt.age=rep(1, (nages-1)), M.age=m[first:(last)],  F.mult=F.Tar.num, sel.age=sel[first:(last)] )   
SPR.Tar.num <- s.per.recr(nages=(nages-1), fec.age=fec[first:(last)], mat.age=mu[first:(last)], M.age= m[first:(last)], F.mult=F.Tar.num, sel.age=sel[first:(last)], spawn.time=spawn.time)/sr0.calc    

Y.Tar.num <- t.ypr*B.Tar.num/SPR.Tar.num

MER.numeric <- matrix(NA, nrow=12, ncol=1)
rownames(MER.numeric) <- c("SPR.MER", "B.MER", "F.MER", "Y.MER",
                           "SPR.Lim", "B.Lim", "F.Lim", "Y.Lim",
                           "SPR.Tar", "B.Tar", "F.Tar", "Y.Tar")
MER.numeric[1] <-  SPR.MER.num
MER.numeric[2] <-  B.MER.num
MER.numeric[3] <-  F.MER.num
MER.numeric[4] <-  Y.MER.num
MER.numeric[5] <-  SPR.Lim.num
MER.numeric[6] <-  B.Lim.num
MER.numeric[7] <-  F.Lim.num
MER.numeric[8] <-  Y.Lim.num
MER.numeric[9] <-  SPR.Tar.num
MER.numeric[10] <-  B.Tar.num
MER.numeric[11] <-  F.Tar.num
MER.numeric[12] <-  Y.Tar.num

MER.comparison <- cbind(MER.analytic[,1], MER.numeric[,1])
colnames(MER.comparison) <- c("Analytic", "Numeric")

write.csv(x=MER.comparison, file=paste(wd,"Outputs/","MER.comparison.csv", sep="")) 


##################################################################################################


# #Not used
# 
# #3. Per recruit analysis
# Fe=seq(0,1, by=0.001) #fishing vector
# SPR.LIM=c(.5,.6,.6)   #limit F%SPR
# 
# A.MAX=A.max(AA[1],AA[2],AA[1])
# 
# 
# 
# 
# 
# 
# #Outputs
# par(mfcol=c(3,1),mai=c(1,.8,.2,.1))
# #YPR
# plot(Fe,FeYield,type='l',ylab="Yield per recruit (kg)",xlab="Fishing mortality")
# arrows(Fmax,fn.yield(Fmax)*.92,Fmax,fn.yield(Fmax),col=2,length=.1)
# text(Fmax,fn.yield(Fmax)*.9,paste("Fmax=",Fmax),col=2)
# 
# 
#     #Potential ratio
#     plot(Fe,SPR,type='l',ylab="Potential ratio",xlab="Fishing mortality",ylim=c(0,1),
#          yaxs="i",xaxs="i",las=1)
#     lines(Fe,EPR,col=2)
#     legend('topright',c("SPR","EPR"),bty='n',col=1:2,lty=1)
#     
# 
#     
#     #F limit
#     arrows(FLim,fn.Lim(FLim),FLim,.04,col=2,length=.1)
#     arrows(0,fn.Lim(FLim),FLim,fn.Lim(FLim),col=2,length=0)
#     text(FLim,.03,FLim,col=2)
#     text(.065,SPR.lim*1.05,paste(SPR.lim*100,"%SPR",sep=""),col=2)
#      
# 
# #4. Compare empirical F (from McAuley et al 2007) with derived Fbrp
# Store.sel=Store.age=list(Whiskery=NA,Dusky=NA,Sandbar=NA)
# fn.get.sel=function(species,age,Lo,Linf,K,alphabeta,alpha,beta)
# {
#   mid.FL.fem=Lo+(Linf-Lo)*(1-exp(-K*age))  #use modified version
#   sel.fem=((mid.FL.fem*10/alphabeta)^alpha)*(exp(alpha-(mid.FL.fem*10/beta)))
#   return(sel.fem=sel.fem)
# }
# 
# for (a in 1:N.sp)
# {  
#   if(!a==3)age=seq(0,max(species.list[[a]]$max.age),by=1)
#   if(a==3) age=seq(0,max(species.list[[a]]$max.age),by=.1)
#   GR=Growth.pars[a,];ab=Sel.pars[a,1];al=Sel.pars[a,2];be=Sel.pars[a,3]
#   Store.sel[[a]]=fn.get.sel(species.names[a],age,Size.birth[a],GR[1,1],GR[1,2],ab,al,be)
#   Store.age[[a]]=age
# }
# 
# FMER=c(0.058,0.017)
# tiff(file="Compare_Fs.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
# 
# par(mfcol=c(2,1),mai=c(1.3,1.3,.1,.1))
# for (a in 2:N.sp)
#   {
#     nn=1:length(Empirc.F[[a-1]])
#     if(!a==3)
#     {
#       plot(as.numeric(names(Empirc.F[[a-1]])),Empirc.F[[a-1]],pch=19,cex=2,col=2,
#            ylab="F at age",xlab="age",cex.lab=1.5,xaxt="n", yaxs="i",xaxs="i")
#       points(FMER[a-1]*Store.sel[[a]][nn],pch=19,col=3,cex=2)
#       legend("right",c("Empirical","FThr"),pch=19,col=2:3,cex=1.25,bty="n")
#       axis(1,at=as.numeric(names(Empirc.F[[a-1]])),
#            labels=paste(as.numeric(names(Empirc.F[[a-1]])),"+",sep=""),tck=-0.04,cex.axis=1.25)
#     }
#     if(a==3)
#     {
#       plot( nn,Empirc.F[[a-1]],pch=19,cex=2,col=2,ylab="F at age",xlab="age",cex.lab=1.5,xaxt="n")
#       rango=list(seq(0,2.9,.1),seq(3,5.9,.1),seq(6,8.9,.1),seq(9,11.9,.1),seq(12,14.9,.1),
#                  seq(15,17.9,.1),seq(18,24,.1),25,26,27,28,29,30)
#       Mean.Sel=length(rango)
#       for(p in 1:length(rango))
#       {
#         id=match(round(rango[[p]],1),round(Store.age[[a]],1))
#         Mean.Sel[p]=mean(Store.sel[[a]][id])        
#       }
#       points(FMER[a-1]*Mean.Sel,pch=19,col=3,cex=2)
#       axis(1,at=nn,labels=names(Empirc.F[[a-1]]),tck=-0.04,cex.axis=1)
#       
#     }
#     legend("topright",species.names[a],cex=1.5,bty="n")    
#   }
# dev.off()
# 

# 


# EqPerRec=function(Fe=0)
# {
#   lx[1]=1    				# Unfished survival at age 0
#   lz[1]=1						# Fished survival at age 0
#   
#   for(i in 2:length(age))			
#   {
#     lx[i]=lx[i-1]*exp(-mm[i])				#Unfished survivorship curve 
#     lz[i]=lz[i-1]*exp(-mm[i]-(Fe*Sel.A[i-1]))		#Fished survivorship curve 
#   }	
#   
#   #Unfished eggs
#   phieggs0=sum(lx*ma*fec)				#total number of eggs per recruit
#   
#   #Fished eggs
#   phieggs=sum(lz*ma*fec)				#total number of eggs per recruit
#   
#   #Spawning potential ratio
#   SPR=phieggs/phieggs0
#   
#   #Objective function
#   epsilon=(Observed-SPR)^2   				#objective function
#   
#   
#   return(list(phieggs0=phieggs0,phieggs=phieggs,epsilon=epsilon,SPR=SPR))
# }

#     BTHR=spr.temp$Dep.MER
#     
#     #4. calculate F reference points for MER
#       #4.1. Set up useful vectors
#         #fecundity  
#     fecundity=spr.temp$fecundity
#     
#         #maturity
#     maturity=spr.temp$maturity
#     
#         #survivorship vectors
#     lx=vector(length=length(age))    #no fishing
#     lz=vector(length=length(age))    #with fishing
#     
#       #4.2 F.THRe  (F.mer)
#     Observed=spr.temp$SPR.mer
#     phieggs0=spr.temp$phi.o
#     fn.yield=function(Fe) EqPerRec(Fe=Fe)$epsilon
#         #numerical approximation
#     Fmer.fit <- nlminb( start=0.1 , objective=fn.yield, lower=0.0001, upper=5.0,
#                         control=list(eval.max=500, iter.max=500  )   )
#     FTHR=Fmer.fit$par
#     FTHR.Conv=Fmer.fit$convergence
#     
#     
#       #4.3 Target reference points
#         #first find SPR from FLim
#     Flim.man.par=1/Prop.FMsy
#     FLim=FTHR*Flim.man.par
#     SPR.lim=EqPerRec(FLim)$SPR
#     
#         #then find B from SPR assuming B-H
#     B.lim=Find.B(spr.temp$alpha,R0=1,spr=SPR.lim,spr0=1)
#     
# #       #4.3 Limit reference points
# #         #first find SPR from BLIMIT
# #     BLIMIT=spr.temp$Dep.MER*spr.temp$p
# #     Observed=BLIMIT
# #     fn.SPR=function(SPR) Find.SPR(SPR=SPR)$epsilon
# #     SPR.fit=nlminb( start=0.6 , objective=fn.SPR, lower=0.001, upper=1.0,
# #                     control=list(eval.max=500, iter.max=500  )   )
# #     SPR.lim.Conv=SPR.fit$convergence
# #     
# #         #then find F from SPR.fit
# #     Observed=SPR.fit$par
# #     Flim.fit=nlminb( start=FTHR , objective=fn.yield, lower=0.0001, upper=5.0,
# #                      control=list(eval.max=500, iter.max=500  )   )
# #     FLim=Flim.fit$par
# #     Flim.Conv=Flim.fit$convergence
#     
#         #4.4 Target reference points
#           #first find SPR from Ftar
#     FTar=FTHR*Prop.FMsy
#     SPR.tar=EqPerRec(FTar)$SPR
#     
#           #then find B from SPR assuming B-H
#     B.tar=Find.B(spr.temp$alpha,R0=1,spr=SPR.tar,spr0=1)
#     
#     
# #     #first find SPR from BTAR
# #     BTAR=spr.temp$Dep.MER*1.2
# #     Observed=BTAR
# #     SPR.fit=nlminb( start=0.8 , objective=fn.SPR, lower=0.001, upper=1.0,
# #                     control=list(eval.max=500, iter.max=500  )   )
# #     SPR.tar.Conv=SPR.fit$convergence
# #     
# #     #then find F from SPR.fit
# #     Observed=SPR.fit$par
# #     Ftar.fin=nlminb( start=FTHR , objective=fn.yield, lower=0.0001, upper=5.0,
# #                      control=list(eval.max=500, iter.max=500  )   )
# #     FTar=Ftar.fin$par
# #     FTar.Conv=Ftar.fin$convergence
# 

# system.time(for (a in 1:N.sp)
# {
#   spec.=species.names[a]
#   WT=Weight.conv[a,]  
#   GROW=Growth.pars[a,]
#   TO=Growth.pars[a,3]
#   SIG=SIGMA[[a]]
#   Lo=Size.birth[a]
#   AA=species.list[[a]]$max.age
#   FF=species.list[[a]]$fec
#   BF=species.list[[a]]$breed.freq
#   b.fem=Length.conv[a,1]
#   a.fem=Length.conv[a,2]
#   alphabeta.S=Sel.pars[a,1]
#   alpha.S=Sel.pars[a,2]
#   beta.S=Sel.pars[a,3]
#   sex.ratio=0.5
#   AMat=species.list[[a]]$age.mat
#   Temper=Temperature[a]
#   r=1
#   spawn.time = 0 # specify time of the year when spawning (or pupping) occurs as a fraction beteween 0 and 1
#   
#   w.g=Wght.G[a]
#   w.ng=Wght.noG[a] 
#   
#   SPR.out=f.out=b.out=Store.M=Store.Fec=Convergence=Store.Sel=Store.len=f.out.30=f.out.40=
#     f.zhou=vector("list", length=N.sim)
#   #names(f.out)=c("threshold","limit","target")
#   Trigged=rep(NA,N.sim)
#   for (j in 1:N.sim)
#   {
#     #1. draw random samples of input parameters
#     A.MAX=A.max(AA[1],AA[2],AA[1])
#     age=0:A.MAX
#     
#     if(Growth.var=="YES") GR=GROWTH(GROW[[1]],GROW[[2]],SIG)
#     if(Growth.var=="NO") GR=matrix(c(GROW[[1]],GROW[[2]]),ncol=2)
#     
#     Rep=REP.PER(BF[2],BF[1])
#     A.MAT=AGE.MAT(AMat[1],AMat[2])    
#     
#     if(!spec.=="gummy")
#     {
#       mid.FL.fem=Lo+(GR[1]-Lo)*(1-exp(-GR[2]*age))  #modified version
#       total.length=b.fem+a.fem*mid.FL.fem
#     }
#     
#     if(spec.=="gummy")
#     {
#       total.length=Lo+(GR[1]-Lo)*(1-exp(-GR[2]*age)) #gummy pars are fn(total length)
#       mid.FL.fem=total.length
#     }
#     
#     Fec=FEC(spec.)
#     
#     #put a cap on gummy Fecundity to avoid predicting beyond data range
#     if(spec.=="gummy")
#     {
#       Fec=ifelse(Fec>Max.fec[1],Max.fec[1],Fec)
#       Fec=ifelse(Fec<0,0,Fec)
#     }
#     
#     
#     mm=M.fun(A.MAX,GR[2],GR[1],Temper,A.MAT,TO,WT[1,1],WT[1,2],w.g,w.ng)
#     Sel.A=Select(alphabeta.S,alpha.S,beta.S)
#     
#     
#     #2. calculate SPR.mer quantities 
#     spr.temp=SPR(A.MAX,mm,Fec,Rep,sex.ratio,A.MAT,Sel.A,0,"Y","ogive")
#     
#     if(length(spr.temp)>1)
#     {
#       #3. Calculate reference points
#       
#       #3.1. Threshold (i.e. MER)    
#       
#       #3.1.1 Biomass (Bmer)
#       BTHR=spr.temp$Dep.MER
#       
#       #3.1.2. F (Fmer)
#       #Set up useful vectors
#       fecundity=spr.temp$fecundity
#       maturity=spr.temp$maturity
#       lx=vector(length=length(age))    #no fishing
#       lz=vector(length=length(age))    #with fishing
#       Observed=spr.temp$SPR.mer
#       phieggs0=spr.temp$phi.o
#       fn.yield=function(Fe) EqPerRec(Fe=Fe)$epsilon
#       
#       #numerical approximation
#       Fstart=0.1
#       Fmer.fit <- nlminb( start=Fstart, objective=fn.yield, lower=0.0001, upper=5.0,
#                           control=list(eval.max=500, iter.max=500  )   )
#       FTHR=Fmer.fit$par
#       FTHR.Conv=Fmer.fit$convergence
#       
#       
#       #3.2. Limits 
#       
#       #3.2.1. Biomass (Blim = p x BTHR)
#       B.lim=Prop.FMsy*BTHR
#       
#       #3.2.2 F
#       obs.biom <- B.lim   
#       b.f <- function(Fstart)
#       {
#         temp.spr = EqPerRec(Fstart)$SPR    
#         temp.SSB =Find.B(spr.temp$alpha,R0=1,spr=temp.spr,spr0=1)  #R0 and spr0 =1 as it's on relative scale
#         yy<-abs(temp.SSB - obs.biom )
#         return(yy)
#       }
#       Flim.fit = nlminb(start=Fstart, objective=b.f, lower=0, upper=3)
#       FLim = Flim.fit$par
#       Flim.Conv = Flim.fit$convergence
#       
#       
#       #3.3. Targets
#       
#       #3.3.1. Biomass (Blim = p x BTHR)
#       B.tar=(1/Prop.FMsy)*BTHR
#       
#       #3.3.2 F
#       obs.biom <- B.tar   
#       FTar.fit = nlminb(start=Fstart, objective=b.f, lower=0, upper=3)
#       FTar = FTar.fit$par
#       FTar.Conv = FTar.fit$convergence
#       
#       
#     }
#     
#     
#     #3. Repeat if nonsense outputs (e.g. h<0.2, negative Biomass RP), i.e. obtain sensible joint priors 
#     Trigger=0
#     if (length(spr.temp)==1|B.tar<0|B.lim<0|BTHR<0)Trigger=1
#     Trigged[j]=Trigger
#     if(Trigger==1)    
#     { repeat 
#     {
#       #3.1 draw random samples
#       A.MAX=A.max(AA[1],AA[2],AA[1])
#       age=0:A.MAX
#       
#       if(Growth.var=="YES") GR=GROWTH(GROW[[1]],GROW[[2]],SIG)
#       if(Growth.var=="NO") GR=matrix(c(GROW[[1]],GROW[[2]]),ncol=2)
#       
#       Rep=REP.PER(BF[2],BF[1])
#       A.MAT=AGE.MAT(AMat[1],AMat[2])    
#       
#       if(!spec.=="gummy")
#       {
#         mid.FL.fem=Lo+(GR[1]-Lo)*(1-exp(-GR[2]*age))  #modified version
#         total.length=b.fem+a.fem*mid.FL.fem
#       }
#       
#       if(spec.=="gummy")
#       {
#         total.length=Lo+(GR[1]-Lo)*(1-exp(-GR[2]*age)) #gummy pars are fn(total length)
#         mid.FL.fem=total.length
#       }
#       
#       Fec=FEC(spec.)
#       
#       #3.2 put a cap on gummy Fecundity to avoid predicting beyond data range
#       if(spec.=="gummy")
#       {
#         Fec=ifelse(Fec>Max.fec[1],Max.fec[1],Fec)
#         Fec=ifelse(Fec<0,0,Fec)
#       }
#       
#       
#       mm=M.fun(A.MAX,GR[2],GR[1],Temper,A.MAT,TO,WT[1,1],WT[1,2],w.g,w.ng)
#       Sel.A=Select(alphabeta.S,alpha.S,beta.S)
#       
#       
#       #3.3 calculate SPR.mer quantities
#       spr.temp=SPR(A.MAX,mm,Fec,Rep,sex.ratio,A.MAT,Sel.A,0,"Y","ogive")
#       
#       if(length(spr.temp)>1)
#       {
#         #3.3.1 Calculate reference points
#         
#         #3.3.1.1 Threshold (i.e. MER)    
#         
#         #Biomass (Bmer)
#         BTHR=spr.temp$Dep.MER
#         
#         #F (Fmer)
#         #Set up useful vectors
#         fecundity=spr.temp$fecundity
#         maturity=spr.temp$maturity
#         lx=vector(length=length(age))    #no fishing
#         lz=vector(length=length(age))    #with fishing
#         Observed=spr.temp$SPR.mer
#         phieggs0=spr.temp$phi.o
#         fn.yield=function(Fe) EqPerRec(Fe=Fe)$epsilon
#         
#         #numerical approximation
#         Fmer.fit <- nlminb( start=0.1 , objective=fn.yield, lower=0.0001, upper=5.0,
#                             control=list(eval.max=500, iter.max=500  )   )
#         FTHR=Fmer.fit$par
#         FTHR.Conv=Fmer.fit$convergence
#         
#         
#         #3.2. Limits 
#         
#         #3.2.1. Biomass (Blim = p x BTHR)
#         B.lim=Prop.FMsy*BTHR
#         
#         #3.2.2 F
#         obs.biom <- B.lim   
#         b.f <- function(Fstart)
#         {
#           temp.spr = EqPerRec(Fstart)$SPR    
#           temp.SSB =Find.B(spr.temp$alpha,R0=1,spr=temp.spr,spr0=1)  #R0 and spr0 =1 as it's on relative scale
#           yy<-abs(temp.SSB - obs.biom )
#           return(yy)
#         }
#         Flim.fit = nlminb(start=Fstart, objective=b.f, lower=0, upper=3)
#         FLim = Flim.fit$par
#         Flim.Conv = Flim.fit$convergence
#         
#         
#         #3.3. Targets
#         
#         #3.3.1. Biomass (Blim = p x BTHR)
#         B.tar=(1/Prop.FMsy)*BTHR
#         
#         #3.3.2 F
#         obs.biom <- B.tar   
#         FTar.fit = nlminb(start=Fstart, objective=b.f, lower=0, upper=3)
#         FTar = FTar.fit$par
#         FTar.Conv = FTar.fit$convergence         
#         
#         
#       }
#       
#       #break if sensible joint prior
#       Trigger=0
#       if (length(spr.temp)==1|B.tar<0|B.lim<0|BTHR<0)Trigger=1
#       
#       
#       if(Trigger==0)break
#     }
#     }
#     
#     
#     #store quantities of interest
#     SPR.out[[j]]=spr.temp
#     Store.M[[j]]=mm
#     Store.Sel[[j]]=Sel.A
#     Store.len[[j]]=mid.FL.fem
#     Store.Fec[[j]]=Fec
#     
#     
#     
#     #4.5 Store all Bs and Fs
#     b.out[[j]]=c(BTHR,B.lim,B.tar)
#     f.out[[j]]=c(FTHR,FLim,FTar)
#     
#     #4.6 Store convergence code
#     Convergence[[j]]=c(FTHR.Conv=FTHR.Conv)
#     
#     
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
#     
#     
#   }
#   Species.SPR.MER[[a]]=SPR.out
#   Species.M[[a]]=Store.M
#   Species.Fec[[a]]=Store.Fec
#   Species.Sel[[a]]=Store.Sel
#   Species.len[[a]]=Store.len
#   Species.B.MER[[a]]=b.out
#   Species.F.MER[[a]]=f.out
#   Species.F.est.Conv[[a]]=do.call(rbind,Convergence)
#   Species.F.40[[a]]=f.out.40
#   Species.F.30[[a]]=f.out.30
#   Species.F.zhou[[a]]=f.zhou
#   Species.Trigged[[a]]=Trigged
# })

# for (a in 1:N.sp)
# {
#   spec.=species.names[a]
#   WT=Weight.conv[a,]  
#   GROW=Growth.pars[a,]
#   TO=Growth.pars[a,3]
#   SIG=SIGMA[[a]]
#   Lo=Size.birth[a]
#   AA=species.list[[a]]$max.age
#   FF=species.list[[a]]$fec
#   BF=species.list[[a]]$breed.freq
#   b.fem=Length.conv[a,1]
#   a.fem=Length.conv[a,2]
#   alphabeta.S=Sel.pars[a,1]
#   alpha.S=Sel.pars[a,2]
#   beta.S=Sel.pars[a,3]
#   sex.ratio=0.5
#   AMat=species.list[[a]]$age.mat
#   Temper=Temperature[a]
#   r=1
#   spawn.time = 0 # specify time of the year when spawning (or pupping) occurs as a fraction beteween 0 and 1
#   
#   w.g=Wght.G[a]
#   w.ng=Wght.noG[a] 
#   
#   SPR.out=f.out=b.out=Store.M=Store.Fec=Convergence=Store.Sel=Store.len=f.out.30=f.out.40=
#     f.zhou=vector("list", length=N.sim)
#   #names(f.out)=c("threshold","limit","target")
#   Trigged=rep(NA,N.sim)
#   for (j in 1:N.sim)
#   {
#     #1. draw random samples of input parameters
#     A.MAX=A.max(AA[1],AA[2],AA[1])
#     age=0:A.MAX
#     
#     if(Growth.var=="YES") GR=GROWTH(GROW[[1]],GROW[[2]],SIG)
#     if(Growth.var=="NO") GR=matrix(c(GROW[[1]],GROW[[2]]),ncol=2)
#     
#     Rep=REP.PER(BF[2],BF[1])
#     A.MAT=AGE.MAT(AMat[1],AMat[2])    
#     
#     if(!spec.=="gummy")
#     {
#       mid.FL.fem=Lo+(GR[1]-Lo)*(1-exp(-GR[2]*age))  #modified version
#       total.length=b.fem+a.fem*mid.FL.fem
#     }
#     
#     if(spec.=="gummy")
#     {
#       total.length=Lo+(GR[1]-Lo)*(1-exp(-GR[2]*age)) #gummy pars are fn(total length)
#       mid.FL.fem=total.length
#     }
#     
#     Fec=FEC(spec.)
#     
#     #put a cap on gummy Fecundity to avoid predicting beyond data range
#     if(spec.=="gummy")
#     {
#       Fec=ifelse(Fec>Max.fec[1],Max.fec[1],Fec)
#       Fec=ifelse(Fec<0,0,Fec)
#     }
#     
#     
#     mm=M.fun(A.MAX,GR[2],GR[1],Temper,A.MAT,TO,WT[1,1],WT[1,2],w.g,w.ng)
#     Sel.A=Select(alphabeta.S,alpha.S,beta.S)
#     
#     
#     #2. calculate SPR.mer quantities 
#     spr.temp=SPR(A.MAX,mm,Fec,Rep,sex.ratio,A.MAT,Sel.A,0,"Y","ogive")
#     if(length(spr.temp)>1)
#     {
#       #3. Calculate reference points
#       
#       #3.1. Threshold (i.e. MER)    
#       
#       #3.1.1 Biomass (Bmer)
#       BTHR=spr.temp$Dep.MER
#       
#       #3.1.2. F (Fmer)
#       #Set up useful vectors
#       fecundity=spr.temp$fecundity
#       maturity=spr.temp$maturity
#       lx=vector(length=length(age))    #no fishing
#       lz=vector(length=length(age))    #with fishing
#       Observed=spr.temp$SPR.mer
#       phieggs0=spr.temp$phi.o
#       fn.yield=function(Fe) EqPerRec(Fe=Fe)$epsilon
#       
#       #numerical approximation
#       Fmer.fit <- nlminb( start=0.1 , objective=fn.yield, lower=0.0001, upper=5.0,
#                           control=list(eval.max=500, iter.max=500  )   )
#       FTHR=Fmer.fit$par
#       FTHR.Conv=Fmer.fit$convergence
#       
#       
#       #3.2. Limits 
#       
#       #3.2.1. F (Flim=(1/p) x Fmer)
#       #first find SPR from FLim
#       Flim.man.par=1/Prop.FMsy
#       FLim=FTHR*Flim.man.par
#       SPR.lim=EqPerRec(FLim)$SPR
#       
#       #3.2.2 Biomass (Blim)
#       #then find B from SPR assuming B-H
#       B.lim=Find.B(spr.temp$alpha,R0=1,spr=SPR.lim,spr0=1)
#       
#       
#       #3.3. Targets
#       
#       #3.3.1. F (Ftar=p x Fmer)
#       #first find SPR from Ftar
#       FTar=FTHR*Prop.FMsy
#       SPR.tar=EqPerRec(FTar)$SPR
#       
#       #3.3.2 Biomass (Btar)
#       #then find B from SPR assuming B-H
#       B.tar=Find.B(spr.temp$alpha,R0=1,spr=SPR.tar,spr0=1)
#       
#     }
#     
#     
#     #3. Repeat if nonsense outputs (e.g. h<0.2, negative Biomass RP), i.e. obtain sensible joint priors 
#     Trigger=0
#     if (length(spr.temp)==1|B.tar<0|B.lim<0|BTHR<0)Trigger=1
#     Trigged[j]=Trigger
#     if(Trigger==1)    
#     { repeat 
#     {
#       #3.1 draw random samples
#       A.MAX=A.max(AA[1],AA[2],AA[1])
#       age=0:A.MAX
#       
#       if(Growth.var=="YES") GR=GROWTH(GROW[[1]],GROW[[2]],SIG)
#       if(Growth.var=="NO") GR=matrix(c(GROW[[1]],GROW[[2]]),ncol=2)
#       
#       Rep=REP.PER(BF[2],BF[1])
#       A.MAT=AGE.MAT(AMat[1],AMat[2])    
#       
#       if(!spec.=="gummy")
#       {
#         mid.FL.fem=Lo+(GR[1]-Lo)*(1-exp(-GR[2]*age))  #modified version
#         total.length=b.fem+a.fem*mid.FL.fem
#       }
#       
#       if(spec.=="gummy")
#       {
#         total.length=Lo+(GR[1]-Lo)*(1-exp(-GR[2]*age)) #gummy pars are fn(total length)
#         mid.FL.fem=total.length
#       }
#       
#       Fec=FEC(spec.)
#       
#       #3.2 put a cap on gummy Fecundity to avoid predicting beyond data range
#       if(spec.=="gummy")
#       {
#         Fec=ifelse(Fec>Max.fec[1],Max.fec[1],Fec)
#         Fec=ifelse(Fec<0,0,Fec)
#       }
#       
#       
#       mm=M.fun(A.MAX,GR[2],GR[1],Temper,A.MAT,TO,WT[1,1],WT[1,2],w.g,w.ng)
#       Sel.A=Select(alphabeta.S,alpha.S,beta.S)
#       
#       
#       #3.3 calculate SPR.mer quantities
#       spr.temp=SPR(A.MAX,mm,Fec,Rep,sex.ratio,A.MAT,Sel.A,0,"Y","ogive")
#       
#       if(length(spr.temp)>1)
#       {
#         #3.3.1 Calculate reference points
#         
#         #3.3.1.1 Threshold (i.e. MER)    
#         
#         #Biomass (Bmer)
#         BTHR=spr.temp$Dep.MER
#         
#         #F (Fmer)
#         #Set up useful vectors
#         fecundity=spr.temp$fecundity
#         maturity=spr.temp$maturity
#         lx=vector(length=length(age))    #no fishing
#         lz=vector(length=length(age))    #with fishing
#         Observed=spr.temp$SPR.mer
#         phieggs0=spr.temp$phi.o
#         fn.yield=function(Fe) EqPerRec(Fe=Fe)$epsilon
#         
#         #numerical approximation
#         Fmer.fit <- nlminb( start=0.1 , objective=fn.yield, lower=0.0001, upper=5.0,
#                             control=list(eval.max=500, iter.max=500  )   )
#         FTHR=Fmer.fit$par
#         FTHR.Conv=Fmer.fit$convergence
#         
#         
#         #3.3.1.2 Limits 
#         #F (FLim)
#         #note: Flim=(1/p) x Fmer
#         #      first find SPR from FLim
#         Flim.man.par=1/Prop.FMsy
#         FLim=FTHR*Flim.man.par
#         SPR.lim=EqPerRec(FLim)$SPR
#         
#         #Biomass (Blim)
#         #then find B from SPR assuming B-H
#         B.lim=Find.B(spr.temp$alpha,R0=1,spr=SPR.lim,spr0=1)
#         
#         
#         #3.3.1.3 Targets
#         #F (Tar)
#         #note:Ftar=p x Fmer)
#         #     first find SPR from Ftar
#         FTar=FTHR*Prop.FMsy
#         SPR.tar=EqPerRec(FTar)$SPR
#         
#         #Biomass (Btar)
#         #then find B from SPR assuming B-H
#         B.tar=Find.B(spr.temp$alpha,R0=1,spr=SPR.tar,spr0=1)          
#         
#       }
#       
#       #break if sensible joint prior
#       Trigger=0
#       if (length(spr.temp)==1|B.tar<0|B.lim<0|BTHR<0)Trigger=1
#       
#       
#       if(Trigger==0)break
#     }
#     }
#     
#     
#     #store quantities of interest
#     SPR.out[[j]]=spr.temp
#     Store.M[[j]]=mm
#     Store.Sel[[j]]=Sel.A
#     Store.len[[j]]=mid.FL.fem
#     Store.Fec[[j]]=Fec
#     
#     
#     
#     #4.5 Store all Bs and Fs
#     b.out[[j]]=c(BTHR,B.lim,B.tar)
#     f.out[[j]]=c(FTHR,FLim,FTar)
#     
#     #4.6 Store convergence code
#     Convergence[[j]]=c(FTHR.Conv=FTHR.Conv)
#     #     Convergence[[j]]=c(FTHR.Conv=FTHR.Conv,SPR.lim.Conv=SPR.lim.Conv,Flim.Conv=Flim.Conv,
#     #                        SPR.tar.Conv=SPR.tar.Conv,FTar.Conv=FTar.Conv)
#     
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
#     
#     
#     #     #7. Calculate F from standard per recruit analysis   #fix, not working well
#     #     age=0:A.MAX
#     #    
#     #       #7.1 YPR
#     #     fn.yield=function(Fe) Per.recruit(Linf=GR[1,1],K=GR[1,2],to=GR[1,3],b=b.fem,a=a.fem,
#     #             M=mm,bwt=WT[1,1],awt=WT[1,2],age.mat=A.MAT,fec=Fec,breed.freq=Rep,fishing=Fe)$Epsilon  
#     #     Fmax.optim=optimize(fn.yield,lower=0.0001,upper=5)          
#     #     Fmax=Fmax.optim$minimum
#     #         
#     #       #7.2 Potential ratio
#     #     fn.yield=function(Fe) Per.recruit(Linf=GR[1,1],K=GR[1,2],to=GR[1,3],b=b.fem,a=a.fem,
#     #             M=mm,bwt=WT[1,1],awt=WT[1,2],age.mat=A.MAT,fec=Fec,breed.freq=Rep,fishing=Fe)$Epsilon.mer  
#     #     Fmer.optim=optimize(fn.yield,lower=0.0001,upper=5)        	
#     #     Fmer=Fmer.optim$minimum
#     
#     
#   }
#   Species.SPR.MER[[a]]=SPR.out
#   Species.M[[a]]=Store.M
#   Species.Fec[[a]]=Store.Fec
#   Species.Sel[[a]]=Store.Sel
#   Species.len[[a]]=Store.len
#   Species.B.MER[[a]]=b.out
#   Species.F.MER[[a]]=f.out
#   Species.F.est.Conv[[a]]=do.call(rbind,Convergence)
#   Species.F.40[[a]]=f.out.40
#   Species.F.30[[a]]=f.out.30
#   Species.F.zhou[[a]]=f.zhou
#   Species.Trigged[[a]]=Trigged
# }



