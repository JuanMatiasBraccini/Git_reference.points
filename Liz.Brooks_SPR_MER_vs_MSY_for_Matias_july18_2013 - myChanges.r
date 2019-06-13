rm(list=ls(all=TRUE))


### Matias' implementation  ####


#1.1 life history
#Dusky parameters as per Brooks et al 2010
#ages
age.min=0
age.max=40
age=age.min:age.max

#age at recruitment
r=1    
first=match(r,age)
last=match(age.max,age)

Numb.pups=7.1     #mean number of pups
cycle.length=3    #reproductive cycle length
sex.ratio=0.5     #embryo sex ratio
FEC= Numb.pups*sex.ratio/cycle.length   #annual female fecundity
fec=c(rep(0,21),rep(FEC,length(age)-21)) #productivity vector (fecundity X sex ratio X breeding cycle)

#proportion mature at age
ma=c(rep(0,11),.001,.002,.006,.013,.03,.062,.124,.226,.37,.535,.687,.803,.881,.929,.958,.975,.985,
     .991,.994,.996,.998,.998,.999,.999,rep(1,6))

#survivorship
surv=c(.78,.865,.875,.883,.889,.894,.899,.903,.907,.91,.913,.915,.918,.919,.922,.923,.925,.926,.927,.928,
       .93,.93,.932,.933,.934,.934,.935,.935,.936,.937,.937,.938,.939,.939,.939,.94,.94,.94,.94,.941,
       .942)

m=-log(surv)



#----PROCEDURE ----

SPR=function(max.age,M,Numb.pups,breed.freq,sex.ratio,age.mat,ma,sel.age,F.mult=0,Plus.gp="Y")
{
  #age
  age=0:max.age
  
  #survivorship
  surv=exp(-M)    
  
  #fecundity  
  fecundity=ifelse(age>=age.mat,Numb.pups*breed.freq*sex.ratio,0)

  #maturity
  maturity=ma   


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
  
  return(list(phi.o=phi.o,alpha=alpha,h=h,SPR.mer=SPR.mer,Dep.MER=Dep.MER,p=p))
  
}

MER.pars=SPR(age.max,m,Numb.pups,1/cycle.length,sex.ratio,21,ma,sel.age,F.mult=0,Plus.gp="Y")
Observed=MER.pars$SPR.mer

                                                     

#Esimate F
  #some useful vecs
#Sel.A=c(1,rep(0,age.max))

##
#dusky                #Simpfendorfer & Unsworth 1998
theta1=130.13
theta2=29237
mesh= 6.5 #(in inches, = 16.5 cm)
alphabeta=theta1*mesh
beta=-0.5*(theta1*mesh-((theta1^2)*(mesh^2)+4*theta2)^0.5)
alpha=alphabeta/beta
mid.FL.fem=374.4*(1-exp(-.0367*((0:age.max)+0.5-(-3.3) )))
Sel.A=((mid.FL.fem*10/alphabeta)^alpha)*(exp(alpha-(mid.FL.fem*10/beta)))

##


lx=vector(length=length(age))  	#survivorship vector with no fishing
lz=vector(length=length(age))		#survivorship vector with fishing


  #Equilibrium Per recruit function
EqPerRec=function(Fe=0)
{
  lx[1]=1						# Unfished survival at age 0
  lz[1]=1						# Fished survival at age 0
  
  for(i in 2:length(age))			
  {
    lx[i]=lx[i-1]*exp(-m[i])				#Unfished survivorship curve 
    lz[i]=lz[i-1]*exp(-m[i]-(Fe*Sel.A[i-1]))		#Fished survivorship curve 
  }	
  
  #Unfished eggs
  phieggs0=sum(lx*ma*fec)				#total number of eggs per recruit
  
  #Fished eggs
  phieggs=sum(lz*ma*fec)				#total number of eggs per recruit
  
  #Spawning potential ratio
  SPR=phieggs/phieggs0
  
  #Objective function
  epsilon=(Observed-SPR)^2   				#objective function
  
  
  return(list(phieggs0=phieggs0,phieggs=phieggs,epsilon=epsilon,SPR=SPR))
}



Fe=seq(0,1, by=0.002)    
fn.yield=function(Fe) EqPerRec(Fe=Fe)$epsilon
fn.SPR=function(Fe) EqPerRec(Fe=Fe)$SPR
FeYield=sapply(Fe,fn.yield)					
Fmer.LP=Fe[which.min(FeYield)]	#likelihood profile 

fit.F=optimize(fn.yield,lower=0.0001,upper=5)   #numerical approximation
Fmer=fit.F$minimum


plot(Fe,sapply(Fe,fn.SPR),type='l',xlab="F",ylab="SPR")
points(Fmer,Observed,col=2,pch=19)


SPR.Fmer=EqPerRec(Fmer)$SPR  #SPR for estimated F==Observed
######








###################################################################



### Liz's implementation  ####



#------------------USER MUST DEFINE THESE---------------------------#


# specify the biological parameters at age

nages <- age.max  # specify the number of age classes
ages <- seq(1,nages) # this will create a vector from 1 to the number of ages

#Natural mortality
M.age=-log(surv)

# mortality of pups (survival from age 0 to age 1)
# note: this should be the value that gives maximal pup survival when stock size is really low
M.pup <-  M.age[1]

# natural mortality for ages 1-40
M.age=M.age[2:length(M.age)]



  # maturity at age (ages 1:40)
mat.age=ma[2:length(age)]

  # pup production at age (number of pups born per female per year, ages 1:40)
pup.prod.age <- fec[2:length(age)]

  # sex ratio (if 50:50 male female, enter 0.5)
sex.ratio <- 0.5
  # specify time of the year when spawning (or pupping) occurs as a fraction beteween 0 and 1
spawn.time = 0


#   # weight at age
Linf.d=374.4
K.d=.0367
to.d=-3.3
FL.fem=Linf.d*(1-exp(-K.d*(ages-to.d))) #in cm
bwt.d=1.2334e-5
awt.d=2.855
wgt.age=bwt.d*FL.fem^awt.d #in kg


#   # fishery selectivity at age
#sel.age=Sel.A[1:nages]
sel.age=Sel.A
 
##-----------------------That's all----------------------------------
##-----------------Do not change anything below this line------------
##-------------------------------------------------------------------


####     FUNCTIONS   ###########################################################
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
#---------------------------------------------------------
#-------Equilibrium SSB (includes Ricker option) ------------------------
ssb.eq<-function(recr.par, R0, spr, spr0, is.steepness=T, SRtype ) {
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
#------------------------------------  
 bev.holt.alpha<-function(S,R0,recr.par,spr0, is.steepness=T){

 if (is.steepness==T)  alpha.BH <- 4*recr.par/(1-recr.par)
 if (is.steepness==F)  alpha.BH <- recr.par


  y=rep(0,length(S))
  y=R0*S*alpha.BH/(R0*spr0 +(alpha.BH-1.0)*S)
  return(y)

}  
#------------------------------------  
##----  Calculate Fmsy, MSY given selectivity and SR relationship  ------------------------------


get.MSY.vector <- function(nages, pup.prod.age, sex.ratio, mat.age, M.age,
                           sel.age, spawn.time, wgt.age, alpha.hat, R.vals=1, srtype=1) {
 
  R0 <- 1
  sr0 <- 1
  sr0.calc <- s.per.recr(nages=nages, fec.age=pup.prod.age, mat.age=mat.age, M.age= M.age, F.mult=0, sel.age=sel.age, spawn.time=spawn.time)  
  ahat.val <- alpha.hat
    
  MSY.soln <- rep(NA, 8)
  names(MSY.soln) <- c("Fmsy", "MSY", "SPRmsy", "Rel.SSBmsy", "Rel.Rmsy", "YPRmsy", 
  "Convergence.code", "steepness")



  F.start=0.12
  get.yield.f.min <- function(F.start) {

      temp.spr = s.per.recr(nages=nages, fec.age=(pup.prod.age*sex.ratio), mat.age=mat.age, M.age= M.age, F.mult=F.start, sel.age=sel.age, spawn.time=spawn.time)/sr0.calc   
      temp.ypr = ypr(nages=nages, wgt.age=wgt.age, M.age=M.age,  F.mult=F.start, sel.age=sel.age )
      temp.SSB = ssb.eq( recr.par=ahat.val, R0=R0, spr=temp.spr, spr0=sr0, is.steepness=F, SRtype=srtype )
      yield = temp.ypr*temp.SSB/temp.spr          
      yy=-1*yield
  return(yy)  
              }   # end get.yield.f.min function

F.nlmin <- nlminb( start=F.start , objective=get.yield.f.min, lower=0.0001, upper=5.0,
              control=list(eval.max=500, iter.max=500  )   )



MSY.soln[1] <- F.nlmin$par  #Fmsy
MSY.soln[2] <- -1*F.nlmin$objective #MSY

# calculate other MSY quantities

spr.nlmin <- s.per.recr(nages=nages, fec.age=pup.prod.age*sex.ratio , mat.age=mat.age , M.age= M.age , F.mult=F.nlmin$par, sel.age=sel.age , spawn.time=spawn.time)/sr0.calc   
MSY.soln[3] <- spr.nlmin/sr0  #SPRmsy

ssb.nlmin <- ssb.eq( recr.par=ahat.val, R0=R0 , spr=spr.nlmin, spr0=sr0 , is.steepness=F, SRtype=srtype )
MSY.soln[4] <-  ssb.nlmin #relative SSBmsy

  if (srtype==1) {
Rmsy.nlmin <- bev.holt.alpha(S=ssb.nlmin,R0=R0 ,recr.par=ahat.val ,spr0=sr0, is.steepness=F )
MSY.soln[5 ] <- Rmsy.nlmin #relative Rmsy
               }


               
ypr.nlmin <- ypr(nages=nages, wgt.age=wgt.age , M.age=M.age ,  F.mult=F.nlmin$par, sel.age=sel.age  )
MSY.soln[6 ] <- ypr.nlmin #YPRmsy


MSY.soln[7] <-  F.nlmin$convergence  #code=0 indicates convergence
MSY.soln[8] <-  steep  #steepness
    
    return(MSY.soln)

   }  # end function
  
  
###-------------------------------------------------------------------------------------

#===============================================================================
#===============================================================================
##     End of Functions
#===============================================================================
#===============================================================================

#unexploited spawners per recruit
s.per.r0 <- s.per.recr(nages=nages, fec.age=pup.prod.age, mat.age=mat.age, M.age=M.age, 
          F.mult=0, sel.age=sel.age, spawn.time=spawn.time )

#maximum reproductive rate (sensu Myers)
alpha.hat <- exp(-M.pup)*s.per.r0

#steepness 
steep <- alpha.hat/(alpha.hat+4)

1/sqrt(alpha.hat)


# solve for MER (using input selectivity vector)
MER_user_sel <- get.MSY.vector(nages=nages, pup.prod.age=pup.prod.age, sex.ratio=sex.ratio, 
                     mat.age=mat.age, M.age=M.age, sel.age=sel.age, spawn.time=spawn.time, 
                     wgt.age=rep(1,nages), alpha.hat=alpha.hat, R.vals=1, srtype=1)


# solve for MSY (using input selectivity vector)
MSY_user_sel <- get.MSY.vector(nages=nages, pup.prod.age=pup.prod.age, sex.ratio=sex.ratio, 
                     mat.age=mat.age, M.age=M.age, sel.age=sel.age, spawn.time=spawn.time, 
                     wgt.age=wgt.age, alpha.hat=alpha.hat, R.vals=1, srtype=1)


MSY_Table <- cbind(MER_user_sel, MSY_user_sel)




#COMPARE APPROACHES
SpawnPerRec.0=c(MER.pars$phi.o,s.per.r0)
ALPHAS=c(MER.pars$alpha,alpha.hat)
FMSYs=c(Fmer,MER_user_sel[1])
names(ALPHAS)=names(SpawnPerRec.0)=names(FMSYs)=c("Matias(DuskyPaper)","Liz")
