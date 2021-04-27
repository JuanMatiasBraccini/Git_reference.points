#  Derive Reference points from life history parameters and selectivity 
#   (for Matias Braccini)
#
#  Created on: 14 June 2013
#  Last updated: 18 July 2013
#
#  Liz Brooks  (Liz.Brooks@noaa.gov)
#  
#  Example assumes 16 age classes with selectivity on first age only (fake data)
#    the age 16 category is assumed to be a plus group
#  Beverton-Holt functional form is assumed for the stock recruit relationship
#
#  ***NOTE: Units for SSB are number of pups per female (not in biomass)
# 
# This code creates a .csv file with estimates for MER and MSY


rm(list=ls(all=TRUE)) # Remove all variables, etc. from the session memory
graphics.off()  # close any open graphics windows

#------------------USER MUST DEFINE THESE---------------------------#

handl_OneDrive=function(x)paste('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias',x,sep='/')

#working directory (copy from Explorer window and double all "\")
wd <-  setwd(handl_OneDrive("Analyses/Reference Points/From Liz Brooks"))    

# specify the biological parameters at age

nages <- 16  # specify the number of age classes
ages <- seq(1,nages) # this will create a vector from 1 to the number of ages

  # natural mortality for ages 1-12
M.age <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
  # mortality of pups (survival from age 0 to age 1)
  # note: this should be the value that gives maximal pup survival when stock size is really low
M.pup <-  0.2
  # maturity at age
mat.age <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0.75, 1, 1)
  # pup production at age (number of pups born per female per year)
pup.prod.age <- c(0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 4, 6, 8, 8)
  # sex ratio (if 50:50 male female, enter 0.5)
sex.ratio <- 0.5
  # specify time of the year when spawning (or pupping) occurs as a fraction beteween 0 and 1
spawn.time = 0
  # weight at age
wgt.age <- c(10.8,  60.9, 147.0, 253.0, 363.0, 469.0, 564.0, 646.0, 716.0, 773.0,
      820.0, 858.0, 888.0, 912.0, 931.0, 946.0)
  # fishery selectivity at age
sel.age <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

 
##-----------------------That's all----------------------------------
##-----------------Do not change anything below this line------------
##-------------------------------------------------------------------

setwd(wd)


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
  sr0.calc <- s.per.recr(nages=nages, fec.age=(pup.prod.age*sex.ratio), mat.age=mat.age, M.age= M.age, F.mult=0, sel.age=sel.age, spawn.time=spawn.time)  
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
s.per.r0 <- s.per.recr(nages=nages, fec.age=pup.prod.age*sex.ratio, mat.age=mat.age, M.age=M.age, 
          F.mult=0, sel.age=sel.age, spawn.time=spawn.time )

#maximum reproductive rate (sensu Myers)
alpha.hat <- exp(-M.pup)*s.per.r0

#steepness 
steep <- alpha.hat/(alpha.hat+4)

1/sqrt(alpha.hat)

# solve for MER (assuming selectivity=1 for all ages)
MER_full_sel <- get.MSY.vector(nages=nages, pup.prod.age=pup.prod.age, sex.ratio=sex.ratio, 
                     mat.age=mat.age, M.age=M.age, sel.age=rep(1,nages),spawn.time=spawn.time, 
                     wgt.age=rep(1,nages), alpha.hat=alpha.hat, R.vals=1, srtype=1)


# solve for MER (using input selectivity vector)
MER_user_sel <- get.MSY.vector(nages=nages, pup.prod.age=pup.prod.age, sex.ratio=sex.ratio, 
                     mat.age=mat.age, M.age=M.age, sel.age=sel.age, spawn.time=spawn.time, 
                     wgt.age=rep(1,nages), alpha.hat=alpha.hat, R.vals=1, srtype=1)



# solve for MSY (assuming selectivity=1 for all ages)
MSY_full_sel <- get.MSY.vector(nages=nages, pup.prod.age=pup.prod.age, sex.ratio=sex.ratio, 
                     mat.age=mat.age, M.age=M.age, sel.age=rep(1,nages), spawn.time=spawn.time, 
                     wgt.age=wgt.age, alpha.hat=alpha.hat, R.vals=1, srtype=1)



# solve for MSY (using input selectivity vector)
MSY_user_sel <- get.MSY.vector(nages=nages, pup.prod.age=pup.prod.age, sex.ratio=sex.ratio, 
                     mat.age=mat.age, M.age=M.age, sel.age=sel.age, spawn.time=spawn.time, 
                     wgt.age=wgt.age, alpha.hat=alpha.hat, R.vals=1, srtype=1)


MSY_Table <- cbind(MER_full_sel, MER_user_sel, MSY_full_sel, MSY_user_sel)


# write out MSY_Table as a .csv file
write.csv(MSY_Table, file="MSY_Table.csv", row.names=T)

