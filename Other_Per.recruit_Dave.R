What=1#whiskery
#What=2#dhufish

if(What==1)
{
  i=1
  Max.age=species.list.Min[[i]]$max.age
  Linf=Growth.pars[i,1]
  K=Growth.pars[i,2]
  to=Growth.pars[i,3]
  b=Length.conv[i,1]
  a=Length.conv[i,2]
  alpha=Sel.pars[i,2]
  beta=Sel.pars[i,3]
  alphabeta=alpha*beta
  M=species.list.Min[[i]]$M
  bwt=Weight.conv[i,1]
  awt=Weight.conv[i,2]
  age.mat=species.list.Min[[i]]$age.mat

  fec=species.list.Min[[i]]$fec

  breed.freq=species.list.Min[[i]]$breed.freq

  sex.ratio=species.list.Min[[i]]$sex.ratio
  
  Selec.50=1
  Selec.95=2
  
}

#Dhufish
if(What==2)
{
  Max.age=100
  M=0.11
  Z=0.51
  fishing=0.4
  
  Linf=929
  K=0.111
  to=(-0.141)
  
  awt=2.856641924
  bwt=-10.07645271
  
  L50.mature=331
  L95.mature=509
  
  Fecundity.1=0.0841
  Fecundity.2=10.432
  
  Selec.50=7.45
  Selec.95=8.75  
}



#Model
age=0:Max.age
Per.rec=function(fishing)
{

  N=length(age)
  
  #length at mid age
  if(What==1)FL=Linf*(1-exp(-K*(age+0.5-to)))
  if(What==2)
    {
      FL=Linf*(1-exp(-K*(age-to)))
      FL[1]=0.01
    }
  if(What==1)TL=b+a*FL
  
  #selectivity at age
  if(What==1)selectivity<<-((FL*10/alphabeta)^alpha)*(exp(alpha-(FL*10/beta)))
  #if(What==2)selectivity<<-1/(1+exp((-log(19))*(age-Selec.50)/(Selec.95-Selec.50)))
  #selectivity<<-1/(1+exp((-log(19))*(age-Selec.50)/(Selec.95-Selec.50)))

  #mortality at age
  Fat.age=fishing*selectivity
  Z.at.age=Fat.age+M
  
  #numbers at age
  Rel.numbers=vector(length=N)
  Rel.numbers[1]=1
  for(i in 2:N)Rel.numbers[i]=Rel.numbers[i-1]*exp(-Z.at.age[i])
  
  
  #harvest fraction
  harvest=(Fat.age/Z.at.age)*Rel.numbers*(1-exp(-Z.at.age))
  
  
  #weigth at age
  if(What==1)Wt=bwt*TL^awt
  if(What==2)Wt=exp(awt*log(FL)+bwt)
  
  #fecundity at age (considering Maturity and breeding freq)
  if(What==1)fecundity=ifelse(age>=age.mat,fec*breed.freq,0)
  if(What==2)
    {
      maturity=1/(1+exp((-log(19))*(FL-L50.mature)/(L95.mature-L50.mature)))
      N.mature=Rel.numbers*maturity
    }
    
  #YPR
  YPR=Wt*harvest
  
  #Spawning stock biomass per recruit                           #Check if need breeding frequency!!!
  if(What==1)SSB=ifelse(age>=age.mat,Wt*Rel.numbers*breed.freq,0)
  if(What==2)SSB=N.mature*Wt
    
  #Eggs per recruit                                             #Check if need sex ratio!!!
  if(What==1)Eggs=fecundity*sex.ratio*Rel.numbers
  if(What==2)
    {
      Eggs=N.mature*(Fecundity.1*FL-Fecundity.2)^3
      Eggs[1:2]=0 
  }
  #Totals per recruit
  N.per.rec=sum(Rel.numbers)
  Yield.per.rec=sum(YPR) #(in kg)
  SSB.per.rec=sum(SSB) #(in kg)
  Eggs.per.rec=sum(Eggs)
  
  return(list(N.per.rec=N.per.rec,Yield.per.rec=Yield.per.rec,
              SSB.per.rec=SSB.per.rec, Eggs.per.rec=Eggs.per.rec))
}

#Apply F vector
Fe=seq(0,1,by=.001)


fn.yield=function(Fe) Per.rec(fishing=Fe)$Yield.per.rec  
FeYield=sapply(Fe,fn.yield)  				
Fmax=Fe[which.max(FeYield)]

#Yield per recruit
plot(Fe,FeYield,type='l',ylab="Yield per recruit (kg)",xlab="Fishing mortality")
arrows(Fmax,fn.yield(Fmax)*.92,Fmax,fn.yield(Fmax),col=2,length=.1)
text(Fmax,fn.yield(Fmax)*.9,paste("Fmax=",Fmax),col=2)




FeYield=sapply(Fe,Per.rec)
N.per.rec=unlist(FeYield[1,])
SSB.per.rec=unlist(FeYield[3,])
Eggs.per.rec=unlist(FeYield[4,])


#Potential ratio
SPR=SSB.per.rec/SSB.per.rec[1]
EPR=Eggs.per.rec/Eggs.per.rec[1]
plot(Fe,SPR,type='l',ylab="Potential ratio",xlab="Fishing mortality",ylim=c(0,1),
     yaxs="i",xaxs="i",las=1)
lines(Fe,EPR,col=2)
legend('topright',c("SPR","EPR"),bty='n',col=1:2,lty=1)

#Get F for limit
SPR.lim=0.500
fn.Lim=function(Fe) Per.rec(fishing=Fe)$SSB.per.rec/Per.rec(fishing=0)$SSB.per.rec  
Fe.Lim=round(sapply(Fe,fn.Lim),2)    			
FLim=mean(Fe[which(Fe.Lim==SPR.lim)])
arrows(FLim,fn.Lim(FLim),FLim,.04,col=2,length=.1)
arrows(0,fn.Lim(FLim),FLim,fn.Lim(FLim),col=2,length=0)
text(FLim,.03,FLim,col=2)
text(.065,SPR.lim*1.05,paste(SPR.lim*100,"%SPR",sep=""),col=2)




#Selectivity
plot(age,selectivity)