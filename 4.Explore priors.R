Estore=vector('list',length=3)
for (a in 1:N.sp)
{
WT=Weight.conv[a,]  
GR=Growth.pars[a,]
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



AA.sim=GR.sim=Fec.sim=Rep.sim=A.mat.sim=vector("list", length=N.sim)
for (j in 1:N.sim)
{
  #1. draw random samples of input parameters
  A.MAX=A.max(AA[1],AA[2],AA[1])
  age=0:A.MAX
  GR=GROWTH(GR[[1]],GR[[2]],SIG)     
  Fec=FEC(FF[1],FF[2],round(mean(c(FF[1],FF[2]))))
  Rep=REP.PER(BF[2],BF[1])
  A.MAT=AGE.MAT(AMat[1],AMat[2])    
  mid.FL.fem=Lo+(GR[1]-Lo)*(1-exp(-GR[2]*age))  #modified version
  total.length=b.fem+a.fem*mid.FL.fem
  mm=M.fun(A.MAX,GR[2],GR[1],Temper,A.MAT,TO,WT[1,1],WT[1,2],w.g,w.ng)
  Sel.A=Select(alphabeta.S,alpha.S,beta.S)
  
  
  #2. calculate SPR.mer quantities 
  spr.temp=SPR(A.MAX,mm,Fec,Rep,sex.ratio,A.MAT,Sel.A,0,"Y","ogive")

  AA.sim[[j]]=A.MAX
  GR.sim[[j]]=GR
  Fec.sim[[j]]=Fec
  Rep.sim[[j]]=Rep
  A.mat.sim[[j]]=A.MAT
}

Estore[[a]]=list(AA.sim,GR.sim,Fec.sim,Rep.sim,A.mat.sim)
}

#Max age
par(mfcol=c(3,1),las=1,mai=c(.6,.5,.075,.15),omi=c(.1,.1,.1,.05),mgp=c(1.5,0.6,0))
Rango=c(15,55)
hist(unlist(Estore[[1]][[1]]),main="whisk",xlab="max age",xlim=Rango)
legend("topright",paste(round( sd(unlist(Estore[[1]][[1]])),4)),bty="n")
hist(unlist(Estore[[2]][[1]]),main="dusk",xlab="max age",xlim=Rango)
legend("topright",paste(round( sd(unlist(Estore[[2]][[1]])),4)),bty="n")
hist(unlist(Estore[[3]][[1]]),main="sand",xlab="max age",xlim=Rango)
legend("topright",paste(round( sd(unlist(Estore[[3]][[1]])),4)),bty="n")



aa=do.call(rbind,Estore[[1]][[2]])
bb=do.call(rbind,Estore[[2]][[2]])
cc=do.call(rbind,Estore[[3]][[2]])
#Linf
Rango=c(122,400)
par(mfcol=c(3,1),las=1,mai=c(.6,.5,.075,.15),omi=c(.1,.1,.1,.05),mgp=c(1.5,0.6,0))
hist(aa[,1],main="whisk",xlab="Linf",xlim=Rango)
legend("topright",paste(round( sd(aa[,1]),3)),bty="n")
hist(bb[,1],main="dusk",xlab="Linf",xlim=Rango)
legend("topright",paste(round( sd(bb[,1]),3)),bty="n")
hist(cc[,1],main="sand",xlab="Linf",xlim=Rango)
legend("topright",paste(round( sd(cc[,1]),3)),bty="n")


#k
Rango=c(0.01,0.4)
par(mfcol=c(3,1),las=1,mai=c(.6,.5,.075,.15),omi=c(.1,.1,.1,.05),mgp=c(1.5,0.6,0))
hist(aa[,2],main="whisk",xlab="k",xlim=Rango)
legend("topright",paste(round( sd(aa[,2]),4)),bty="n")
hist(bb[,2],main="dusk",xlab="k",xlim=Rango)
legend("topright",paste(round( sd(bb[,2]),4)),bty="n")
hist(cc[,2],main="sand",xlab="k",xlim=Rango)
legend("topright",paste(round( sd(cc[,2]),4)),bty="n")


#Fec
Rango=c(2,28)
par(mfcol=c(3,1),las=1,mai=c(.6,.5,.075,.15),omi=c(.1,.1,.1,.05),mgp=c(1.5,0.6,0))
hist(unlist(Estore[[1]][[3]]),main="whisk",xlab="Fec",xlim=Rango)
legend("topright",paste(round( sd(unlist(Estore[[1]][[3]])),4)),bty="n")
hist(unlist(Estore[[2]][[3]]),main="dusk",xlab="Fec",xlim=Rango)
legend("topright",paste(round( sd(unlist(Estore[[2]][[3]])),4)),bty="n")
hist(unlist(Estore[[3]][[3]]),main="sand",xlab="Fec",xlim=Rango)
legend("topright",paste(round( sd(unlist(Estore[[3]][[3]])),4)),bty="n")



#Rep
Rango=c(1,3)
par(mfcol=c(3,1),las=1,mai=c(.6,.5,.075,.15),omi=c(.1,.1,.1,.05),mgp=c(1.5,0.6,0))
hist(unlist(Estore[[1]][[4]]),main="whisk",xlab="Rep",xlim=Rango)
legend("topright",paste(round( sd(unlist(Estore[[1]][[4]])),4)),bty="n")
hist(unlist(Estore[[2]][[4]]),main="dusk",xlab="Rep",xlim=Rango)
legend("topright",paste(round( sd(unlist(Estore[[2]][[4]])),4)),bty="n")
hist(unlist(Estore[[3]][[4]]),main="sand",xlab="Rep",xlim=Rango)
legend("topright",paste(round( sd(unlist(Estore[[3]][[4]])),4)),bty="n")


#Maturity
Rango=c(5,35)
par(mfcol=c(3,1),las=1,mai=c(.6,.5,.075,.15),omi=c(.1,.1,.1,.05),mgp=c(1.5,0.6,0))
hist(unlist(Estore[[1]][[5]]),main="whisk",xlab="maturity",xlim=Rango)
legend("topright",paste(round( sd(unlist(Estore[[1]][[5]])),4)),bty="n")
hist(unlist(Estore[[2]][[5]]),main="dusk",xlab="maturity",xlim=Rango)
legend("topright",paste(round( sd(unlist(Estore[[2]][[5]])),4)),bty="n")
hist(unlist(Estore[[3]][[5]]),main="sand",xlab="maturity",xlim=Rango)
legend("topright",paste(round( sd(unlist(Estore[[3]][[5]])),4)),bty="n")
