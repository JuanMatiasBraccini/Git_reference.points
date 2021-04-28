##################### Catch-based MSY estimation #############################

#note: implementation of Martell & Froese 2012.
#       total catches must be used (i.e. all sources of F)


rm(list=ls(all=TRUE))

set.seed(999)  ## for same random sequence



#---DATA SECTION-----
if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')

#TDGDLF
setwd(handl_OneDrive("Analyses/Reference Points"))
Dus.catch=read.csv("Tot.c.dusky.csv")
Gum.catch=read.csv("Tot.c.gummy.csv")
San.catch=read.csv("Tot.c.sandbar.csv")
Whi.catch=read.csv("Tot.c.whiskery.csv")

#Other Fisheries
Dus.catch.other=read.csv("Tot.c.dusky.other.csv")
Dus.catch.other=Dus.catch.other[,1:2]
San.catch.other=read.csv("Tot.c.sandbar.NSF.csv")

Dus.catch[,2]=Dus.catch[,2]+Dus.catch.other[,2]
San.catch[,2]=San.catch[,2]+San.catch.other[,2]



yr=as.numeric(sapply(strsplit(as.character(Dus.catch$finyear),"-"), "[", 1))
FinYear=as.character(Dus.catch$finyear)

Dus.catch$finyear=Gum.catch$finyear=San.catch$finyear=Whi.catch$finyear=yr

names(Dus.catch)[1]="yr"  
names(Gum.catch)[1]="yr"
names(San.catch)[1]="yr"
names(Whi.catch)[1]="yr"

names(Dus.catch)[2]="ct"  
names(Gum.catch)[2]="ct"
names(San.catch)[2]="ct"
names(Whi.catch)[2]="ct"

Sp=c("Dusky","Gummy","Sandbar","Whiskery")

Dus.catch$stock=Sp[1]
Gum.catch$stock=Sp[2]
San.catch$stock=Sp[3]
Whi.catch$stock=Sp[4]

Dus.catch$res="Very low"
Gum.catch$res="Low"
San.catch$res="Very low"
Whi.catch$res="Low"

Conv.fac=1e6    #convert kg to 1000s tons
Dus.catch$ct=Dus.catch$ct/Conv.fac  
Gum.catch$ct=Gum.catch$ct/Conv.fac
San.catch$ct=San.catch$ct/Conv.fac
Whi.catch$ct=Whi.catch$ct/Conv.fac

San.catch=subset(San.catch,!is.na(ct)) #remove NA catches

cdat=rbind(Dus.catch,Gum.catch,San.catch,Whi.catch)




#---PARAMETERS SECTION-----

#PRIORS

# 1. r (set boundaries of lognormal prior) 

# for Dusky and Sandbar use McAuley et al 2007
# for Gummy and Whiskery I did my own demographic analysis (see Demographic_analysis.R in R scripts)
# Xiao & Walker used wrong life history pars

r.list <- list(Dusky=c(0.001,0.05), Gummy=c(0.1,0.47), Sandbar=c(0.001,0.05),Whiskery=c(0.03,0.20))  
CV=0.3  #SD used for r lognormal prior


# 2. k (in tonnes)

K.upper=50
k.list <- list(Dusky=c(max(Dus.catch$ct),K.upper*max(Dus.catch$ct)),     #MISSING (discuss with managers)
               Gummy=c(max(Gum.catch$ct),K.upper*max(Gum.catch$ct)),
               Sandbar=c(max(San.catch$ct),K.upper*max(San.catch$ct)),
               Whiskery=c(max(Whi.catch$ct),K.upper*max(Whi.catch$ct))) 


#3. boundaries of initial and final depletion ranges

  #3.1 User defined             #(discuss with managers)
    #3.1.1. biomass range at start of time series (as fraction of k)    
startbio.list   <- list(Dusky=c(0.8,.95), Gummy=c(0.8,.95), Sandbar=c(0.25,0.75),Whiskery=c(0.8,.95))

    #3.1.2. biomass range after last catch (as fraction of k)
finalbio.list   <- list(Dusky=c(0.3, 0.7), Gummy=c(0.3, 0.7), Sandbar=c(0.3, 0.7),Whiskery=c(0.3, 0.7))


  #3.2 Default
#note: this assumes C trend is Abundance trend proxy, it ignores targeting changes, etc
    #3.2.1. biomass range at start of time series
B.i.def.fn=function(DAT)if(DAT$ct[nrow(DAT)]/max(DAT$ct)<0.5)B.i=c(.5,.9)else{B.i=c(.3,.6)}
#startbio.list   <- list(Dusky=B.i.def.fn(Dus.catch), Gummy=B.i.def.fn(Gum.catch),
#                   Sandbar=B.i.def.fn(San.catch),Whiskery=B.i.def.fn(Whi.catch))

    #3.2.2. biomass range at after last catch 
B.f.def.fn=function(DAT)if(DAT$ct[nrow(DAT)]/max(DAT$ct)>0.5)B.f=c(.3,.7)else{B.f=c(.01,.4)}
# finalbio.list   <- list(Dusky=B.f.def.fn(Dus.catch), Gummy=B.f.def.fn(Gum.catch),
#                   Sandbar=B.f.def.fn(San.catch),Whiskery=B.f.def.fn(Whi.catch))



#NUMBER OF ITERATIONS (may no need as high for gummy and whiskery)
n           <- 1000000  


#PROCESS ERROR; 0 if deterministic model; 0.05 reasonable value? 0.2 is too high
sigR        <- 0.05   #source Eval Lai      


#REFERENCE POINTS
  #Quantiles approach
Threshold=.8
Th.low=0+((1-Threshold)/2)
Th.up=1-((1-Threshold)/2)
Limit=.9
Lim.low=0+((1-Limit)/2)
Lim.up=1-((1-Limit)/2)
Target=.5


  #Consistence with Bbrp
Ref.p=0.75



#---PROCEDURE SECTION-----

## FUNCTIONS
#Surplus production model for calculating biomass
.schaefer  <- function(theta)
{
  with(as.list(theta), {  ## for all combinations of ri & ki
    bt=vector()
    ell = 0  ## initialize ell
    for (j in startbt)
    {
      if(ell == 0) 
      {
        bt[1]=j*k*exp(rnorm(1,0, sigR))  ## set biomass in first year
        for(i in 1:nyr) ## for all years in the time series
        {
          xt=rnorm(1,0, sigR)
          bt[i+1]=(bt[i]+r*bt[i]*(1-bt[i]/k)-ct[i])*exp(xt) ## calculate biomass as function of previous year's biomass plus net production minus catch
        }
        
        #Bernoulli likelihood, assign 0 or 1 to each combination of r and k
        ell = 0
        if(bt[nyr+1]/k>=lam1 && bt[nyr+1]/k <=lam2 && min(bt) > 0 && max(bt) <=k && bt[which(yr==interyr)]/k>=interbio[1] && bt[which(yr==interyr)]/k<=interbio[2]) 
          ell = 1
      }	
    }
    return(list(ell=ell))
    
    
  })
}


#This function conducts the stock reduction analysis for N trials
#args:
#  theta - a list object containing:
#		r (lower and upper bounds for r)
#		k (lower and upper bounds for k)
#		lambda (limits for current depletion)
sraMSY	<-function(theta, N)
{
  with(as.list(theta), 
{
  ## get N values of r, assign to ri
      #lognormal r  
  if(user=="Yes")
    {
      logSD=sqrt(log(1+((CV*mean(r))/mean(r))^2))
      ri = rlnorm(N, log(mean(r)), logSD)
    } 
      #uniform r   
  if(user=="No") ri = exp(runif(N, log(r[1]), log(r[2])))  
  
  ## get N values between k[1] and k[2], assing to ki
  ki = exp(runif(N, log(k[1]), log(k[2])))  
  
  ## assign ri, ki, and final biomass range to itheta
  itheta=cbind(r=ri,k=ki, lam1=lambda[1],lam2=lambda[2], sigR=sigR) 
  
  ## call Schaefer function with parameters in itheta
  M = apply(itheta,1,.schaefer) 
  i=1:N
  ## prototype objective function
  get.ell=function(i) M[[i]]$ell
  ell = sapply(i, get.ell) 
  return(list(r=ri,k=ki, ell=ell))	
})
}



#---MAIN SECTION-----

  #Define if using specific priors(Yes) or default priors (No)
user="Yes"
#user="No"

setwd(handl_OneDrive("Analyses/Reference Points/Outputs/MSY_Catch"))

## Loop through stocks
stock_id <- unique(as.character(cdat$stock))
system.time(for(stock in stock_id)
  {
  p=which(stock==stock_id)
  
  yr   <- cdat$yr[as.character(cdat$stock)==stock]
  ct   <- as.numeric(cdat$ct[as.character(cdat$stock)==stock])
  res  <- unique(as.character(cdat$res[as.character(cdat$stock)==stock])) ## resilience 
  nyr  <- length(yr)    ## number of years in the time series
  
#  cat("\n","Stock",stock,"\n")
#  flush.console()
  
  ## PARAM SECTION
  #user defined pars
  if(user=="Yes")
  {
    start_r=r.list[[p]]
    start_k=k.list[[p]]
    startbio=startbio.list[[p]]
    finalbio=finalbio.list[[p]]
    
  }else
    #defaults
  {
    start_r  <- if(res == "Very low")             ## for unknown resilience 
    {c(0.015, 0.1)}else if
    (res == "Low")
    {c(0.05,0.5)}else if
    (res == "High")
    {c(0.6,1.5)}else
    {c(0.2,1)} 
    
    start_k     <- c(max(ct),100*max(ct))          ## for unknown upper K 
    
    startbio    <- if(ct[1]/max(ct) < 0.5) {c(0.5,0.9)} else {c(0.3,0.6)} # for unknown startbio
    finalbio    <- if(ct[nyr]/max(ct) > 0.5) {c(0.3,0.7)} else {c(0.01,0.4)} # for unknown finalbio 
    
  }

  interyr 	<- yr[2]   ## interim year within time series for which biomass estimate is available;
                       ## set to yr[2] if no estimates are available
  interbio 	<- c(0, 1) ## biomass range for interim year, as fraction of k;
                       ## set to 0 and 1 if not available
  
  startbt     <- seq(startbio[1], startbio[2], by = 0.05) ## apply range of start biomass in steps of 0.05	
  parbound <- list(r = start_r, k = start_k, lambda = finalbio, sigR)
  
  
  #print out some data and assumption (checks)
  cat("Last year =",max(yr),", last catch (tons)=",(Conv.fac/1000)*ct[nyr],"\n")
  cat("Resilience =",res,"\n")
  cat("Process error =", sigR,"\n")
  cat("Assumed initial biomass (B/k) =", startbio[1],"-", startbio[2], " k","\n")
  cat("Assumed intermediate biomass (B/k) in", interyr, " =", interbio[1],"-",interbio[2]," k","\n")
  cat("Assumed final biomass (B/k) =", parbound$lambda[1],"-",parbound$lambda[2]," k","\n")
  cat("Initial bounds for r =", parbound$r[1], "-", parbound$r[2],"\n")
  cat("Initial bounds for k (tons)=", format((Conv.fac/1000)*parbound$k[1], digits=3), "-", format((Conv.fac/1000)*parbound$k[2],digits=3),"\n")
  
  #flush.console()
  
  

  
  ## MAIN
  R1 = sraMSY(parbound, n)  
  
  
  ## Get statistics on r, k, MSY and determine new bounds for r and k
  r1 	<- R1$r[R1$ell==1]
  k1 	<- R1$k[R1$ell==1]
  msy1  <- r1*k1/4
  mean_msy1 <- exp(mean(log(msy1))) 
  
  if(user=="Yes") max_k1 <- max(k1[r1*k1/4<mean_msy1])                            #REVIEW!!!!!
  
  if(!user=="Yes")
  {
    max_k1a  <- min(k1[r1<1.1*parbound$r[1]]) ## smallest k1 near initial lower bound of r
    max_k1b  <- max(k1[r1*k1/4<mean_msy1]) ## largest k1 that gives mean MSY
    max_k1 <- if(max_k1a < max_k1b) {max_k1a} else {max_k1b}
  }
  
  if(length(r1)<10)
    {
      cat("Too few (", length(r1), ") possible r-k combinations, check input parameters","\n")
 #   flush.console()
    }
  
  if(length(r1)>=10)
    {
        ## set new upper bound of r to 1.2 max r1
    parbound$r[2] <- 1.2*max(r1)
    ## set new lower bound for k to 0.9 min k1 and upper bound to max_k1 
    parbound$k 	  <- c(0.9 * min(k1), max_k1)
    
    
    cat("First MSY =", format((Conv.fac/1000)*mean_msy1, digits=3),"\n")
    cat("First r =", format(exp(mean(log(r1))), digits=3),"\n")
    cat("New upper bound for r =", format(parbound$r[2],digits=2),"\n")	
    cat("New range for k (tons)=", format((Conv.fac/1000)*parbound$k[1], digits=3), "-", format((Conv.fac/1000)*parbound$k[2],digits=3),"\n")
    
    
    ## Repeat analysis with new r-k bounds
    R1 = sraMSY(parbound, n)
    
    ## Get statistics on r, k and msy
    r = R1$r[R1$ell==1]
    k = R1$k[R1$ell==1]
    k = k*(Conv.fac/1000) #convert to tonnes
    msy = r * k / 4
    mean_ln_msy = mean(log(msy))
    
    ct1=ct*(Conv.fac/1000) #convert back to tonnes
    
    ## plotting
    Mean.MSY=exp(mean(log(msy)))
    Limit.low.MSY=quantile(msy, Lim.low)
    Threshold.low.MSY=quantile(msy, Th.low)
    Limit.up.MSY=quantile(msy, Lim.up)
    Threshold.up.MSY=quantile(msy, Th.up)
    
    Thr.catch=Mean.MSY
    Lim.catch=Thr.catch*1/0.75
    Tar.catch=Thr.catch*0.75
    
    #outputs for plotting Figure 5
    if(p==1)  Store.Dusky=c(Mean.MSY,Limit.low.MSY,Threshold.low.MSY,Limit.up.MSY,Threshold.up.MSY)
    if(p==2)  Store.Gummy=c(Mean.MSY,Limit.low.MSY,Threshold.low.MSY,Limit.up.MSY,Threshold.up.MSY)
    if(p==3)  Store.Sandbar=c(Mean.MSY,Limit.low.MSY,Threshold.low.MSY,Limit.up.MSY,Threshold.up.MSY)
    if(p==4)  Store.Whiskery=c(Mean.MSY,Limit.low.MSY,Threshold.low.MSY,Limit.up.MSY,Threshold.up.MSY)

    
    

    tiff(file=paste("CatchMSY_Plots.",Sp[p],".tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
    par(mai=c(1,1.1,.1,.1),las=1,mgp=c(3,.5,0))
    plot(yr, ct1, type="l", ylim = c(0, max(ct1)), xlab = "Year", ylab = "Catch (t)", main = "",lwd=2,
         cex.lab=2,cex.axis=1.5)
    
    abline(h=Tar.catch,col="green", lwd=2.5)
    abline(h=Thr.catch,col="orange", lwd=2.5)
    abline(h=Lim.catch,col="red", lwd=2.5)
    
#     abline(h=exp(mean(log(msy))),col="green", lwd=2)
#     abline(h=Limit.low.MSY,col="red", lwd=2)
#     abline(h=Threshold.low.MSY,col="orange", lwd=2)
#     abline(h=Limit.up.MSY,col="red", lwd=2)
#     abline(h=Threshold.up.MSY,col="orange", lwd=2)
    legend("topright",c("","target","limit","threshold","",""),bty='n',col=c(NA,"green","red","orange",NA,NA),
           lty=1,lwd=2,cex=1.25)
    dev.off()

    tiff(file=paste("CatchMSY_Multi.",Sp[p],".tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
    par(mfcol=c(2,3))    
    plot(yr, ct1, type="l", ylim = c(0, max(ct1)), xlab = "Year", ylab = "Catch (t)", main = stock,lwd=2)
    abline(h=exp(mean(log(msy))),col="black", lwd=2)
    abline(h=exp(mean_ln_msy - 2 * sd(log(msy))),col="red")
    abline(h=exp(mean_ln_msy + 2 * sd(log(msy))),col="red")

    
    hist(r, freq=F, xlim=c(0, 1.2 * max(r)), main = "")
    abline(v=exp(mean(log(r))),col="black",lwd=2)
    abline(v=exp(mean(log(r))-2*sd(log(r))),col="red")
    abline(v=exp(mean(log(r))+2*sd(log(r))),col="red")
    
    plot(r1, k1*(Conv.fac/1000), xlim = start_r, ylim = start_k*(Conv.fac/1000), xlab="r", ylab="k (t)")
    
    hist(k, freq=F, xlim=c(0, 1.2 * max(k)), xlab="k (t)", main = "")
    abline(v=exp(mean(log(k))),col="black", lwd=2)	
    abline(v=exp(mean(log(k))-2*sd(log(k))),col="red")
    abline(v=exp(mean(log(k))+2*sd(log(k))),col="red")
    
    
    plot(log(r), log(k),xlab="ln(r)",ylab="ln(k)")
    abline(v=mean(log(r)))
    abline(h=mean(log(k)))
    abline(mean(log(msy))+log(4),-1, col="black",lwd=2)
    abline(mean(log(msy))-2*sd(log(msy))+log(4),-1, col="red")
    abline(mean(log(msy))+2*sd(log(msy))+log(4),-1, col="red")
    
    hist(msy, freq=F, xlim=c(0, 1.2 * max(msy)), xlab="MSY (t)",main = "")
    abline(v=exp(mean(log(msy))),col="black", lwd=2)
    abline(v=exp(mean_ln_msy - 2 * sd(log(msy))),col="red")
    abline(v=exp(mean_ln_msy + 2 * sd(log(msy))),col="red")
    dev.off()
    
    cat("Possible combinations = ", length(r),"\n")
    cat("geom. mean r =", format(exp(mean(log(r))),digits=3), "\n")
    cat("r +/- 2 SD =", format(exp(mean(log(r))-2*sd(log(r))),digits=3),"-",format(exp(mean(log(r))+2*sd(log(r))),digits=3), "\n")
    cat("geom. mean k (tons)=", format(exp(mean(log(k))),digits=3), "\n")
    cat("k +/- 2 SD (tons)=", format(exp(mean(log(k))-2*sd(log(k))),digits=3),"-",format(exp(mean(log(k))+2*sd(log(k))),digits=3), "\n")
    cat("geom. mean MSY (tons)=", format(exp(mean(log(msy))),digits=3),"\n")
    cat("MSY +/- 2 SD (tons)=", format(exp(mean_ln_msy - 2 * sd(log(msy))),digits=3), "-", format(exp(mean_ln_msy + 2 * sd(log(msy))),digits=3), "\n")
    
    
    ## Export results
    outfile  <- paste("CatchMSY_Output.",Sp[p],".csv",sep="")
    output = data.frame(Species=stock, Proc.Err=sigR, B.init.low=startbio[1], B.init.high=startbio[2], 
            B.2.low=interbio[1], B.2.high=interbio[2], 
            B.fin.low=finalbio[1], B.fin.high=finalbio[2], Min.yr=min(yr), Max.yr=max(yr), 
            Resilience=res, Catch.max=max(ct), Catch.init=ct[1], Catch.fin=ct[nyr], 
            Lenth.r=length(r), Mean.r=exp(mean(log(r))), SD.r=sd(log(r)), Min.r=min(r), CI.95.low.r=quantile(r,0.025), 
            Median.r=median(r), CI.95.up.r=quantile(r,0.975), Max.r=max(r), 
            Mean.k=exp(mean(log(k))), SD.k=sd(log(k)), Min.k=min(k), CI.95.low.k=quantile(k, 0.025), 
            Median.k=median(k), CI.95.up.k=quantile(k, 0.975),Max.k=max(k),
            Median.MSY=median(msy),Mean.MSY=Mean.MSY,SD.MSY=sd(log(msy)),Min.MSY=min(msy),Max.MSY=max(msy),
            CI.95.low.MSY=quantile(msy,0.025),CI.95.up.MSY=quantile(msy,0.975)) 
    
    
        #add simple catch-based assessment as per Froese et al 2012
    Status="whatever"
    if(ct[nyr]/max(ct) >.5) Status="Fully exploited"
    if(ct[nyr]/max(ct)>=0.1 & ct[nyr]/max(ct)<=0.5) Status="Over exploited"
    if(ct[nyr]/max(ct)<0.1) Status="Collapsed"
    output$Froese.Status=Status
    
    
    write.table(output, file = outfile, sep = ",",row.names = FALSE, col.names =T)
    
  }
})  



#Paper Figure 5
#Catch-based

fn.plotFig.5=function(OUT,DAT,SPEC)
{
  Thr.catch=OUT[1] 
  Lim.catch=Thr.catch*1/Ref.p
  Tar.catch=Thr.catch*Ref.p
  N=length(DAT$ct)
  yy=match(DAT[1,1],yr)
  
  max1=max(DAT$ct*Conv.fac/1000,na.rm=T)*1.05
  plot(1:N,(DAT$ct*Conv.fac/1000),ylab="",xlab="",xaxt="n",las=2,type="o",pch=19,cex=1.25,ylim=c(0,max1)
       ,cex.axis=1.2,lwd=1.75)
  axis(1,at=1:N,labels=F,tck=-0.015)
  axis(1,at=seq(1,N,5),labels=FinYear[seq(yy,length(yr),5)],tck=-0.03,cex.axis=1.25)
  abline(h=Tar.catch,col=Col.Tar, lwd=2.5,lty=5)
  abline(h=Thr.catch,col=Col.Thr, lwd=2.5)
  abline(h=Lim.catch,col=Col.Lim, lwd=2.5)
  if(SPEC%in%c("whiskery shark","dusky shark"))legend('topright',SPEC,bty='n',cex=1.5)
  if(SPEC%in%c("gummy shark","sandbar shark"))legend('topleft',SPEC,bty='n',cex=1.5)
  
}

tiff(file=handl_OneDrive("Analyses/Reference Points/Outputs/Figure5.tiff"),width = 2400,
     height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff

par(mfcol=c(2,2),mai=c(.6,.6,.1,.1),oma=c(.01,.1,.01,.1),mgp=c(1,.6,0))
Col.Tar="grey40"
Col.Thr="grey60"
Col.Lim="black"
fn.plotFig.5(Store.Gummy,Gum.catch,"gummy shark")
legend("bottomright",c("","Limit","Threshold","Target"),bty='n',col=c(NA,Col.Lim,Col.Thr,Col.Tar),
       lty=c(1,1,1,5),lwd=2.5,cex=1.4)
fn.plotFig.5(Store.Whiskery,Whi.catch,"whiskery shark")
fn.plotFig.5(Store.Sandbar,San.catch,"sandbar shark")
fn.plotFig.5(Store.Dusky,Dus.catch,"dusky shark")

mtext("Financial year",1,-1.5,outer=T,cex=1.45)
mtext("Catch (t)",2,-1.5,outer=T,las=3,cex=1.45)
dev.off()




#not used
tiff(file=handl_OneDrive("Analyses/Reference Points/Outputs/Figure5.tiff"),width = 2400,
     height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff

par(mfcol=c(2,1),mai=c(.8,.8,.01,.1),oma=c(.1,.1,.1,.1),mgp=c(1,.6,0))

# Col.Tar="green"
# Col.Thr="orange"
# Col.Lim="red"
Col.Tar="grey40"
Col.Thr="grey60"
Col.Lim="black"


#5.a CPUE-based
Whisk.cpue=read.csv(handl_OneDrive("Analyses/Reference Points/CPUE.whiskery.csv"))

Binit=0.95
Tar.cpue=0.461*Whisk.cpue[1,2]/Binit    #update from Reference Points
Lim.cpue=0.242*Whisk.cpue[1,2]/Binit
Thr.cpue=0.354*Whisk.cpue[1,2]/Binit

max1=max(Whisk.cpue[,2],na.rm=T)*1.05
FInYEAR=as.character(unique(Whisk.cpue$FINYEAR))
N=length(FInYEAR)
plot(1:N,Whisk.cpue[,2],ylab="",xlab="",xaxt="n",las=2,type="o",pch=19,cex=1.25,ylim=c(0,max1)
     ,cex.axis=1.2,lwd=1.75)
axis(1,at=1:N,labels=F,tck=-0.015)
axis(1,at=seq(1,N,5),labels=F,tck=-0.03,cex.axis=1.1)
#mtext("Year",1,-2,outer=T,cex=1.35)
mtext("CPUE (kg/km gn.d)",2,2.5,outer=F,las=3,cex=1.35)
abline(h=Tar.cpue,col=Col.Tar, lwd=2.5, lty=5)
abline(h=Lim.cpue,col=Col.Lim, lwd=2.5)
abline(h=Thr.cpue,col=Col.Thr, lwd=2.5)
legend("topright",c("","Limit","Threshold","Target"),bty='n',col=c(NA,Col.Lim,Col.Thr,Col.Tar),
       lty=c(1,1,1,5),lwd=2.5,cex=1.15)
legend("topleft",'(a)',bty='n',cex=1.25)

#5.b Catch-based
Thr.catch=Store.Dusky[1] 
# Low.Lim.catch=Store.Dusky[2]
# Up.Lim.catch=Store.Dusky[4]
# Low.Thr.catch=Store.Dusky[3]
# Up.Thr.catch=Store.Dusky[5]
Lim.catch=Thr.catch*1/Ref.p
Tar.catch=Thr.catch*Ref.p

max1=max(Dus.catch$ct*Conv.fac/1000,na.rm=T)*1.05
plot(1:N,(Dus.catch$ct*Conv.fac/1000),ylab="",xlab="",xaxt="n",las=2,type="o",pch=19,cex=1.25,ylim=c(0,max1)
     ,cex.axis=1.2,lwd=1.75)
axis(1,at=1:N,labels=F,tck=-0.015)
axis(1,at=seq(1,N,5),labels=FInYEAR[seq(1,N,5)],tck=-0.03,cex.axis=1.1)
mtext("Financial year",1,-2,outer=T,cex=1.35)
mtext("Catch (t)",2,2.5,outer=F,las=3,cex=1.35)
abline(h=Tar.catch,col=Col.Tar, lwd=2.5,lty=5)
abline(h=Thr.catch,col=Col.Thr, lwd=2.5)
abline(h=Lim.catch,col=Col.Lim, lwd=2.5)
#abline(h=Low.Lim.catch,col=Col.Lim, lwd=2.5)
#abline(h=Up.Lim.catch,col=Col.Lim, lwd=2.5)
#abline(h=Low.Thr.catch,col=Col.Thr, lwd=2.5)
#abline(h=Up.Thr.catch,col=Col.Thr, lwd=2.5)
#legend("topright",c("","target","limit","threshold","",""),bty='n',col=c(NA,Col.Tar,Col.Lim,Col.Thr,NA,NA),
#       lty=1,lwd=2,cex=1.25)
legend('topleft','(b)',bty='n',cex=1.25)
dev.off()
