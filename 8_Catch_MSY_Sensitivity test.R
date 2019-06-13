#1. Test effects of sigma, prior definitions and K.upper level


SIGR=c(0,0.05,0.1)
USER=c("No","Yes")
K.UPPER=c(10,50,100)



for (w in 1: length(USER))
{
  user=USER[w]
  for (qq in 1:length(SIGR))
  {
    sigR=SIGR[qq]
    for (x in 1:length(K.UPPER))
    {
      K.upper=K.UPPER[x]
      setwd(paste("C:/Matias/Analyses/Reference Points/Outputs/MSY_Catch/sensitivity tests/NEW/",
                  K.upper,"/",sigR,"/",user,sep=""))
      for(stock in stock_id)
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
        
        interyr   <- yr[2]   ## interim year within time series for which biomass estimate is available;
        ## set to yr[2] if no estimates are available
        interbio   <- c(0, 1) ## biomass range for interim year, as fraction of k;
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
          
          #output Dusky for plotting e.g.
          if(p==1)  Store.Dusky=c(Mean.MSY,Limit.low.MSY,Threshold.low.MSY,Limit.up.MSY,Threshold.up.MSY)
          
          
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
          
          
          
          
        }
      }
    }
  }
}



#2. Test effect of catch time series
#Dusky shark (one default, one swapped to gummy catch series )
SIGR=c(0)
USER=c("No","Yes")
K.UPPER=50
Sp.sens=c("Dusky.dusky.ct","Dusky.gummy.ct")

Sens.dat=cdat[1:(37+37),] #select dusky and gummy data
stock_id.sens=c("Dusky" ,   "Gummy")
Sens.dat$res="Very low"

for (w in 1: length(USER))
{
  user=USER[w]
  for (qq in 1:length(SIGR))
  {
    sigR=SIGR[qq]
    for (x in 1:length(K.UPPER))
    {
      K.upper=K.UPPER[x]
      setwd(paste("C:/Matias/Analyses/Reference Points/Outputs/MSY_Catch/sensitivity tests/NEW/Catch.effect/dusky/k.50/sigma.0/",
                  user,sep=""))
      for(stock in stock_id.sens)
      {
        p=which(stock==stock_id.sens)
        
        yr   <- Sens.dat$yr[as.character(Sens.dat$stock)==stock]
        ct   <- as.numeric(Sens.dat$ct[as.character(Sens.dat$stock)==stock])
        res  <- unique(as.character(Sens.dat$res[as.character(Sens.dat$stock)==stock])) ## resilience 
        nyr  <- length(yr)    ## number of years in the time series
        
        #  cat("\n","Stock",stock,"\n")
        #  flush.console()
        
        ## PARAM SECTION
        #user defined pars
        if(user=="Yes")
        {
          start_r=r.list[[1]]
          start_k=k.list[[p]]
          startbio=startbio.list[[1]]
          finalbio=finalbio.list[[1]]
          
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
        
        interyr   <- yr[2]   ## interim year within time series for which biomass estimate is available;
        ## set to yr[2] if no estimates are available
        interbio   <- c(0, 1) ## biomass range for interim year, as fraction of k;
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
          
          #output Dusky for plotting e.g.
          if(p==1)  Store.Dusky=c(Mean.MSY,Limit.low.MSY,Threshold.low.MSY,Limit.up.MSY,Threshold.up.MSY)
          
          
          tiff(file=paste("CatchMSY_Plots.",Sp.sens[p],".tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
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
          
          tiff(file=paste("CatchMSY_Multi.",Sp.sens[p],".tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
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
          outfile  <- paste("CatchMSY_Output.",Sp.sens[p],".csv",sep="")
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
          
          
          
          
        }
      }
    }
  }
}




#Gummy shark shark (one default, one series increasing, one decreasing)
SIGR=c(0)
USER=c("No","Yes")
K.UPPER=50
Sp.sens=c("Gummy.gummy.ct","Gummy.increase.ct","Gummy.decrease.ct")

#ACA
Sens.dat=cdat[38:(37+37),] #select dusky and gummy data
MIN.gum=min(Sens.dat$ct)
MAX.gum=max(Sens.dat$ct)

Sens.dat=rbind(Sens.dat,Sens.dat,Sens.dat)
Sens.dat$stock=c(rep("Gummy",37),rep("Gummy.inc",37),rep("Gummy.dec",37))
Sens.dat$ct[38:74]=jitter(seq(MIN.gum,MAX.gum,length.out=37),10)
Sens.dat$ct[75:111]=jitter(seq(MAX.gum,MIN.gum,length.out=37),10)


stock_id.sens=c("Gummy" ,   "Gummy.inc","Gummy.dec")



for (w in 1: length(USER))
{
  user=USER[w]
  for (qq in 1:length(SIGR))
  {
    sigR=SIGR[qq]
    for (x in 1:length(K.UPPER))
    {
      K.upper=K.UPPER[x]
      setwd(paste("C:/Matias/Analyses/Reference Points/Outputs/MSY_Catch/sensitivity tests/NEW/Catch.effect/gummy/k.50/sigma.0/",
                  user,sep=""))
      for(stock in stock_id.sens)
      {
        p=which(stock==stock_id.sens)
        
        yr   <- Sens.dat$yr[as.character(Sens.dat$stock)==stock]
        ct   <- as.numeric(Sens.dat$ct[as.character(Sens.dat$stock)==stock])
        res  <- unique(as.character(Sens.dat$res[as.character(Sens.dat$stock)==stock])) ## resilience 
        nyr  <- length(yr)    ## number of years in the time series
        
        #  cat("\n","Stock",stock,"\n")
        #  flush.console()
        
        ## PARAM SECTION
        #user defined pars
        if(user=="Yes")
        {
          start_r=r.list[[2]]
          start_k=k.list[[2]]
          startbio=startbio.list[[2]]
          finalbio=finalbio.list[[2]]
          
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
        
        interyr   <- yr[2]   ## interim year within time series for which biomass estimate is available;
        ## set to yr[2] if no estimates are available
        interbio   <- c(0, 1) ## biomass range for interim year, as fraction of k;
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
        r1   <- R1$r[R1$ell==1]
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
          
          #output Dusky for plotting e.g.
          if(p==1)  Store.Dusky=c(Mean.MSY,Limit.low.MSY,Threshold.low.MSY,Limit.up.MSY,Threshold.up.MSY)
          
          
          tiff(file=paste("CatchMSY_Plots.",Sp.sens[p],".tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
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
          
          tiff(file=paste("CatchMSY_Multi.",Sp.sens[p],".tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
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
          outfile  <- paste("CatchMSY_Output.",Sp.sens[p],".csv",sep="")
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
          
          
          
          
        }
      }
    }
  }
}