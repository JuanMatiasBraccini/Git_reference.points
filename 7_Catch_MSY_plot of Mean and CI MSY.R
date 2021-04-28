
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


Mean.MSY=c(280.8069921,475.4528215,226.8289514,253.5824548)
CI.95.low.MSY=c(149.6962406,365.4134303,72.84143994,187.0280011)
CI.95.up.MSY=c(538.8756541,630.9149086,732.8367652,392.7608434)

names(Mean.MSY)=names(CI.95.low.MSY)=names(CI.95.up.MSY)=c("dusky shark","gummy shark","sandbar shark","whiskery shark")



yr=as.numeric(sapply(strsplit(as.character(Dus.catch$finyear),"-"), "[", 1))
FinYear=as.character(Dus.catch$finyear)

Dus.catch$finyear=Gum.catch$finyear=San.catch$finyear=Whi.catch$finyear=yr

Conv.fac=1e6    #convert kg to 1000s tons

fn.plotFig.5=function(OUT,OUT.L,OUT.U,DAT,SPEC)
{
  Thr.catch=OUT 
  Lim.catch=OUT.L
  Tar.catch=OUT.U
  N=length(DAT$Total.catch)
  yy=match(DAT[1,1],yr)
  
  max1=max(DAT$Total.catch/1000,na.rm=T)*1.05
  plot(1:N,(DAT$Total.catch/1000),ylab="",xlab="",xaxt="n",las=2,type="o",pch=19,cex=1.25,ylim=c(0,max1)
       ,cex.axis=1.2,lwd=1.75)
  axis(1,at=1:N,labels=F,tck=-0.015)
  axis(1,at=seq(1,N,5),labels=FinYear[seq(yy,length(yr),5)],tck=-0.03,cex.axis=1.25)
  abline(h=Tar.catch,col=2, lwd=2.5)
  abline(h=Thr.catch,col=1, lwd=2.5)
  abline(h=Lim.catch,col=2, lwd=2.5)
  if(SPEC%in%c("whiskery shark","dusky shark"))legend('topright',SPEC,bty='n',cex=1.5)
  if(SPEC%in%c("gummy shark","sandbar shark"))legend('topleft',SPEC,bty='n',cex=1.5)
  
}

tiff(file=handl_OneDrive("Analyses/Reference Points/Outputs/Figure5_Mean.CI.MSY.tiff"),width = 2400,
     height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff

par(mfcol=c(2,2),mai=c(.6,.6,.1,.1),oma=c(.01,.1,.01,.1),mgp=c(1,.6,0))
Col.Tar="grey40"
Col.Thr="grey60"
Col.Lim="black"
fn.plotFig.5(Mean.MSY[2],CI.95.low.MSY[2],CI.95.up.MSY[2],Gum.catch,"gummy shark")
fn.plotFig.5(Mean.MSY[4],CI.95.low.MSY[4],CI.95.up.MSY[4],Whi.catch,"whiskery shark")
fn.plotFig.5(Mean.MSY[3],CI.95.low.MSY[3],CI.95.up.MSY[3],San.catch,"sandbar shark")
fn.plotFig.5(Mean.MSY[1],CI.95.low.MSY[1],CI.95.up.MSY[1],Dus.catch,"dusky shark")

mtext("Financial year",1,-1.5,outer=T,cex=1.45)
mtext("Catch (t)",2,-1.5,outer=T,las=3,cex=1.45)
dev.off()


