#HARVEST STRATEGIES AND REFERENCE POINTS EXPLORATION

library(lattice)

#Compare Blim DoF and NMFS

limit.Dof=.3

Min.M=exp(1.44-0.982*log(7))
Max.M=exp(1.44-0.982*log(80))

M=seq(Max.M,Min.M,by=.05)
BMSY.conv=c(.45,.5,.6)

Mat=expand.grid(M,BMSY.conv)
Mat=as.data.frame(Mat)
names(Mat)=c("M","BMSY.conv")
B.unfished=1

Mat$Dof.limit=B.unfished*limit.Dof
Mat$NMFS.limit=(1-Mat$M)*(B.unfished*Mat$BMSY.conv)
Mat=Mat[order(Mat$BMSY.conv),]

plot(1:nrow(Mat),Mat$Dof.limit,ylim=c(0,.6),xlab="Case #", ylab="Limit ref level (Bunfished=1)",las=1)
points(1:nrow(Mat),Mat$NMFS.limit,col=2,pch=21)
points(1:nrow(Mat),Mat$NMFS.limit,col=2,pch=Mat$BMSY.conv)



