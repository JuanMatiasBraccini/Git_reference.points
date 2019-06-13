# #Selectivity
# #whiskery                #Simpfendorfer & Unsworth 1998
# alpha=64.01339
# beta=18.53164
# alphabeta=alpha.w*beta.w
# age=0:20
# ext=seq(130,200,10)
# Linf=120.7
# K=0.369
# to=-0.544
# 
# 
# 
# #Dusky
# theta1=130.13
# theta2=29237
# mesh= 6.5 #(in inches, = 16.5 cm)
# alphabeta=theta1*mesh
# beta=-0.5*(theta1*mesh-((theta1^2)*(mesh^2)+4*theta2)^0.5)
# alpha=alphabeta/beta
# 
# age=0:55
# ext=seq(130,200,10)
# Linf=374.4
# K=.0367
# to=-3.3


#sandbar
theta1=135.58
theta2=117001
mesh= 6.5
alphabeta=theta1*mesh
beta=-0.5*(theta1*mesh-((theta1^2)*(mesh^2)+4*theta2)^0.5)
alpha=alphabeta/beta
age=0:39
ext=seq(130,200,10)
Linf=244.2
K=.040
to=-4.8

Lo=42.5 #(McAuley e tal 2005)

mid.FL.fem=Linf*(1-exp(-K*(age+0.5-(-to))))
sel.fem=((mid.FL.fem*10/alphabeta)^alpha)*(exp(alpha-(mid.FL.fem*10/beta)))

#Modified vonB (Simpfendorfer et al 2000)
mid.FL.fem.mod=Lo+(Linf-Lo)*(1-exp(-K*age))

plot(age,mid.FL.fem,pch=19,cex=1.5)
plot(mid.FL.fem,sel.fem,xlim=c(0,200),pch=19,cex=1.5)
#if species ==whiskery
#sel.fem.ext=((c(mid.FL.fem,ext)*10/alphabeta)^alpha)*(exp(alpha-(c(mid.FL.fem,ext)*10/beta)))
#lines(c(mid.FL.fem,ext),sel.fem.ext,col=2,lwd=2)
plot(age,sel.fem,pch=19,cex=1.5)