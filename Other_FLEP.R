  ######################### FRACTION LIFETIME EGG PRODUCTION ########################
  
#note: This script calculates the Fraction of Lifetime Egg Production (FLEP, Chapple and Botsford 2012) as
  #     a biological reference point when compared to F60%
  #     It also calculates the necessary Fishing mortality for maintaining FLEP at F60%   (MISSING OPTIMISATION 
  #     FOR FINDING F60 AND F FRACTION OF M, SEE EXCEL!!!!)
  
  
#assumptions: knife edge selectivity starting at age at 50% maturity

  
  
  

#---PARAMETERS SECTION----
species.names <- c("Whiskery", "Dusky")  
species.list <- vector("list", length(species.names))
names(species.list) <- species.names
  
pars.names <- c("max.age", "age.mat")
pars.list <- vector("list", length(pars.names))
names(pars.list) <- pars.names
  
#Whiskery shark
species.list$Whiskery=pars.list
species.list$Whiskery$max.age=15      
species.list$Whiskery$age.mat=6


#Dusky shark
species.list$Dusky=pars.list
species.list$Dusky$max.age=40      
species.list$Dusky$age.mat=21
  
  
#Reference Point
F60=0.6
  
  
  
#---PROCEDURE SECTION----
  
#1. Spawning potential ratio
#dummies
max.age=age.mat=1

theta=c(max.age,age.mat)

FLEP.fn=function(max.age,age.mat)
{
  #Hoenig's natural mortality
  M=exp(1.44-0.982*log(max.age))
  
  f=1*M
  
  #FLEP
  flep=(((exp(-(M+f)*max.age))-(exp(-(M+f)*age.mat)))/(M+f))/((exp(-M*max.age)-exp(-M*age.mat))/(M))
  
  #F60%
  Status=ifelse(flep>=F60,"OK","BAD")
  
  delta=abs(flep-F60)

  return(list(flep=flep,Status=Status,delta=delta))
}
  
  
#---MAIN SECTION----  
  
Whiskery=FLEP.fn(species.list$Whiskery$max.age,species.list$Whiskery$age.mat)
  
Dusky=FLEP.fn(species.list$Dusky$max.age,species.list$Dusky$age.mat)  

  
  