setwd("~/Activité Professionnelle/LIMNO 2019-2023/Modeling/Indirect Evolutionary Processes")

rm(list=ls())

library(actuar)
library(cowplot)
library(data.table)
library(deSolve)
library(directlabels)
library(dplyr)
library(dynlm)
library(emmeans)
library(ggalt)
library(ggpattern)
library(ggplot2)
library(ggrepel)
library(grid)
library(gridExtra)
library(lme4)
library(magrittr)
library(nlme)
library(plotly)
library(plyr)
library(pracma)
library(quantmod)
library(rlist)
library(reshape2)
library(scales)
library(tidyr)
library(zoo)

#############################################################################
#############################################################################
##### PREDATOR INDIRECT EVOLUTIONARY PROCESSES FOR PREDATOR-PREY SYSTEM #####
#############################################################################
#############################################################################

# Import the datasets
DataC=read.table("~/Activité Professionnelle/LIMNO 2019-2023/Modeling/Data_CEB.txt", h=T, dec=",")
summary(DataC)
names(DataC)

# Specify the variables as numeric or factor
DataC[,c(2:3)] %<>% mutate_if(is.character,as.numeric)

# Match conversion efficiency
Bead=as.list(round(seq(0.0,0.05,length.out=11),3))
Conv=lapply(Bead, function(x) {DataC %>% group_by(Strain) %>% dplyr::slice(which.min(abs(DensBP-x)))})
DataC2=as.data.frame(do.call("rbind",(Conv)))
DataC2$DensBM=as.character(rep(Bead, each=6))

# Extract conversion efficiency
Conv=as.numeric(round(subset(DataC2, Strain=="CR6")[,3],2))

# Parameters
a=0.031         # Predator attack rate on species
h=35.99         # Predator handling time on species
Ca=0.50         # Predator fraction of time for attacking
Xa=0.0027       # Prey conversion efficiency
Xp=170          # Predator conversion efficiency

Param=list(a,h,Ca,Xa,Xp)

# Simulation time
Time=seq(0,500,by=1)

# Range of nitrogen concentration
Nutri=200

# Range of dilution rate
Dilu=0.8

# Range of ingestion defense values
Ing1=seq(0.1,1.0,length.out=19)

# Range of competitiveness cost values
Comp1=round(seq(5.0,1.0,length.out=19),2)

# Range of conversion efficiency values
Conv=round(Conv,0)

# Range of beads densities
Bead=round(seq(0.0,0.05,length.out=11),3)

# Beads interference coefficient 
c=9.98

# Plot of trade-off shape
DataTO=data.frame(Ing1,Comp1)
ggplot(DataTO, aes(Comp1, Ing1)) +
  geom_smooth(method="loess", color="black", size=1, span=0.5, se=F) +
  ylab(expression('A'[1]~'ingestion probability'~'('*p[A~1]*')')) + 
  xlab(expression('A'[1]~'half-saturation constant'~'('*K[A~1]*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(y) sprintf("%.1f", y), breaks=seq(0.1,1.0,by=0.1), limits=c(0.1,1.0)) +
  scale_x_continuous(labels=function(y) sprintf("%.1f", y), breaks=seq(1.0,5.0,by=0.5), limits=c(1.0,5.0)) +
  theme(panel.background=element_blank(), panel.border=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.line=element_line(colour="black", linetype="solid", size=0.5))


#################################################
### ONE VARIABLE PRE-DEFENDED PREY WITH BEADS ###
#################################################

# Number of simulations
Nsim=209

# Starting values
Start=c(N=200, A1=0.2, P=1.0)

#######################################
### Ordinary differential equations ###
#######################################

# Iterate the function for each combination
System=function (Time, Start, Param) {
  with(as.list(c(Start,Param)),{
    
    if (Time >= 100) {B = B} 
    else if (Time < 100) {B = 0}
    
    if (Time >= 100) {Xp = Xp} 
    else if (Time < 100) {Xp = 143}
    
    Q1=1.0
    B1=1.5
    b=1.05
    
    G1 = B1*(N/(K1 + N))                                             
    F1 = (a*P1*Q1*A1^b)/(1 + a*P1*(Ca*h + Q1*(1 - Ca)*h)*A1^b + c*B)
    
    dN = D*(Ni - N) - (1/Xa)*G1*A1
    dA1 = G1*A1 - F1*P - D*A1
    dP = Xp*F1*P - D*P
    return(list(c(dN,dA1,dP)))})
}

# Create saving objects
ListODE=list()
Save=data.frame(Time=integer(0), N=integer(0), A1=integer(0), P=integer(0), Ni=integer(0), D=integer(0), P1=integer(0), B=integer(0))
SaveODE1=data.frame(Time=integer(0), N=integer(0), A1=integer(0), P=integer(0), Ni=integer(0), D=integer(0), P1=integer(0), B=integer(0))

# Calculate the dynamics for each combination
for(i in 1:length(Nutri)) {
  for(j in 1:length(Dilu)) {
    for(k in 1:length(Ing1)) {
      for(l in 1:length(Bead)) {
        P1=Ing1[k]
        K1=Comp1[k]
        B=Bead[l]
        Xp=Conv[l]
        if(Ing1[k]) Comp1[k]
        if(Bead[l]) Conv[l]
        SolveODE=function(time, y=c(N=Nutri[i], A1=0.2, P=1.0), parms=c(Ni=Nutri[i], D=Dilu[j], a=a, h=h, Ca=Ca, Xa=Xa)) {with(as.list(y), {y[which(y < 10^-20)]=0; return(y)})}
        ListODE=ode(y=c(N=Nutri[i], A1=0.2, P=1.0), Time, System, parms=c(Ni=Nutri[i], D=Dilu[j], a=a, h=h, Ca=Ca, Xa=Xa), method="lsoda", events=list(func=SolveODE, time=Time))
        ListODE[ListODE<0]=0
        Save=data.frame(Time=ListODE[,1],
                        N=ListODE[,2],
                        A1=ListODE[,3],
                        P=ListODE[,4],
                        Ni=rep(Nutri[i],length(Time)),
                        D=rep(Dilu[j],length(Time)),
                        P1=rep(Ing1[k],length(Time)),
                        B=rep(Bead[l],length(Time)))
        SaveODE1=rbind(SaveODE1,Save)
      }
    } 
  }
}


##################################
### Population dynamics panels ###
##################################

# Name columns according to combinations
Comb=unique(SaveODE1[c("P1","B")])
colnames(Comb)=c("Var1","Var2")
Comb=cbind(Comb,Var3=rep("P1",Nsim),Var4=rep("B",Nsim))
Names=paste(paste(Comb$Var3,Comb$Var1,sep="-"),paste(Comb$Var4,Comb$Var2,sep="-"),sep=" / ")

# Create a dataset
MeltODE1=SaveODE1[,-c(2)]
MeltODE1$A1=MeltODE1$A1*10^6
MeltODE1=melt(MeltODE1, id.vars=c("Time","Ni","D","P1","B"), variable.name="Organism", value.name="Density")
MeltODE1=MeltODE1[order(-MeltODE1$P1),]
MeltODE1$P1=factor(MeltODE1$P1, levels=unique(MeltODE1$P1))

# Rescale densities
Plot1ODE1=MeltODE1[,c(1:7)]
Plot1ODE1[Plot1ODE1$Organism=="A1",]$Density=log(MeltODE1[MeltODE1$Organism=="A1",]$Density+1)/3
Plot1ODE1[Plot1ODE1$Organism=="P",]$Density=log(MeltODE1[MeltODE1$Organism=="P",]$Density+1)

# Calculate effective densities
Plot2ODE1=MeltODE1[,c(1:7)]
Plot2ODE1$Organism=gsub("A1", "A", Plot2ODE1$Organism)
DensityA1=MeltODE1[MeltODE1$Organism=="A1",]$Density
DefenseA1=as.numeric(MeltODE1[MeltODE1$Organism=="A1",]$P1)
Plot2ODE1[Plot2ODE1$Organism=="A",]$Density=(log(DensityA1+1)*DefenseA1)/125
Plot2ODE1[Plot2ODE1$Organism=="P",]$Density=log(MeltODE1[MeltODE1$Organism=="P",]$Density+1)

# Represent the population dynamics
tiff('Panel Population Dynamics No Evolution 1.tiff', units="in", width=12, height=20, res=1000)
ggplot(Plot1ODE1, aes(Time, Density, group=Organism)) +
  geom_line(aes(color=Organism, linetype=Organism), size=0.5) +
  ylab(expression('Predator density'~'('*ln~individuals~mL^-1*')')) + xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=14)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=14)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=14)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=14)) +
  theme(axis.title.y.left=element_text(angle=90), axis.title.y.right=element_text(angle=90)) +
  scale_y_continuous(labels=function(y) sprintf("%.0f", y), breaks=seq(0,6,by=3), limits=c(0,6),
  sec.axis=sec_axis(~./5*5, expression('Prey density'~'('*ln~cells~mL^-1*')'), labels=sprintf(seq(0,14,by=7), fmt="%.0f"), breaks=seq(0,6,by=3))) +
  scale_x_continuous(labels=function(y) sprintf("%.0f", y), breaks=seq(0,200,by=100), limits=c(0,200)) +
  theme(panel.background=element_blank(), panel.border=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.line=element_line(colour="black", linetype="solid", size=0.5)) +
  scale_color_manual(values=alpha(c("A1"="royalblue2","P"="firebrick3"),0.8)) +
  scale_linetype_manual(values=c("A1"="solid","P"="solid")) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color="black", size=0.9) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color="black", size=0.9) +
  facet_wrap(P1~B, ncol=11, nrow=19) +
  theme(legend.position="none")
dev.off()

tiff('Panel Population Dynamics No Evolution 2.tiff', units="in", width=15, height=25, res=1000)
ggplot(Plot2ODE1, aes(Time, Density, group=Organism)) +
  geom_line(aes(color=Organism, linetype=Organism), size=0.5) +
  ylab(expression('Predator density'~'('*ln~individuals~mL^-1*')')) + xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=14)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=14)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=14)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=14)) +
  theme(axis.title.y.left=element_text(angle=90), axis.title.y.right=element_text(angle=90)) +
  scale_y_continuous(labels=function(y) sprintf("%.0f", y), breaks=seq(0,4,by=1), limits=c(0,4),
  sec.axis=sec_axis(~.*5/4, expression(italic('C. reinhardtii')~'effective density'~'('*10^3*')'), labels=function(y) sprintf("%.0f", y), breaks=seq(0,5,by=1))) +
  scale_x_continuous(labels=function(y) sprintf("%.0f", y), breaks=seq(0,200,by=100), limits=c(0,200)) +
  theme(panel.background=element_blank(), panel.border=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.line=element_line(colour="black", linetype="solid", size=0.5)) +
  scale_color_manual(values=alpha(c("A"="royalblue2","P"="firebrick3"),0.8)) +
  scale_linetype_manual(values=c("A"="solid","R"="solid")) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color="black", size=0.9) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color="black", size=0.9) +
  facet_wrap(P1~B, ncol=11, nrow=19) +
  theme(legend.position="none")
dev.off()


###############################################
### Eco-evolutionary states of the dynamics ###
###############################################

# Define phases for beads
SaveODE1$Phase=ifelse(SaveODE1$Time < 100, "Pre", "Post")
SaveODE1$Phase=factor(SaveODE1$Phase, levels=unique(SaveODE1$Phase))

# Calculate the difference in ingestion probabilities
SaveODE1$P0=as.character(round(1-as.numeric(SaveODE1$P1),2))

# Calculate density boundaries
CutODE1=SaveODE1[SaveODE1$Time %in% c(90:99,490:500),]
DataODE1=setDT(CutODE1)[, .(MinA1=min(A1), MaxA1=max(A1), MinP=min(P), MaxP=max(P)), by=list(Ni,D,P1,P0,B,Phase)]
DataODE1=as.data.frame(DataODE1)

# Calculate the difference between minimum and maximum values
DeltaA1=rep(NA,Nsim*2)
DeltaP=rep(NA,Nsim*2)
for (j in 1:length(DataODE1[,1])) {
  DeltaA1[j]=DataODE1[j,8] - DataODE1[j,7]
  DeltaP[j]=DataODE1[j,10] - DataODE1[j,9]
}

# Create a dataset with columns of interest
Data1=data.frame(DataODE1[,c(1:6,7,10)])
Data1[,c(1:5)]=sapply(Data1[,c(1:5)], as.factor)
Data1$DeltaA1=DeltaA1
Data1$DeltaP=DeltaP

# Determine the type of dynamic for each combination
State=rep(NA,Nsim*2)
for (j in 1:length(Data1[,1])) {
  if (Data1[j,7] > 10^-3 & Data1[j,8] > 1 & Data1[j,9] < 10^-1 & Data1[j,10] < 5)     
    State[j]="SS"
  else if (Data1[j,7] > 10^-3 & Data1[j,8] < 1)
    State[j]="EX"
  else if (Data1[j,7] < 10^-3 & Data1[j,8] < 1)
    State[j]="EX"
  else (State[j]="CRC")
}
Data1$State=State

System=rep(NA,Nsim*2)
for (j in 1:length(Data1[,1])) {
  if (Data1[j,7] > 10^-3 & Data1[j,8] > 1 & Data1[j,9] < 10^-1 & Data1[j,10] < 5)     
    System[j]="Yes"
  else if (Data1[j,7] > 10^-3 & Data1[j,8] < 1)
    System[j]="No"
  else if (Data1[j,7] < 10^-3 & Data1[j,8] < 1)
    System[j]="No"
  else (System[j]="No")
}
Data1$System=System


#######################################
### Plot of eco-evolutionary states ###
#######################################

States1=ggplot(Data1, aes(B, P0, group=State, width=1)) +
  geom_tile(aes(fill=State), color="grey90") + coord_cartesian(ylim=c(1.10,18.90), xlim=c(1.10,10.90)) +
  ylab(expression('A'~'ingestion probability difference'~'('*p[Delta]*')')) + xlab(expression('Microplastic density'~'('*10^4~particles~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_blank()) + theme(axis.title.x=element_blank()) +
  scale_y_discrete(breaks=seq(0.9,0.1,length.out=5), labels=sprintf(seq(0.9,0.1,length.out=5), fmt="%.1f")) +
  scale_x_discrete(breaks=seq(0,0.05,length.out=6), labels=sprintf(seq(0,5.0,length.out=6), fmt="%.1f")) +
  theme(panel.background=element_blank(), panel.border=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(plot.title=element_text(face="plain", colour="black", size=22, hjust=0.5)) +
  theme(axis.line=element_line(colour="black", linetype="solid", size=0.5)) +
  scale_fill_manual(values=alpha(c("SS"="royalblue2","CRC"="forestgreen","EX"="firebrick3"),0.8)) +
  theme(panel.spacing.y=unit(1,'lines'), panel.spacing.x=unit(0,'lines')) + 
  theme(strip.background=element_blank(), strip.text.y=element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=Inf, yend=Inf, color="black", size=0.9) +
  facet_wrap(~Phase, ncol=1, nrow=2, strip.position="right") +
  theme(legend.position="none")


##########################################
### Predator densities in the dynamics ###
##########################################

# Calculate density boundaries
CutODE1=SaveODE1[SaveODE1$Time %in% c(90:99,490:500),]
DataODE1=setDT(CutODE1)[, .(MinA1=min(A1), MaxA1=max(A1), MinP=min(P), MaxP=max(P)), by=list(Ni,D,P1,P0,B,Phase)]
DataODE1=as.data.frame(DataODE1)

# Calculate the mean predator final densities
DataPEnd=setDT(CutODE1)[, .(PEnd=mean(P)), by=list(Ni,D,P1,P0,B,Phase)]
DataPEnd=as.data.frame(DataPEnd)

# Select columns of interest
Data2=Data1[,c(1:6,11,12)]
Data1$PEnd=round(DataPEnd[,7],0)
Data2$PEnd=round(DataPEnd[,7],0)
Data1$PEnd[Data1$PEnd < 1]=0
Data2$PEnd[Data2$PEnd < 1]=0

# Filter the dataset
Data2$PEnd[Data2$System=="No"]=NA

Predator1=ggplot(Data2, aes(B, P0, group=PEnd, width=1)) +
  geom_tile(aes(fill=PEnd), color="grey90") + coord_cartesian(ylim=c(1.10,18.90), xlim=c(1.10,10.90)) +
  geom_tile_pattern(aes(pattern=State), fill=NA, pattern_fill="black", pattern_color="black", pattern_density=0.1, pattern_spacing=0.02) +
  ylab(expression('A'~'ingestion probability difference'~'('*p[Delta]*')')) + xlab(expression('Microplastic density'~'('*10^4~particles~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_blank()) + theme(axis.title.x=element_blank()) +
  scale_y_discrete(breaks=seq(0.9,0.1,length.out=5), labels=sprintf(seq(0.9,0.1,length.out=5), fmt="%.1f")) +
  scale_x_discrete(breaks=seq(0,0.05,length.out=6), labels=sprintf(seq(0,5.0,length.out=6), fmt="%.1f")) +
  theme(panel.background=element_blank(), panel.border=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(plot.title=element_text(face="plain", colour="black", size=22, hjust=0.5)) +
  theme(axis.line=element_line(colour="black", linetype="solid", size=0.5)) +
  scale_fill_gradient2(name=expression(atop('Predator density', paste('('*individuals~mL^-1*')'))), 
  low="white", high="royalblue2", na.value="white", breaks=c(2,32), limits=c(0,34)) +
  scale_pattern_manual(values=c("EX"="circle","CRC"="none","EC"="none","SS"="none")) +
  theme(panel.spacing.y=unit(1,'lines'), panel.spacing.x=unit(0,'lines')) + 
  theme(strip.background=element_blank(), strip.text.y=element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=Inf, yend=Inf, color="black", size=0.9) +
  facet_wrap(~Phase, ncol=1, nrow=2, strip.position="right") +
  theme(legend.position="none")


#######################################
### Preys densities in the dynamics ###
#######################################

# Calculate the mean prey final densities
DataAEnd=setDT(CutODE1)[, .(AEnd=mean(A1)), by=list(Ni,D,P1,P0,B,Phase)]
DataAEnd=as.data.frame(DataAEnd)

# Select columns of interest
Data2$AEnd=round(DataAEnd[,7],2)
Data2$AEnd[Data2$AEnd < 10^-3]=0

# Filter the dataset
Data2$AEnd[Data2$System=="No"]=NA

Prey1=ggplot(Data2, aes(B, P0, group=AEnd, width=1)) +
  geom_tile(aes(fill=AEnd), color="grey90") + coord_cartesian(ylim=c(1.10,18.90), xlim=c(1.10,10.90)) +
  geom_tile_pattern(aes(pattern=State), fill=NA, pattern_fill="black", pattern_color="black", pattern_density=0.1, pattern_spacing=0.02) +
  ylab(expression('A'~'ingestion probability difference'~'('*p[Delta]*')')) + xlab(expression('Microplastic density'~'('*10^4~particles~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_blank()) + theme(axis.title.x=element_blank()) +
  scale_y_discrete(breaks=seq(0.9,0.1,length.out=5), labels=sprintf(seq(0.9,0.1,length.out=5), fmt="%.1f")) +
  scale_x_discrete(breaks=seq(0,0.05,length.out=6), labels=sprintf(seq(0,5.0,length.out=6), fmt="%.1f")) +
  theme(panel.background=element_blank(), panel.border=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(plot.title=element_text(face="plain", colour="black", size=22, hjust=0.5)) +
  theme(axis.line=element_line(colour="black", linetype="solid", size=0.5)) +
  scale_fill_gradient2(name=expression(atop('Prey density', paste('('*10^6~cells~mL^-1*')'))), 
  low="white", high="royalblue2", na.value="white", breaks=c(0.25,0.50), limits=c(0.20,0.55)) +
  scale_pattern_manual(values=c("EX"="circle","CRC"="none","EC"="none","SS"="none")) +
  theme(panel.spacing.y=unit(1,'lines'), panel.spacing.x=unit(0,'lines')) + 
  theme(strip.background=element_blank(), strip.text.y=element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=Inf, yend=Inf, color="black", size=0.9) +
  facet_wrap(~Phase, ncol=1, nrow=2, strip.position="right") +
  theme(legend.position="none")




#####################################################################
### ONE UNDEFENDED PREY AND ONE VARIABLE DEFENDED PREY WITH BEADS ###
#####################################################################

# Number of simulations
Nsim=209

# Starting values
Start=c(N=200, A1=0.1, A2=0.1, P=1.0)

#######################################
### Ordinary differential equations ###
#######################################

# Iterate the function for each combination
System=function (Time, Start, Param) {
  with(as.list(c(Start,Param)),{
    
    if (Time >= 100) {B = B} 
    else if (Time < 100) {B = 0}
    
    if (Time >= 100) {Xp = Xp} 
    else if (Time < 100) {Xp = 143}
    
    P2=1.0
    Q1=1.0
    Q2=1.0
    K2=1.0
    B1=1.5
    B2=1.5
    b=1.05
    
    G1 = B1*(N/(K1 + N))                                                                             
    G2 = B2*(N/(K2 + N))                                                                                   
    F1 = (a*P1*Q1*A1^b)/(1 + a*P1*(Ca*h + Q1*(1 - Ca)*h)*A1^b + a*P2*(Ca*h + Q2*(1 - Ca)*h)*A2^b + c*B)
    F2 = (a*P2*Q2*A2^b)/(1 + a*P1*(Ca*h + Q1*(1 - Ca)*h)*A1^b + a*P2*(Ca*h + Q2*(1 - Ca)*h)*A2^b + c*B)
    
    dN = D*(Ni - N) - (1/Xa)*G1*A1 - (1/Xa)*G2*A2
    dA1 = G1*A1 - F1*P - D*A1
    dA2 = G2*A2 - F2*P - D*A2
    dP = Xp*F1*P + Xp*F2*P - D*P
    return(list(c(dN,dA1,dA2,dP)))})
}

# Create saving objects
ListODE=list()
Save=data.frame(Time=integer(0), N=integer(0), A1=integer(0), A2=integer(0), A=integer(0), P=integer(0), Ni=integer(0), D=integer(0), P1=integer(0), B=integer(0))
SaveODE2=data.frame(Time=integer(0), N=integer(0), A1=integer(0), A2=integer(0), A=integer(0), P=integer(0), Ni=integer(0), D=integer(0), P1=integer(0), B=integer(0))

# Calculate the dynamics for each combination
for(i in 1:length(Nutri)) {
  for(j in 1:length(Dilu)) {
    for(k in 1:length(Ing1)) {
      for(l in 1:length(Bead)) {
        P1=Ing1[k]
        K1=Comp1[k]
        B=Bead[l]
        Xp=Conv[l]
        if(Ing1[k]) Comp1[k]
        if(Bead[l]) Conv[l]
        SolveODE=function(time, y=c(N=Nutri[i], A1=0.1, A2=0.1, P=1.0), parms=c(Ni=Nutri[i], D=Dilu[j], a=a, h=h, Ca=Ca, Xa=Xa)) {with(as.list(y), {y[which(y < 10^-20)]=0; return(y)})}
        ListODE=ode(y=c(N=Nutri[i], A1=0.1, A2=0.1, P=1.0), Time, System, parms=c(Ni=Nutri[i], D=Dilu[j], a=a, h=h, Ca=Ca, Xa=Xa), method="lsoda", events=list(func=SolveODE, time=Time))
        ListODE[ListODE<0]=0
        Save=data.frame(Time=ListODE[,1],
                        N=ListODE[,2],
                        A1=ListODE[,3],
                        A2=ListODE[,4],
                        A=ListODE[,3]+ListODE[,4],
                        P=ListODE[,5],
                        Ni=rep(Nutri[i],length(Time)),
                        D=rep(Dilu[j],length(Time)),
                        P1=rep(Ing1[k],length(Time)),
                        B=rep(Bead[l],length(Time)))
        SaveODE2=rbind(SaveODE2,Save)
      }
    }
  } 
}


##################################
### Population dynamics panels ###
##################################

# Name columns according to combinations
Comb=unique(SaveODE2[c("P1","B")])
colnames(Comb)=c("Var1","Var2")
Comb=cbind(Comb,Var3=rep("P1",Nsim),Var4=rep("B",Nsim))
Names=paste(paste(Comb$Var3,Comb$Var1,sep="-"),paste(Comb$Var4,Comb$Var2,sep="-"),sep=" / ")

# Create a dataset
MeltODE2=SaveODE2[,-c(2,5,12,13)] 
MeltODE2$A1=MeltODE2$A1*10^6
MeltODE2$A2=MeltODE2$A2*10^6
MeltODE2=melt(MeltODE2, id.vars=c("Time","Ni","D","P1","B"), variable.name="Organism", value.name="Density")
MeltODE2=MeltODE2[order(-MeltODE2$P1),]
MeltODE2$P1=factor(MeltODE2$P1, levels=unique(MeltODE2$P1))

# Rescale densities
Plot1ODE2=MeltODE2[,c(1:7)]
Plot1ODE2[Plot1ODE2$Organism=="A1",]$Density=log(MeltODE2[MeltODE2$Organism=="A1",]$Density+1)/3
Plot1ODE2[Plot1ODE2$Organism=="A2",]$Density=log(MeltODE2[MeltODE2$Organism=="A2",]$Density+1)/3
Plot1ODE2[Plot1ODE2$Organism=="P",]$Density=log(MeltODE2[MeltODE2$Organism=="P",]$Density+1)

# Calculate effective densities
Plot2ODE2=MeltODE1[,c(1:7)]
Plot2ODE2$Organism=gsub("A1", "A", Plot2ODE2$Organism)
DensityA1=MeltODE2[MeltODE2$Organism=="A1",]$Density
DensityA2=MeltODE2[MeltODE2$Organism=="A2",]$Density
DefenseA1=as.numeric(MeltODE2[MeltODE2$Organism=="A1",]$P1)
DefenseA2=as.numeric(MeltODE2[MeltODE2$Organism=="A2",]$P1)
Plot2ODE2[Plot2ODE2$Organism=="A",]$Density=(log(DensityA1+1)*DefenseA1+log(DensityA2+1)*DefenseA2)/125
Plot2ODE2[Plot2ODE2$Organism=="P",]$Density=log(MeltODE2[MeltODE2$Organism=="P",]$Density+1)

# Represent the population dynamics
tiff('Panel Population Dynamics Evolution 1.tiff', units="in", width=12, height=20, res=1000)
ggplot(Plot1ODE2, aes(Time, Density, group=Organism)) +
  geom_line(aes(color=Organism, linetype=Organism), size=0.5) +
  ylab(expression('Predator density'~'('*ln~individuals~mL^-1*')')) + xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=14)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=14)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=14)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=14)) +
  theme(axis.title.y.left=element_text(angle=90), axis.title.y.right=element_text(angle=90)) +
  scale_y_continuous(labels=function(y) sprintf("%.0f", y), breaks=seq(0,6,by=3), limits=c(0,6),
  sec.axis=sec_axis(~./5*5, expression('Prey density'~'('*ln~cells~mL^-1*')'), labels=sprintf(seq(0,14,by=7), fmt="%.0f"), breaks=seq(0,6,by=3))) +
  scale_x_continuous(labels=function(y) sprintf("%.0f", y), breaks=seq(0,200,by=200), limits=c(0,200)) +
  theme(panel.background=element_blank(), panel.border=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.line=element_line(colour="black", linetype="solid", size=0.5)) +
  scale_color_manual(values=alpha(c("A1"="royalblue2","A2"="royalblue2","P"="firebrick3"),0.8)) +
  scale_linetype_manual(values=c("A1"="solid","A2"="11","P"="solid")) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color="black", size=0.9) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color="black", size=0.9) +
  facet_wrap(P1~B, ncol=11, nrow=19) +
  theme(legend.position="none")
dev.off()

tiff('Panel Population Dynamics Evolution 2.tiff', units="in", width=12, height=20, res=1000)
ggplot(Plot2ODE2, aes(Time, Density, group=Organism)) +
  geom_line(aes(color=Organism, linetype=Organism), size=0.5) +
  ylab(expression('Predator density'~'('*ln~individuals~mL^-1*')')) + xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=14)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=14)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=14)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=14)) +
  theme(axis.title.y.left=element_text(angle=90), axis.title.y.right=element_text(angle=90)) +
  scale_y_continuous(labels=function(y) sprintf("%.0f", y), breaks=seq(0,4,by=1), limits=c(0,4),
  sec.axis=sec_axis(~.*5/4, expression(italic('C. reinhardtii')~'effective density'~'('*10^3*')'), labels=function(y) sprintf("%.0f", y), breaks=seq(0,5,by=1))) +
  scale_x_continuous(labels=function(y) sprintf("%.0f", y), breaks=seq(0,200,by=200), limits=c(0,200)) +
  theme(panel.background=element_blank(), panel.border=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.line=element_line(colour="black", linetype="solid", size=0.5)) +
  scale_color_manual(values=alpha(c("A"="royalblue2","P"="firebrick3"),0.8)) +
  scale_linetype_manual(values=c("A"="solid","P"="solid")) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color="black", size=0.9) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color="black", size=0.9) +
  facet_wrap(P1~B, ncol=11, nrow=19) +
  theme(legend.position="none")
dev.off()


###############################################
### Eco-evolutionary states of the dynamics ###
###############################################

# Define phases for beads
SaveODE2$Phase=ifelse(SaveODE2$Time < 100, "Pre", "Post")
SaveODE2$Phase=factor(SaveODE2$Phase, levels=unique(SaveODE1$Phase))

# Calculate the difference in ingestion probabilities
SaveODE2$P0=as.character(round(1-as.numeric(SaveODE2$P1),2))

# Calculate density boundaries
CutODE2=SaveODE2[SaveODE2$Time %in% c(90:99,491:500),]
DataODE2=setDT(CutODE2)[, .(MinA1=min(A1), MaxA1=max(A1), MinA2=min(A2), MaxA2=max(A2), MinP=min(P), MaxP=max(P)), by=list(Ni,D,P1,P0,B,Phase)]
DataODE2=as.data.frame(DataODE2)

# Calculate the difference between minimum and maximum values
DeltaA1=rep(NA,Nsim*2)
DeltaA2=rep(NA,Nsim*2)
DeltaP=rep(NA,Nsim*2)
for (j in 1:length(DataODE2[,1])) {
  DeltaA1[j]=DataODE2[j,8] - DataODE2[j,7]
  DeltaA2[j]=DataODE2[j,10] - DataODE2[j,9]
  DeltaP[j]=DataODE2[j,12] - DataODE2[j,11]
}

# Create a dataset with columns of interest
Data3=data.frame(DataODE2[,c(1:6,7,9,11)])
Data3[,c(1:5)]=sapply(Data3[,c(1:5)], as.factor)
Data3$DeltaA1=DeltaA1
Data3$DeltaA2=DeltaA2
Data3$DeltaP=DeltaP

# Determine the type of dynamic for each combination
State=rep(NA,Nsim)
for (j in 1:length(Data3[,1])) {
  if (Data3[j,7] > 10^-3 & Data3[j,8] > 10^-3 & Data3[j,9] > 1 & Data3[j,10] > 10^-1 & Data3[j,11] > 10^-1 & Data3[j,12] > 5)     
    State[j]="EC"
  else if (Data3[j,7] > 10^-3 & Data3[j,8] > 10^-3 & Data3[j,9] > 1 & Data3[j,10] < 10^-1 & Data3[j,11] < 10^-1 & Data3[j,12] < 5)
    State[j]="SS"
  else if (Data3[j,7] > 10^-3 & Data3[j,8] < 10^-3 & Data3[j,9] > 1 & Data3[j,10] < 10^-1 & Data3[j,12] < 5)
    State[j]="CRC"
  else if (Data3[j,7] < 10^-3 & Data3[j,8] > 10^-3 & Data3[j,9] > 1 & Data3[j,11] < 10^-1 & Data3[j,12] < 5)
    State[j]="CRC"
  else if (Data3[j,7] < 10^-3 & Data3[j,8] > 10^-3 & Data3[j,9] < 1)     
    State[j]="EX"
  else if (Data3[j,7] > 10^-3 & Data3[j,8] < 10^-3 & Data3[j,9] < 1)
    State[j]="EX"
  else if (Data3[j,7] < 10^-3 & Data3[j,8] < 10^-3 & Data3[j,9] < 1)
    State[j]="EX"
  else (State[j]="CRC")
}
Data3$State=State

System=rep(NA,Nsim*2)
for (j in 1:length(Data3[,1])) {
  if (Data3[j,7] > 10^-3 & Data3[j,8] > 10^-3 & Data3[j,9] > 1 & Data3[j,10] > 10^-1 & Data3[j,11] > 10^-1 & Data3[j,12] > 5)     
    System[j]="No"
  else if (Data3[j,7] > 10^-3 & Data3[j,8] > 10^-3 & Data3[j,9] > 1 & Data3[j,10] < 10^-1 & Data3[j,11] < 10^-1 & Data3[j,12] < 5)
    System[j]="Yes"
  else if (Data3[j,7] > 10^-3 & Data3[j,8] < 10^-3 & Data3[j,9] > 1 & Data3[j,10] < 10^-1 & Data3[j,12] < 5)
    System[j]="Yes"
  else if (Data3[j,7] < 10^-3 & Data3[j,8] > 10^-3 & Data3[j,9] > 1 & Data3[j,11] < 10^-1 & Data3[j,12] < 5)
    System[j]="Yes"
  else if (Data3[j,7] < 10^-3 & Data3[j,8] > 10^-3 & Data3[j,9] < 1)     
    System[j]="No"
  else if (Data3[j,7] > 10^-3 & Data3[j,8] < 10^-3 & Data3[j,9] < 1)
    System[j]="No"
  else if (Data3[j,7] < 10^-3 & Data3[j,8] < 10^-3 & Data3[j,9] < 1)
    System[j]="No"
  else (System[j]="No")
}
Data3$System=System


#######################################
### Plot of eco-evolutionary states ###
#######################################

# Correct for dynamics with same preys
Data3$State=ifelse(Data3$P1=="1" & Data3$State=="EC", "CRC", Data3$State)

States2=ggplot(Data3, aes(B, P0, group=State, width=1)) +
  geom_tile(aes(fill=State), color="grey90") + coord_cartesian(ylim=c(1.10,18.90), xlim=c(1.10,10.90)) +
  ylab(expression('A'~'ingestion probability difference'~'('*p[Delta]*')')) + xlab(expression('Microplastic density'~'('*10^4~particles~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_blank()) + theme(axis.title.x=element_blank()) +
  scale_y_discrete(breaks=seq(0.9,0.1,length.out=5), labels=sprintf(seq(0.9,0.1,length.out=5), fmt="%.1f")) +
  scale_x_discrete(breaks=seq(0,0.05,length.out=6), labels=sprintf(seq(0,5.0,length.out=6), fmt="%.1f")) +
  theme(panel.background=element_blank(), panel.border=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(plot.title=element_text(face="plain", colour="black", size=22, hjust=0.5)) +
  theme(axis.line=element_line(colour="black", linetype="solid", size=0.5)) +
  scale_fill_manual(values=alpha(c("SS"="royalblue2","CRC"="darkorange3","EC"="darkorange3","EX"="firebrick3"),0.8)) +
  theme(panel.spacing.y=unit(1,'lines'), panel.spacing.x=unit(0,'lines')) + 
  theme(strip.background=element_blank(), strip.text.y=element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=Inf, yend=Inf, color="black", size=0.9) +
  facet_wrap(~Phase, ncol=1, nrow=2, strip.position="right") +
  theme(legend.position="none")


##########################################
### Predator densities in the dynamics ###
##########################################

# Calculate density boundaries
CutODE2=SaveODE2[SaveODE2$Time %in% c(90:99,491:500),]
DataODE2=setDT(CutODE2)[, .(MinA1=min(A1), MaxA1=max(A1), MinA2=min(A2), MaxA2=max(A2), MinP=min(P), MaxP=max(P)), by=list(Ni,D,P1,P0,B,Phase)]
DataODE2=as.data.frame(DataODE2)

# Calculate the mean predator final densities
DataPEnd=setDT(CutODE2)[, .(PEnd=mean(P)), by=list(Ni,D,P1,P0,B,Phase)]
DataPEnd=as.data.frame(DataPEnd)

# Select columns of interest
Data4=Data3[,c(1:6,13,14)]
Data3$PEnd=round(DataPEnd[,7],0)
Data4$PEnd=round(DataPEnd[,7],0)
Data3$PEnd[Data4$PEnd < 1]=0
Data4$PEnd[Data4$PEnd < 1]=0

# Filter the dataset
Data4$PEnd[Data4$System=="No"]=NA

Predator2=ggplot(Data4, aes(B, P0, group=PEnd, width=1)) +
  geom_tile(aes(fill=PEnd), color="grey90") + coord_cartesian(ylim=c(1.10,18.90), xlim=c(1.10,10.90)) +
  geom_tile_pattern(aes(pattern=State), fill=NA, pattern_fill="black", pattern_color="black", pattern_density=0.1, pattern_spacing=0.02) +
  ylab(expression('A'~'ingestion probability difference'~'('*p[Delta]*')')) + xlab(expression('Microplastic density'~'('*10^4~particles~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_blank()) + theme(axis.title.x=element_blank()) +
  scale_y_discrete(breaks=seq(0.9,0.1,length.out=5), labels=sprintf(seq(0.9,0.1,length.out=5), fmt="%.1f")) +
  scale_x_discrete(breaks=seq(0,0.05,length.out=6), labels=sprintf(seq(0,5.0,length.out=6), fmt="%.1f")) +
  theme(panel.background=element_blank(), panel.border=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(plot.title=element_text(face="plain", colour="black", size=22, hjust=0.5)) +
  theme(axis.line=element_line(colour="black", linetype="solid", size=0.5)) +
  scale_fill_gradient2(name=expression(atop('Predator density', paste('('*individuals~mL^-1*')'))), 
  low="white", high="royalblue2", na.value="white", breaks=c(2,32), limits=c(0,34)) +
  scale_pattern_manual(values=c("EX"="circle","CRC"="none","EC"="none","SS"="none")) +
  theme(panel.spacing.y=unit(1,'lines'), panel.spacing.x=unit(0,'lines')) + 
  theme(strip.background=element_blank(), strip.text.y=element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=Inf, yend=Inf, color="black", size=0.9) +
  facet_wrap(~Phase, ncol=1, nrow=2, strip.position="right") +
  theme(legend.position="none")


######################################
### Prey densities in the dynamics ###
######################################

# Calculate the prey frequencies
SaveODE2$FA1=SaveODE2$A1/(SaveODE2$A1+SaveODE2$A2)
SaveODE2$FA2=SaveODE2$A2/(SaveODE2$A1+SaveODE2$A2)

# Calculate the mean final prey frequencies
CutODE3=SaveODE2[SaveODE2$Time %in% c(90:99,491:500),]
Data4$FA1End=as.data.frame(setDT(CutODE3)[, .(FA1End=mean(FA1)*100), by=list(P1,P0,B,Phase)])[,5]
Data4$FA2End=as.data.frame(setDT(CutODE3)[, .(FA2End=mean(FA2)*100), by=list(P1,P0,B,Phase)])[,5]

# Select columns of interest
Data4$FA1End=round(Data4$FA1End,2)
Data4$FA2End=round(Data4$FA2End,2)
Data4$FA1End[Data4$FA1End < 10^-2]=0
Data4$FA2End[Data4$FA2End < 10^-2]=0

# Calculate the mean prey final densities
DataAEnd=setDT(CutODE2)[, .(AEnd=mean(A1+A2)), by=list(Ni,D,P1,P0,B,Phase)]
DataAEnd=as.data.frame(DataAEnd)

# Select columns of interest
Data4$AEnd=round(DataAEnd[,7],2)
Data4$AEnd[Data4$AEnd < 10^-3]=0

# Filter the dataset
Data4$FA1End[Data4$System=="No"]=NA
Data4$FA2End[Data4$System=="No"]=NA
Data4$AEnd[Data4$System=="No"]=NA

Prey2=ggplot(Data4, aes(B, P0, group=AEnd, width=1)) +
  geom_tile(aes(fill=AEnd), color="grey90") + coord_cartesian(ylim=c(1.10,18.90), xlim=c(1.10,10.90)) +
  geom_tile_pattern(aes(pattern=State), fill=NA, pattern_fill="black", pattern_color="black", pattern_density=0.1, pattern_spacing=0.02) +
  ylab(expression('A'~'ingestion probability difference'~'('*p[Delta]*')')) + xlab(expression('Microplastic density'~'('*10^4~particles~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_blank()) + theme(axis.title.x=element_blank()) +
  scale_y_discrete(breaks=seq(0.9,0.1,length.out=5), labels=sprintf(seq(0.9,0.1,length.out=5), fmt="%.1f")) +
  scale_x_discrete(breaks=seq(0,0.05,length.out=6), labels=sprintf(seq(0,5.0,length.out=6), fmt="%.1f")) +
  theme(panel.background=element_blank(), panel.border=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(plot.title=element_text(face="plain", colour="black", size=22, hjust=0.5)) +
  theme(axis.line=element_line(colour="black", linetype="solid", size=0.5)) +
  scale_fill_gradient2(name=expression(atop('Prey density', paste('('*10^6~cells~mL^-1*')'))), 
  low="white", high="royalblue2", na.value="white", breaks=c(0.25,0.50), limits=c(0.20,0.55)) +
  scale_pattern_manual(values=c("EX"="circle","CRC"="none","EC"="none","SS"="none")) +
  theme(panel.spacing.y=unit(1,'lines'), panel.spacing.x=unit(0,'lines')) + 
  theme(strip.background=element_blank(), strip.text.y=element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=Inf, yend=Inf, color="black", size=0.9) +
  facet_wrap(~Phase, ncol=1, nrow=2, strip.position="right") +
  theme(legend.position="none")




##################################################
##################################################
##### INDIRECT ECO-EVOLUTIONARY FACILITATION #####
##################################################
##################################################

###############################
### Eco-evolutionary states ###
###############################

# Panel of evolutionary states
tiff('Evolutionary States.tiff', units="in", width=15, height=15, res=1000)
States1=States1 + theme(strip.text.y=element_blank())
Panel=plot_grid(States1, States2, align="h", vjust=1, nrow=1, ncol=2)
Yaxis=textGrob(expression('A'~'ingestion probability difference'~'('*p[Delta]*')'), gp=gpar(fontface="bold", fontsize=18), rot=90)
Xaxis=textGrob(expression('Microplastic density'~'('*10^4~particles~mL^-1*')'), gp=gpar(fontface="bold", fontsize=18))
grid.arrange(arrangeGrob(Panel, left=Yaxis, bottom=Xaxis))
dev.off()

# Panel of prey final densities
Prey1=Prey1 + theme(strip.text.y=element_blank())
Panel=plot_grid(Prey1, Prey2, align="h", vjust=1, nrow=1, ncol=2)
Yaxis=textGrob(expression('A'~'ingestion probability difference'~'('*p[Delta]*')'), gp=gpar(fontface="bold", fontsize=18), rot=90)
Xaxis=textGrob(expression('Microplastic density'~'('*10^4~particles~mL^-1*')'), gp=gpar(fontface="bold", fontsize=18))
grid.arrange(arrangeGrob(Panel, left=Yaxis, bottom=Xaxis))

# Panel of predator final densities
Predator1=Predator1 + theme(strip.text.y=element_blank())
Panel=plot_grid(Predator1, Predator2, align="h", vjust=1, nrow=1, ncol=2)
Yaxis=textGrob(expression('A'~'ingestion probability difference'~'('*p[Delta]*')'), gp=gpar(fontface="bold", fontsize=18), rot=90)
Xaxis=textGrob(expression('Microplastic density'~'('*10^4~particles~mL^-1*')'), gp=gpar(fontface="bold", fontsize=18))
grid.arrange(arrangeGrob(Panel, left=Yaxis, bottom=Xaxis))


###################################################################
### Eco-evolutionary rescue or facilitation: predator densities ###
###################################################################

# Create the datasets
Data5=Data4[,c(1:8)]
Data6=Data4[,c(1:8)]
Data6=Data6[order(Data6$Phase),][,-6]
Data6$Organism=c(rep("1P",209),rep("2P",209))

# Calculate the mean predator final densities
Data2PEnd=setDT(CutODE1)[, .(PEnd=mean(P)), by=list(Ni,D,P1,P0,B,Phase)]
Data4PEnd=setDT(CutODE2)[, .(PEnd=mean(P)), by=list(Ni,D,P1,P0,B,Phase)]
Data2PEnd=as.data.frame(Data2PEnd)
Data4PEnd=as.data.frame(Data4PEnd)

# Calculate the difference of predator final densities between systems
Data5$PEnd1=Data2PEnd$PEnd
Data5$PEnd2=Data4PEnd$PEnd
Data5$DeltaP=Data5$PEnd2-Data5$PEnd1

# Calculate the difference of predator final densities between phases
DiffeP1=subset(Data2PEnd, Phase=="Post")[,7]-subset(Data2PEnd, Phase=="Pre")[,7]
DiffeP2=subset(Data4PEnd, Phase=="Post")[,7]-subset(Data4PEnd, Phase=="Pre")[,7]
Data6$DiffeP=c(DiffeP1,DiffeP2)

# Replace states and systems
Data5$State=Data4$State; Data5$System=Data4$System
Data6$State=c(Data2$State[Data2$Phase=="Post"], Data4$State[Data4$Phase=="Post"])
Data6$System=c(Data2$System[Data2$Phase=="Post"], Data4$System[Data4$Phase=="Post"])

# Define evolutionary rescue and facilitation
Mechanism=rep(NA,Nsim)
DataP1=subset(Data5, Phase=="Pre")
DataP2=subset(Data5, Phase=="Post")
for (i in 1:length(DataP1[,1])) {
  if (DataP1[i,9] < 10^-1 & DataP2[i,10] > 10^-1)
    Mechanism[i]="IEP"
  else if (DataP2[i,9] < 10^-1 & DataP2[i,10] > 10^-1)
    Mechanism[i]="IER"
  else if (DataP2[i,9] < DataP2[i,10] & DataP2[i,11] > 1)
    Mechanism[i]="IEF"
  else (Mechanism[i]="No")
}
Data6$Mechanism=rep(Mechanism,2)

tiff('Predator Densities Phase No Evolution.tiff', units="in", width=8, height=8, res=1000)
ggplot(subset(Data6, Organism=="1P"), aes(B, P0, width=1)) +
  geom_tile(aes(fill=DiffeP), color="grey90") + coord_cartesian(ylim=c(1.10,18.90), xlim=c(1.10,10.90)) +
  geom_tile_pattern(aes(pattern=State), fill=NA, pattern_fill="black", pattern_color="black", pattern_density=0.1, pattern_spacing=0.02) +
  ylab(expression('A'~'ingestion probability difference'~'('*p[A[Delta]]*')')) + xlab(expression('Microplastic density'~'('*10^4~particles~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_discrete(breaks=seq(0.9,0.1,length.out=5), labels=sprintf(seq(0.9,0.1,length.out=5), fmt="%.1f")) +
  scale_x_discrete(breaks=seq(0,0.05,length.out=6), labels=sprintf(seq(0,5.0,length.out=6), fmt="%.1f")) +
  theme(panel.background=element_blank(), panel.border=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.line=element_line(colour="black", linetype="solid", size=0.5)) +
  scale_fill_gradientn(name=expression(Delta~'Predator density'~'('*individuals~mL^-1*')'), 
  colours=c("firebrick3","white","royalblue2"), values=rescale(x=c(0,32/40,1.0),from=c(0,1.0)), na.value="white", breaks=c(-30,6), limits=c(-32,8)) +
  scale_pattern_manual(values=c("EX"="circle","CRC"="none","EC"="none","SS"="none")) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  theme(plot.title=element_text(face="plain", colour="black", size=22, vjust=2, hjust=0.5)) +
  annotate("segment", x=-Inf, xend=Inf, y=Inf, yend=Inf, color="black", size=0.9) +
  guides(fill=guide_colorbar(title.position="top", title.hjust=0.5, ticks.colour=NA, barwidth=15), pattern="none") +
  theme(legend.title=element_text(size=14), legend.text=element_text(size=14)) +
  theme(legend.key.size=unit(1.5,"line"), legend.position="top") +
  facet_wrap(~Organism, ncol=2, nrow=1) +
  theme(legend.position="none")
dev.off()

tiff('Predator Densities Phase Evolution.tiff', units="in", width=8, height=8, res=1000)
ggplot(subset(Data6, Organism=="2P"), aes(B, P0, width=1)) +
  geom_tile(aes(fill=DiffeP), color="grey90") + coord_cartesian(ylim=c(1.10,18.90), xlim=c(1.10,10.90)) +
  geom_tile_pattern(aes(pattern=State), fill=NA, pattern_fill="black", pattern_color="black", pattern_density=0.1, pattern_spacing=0.02) +
  ylab(expression('A'~'ingestion probability difference'~'('*p[A[Delta]]*')')) + xlab(expression('Microplastic density'~'('*10^4~particles~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_discrete(breaks=seq(0.9,0.1,length.out=5), labels=sprintf(seq(0.9,0.1,length.out=5), fmt="%.1f")) +
  scale_x_discrete(breaks=seq(0,0.05,length.out=6), labels=sprintf(seq(0,5.0,length.out=6), fmt="%.1f")) +
  theme(panel.background=element_blank(), panel.border=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.line=element_line(colour="black", linetype="solid", size=0.5)) +
  scale_fill_gradientn(name=expression(Delta~'Predator density'~'('*individuals~mL^-1*')'), 
  colours=c("firebrick3","white","royalblue2"), values=rescale(x=c(0,32/40,1.0),from=c(0,1.0)), na.value="white", breaks=c(-30,6), limits=c(-32,8)) +
  scale_pattern_manual(values=c("EX"="circle","CRC"="none","EC"="none","SS"="none")) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  theme(plot.title=element_text(face="plain", colour="black", size=22, vjust=2, hjust=0.5)) +
  annotate("segment", x=-Inf, xend=Inf, y=Inf, yend=Inf, color="black", size=0.9) +
  guides(fill=guide_colorbar(title.position="top", title.hjust=0.5, ticks.colour=NA, barwidth=15), pattern="none") +
  theme(legend.title=element_text(size=14), legend.text=element_text(size=14)) +
  theme(legend.key.size=unit(1.5,"line"), legend.position="top") +
  facet_wrap(~Organism, ncol=2, nrow=1) +
  theme(legend.position="none")
dev.off()


###############################################################
### Eco-evolutionary rescue or facilitation: prey densities ###
###############################################################

# Calculate the mean prey final densities
Data2AEnd=setDT(CutODE1)[, .(AEnd=mean(A1)), by=list(Ni,D,P1,P0,B,Phase)]
Data4AEnd=setDT(CutODE2)[, .(AEnd=mean(A1+A2)), by=list(Ni,D,P1,P0,B,Phase)]
Data2AEnd=as.data.frame(Data2AEnd)
Data4AEnd=as.data.frame(Data4AEnd)

# Calculate the difference of prey final densities between systems
Data5$AEnd1=Data2AEnd$AEnd
Data5$AEnd2=Data4AEnd$AEnd
Data5$DeltaA=Data5$AEnd2-Data5$AEnd1

# Calculate the difference of prey final densities between phases
DiffeA1=subset(Data2AEnd, Phase=="Post")[,7]-subset(Data2AEnd, Phase=="Pre")[,7] 
DiffeA2=subset(Data4AEnd, Phase=="Post")[,7]-subset(Data4AEnd, Phase=="Pre")[,7]
Data6$DiffeA=c(DiffeA1,DiffeA2)

# Replace states and systems
Data5$State=Data4$State; Data5$System=Data4$System
Data6$State=c(Data2$State[Data2$Phase=="Post"], Data4$State[Data4$Phase=="Post"])
Data6$System=c(Data2$System[Data2$Phase=="Post"], Data4$System[Data4$Phase=="Post"])

# Define evolutionary rescue and facilitation
Mechanism=rep(NA,Nsim)
DataP1=subset(Data5, Phase=="Pre")
DataP2=subset(Data5, Phase=="Post")
for (i in 1:length(DataP1[,1])) {
  if (DataP1[i,9] < 10^-1 & DataP2[i,10] > 10^-1)
    Mechanism[i]="IEP"
  else if (DataP2[i,9] < 10^-1 & DataP2[i,10] > 10^-1)
    Mechanism[i]="IER"
  else if (DataP2[i,9] < DataP2[i,10] & DataP2[i,11] > 1)
    Mechanism[i]="IEF"
  else (Mechanism[i]="No")
}
Data6$Mechanism=rep(Mechanism,2)

tiff('Prey Densities Phase Evolution.tiff', units="in", width=8, height=8, res=1000)
ggplot(subset(Data6, Organism=="2P"), aes(B, P0, width=1)) +
  geom_tile(aes(fill=DiffeA), color="grey90") + coord_cartesian(ylim=c(1.10,18.90), xlim=c(1.10,10.90)) +
  geom_tile_pattern(aes(pattern=State), fill=NA, pattern_fill="black", pattern_color="black", pattern_density=0.1, pattern_spacing=0.02) +
  ylab(expression('A'~'ingestion probability difference'~'('*p[A[Delta]]*')')) + xlab(expression('Microplastic density'~'('*10^4~particles~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_discrete(breaks=seq(0.9,0.1,length.out=5), labels=sprintf(seq(0.9,0.1,length.out=5), fmt="%.1f")) +
  scale_x_discrete(breaks=seq(0,0.05,length.out=6), labels=sprintf(seq(0,5.0,length.out=6), fmt="%.1f")) +
  theme(panel.background=element_blank(), panel.border=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.line=element_line(colour="black", linetype="solid", size=0.5)) +
  scale_fill_gradient2(name=expression('Prey density'~'('*10^6~cells~mL^-1*')'), 
  low="firebrick3", high="royalblue2", na.value="white", breaks=c(0.00,0.30), limits=c(-0.05,0.35)) +
  scale_pattern_manual(values=c("EX"="circle","CRC"="none","EC"="none","SS"="none")) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  theme(plot.title=element_text(face="plain", colour="black", size=22, vjust=2, hjust=0.5)) +
  annotate("segment", x=-Inf, xend=Inf, y=Inf, yend=Inf, color="black", size=0.9) +
  guides(fill=guide_colorbar(title.position="top", title.hjust=0.5, ticks.colour=NA, barwidth=15), pattern="none") +
  theme(legend.title=element_text(size=14), legend.text=element_text(size=14)) +
  theme(legend.key.size=unit(1.5,"line"), legend.position="top") +
  facet_wrap(~Organism, ncol=2, nrow=1) +
  theme(legend.position="none")
dev.off()

# Add prey final frequencies
Data5$FA1End=Data4$FA1End
Data5$FA2End=Data4$FA2End

# Remove scenarios without prey diversity
Data5$FA1End=ifelse(Data5$P1==1, NA, Data5$FA1End)
Data5$FA2End=ifelse(Data5$P1==1, NA, Data5$FA2End)
Data5$FA1End=ifelse(Data5$FA1End==0, NA, Data5$FA1End)
Data5$FA2End=ifelse(Data5$FA1End==0, NA, Data5$FA2End)

tiff('Prey Frequencies Phase Evolution.tiff', units="in", width=8, height=8, res=1000)
ggplot(subset(Data5, Phase=="Post"), aes(B, P0, width=1)) +
  geom_tile(aes(fill=FA2End), color="grey90") + coord_cartesian(ylim=c(1.10,18.90), xlim=c(1.10,10.90)) +
  geom_tile_pattern(aes(pattern=State), fill=NA, pattern_fill="black", pattern_color="black", pattern_density=0.1, pattern_spacing=0.02) +
  ylab(expression('A'~'ingestion probability difference'~'('*p[A[Delta]]*')')) + xlab(expression('Microplastic density'~'('*10^4~particles~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_discrete(breaks=seq(0.9,0.1,length.out=5), labels=sprintf(seq(0.9,0.1,length.out=5), fmt="%.1f")) +
  scale_x_discrete(breaks=seq(0,0.05,length.out=6), labels=sprintf(seq(0,5.0,length.out=6), fmt="%.1f")) +
  theme(panel.background=element_blank(), panel.border=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.line=element_line(colour="black", linetype="solid", size=0.5)) +
  scale_fill_gradient2(name=expression('Prey'~'A'[2]~'frequency'~'(%)'), 
  low="white", mid=alpha("royalblue2",0.3), high=alpha("royalblue2",0.8), na.value="white", breaks=c(5,95), limits=c(0,100)) +
  scale_pattern_manual(values=c("EX"="circle","CRC"="none","EC"="none","SS"="none")) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  theme(plot.title=element_text(face="plain", colour="black", size=22, vjust=2, hjust=0.5)) +
  annotate("segment", x=-Inf, xend=Inf, y=Inf, yend=Inf, color="black", size=0.9) +
  guides(fill=guide_colorbar(title.position="top", title.hjust=0.5, ticks.colour=NA, barwidth=15), pattern="none") +
  theme(legend.title=element_text(size=14), legend.text=element_text(size=14)) +
  theme(legend.key.size=unit(1.5,"line"), legend.position="top") +
  facet_wrap(~Phase, ncol=2, nrow=1) +
  theme(legend.position="none")
dev.off()


#####################################################################
### Eco-evolutionary rescue or facilitation: predator differences ###
#####################################################################

tiff('Predator Densities System.tiff', units="in", width=8, height=8, res=1000)
ggplot(subset(Data5, Phase=="Post"), aes(B, P0, width=1)) +
  geom_tile(aes(fill=DeltaP), color="grey90") + coord_cartesian(ylim=c(1.10,18.90), xlim=c(1.10,10.90)) +
  geom_tile_pattern(aes(pattern=State), fill=NA, pattern_fill="black", pattern_color="black", pattern_density=0.1, pattern_spacing=0.02) +
  ylab(expression('A'~'ingestion probability difference'~'('*p[A[Delta]]*')')) + xlab(expression('Microplastic density'~'('*10^4~particles~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_discrete(breaks=seq(0.9,0.1,length.out=5), labels=sprintf(seq(0.9,0.1,length.out=5), fmt="%.1f")) +
  scale_x_discrete(breaks=seq(0,0.05,length.out=6), labels=sprintf(seq(0,5.0,length.out=6), fmt="%.1f")) +
  theme(panel.background=element_blank(), panel.border=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.line=element_line(colour="black", linetype="solid", size=0.5)) +
  scale_fill_gradientn(name=expression(Delta~'Predator density'~'('*individuals~mL^-1*')'), 
  colours=c("firebrick3","white","royalblue2"), values=rescale(x=c(0,4/34,1.0),from=c(0,1.0)), na.value="white", breaks=c(-2,28), limits=c(-4,30)) +
  scale_pattern_manual(values=c("EX"="circle","CRC"="none","EC"="none","SS"="none")) +
  theme(plot.title=element_text(face="plain", colour="black", size=22, vjust=2, hjust=0.5)) +
  annotate("segment", x=-Inf, xend=Inf, y=Inf, yend=Inf, color="black", size=0.9) +
  guides(fill=guide_colorbar(title.position="top", title.hjust=0.5, ticks.colour=NA, barwidth=15), pattern="none") +
  theme(legend.title=element_text(size=14), legend.text=element_text(size=14)) +
  theme(legend.key.size=unit(1.5,"line"), legend.position="top") +
  theme(legend.position="none")
dev.off()


##################################################################
### Eco-evolutionary rescue or facilitation: predator decrease ###
##################################################################

# Cut the datasets
DataODE1A=subset(SaveODE1, Time %in% c(100:150))
DataODE1B=subset(SaveODE1, Time==100)
OutODE1A=setDT(DataODE1A)[, .SD[which.min(P)], by=list(P1,B)]
TimeODE1A=as.data.frame(OutODE1A)

DataODE2A=subset(SaveODE2, Time %in% c(100:150))
DataODE2B=subset(SaveODE2, Time==100)
OutODE2A=setDT(DataODE2A)[, .SD[which.min(P)], by=list(P1,B)]
TimeODE2A=as.data.frame(OutODE2A)

# Matching the times for low predator densities
Data7=rbind(TimeODE1A[,c(3,7:8,1,10,2,6)],TimeODE2A[,c(3,9:10,1,12,2,8)])
Data7$State=Data6$State; Data7$System=Data6$System
Data7$Organism=c(rep("1P",209),rep("2P",209))

# Calculate predator density change
Data7$DiffeP=c(TimeODE1A$P-DataODE1B$P,TimeODE2A$P-DataODE2B$P)
Data7$P0=as.character(Data7$P0); Data7$B=as.character(Data7$B)

tiff('Predator Decrease Phase.tiff', units="in", width=8, height=8.5, res=1000)
ggplot(Data7, aes(B, P0, width=1)) +
  geom_tile(aes(fill=DiffeP), color="grey90") + coord_cartesian(ylim=c(1.10,18.90), xlim=c(1.10,10.90)) +
  geom_tile_pattern(aes(pattern=State), fill=NA, pattern_fill="black", pattern_color="black", pattern_density=0.1, pattern_spacing=0.02) +
  ylab(expression('A'~'ingestion probability difference'~'('*p[A[Delta]]*')')) + xlab(expression('Microplastic density'~'('*10^4~particles~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_discrete(breaks=seq(0.9,0.1,length.out=5), labels=sprintf(seq(0.9,0.1,length.out=5), fmt="%.1f")) +
  scale_x_discrete(breaks=seq(0,0.05,length.out=6), labels=sprintf(seq(0,5.0,length.out=6), fmt="%.1f")) +
  theme(panel.background=element_blank(), panel.border=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.line=element_line(colour="black", linetype="solid", size=0.5)) +
  scale_fill_gradient2(name=expression('Predator density'~'('*individuals~mL^-1*')'),
  low="firebrick3", high="royalblue2", na.value="white", breaks=c(-28,-2), limits=c(-30,0)) +
  scale_pattern_manual(values=c("EX"="circle","EC"="none","CRC"="none","EC"="none","SS"="none")) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  theme(plot.title=element_text(face="plain", colour="black", size=22, vjust=2, hjust=0.5)) +
  annotate("segment", x=-Inf, xend=Inf, y=Inf, yend=Inf, color="black", size=0.9) +
  guides(fill=guide_colorbar(title.position="top", title.hjust=0.5, ticks.colour=NA, barwidth=15), pattern="none") +
  theme(legend.title=element_text(size=14), legend.text=element_text(size=14)) +
  theme(legend.key.size=unit(1.5,"line"), legend.position="top") +
  facet_wrap(~Organism, ncol=2, nrow=1) +
  theme(legend.position="none")
dev.off()


###########################################################
### Eco-evolutionary rescue or facilitation: mechanisms ###
###########################################################

tiff('Evolutionary Mechanisms.tiff', units="in", width=8, height=8, res=1000)
ggplot(subset(Data6, Organism=="2P"), aes(B, P0, width=1)) +
  geom_tile(aes(fill=Mechanism), color="grey90") + coord_cartesian(ylim=c(1.10,18.90), xlim=c(1.10,10.90)) +
  geom_tile_pattern(aes(pattern=State), fill=NA, pattern_fill="black", pattern_color="black", pattern_density=0.1, pattern_spacing=0.02) +
  ylab(expression('A'~'ingestion probability difference'~'('*p[A[Delta]]*')')) + xlab(expression('Microplastic density'~'('*10^4~particles~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_discrete(breaks=seq(0.9,0.1,length.out=5), labels=sprintf(seq(0.9,0.1,length.out=5), fmt="%.1f")) +
  scale_x_discrete(breaks=seq(0,0.05,length.out=6), labels=sprintf(seq(0,5.0,length.out=6), fmt="%.1f")) +
  theme(panel.background=element_blank(), panel.border=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.line=element_line(colour="black", linetype="solid", size=0.5)) +
  scale_fill_manual(values=alpha(c("IEP"="firebrick3","IER"="royalblue2","IEF"="forestgreen","No"="white"),0.6)) +
  scale_pattern_manual(values=c("EX"="circle","EC"="none","CRC"="none","EC"="none","SS"="none")) +
  theme(plot.title=element_text(face="plain", colour="black", size=22, vjust=2, hjust=0.5)) +
  annotate("segment", x=-Inf, xend=Inf, y=Inf, yend=Inf, color="black", size=0.9) +
  guides(fill=guide_colorbar(title.position="top", title.hjust=0.5, ticks.colour=NA, barwidth=15), pattern="none") +
  theme(legend.title=element_text(size=14), legend.text=element_text(size=14)) +
  theme(legend.key.size=unit(1.5,"line"), legend.position="top") +
  theme(legend.position="none")
dev.off()


#####################################
### Extract equilibrium densities ###
#####################################

# Cut the datasets before the beads
EquiODE1=subset(SaveODE1, Time==99 & Phase=="Pre")
EquiODE2=subset(SaveODE2, Time==99 & Phase=="Pre")
EquiODE1=cbind(EquiODE1[,c(1:3)], A2=rep(0,209), EquiODE1[,c(4:8)])
EquiODE2=cbind(EquiODE2[,c(1:4)], EquiODE2[,c(6:10)])

# Replace zero values
EquiODE1$A1[EquiODE1$A1 == 0]=0.2
EquiODE1$P[EquiODE1$P == 0]=1.0
EquiODE2$A1[EquiODE2$A1 == 0]=0.1
EquiODE2$A2[EquiODE2$A2 == 0]=0.1
EquiODE2$P[EquiODE2$P == 0]=1.0

# Bind the datasets
DataE=rbind(EquiODE1,EquiODE2)
DataE$Organism=c(rep("1P",209),rep("2P",209))

# Extract equilibrium densities
DataE[,c(2:5)]=round(DataE[,c(2:5)],4)
write.table(DataE, file="~/Activité Professionnelle/LIMNO 2019-2023/Modeling/Data_ED.txt", sep="\t", row.names=F)

# Extract the datasets
Data5[,c(9:14)]=round(Data5[,c(9:14)],4)
Data6[,c(9,11)]=round(Data6[,c(9,11)],4)
write.table(Data5, file="~/Activité Professionnelle/LIMNO 2019-2023/Modeling/Data_MOD1.txt", sep="\t", row.names=F)
write.table(Data6, file="~/Activité Professionnelle/LIMNO 2019-2023/Modeling/Data_MOD2.txt", sep="\t", row.names=F)


#################################################
### Population dynamics for defense and beads ###
#################################################

# Create a dataset
Data8=data.frame(SaveODE1[,c(1:3)], A2=0, SaveODE1[,c(4:10)])
Data9=data.frame(SaveODE2[,c(1:4)], SaveODE2[,c(6:12)])
Data10=rbind(Data8[,c(1:11)],Data9[,c(1:11)])
Data10$System=c(rep("1P",104709),rep("2P",104709))

# Subset the dataset
Data10=subset(Data10, Time %in% c(0:200))

# Rescale densities
Data10$A1=Data10$A1*10^6
Data10$A2=Data10$A2*10^6

# Calculate effective densities
Data10$A=(log(Data10$A1+1)*Data10$P1+log(Data10$A2+1)*1.0)/5

# Calculate the mean defense
Data10$FA1=Data10$A1/(Data10$A1+Data10$A2)
Data10$FA2=Data10$A2/(Data10$A1+Data10$A2)
Data10$Defense=c(Data10[1:42009,]$P1,(Data10[42010:84018,]$FA1*Data10[42010:84018,]$P1+Data10[42010:84018,]$FA2*1.0))

# Melt the datasets
Data10=melt(Data10[,-c(2,13:14)], id.vars=c("Time","Ni","D","P1","P0","B","Defense","Phase","System"), variable.name="Organism", value.name="Density")
Data10$Density=ifelse(Data10$System=="1P" & Data10$Organism=="A2", NA, Data10$Density)

# Extract the dataset
Data10[,c(11)]=round(Data10[,c(11)],4)
write.table(Data10, file="~/Activité Professionnelle/LIMNO 2019-2023/Modeling/Data_MOD3.txt", sep="\t", row.names=F)

# Create the datasets
Data10=subset(Data10, P1 %in% c("0.2","0.4","0.6") & B %in% c("0.02","0.04"))

# Rescale densities
Data10[Data10$Organism=="A1",]$Density=log((Data10[Data10$Organism=="A1",]$Density)+1)/3
Data10[Data10$Organism=="A2",]$Density=log((Data10[Data10$Organism=="A2",]$Density)+1)/3
Data10[Data10$Organism=="P",]$Density=log((Data10[Data10$Organism=="P",]$Density)+1)

# Population dynamic plots
tiff('Population Dynamics No Evolution.tiff', units="in", width=20, height=14, res=1000)
ggplot(subset(Data10, System=="1P" &!Organism=="A"), aes(Time, Density)) +
  geom_line(aes(color=Organism, linetype=Organism), size=2) +
  geom_segment(aes(x=100, y=5.0, xend=100, yend=4.7), color="black", arrow=arrow(length=unit(0.3,"cm"), type="closed"), size=1.5, lineend="round", linejoin="round") +
  ylab(expression('Predator density'~'('*ln~individuals~mL^-1*')')) + xlab(expression('Time (days)')) +    
  theme(axis.text.y=element_text(face="plain", colour="black", size=24)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=24)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=24)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=24)) +
  theme(axis.title.y.left=element_text(face="plain", colour="black", angle=90)) +
  theme(axis.title.y.right=element_text(face="plain", colour="black", angle=90)) +
  scale_y_continuous(labels=function(y) sprintf("%.0f", y), breaks=seq(0,5,by=1), limits=c(0,5),
  sec.axis=sec_axis(~./5*5, expression('Prey density'~'('*ln~cells~mL^-1*')'), labels=sprintf(seq(0,15,by=3), fmt="%.0f"), breaks=seq(0,5,by=1))) +
  scale_x_continuous(labels=function(y) sprintf("%.0f", y), breaks=seq(0,200,by=50), limits=c(0,200)) +
  theme(panel.background=element_blank(), panel.border=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.line=element_line(colour="black", linetype="solid", size=0.5)) +
  scale_color_manual(values=alpha(c("A1"="royalblue2","P"="firebrick3"),0.8)) +
  scale_linetype_manual(values=c("A1"="solid","P"="solid")) +
  theme(panel.spacing.y=unit(1,'lines'), panel.spacing.x=unit(1,'lines')) + 
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color="black", size=0.9) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color="black", size=0.9) +
  facet_wrap(B~P1, ncol=3, nrow=3) +
  theme(legend.position="non")
dev.off()

tiff('Population Dynamics Evolution.tiff', units="in", width=20, height=14, res=1000)
ggplot(subset(Data10, System=="2P" &!Organism=="A"), aes(Time, Density)) +
  geom_line(aes(color=Organism, linetype=Organism), size=2) +
  geom_segment(aes(x=100, y=5.0, xend=100, yend=4.7), color="black", arrow=arrow(length=unit(0.3,"cm"), type="closed"), size=1.5, lineend="round", linejoin="round") +
  ylab(expression('Predator density'~'('*ln~individuals~mL^-1*')')) + xlab(expression('Time (days)')) +    
  theme(axis.text.y=element_text(face="plain", colour="black", size=24)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=24)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=24)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=24)) +
  theme(axis.title.y.left=element_text(face="plain", colour="black", angle=90)) +
  theme(axis.title.y.right=element_text(face="plain", colour="black", angle=90)) +
  scale_y_continuous(labels=function(y) sprintf("%.0f", y), breaks=seq(0,5,by=1), limits=c(0,5),
  sec.axis=sec_axis(~./5*5, expression('Prey density'~'('*ln~cells~mL^-1*')'), labels=sprintf(seq(0,15,by=3), fmt="%.0f"), breaks=seq(0,5,by=1))) +
  scale_x_continuous(labels=function(y) sprintf("%.0f", y), breaks=seq(0,200,by=50), limits=c(0,200)) +
  theme(panel.background=element_blank(), panel.border=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.line=element_line(colour="black", linetype="solid", size=0.5)) +
  scale_color_manual(values=alpha(c("A1"="royalblue2","A2"="royalblue2","P"="firebrick3"),0.8)) +
  scale_linetype_manual(values=c("A1"="solid","A2"="11","P"="solid")) +
  theme(panel.spacing.y=unit(1,'lines'), panel.spacing.x=unit(1,'lines')) + 
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color="black", size=0.9) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color="black", size=0.9) +
  facet_wrap(B~P1, ncol=3, nrow=3) +
  theme(legend.position="none")
dev.off()

# Create a dataset
Data11=rbind(Data8[,c(1:10)],Data9[,c(1:10)])
Data11$A=(Data11$A1+Data11$A2)
Data11$System=c(rep("1P",1809),rep("2P",1809))

# Create elliptic lines
SplitData11=split(Data11, list(Data11$P1,Data11$B,Data11$System))
FuncSpline=function(x) {xspline(x=x$P, y=x$A, shape=0.5, draw=F)}
Data11=as.data.frame(bind_rows(lapply(SplitData11, FuncSpline)))
colnames(Data11)=c("P","A")

# Add columns for identities
SplineDataR=lapply(SplitData11, FuncSpline)
Nrows=unlist(lapply(SplineDataR, function(x) length(x[[1]])))
Data11$P1=c(rep("0.1",Nrows[1]),rep("0.4",Nrows[2]),rep("0.9",Nrows[3]),rep("0.1",Nrows[4]),rep("0.4",Nrows[5]),rep("0.9",Nrows[6]),rep("0.1",Nrows[7]),rep("0.4",Nrows[8]),rep("0.9",Nrows[9]),rep("0.1",Nrows[10]),rep("0.4",Nrows[11]),rep("0.9",Nrows[12]),rep("0.1",Nrows[13]),rep("0.4",Nrows[14]),rep("0.9",Nrows[15]),rep("0.1",Nrows[16]),rep("0.4",Nrows[17]),rep("0.9",Nrows[18]))
Data11$B=c(rep("0.1",Nrows[1]),rep("0.1",Nrows[2]),rep("0.1",Nrows[3]),rep("0.3",Nrows[4]),rep("0.3",Nrows[5]),rep("0.3",Nrows[6]),rep("0.5",Nrows[7]),rep("0.5",Nrows[8]),rep("0.5",Nrows[9]),rep("0.1",Nrows[10]),rep("0.1",Nrows[11]),rep("0.1",Nrows[12]),rep("0.3",Nrows[13]),rep("0.3",Nrows[14]),rep("0.3",Nrows[15]),rep("0.5",Nrows[16]),rep("0.5",Nrows[17]),rep("0.5",Nrows[18]))
Data11$System=c(rep("1P",Nrows[1]),rep("1P",Nrows[2]),rep("1P",Nrows[3]),rep("1P",Nrows[4]),rep("1P",Nrows[5]),rep("1P",Nrows[6]),rep("1P",Nrows[7]),rep("1P",Nrows[8]),rep("1P",Nrows[9]),rep("2P",Nrows[10]),rep("2P",Nrows[11]),rep("2P",Nrows[12]),rep("2P",Nrows[13]),rep("2P",Nrows[14]),rep("2P",Nrows[15]),rep("2P",Nrows[16]),rep("2P",Nrows[17]),rep("2P",Nrows[18]))

# Phase plane plots
ggplot(subset(Data11, System=="1P"), aes(P, A)) +
  geom_path(color="darkorange2", linetype="solid", size=2) +
  ylab(expression('Prey density'~'('*ln~cells~mL^-1*')')) + 
  xlab(expression('Predator density'~'('*ln~individuals~mL^-1*')')) +  
  theme(axis.text.y=element_text(face="plain", colour="black", size=24)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=24)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=24)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=24)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,50,by=10), limits=c(0,50)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,60,by=10), limits=c(0,60)) +
  theme(panel.background=element_blank(), panel.border=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.line=element_line(colour="black", linetype="solid", size=0.5)) +
  theme(panel.spacing.y=unit(1,'lines'), panel.spacing.x=unit(1,'lines')) + 
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color="black", size=0.9) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color="black", size=0.9) +
  facet_wrap(B~P1, ncol=3, nrow=3) +
  theme(legend.position="none")



#####################################################
### Comparing ecology and evolution contributions ###
#####################################################

# Create a dataset
Data12=subset(SaveODE2, Time %in% c(0:200))

# Rescale densities
Data12$A=Data12$A*10^6
Data12$A1=Data12$A1*10^6
Data12$A2=Data12$A2*10^6

# Calculate the mean defense
Data12$FA1=Data12$A1/(Data12$A1+Data12$A2)
Data12$FA2=Data12$A2/(Data12$A1+Data12$A2)
Data12$Defense=c(Data12[1:42009,]$FA1*Data12[1:42009,]$P1+Data12[1:42009,]$FA2*1.0)

# Melt the datasets
Data12a=melt(Data12[,-c(2,11)], id.vars=c("Time","Ni","D","P1","P0","B","A","P","FA1","FA2","Defense"), variable.name="Organism", value.name="Density")
Data12b=melt(Data12[,-c(2,11)], id.vars=c("Time","Ni","D","P1","P0","B","A1","A2","A","P","Defense"), variable.name="Organism", value.name="Frequency")

# Combine the datasets
Data13=data.frame(Data12a[,c(1:6,12)], A=Data12a[,13], R=Data12a[,8], Defense=Data12a[,11], Frequency=Data12b[,13])
Data13[,c(8:9)]=round(Data13[,c(8:9)],0)

# Include prey mean defense and frequency
colnames(Data13)[8:9]=c("DensA","DensR")
Data13[,c(10:11)]=round(Data13[,c(10:11)],2)

# Create a dataset
Data14=setDT(Data13)[, .(DensA=round(mean(DensA),0), DensR=round(mean(DensR),0)), by=list(Time,Ni,D,P1,P0,B,Organism,Defense,Frequency)]
Data14=as.data.frame(Data14)

# Unmelt the dataset
Data15=subset(Data14, Organism=="A")
Data15a=subset(Data14, Organism=="A1")
Data15b=subset(Data14, Organism=="A2")

# Rename column names
colnames(Data15a)[c(8:10)]=paste(colnames(Data15a)[c(8:10)], "DEF", sep="")
colnames(Data15b)[c(8:10)]=paste(colnames(Data15b)[c(8:10)], "UND", sep="")

# Create a dataset
Data16=cbind(Data15a[,c(1:6,8:10)],Data15b[,c(8:11)])

# Arrange the datasets
Data16=Data16 %>% arrange(factor(P1, levels=unique(P1)), factor(B, levels=unique(B)), factor(Time, levels=seq(0,200,by=1)))

# Split the dataset
SplitData16=split(Data16, list(Data16$P1,Data16$B))
SplitData16=SplitData16[sapply(SplitData16, function(x) dim(x)[1]) > 0]

# Subset important columns
SplitData16=lapply(SplitData16, "[", c("Time","DensR","DefenseDEF"))

# Instantaneous predator growth rate function
FuncGR=function(x) {
  FitGR=dynlm(formula=log(DensR+1)~L(Time), data=as.data.frame(x))
  OutGR=c(RateR=summary(FitGR)$coef[2,1])}

# Calculate instantaneous predator growth rates
ModelGR=function(y) {rollapply(y, width=2, FUN=FuncGR, by.column=F)}

# Calculate consecutive predator growth rates
OutGR=round(as.data.frame(do.call("cbind",lapply(SplitData16, ModelGR))),4)
GrowR=data.frame(OutGR=unlist(OutGR,use.names=F))[,1]

# Include predator growth rates
Data16$GrowR=c(t(cbind(rep(NA,209), matrix(GrowR, ncol=200, byrow=T))))
Data16$CombA=Data16$FrequencyDEF*(Data16$DensADEF+Data16$DensAUND)

# Split the dataset
SplitData16=split(Data16, list(Data16$P1,Data16$B))
SplitData16=SplitData16[sapply(SplitData16, function(x) dim(x)[1]) > 0]

# Eco-evolutionary models
ModelEE1=function(x) {lm(x$GrowR~x$CombA)}
OutEE1=lapply(SplitData16, ModelEE1)

# Extract coefficients
CoefEE=lapply(OutEE1,function (x) {c(Coeff1=coef(summary(x))[1,2], Coeff2=coef(summary(x))[2,2])})
CoefEE=as.data.frame(do.call("rbind",CoefEE))
Data16$Coeff1=rep(CoefEE[,1], each=201)
Data16$Coeff2=rep(CoefEE[,2], each=201)

# Split the dataset
SplitData16=split(Data16, list(Data16$P1,Data16$B))
SplitData16=SplitData16[sapply(SplitData16, function(x) dim(x)[1]) > 0]

# Extract combinations of names
Names=unique(Data16[,c("P1","P0","B","Time")])
Names=Names[order(Names$B,Names$P1),]
P1=Names$P1; P0=Names$P0; B=Names$B

# Calculate the contributions
ModelEE2=function(x) {
  a=x[,15] + x[,16] * x[,8] * (x[,9] + x[,12])
  b=x[,15] + x[,16] * lead(x[,8]) * (x[,9] + x[,12])
  c=x[,15] + x[,16] * x[,8] * (lead(x[,9]) + lead(x[,12]))
  d=x[,15] + x[,16] * lead(x[,8]) * (lead(x[,9]) + lead(x[,12]))
  Eco=((b-a)+(d-c))/2
  Evo=((c-a)+(d-b))/2
  Coef=data.frame(Eco=Eco, Evo=Evo, EcoEvo=Eco/Evo)
}
OutEE2=lapply(SplitData16, ModelEE2)

RateEE=round(as.data.frame(do.call("rbind",OutEE2)),4)
RateEE=cbind(P1=Data16[,4],P0=Data16[,5],B=Data16[,6],Time=Data16[,1],RateEE)
rownames(RateEE)=c()

# Include the contributions
Data16$Eco=RateEE[,5]
Data16$Evo=RateEE[,6]
Data16$EcoEvo=RateEE[,7]

# Correct the contributions
Data16$EcoEvo=log(abs(Data16$EcoEvo))

# Replace infinite values
Data16$Eco=ifelse(is.infinite(Data16$Eco), 0, Data16$Eco)
Data16$Evo=ifelse(is.infinite(Data16$Evo), 0, Data16$Evo)
Data16$EcoEvo=ifelse(is.infinite(Data16$EcoEvo), 0, Data16$EcoEvo)

# Include weight for lines
Data16$Weight=ifelse(Data16$Time==0, 1000, 1)

# Specify the variables as numeric or factor
Data16[,c(2:6)] %<>% mutate_if(is.integer,as.factor)

tiff('Contribution Ecology Evolution.tiff', units="in", width=12, height=20, res=1000)
ggplot(Data16, aes(Time, EcoEvo)) +
  geom_hline(yintercept=0, color="black", linetype="11", size=1.5) + 
  geom_smooth(aes(weight=Weight), method="loess", color=alpha("royalblue2",0.8), linetype="solid", size=2, span=0.4, se=F) +
  ylab(expression('Ln Ecology / Evolution')) + xlab(expression('Time (days)')) +    
  theme(axis.text.y=element_text(face="plain", colour="black", size=12)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=12)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=12)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=12)) +
  theme(axis.title.y.left=element_text(angle=90), axis.title.y.right=element_text(angle=90)) +
  scale_y_continuous(labels=function(y) sprintf("%.1f", y), breaks=seq(-3.0,3.0,by=1.5), limits=c(-3.0,3.0)) +
  scale_x_continuous(labels=function(y) sprintf("%.0f", y), breaks=seq(0,50,by=10), limits=c(0,50)) +
  theme(panel.background=element_blank(), panel.border=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.line=element_line(colour="black", linetype="solid", size=0.5)) +
  theme(panel.spacing.y=unit(1,'lines'), panel.spacing.x=unit(1,'lines')) + 
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color="black", size=0.9) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color="black", size=0.9) +
  facet_wrap(P1~B, ncol=11, nrow=19) +
  theme(legend.position="none")
dev.off()
