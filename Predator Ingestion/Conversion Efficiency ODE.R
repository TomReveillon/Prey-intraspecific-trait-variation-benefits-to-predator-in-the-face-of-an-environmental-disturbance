setwd("~/LIMNO 2019-2023/Experiments/Predator Ingestion Beads")

rm(list=ls())

library(bbmle)
library(cowplot)
library(data.table)
library(deSolve)
library(directlabels)
library(dplyr)
library(dynlm)
library(foreach)
library(ggplot2)
library(ggrepel)
library(grid)
library(gridExtra)
library(lme4)
library(lmtest)
library(magrittr)
library(nlme)
library(odeintr)
library(plotly)
library(plyr)
library(propagate)
library(reshape2)
library(rootSolve)

#####################################################################################
#####################################################################################
##### CONVERSION EFFICIENCY FOR PREDATOR-PREY SYSTEM WITH DIFFERENT CLONE TYPES #####
#####################################################################################
#####################################################################################

############################
### Ingestion parameters ###
############################

# Import the dataset
DataFRB=read.table("~/Activité Professionnelle/LIMNO 2019-2023/Experiments/Predator Ingestion Beads/Data_FRBPODE.txt", h=T, dec=",")
summary(DataFRB)
names(DataFRB)

# Specify the variables as numeric or factor
DataFRB[,c(3:6)] %<>% mutate_if(is.character,as.numeric)

# Arrange the dataset
DataFRB$Strain=factor(DataFRB$Strain, levels=unique(DataFRB$Strain))
DataFRB$Bead=factor(DataFRB$Bead, levels=unique(DataFRB$Bead))

# Split the dataset
SplitDataFRB=split(DataFRB, list(DataFRB$Strain,DataFRB$Bead))

# Extract combinations of names
Strain=unique(DataFRB[,c("Strain","Bead")])[,1]
Bead=unique(DataFRB[,c("Strain","Bead")])[,2]

# Functional response model
FuncFRB=function(x) {
  ModFRB=nls(IngesP ~ (a * IDensP) / (1 + a * h * IDensP), start=c(a=1.0, h=0.1), data=x)
  ModFRBL=nls(IngesPLSD ~ (a * IDensP) / (1 + a * h * IDensP), start=c(a=1.0, h=0.1), data=x)
  ModFRBU=nls(IngesPUSD ~ (a * IDensP) / (1 + a * h * IDensP), start=c(a=1.0, h=0.1), data=x)
  AttackB=data.frame(AttackB=coef(ModFRB)[1], AttackBL=coef(ModFRBL)[1], AttackBU=coef(ModFRBU)[1])
  HandlingB=data.frame(HandlingB=coef(ModFRB)[2], HandlingBL=coef(ModFRBU)[2], HandlingBU=coef(ModFRBL)[2])
  Parameters=list(AttackB=AttackB, HandlingB=HandlingB)
}
OutFRB=lapply(SplitDataFRB, FuncFRB)

# Calculate the attack rates
AttacksB=bind_rows(lapply(OutFRB, function (x) x[c("AttackB")]))
AttacksB=round(as.data.frame(do.call("rbind",AttacksB)),4)
AttacksB=cbind(Strain=Strain,Bead=Bead,AttacksB)
rownames(AttacksB)=c()

# Calculate the handling times
HandlingsB=bind_rows(lapply(OutFRB, function (x) x[c("HandlingB")]))
HandlingsB=round(as.data.frame(do.call("rbind",HandlingsB)),4)
HandlingsB=cbind(Strain=Strain,Bead=Bead,HandlingsB)
rownames(HandlingsB)=c()


#########################
### Growth parameters ###
#########################

# Import the dataset
DataRGB=read.table("~/Activité Professionnelle/LIMNO 2019-2023/Experiments/Predator Growth Beads/Data_RGB.txt", h=T, dec=",")
summary(DataRGB)
names(DataRGB)

# Specify the variables as numeric or factor
DataRGB[,c(3:6)] %<>% mutate_if(is.character,as.numeric)

# Arrange the dataset
DataRGB=DataRGB %>% arrange(factor(Bead, levels=c("C","L","M","H")), factor(Strain, levels=c("CR1","CR2","CR3","CR4","CR5","CR6")))
DataRGB$Strain=factor(DataRGB$Strain, levels=unique(DataRGB$Strain))
DataRGB$Bead=factor(DataRGB$Bead, levels=unique(DataRGB$Bead))

# Calculate difference between initial and final densities
DataRGB2=setDT(DataRGB)[, .(Diffe=as.numeric(tail(Count,1) - head(Count,1))), by=list(Strain,Bead,Trial)]
DataRGB2=as.data.frame(DataRGB2)


####################################
### Conversion efficiency models ###
####################################

# Functional response parameters
DensB=c(0,0.03125,0.0625,0.1250)
C=c(11.07,11.07,11.07,11.07,11.07,11.07)

# Convert functional response parameters
A=round((AttacksB[,3]*60*60*24)/(10^6),6)
H=round((HandlingsB[,3]/60/60/24)*(10^6),6)

# Create a dataset
Data=subset(DataRGB, !Day==0)
Data2=subset(DataRGB, !Day==5)
Data3=data.frame(Strain=Data$Strain, Bead=Data$Bead, Day=Data$Day, Trial=Data$Trial, DensR0=Data2$Count, DensR=Data$Count, DensA0=rep(c(rep(0.5,5),rep(NA,20)),6*4), DensB=rep(DensB,each=150), a=rep(A,each=25), h=rep(H,each=25), c=rep(C,each=25*4), m=rep(0.2,each=25*4))
Data3$Inges=rep(NA,600)
Data3$DDensA=rep(NA,600)

# Split the dataset
SplitData3=split(Data3, list(Data3$Strain,Data3$Bead,Data3$Trial))

# Calculate ingestion rate
FuncIR=function(x) {for (i in 1:(nrow(x))) {  
  x$Inges[i]=(x$a[i] * x$DensA0[i])/(1 + x$a[i] * x$h[i] * x$DensA0[i] + x$c[i] * x$DensB[i])
  x$DDensA[i]=x$DensA0[i] - x$Inges[i] * x$DensR0[i]
  x[i+1,7]=x$DDensA[i]}
  return(as.data.frame(x))}

# Bind the lists
Data4=bind_rows(lapply(SplitData3, FuncIR))
Data4=Data4[complete.cases(Data4),]
Data4[,c(7,13:14)]=round(Data4[,c(7,13:14)],6)
rownames(Data4)=c()

# Order the dataset
Data4=Data4 %>% arrange(factor(Strain, levels=c("CR1","CR2","CR3","CR4","CR5","CR6")))
Data4=Data4 %>% arrange(factor(Bead, levels=c("C","L","M","H")))
Data4$Strain=factor(Data4$Strain, levels=unique(Data4$Strain))
Data4$Bead=factor(Data4$Bead, levels=unique(Data4$Bead))

# Include decreasing prey densities
Data4$DensA=as.data.frame(subset(Data4, !Day==1) %>% group_by(Strain,Bead,Trial) %>% group_modify(~add_row(.x,.after=4)))[,7]
Data4=Data4[,c(1:7,15,8:14)]

# Create a dataset
Data5=data.frame(Strain=Data4[,1], Bead=Data4[,2], Trial=Data4[,4], Time=rep(5,600), a=Data4[,10], h=Data4[,11], IDensA=Data4[,7], FDensA=Data4[,8], IDensR=rep(5,600), FDensR=Data4[,6])
Data5=na.omit(Data5)

# Split the dataset
SplitData=split(Data5, list(Data5$Strain,Data5$Bead))
SplitData=SplitData[sapply(SplitData, function(x) dim(x)[1]) > 0]

# Extract combinations of names
Strain=lapply(SplitData, function(x) {unique(x$Strain)})
Bead=lapply(SplitData, function(x) {unique(x$Bead)})
Strain=as.character(unlist(Strain))
Bead=as.character(unlist(Bead))
Param=unique(c("Xp","sigma"))

# Create an integrating formula
Repro='dxdt[0] = (-a * x[0]) / (1 + a * h * x[0]) * x[1];
       dxdt[1] = x[1] * Xp * (a * x[0]) / (1 + a * h * x[0]) - x[1] * x[3];'
compile_sys("Repro", Repro, pars=c("a","h","Xp"), method="rk54")

# Functional response function
Inges=function(Start, a, h, Xp, Time, Time_Length){
  Repro_set_params(a=a, h=h, Xp=Xp)
  Repro(Start, Time, Time_Length)
}

# Densities depletion function
DDensA=c(); NDensR=c()
DensEaten=function(IDensA, a, h, Xp, m, Time, P, steps=100){
  for(i in 1:length(IDensA)){
    DEaten=Inges(Start=c(IDensA[i], P[i]), a=a[i], h=h[i], Xp=Xp, Time=Time[i], Time_Length=Time[i]/steps)
    DDensA[i]=IDensA[i] - DEaten[steps+1,2]
    NDensR[i]=DEaten[steps+1,3] - P[i] - (P[i] * m[i])
  }
  Out=list(DDensA, NDensR, DEaten)
  return(Out)
}

# Maximum likelihood function
Likelihood=function(NDeadA, IDensA, a, h, Xp, m, Time, P, PEnd, sigma, steps=100){
  if(Xp <= 0 || sigma <= 0) return(Inf)
  DEaten2=DensEaten(IDensA=IDensA, a=a, h=h, Xp=Xp, m=m, Time=Time, P=P, steps=steps)
  LR1=-1*sum(dnorm(x=log((IDensA-NDeadA)+1), mean=log((IDensA-DEaten2[[1]])+1), sd=sigma, log=T))
  LR2=-1*sum(dpois(x=PEnd, lambda=P+DEaten2[[2]], log=T))
  LR=LR1+LR2
  return(LR)
}

# Calculate number of decreasing predator densities
SplitDataS=lapply(SplitData, function(x) {x[x$FDensR < 5,]})
Missing=lapply(SplitDataS, function(x) {length(x[,1])})
Missing=data.frame(Strain=Strain,Bead=Bead,Count=do.call("rbind",Missing))

# Remove decreasing predator densities
SplitData=lapply(SplitData, function(x) {x[x$FDensR > 4,]})

# Fitting the model
FuncHII1=function(x) {ModHII=mle2(Likelihood, start=list(Xp=200, sigma=1), control=list(maxit=1000), data=list(a=x$a, h=x$h, m=x$m, IDensA=x$IDensA, NDeadA=x$IDensA-x$FDensA, P=x$IDensR, PEnd=x$FDensR, Time=x$Time))}
OutHII1=lapply(SplitData[c(1:12)], FuncHII1)

FuncHII2=function(x) {ModHII=mle2(Likelihood, start=list(Xp=20, sigma=1), control=list(maxit=1000), data=list(a=x$a, h=x$h, m=x$m, IDensA=x$IDensA, NDeadA=x$IDensA-x$FDensA, P=x$IDensR, PEnd=x$FDensR, Time=x$Time))}
OutHII2=lapply(SplitData[c(13:18)], FuncHII2)

CoefHII=lapply(c(OutHII1,OutHII2), summary)
CoefHII=lapply(CoefHII, coef)
CoefHII=round(as.data.frame(do.call("rbind",CoefHII)),4)
CoefHII=data.frame(Strain=rep(Strain[c(1:12,13:18,19:24)],each=2),Bead=rep(Bead[c(1:12,13:18,19:24)],each=2),Param=rep(Param,6*4),Value=c(CoefHII[,c(1)],rep(0,6*2)),ValueSD=c(CoefHII[,c(2)],rep(0,6*2)))
CoefHII$ValueLSD=CoefHII[,4]-CoefHII[,5]
CoefHII$ValueUSD=CoefHII[,4]+CoefHII[,5]
rownames(CoefHII)=c()

# Specify the variables as numeric or factor
CoefHII[,c(4:7)] %<>% mutate_if(is.character,as.numeric)
CoefHII=subset(CoefHII, Param=="Xp")[,c(1:7)]
colnames(CoefHII)[4:7]=c("Xp","XpSD","XpLSD","XpUSD")

# Pondering conversion efficiencies
CoefHII=CoefHII[order(CoefHII$Bead),]
Missing=Missing[order(Missing$Bead),]
CoefHII$Absent=Missing$Count
CoefHII$Present=20-Missing$Count
CoefHII$XpC=(CoefHII[,4]*CoefHII[,9]+0*CoefHII[,8])/(CoefHII[,8]+CoefHII[,9])

# Create a dataset
Data6=data.frame(Strain=CoefHII[,1], Bead=CoefHII[,2], DensB=rep(DensB[c(1,4,2,3)], each=6), ConvP=CoefHII[,10], ConvPLSD=CoefHII[,10]-CoefHII[,5], ConvPUSD=CoefHII[,10]+CoefHII[,5])
Data6[,c(3:6)]=round(Data6[,c(3:6)],4)
Data6[,c(3:6)][Data6[,c(3:6)]<0]=0

# Arrange the dataset
Data6=Data6 %>% arrange(factor(Bead, levels=c("C","L","M","H")))


###############################
### Fitting elliptic models ###
###############################

# Split the dataset
SplitData6=split(Data6, list(Data6$Strain))

# Elliptic model
plot(Data6$ConvP~Data6$DensB, ylab=c("Conversion efficiency"), xlab=c("Particle concentration"), pch=16)
ModelCE=function(x) {xspline(x=x$DensB, y=x$ConvP, shape=0.8, draw=F)}
ModelCELSD=function(x) {xspline(x=x$DensB, y=x$ConvPLSD, shape=0.8, draw=F)}
ModelCEUSD=function(x) {xspline(x=x$DensB, y=x$ConvPUSD, shape=0.8, draw=F)}

# Calculate the ellipses
SplitData7=lapply(SplitData6, ModelCE)
SplitData8=lapply(SplitData6, ModelCELSD)
SplitData9=lapply(SplitData6, ModelCEUSD)

# Extract lengths of ellipses
Time1=do.call("rbind",lapply(lapply(SplitData7, function (x) x[[c("y")]]), function (x) length(x)))
Time2=do.call("rbind",lapply(lapply(SplitData8, function (x) x[[c("y")]]), function (x) length(x)))
Time3=do.call("rbind",lapply(lapply(SplitData9, function (x) x[[c("y")]]), function (x) length(x)))


################################################
### Plotting predicted conversion efficiency ###
################################################

# Create the datasets
Strain1=c(rep("CR1",Time1[1]),rep("CR2",Time1[2]),rep("CR3",Time1[3]),rep("CR4",Time1[4]),rep("CR5",Time1[5]),rep("CR6",Time1[6]))
Strain2=c(rep("CR1",Time2[1]),rep("CR2",Time2[2]),rep("CR3",Time2[3]),rep("CR4",Time2[4]),rep("CR5",Time2[5]),rep("CR6",Time2[6]))
Strain3=c(rep("CR1",Time3[1]),rep("CR2",Time3[2]),rep("CR3",Time3[3]),rep("CR4",Time3[4]),rep("CR5",Time3[5]),rep("CR6",Time3[6]))
Data7=data.frame(Strain=Strain1, DensBP=bind_rows(SplitData7)[[1]], ConvP=bind_rows(SplitData7)[[2]])
Data8=data.frame(Strain=Strain2, DensBP=bind_rows(SplitData8)[[1]], ConvP=bind_rows(SplitData8)[[2]])
Data9=data.frame(Strain=Strain3, DensBP=bind_rows(SplitData9)[[1]], ConvP=bind_rows(SplitData9)[[2]])

# Split the datasets
SplitData7=split(Data7, list(Data7$Strain))
SplitData8=split(Data8, list(Data8$Strain))
SplitData9=split(Data9, list(Data9$Strain))

# Rescale the datasets
FuncScale=function(x) {x=data.frame(Strain=rep(unique(x$Strain),125), DensBP=approx(x$DensBP,n=125)$y, ConvP=approx(x$ConvP,n=125)$y)}
SplitData7=lapply(SplitData7, FuncScale)
SplitData8=lapply(SplitData8, FuncScale)
SplitData9=lapply(SplitData9, FuncScale)

# Combine the datasets
Data10=data.frame(do.call("rbind",SplitData7)[,c(1:3)], ConvPLSD=do.call("rbind",SplitData8)[,3], ConvPUSD=do.call("rbind",SplitData9)[,3])
Data10[,c(2:5)][is.na(Data10[,c(2:5)])]=0
Data10[,c(2:5)]=round(Data10[,c(2:5)],6)
Data10[,c(2:5)][Data10[,c(2:5)]<0]=0
rownames(Data10)=c()

# Correct standard errors and confidence intervals
Data10$ConvPLSD=ifelse(Data10$ConvPLSD==0, Data10$ConvP-Data10$ConvP*((lead(Data10$ConvP)-lead(Data10$ConvPLSD))/lead(Data10$ConvP)), Data10$ConvPLSD)
Data10$ConvPUSD=ifelse(Data10$ConvPUSD==0, Data10$ConvP+Data10$ConvP*((lead(Data10$ConvPUSD)-lead(Data10$ConvP))/lead(Data10$ConvP)), Data10$ConvPUSD)

# Correct standard errors and confidence intervals
Data10$ConvPLSD=ifelse(Data10$ConvPLSD < Data10$ConvP*0.80, Data10$ConvP*0.80, Data10$ConvPLSD)
Data10$ConvPUSD=ifelse(Data10$ConvPUSD > Data10$ConvP*1.20, Data10$ConvP*1.20, Data10$ConvPUSD)

# Export the dataset
Data10[,c(3:5)]=replace(Data10[,c(3:5)],Data10[,c(3:5)]<0,0)
write.table(Data10, file="Data_CEB.txt", sep="\t", row.names=F)

tiff('Conversion Efficiency Beads.tiff', units="in", width=8, height=8, res=1000)
ggplot(Data10, aes(DensBP*10, ConvP, group=Strain)) +
  geom_line(aes(color=Strain), linetype="solid", size=1) +
  geom_point(data=Data6, aes(DensB*10, ConvP, color=Strain), size=3, pch=16) +
  ylab(expression(italic('B. calyciflorus')~'conversion efficiency'~'('*rotifers~10^-6~cells*')')) +
  xlab(expression('Microplastic density'~'('*10^5~particles~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,160,by=40), limits=c(0,160)) +
  scale_x_continuous(labels=function(x) sprintf("%.2f", x), breaks=seq(0,1.25,by=0.25), limits=c(0,1.25)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="steelblue2","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
  theme(legend.position="none")
dev.off()

tiff('Conversion Efficiency Beads Intervals.tiff', units="in", width=15, height=8, res=1000)
ggplot(Data10, aes(DensBP*10, ConvP, group=Strain)) +
  geom_line(aes(color=Strain), linetype="solid", size=1) +
  geom_point(data=Data6, aes(DensB*10, ConvP, color=Strain), size=3, pch=16) +
  geom_ribbon(aes(ymin=ConvPLSD, ymax=ConvPUSD, fill=Strain, color=Strain), linetype="solid", size=0.5, alpha=0.3) +
  ylab(expression(italic('B. calyciflorus')~'conversion efficiency'~'('*rotifers~10^-6~cells*')')) +
  xlab(expression('Microplastic density'~'('*10^5~particles~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,160,by=40), limits=c(0,160)) +
  scale_x_continuous(labels=function(x) sprintf("%.2f", x), breaks=seq(0,1.25,by=0.25), limits=c(0,1.25)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_fill_manual(values=c("CR1"="mediumpurple3","CR2"="steelblue2","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
  scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="steelblue2","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
  theme(strip.text.x=element_blank()) +
  facet_wrap(~Strain, scales="free", ncol=3, nrow=2) +
  theme(legend.position="none")
dev.off()
