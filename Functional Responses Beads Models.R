setwd("~/LIMNO 2019-2022/Experiments/Predator Ingestion Beads")

rm(list=ls())

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
library(plotly)
library(plyr)
library(propagate)
library(reshape2)
library(scales)

##########################################################################################
##########################################################################################
##### FUNCTIONAL RESPONSE MODELS FOR PREDATOR-PREY SYSTEM WITH DIFFERENT CLONE TYPES #####
##########################################################################################
##########################################################################################

# Import the dataset
DataI=read.table("Data_FRBI.txt", h=T, dec=",")
names(DataI)
summary(DataI)

# Specify the variables as numeric numbers
DataI[,c(3,8)] %<>% mutate_if(is.numeric,as.character())
DataI[,c(4:7,9)] %<>% mutate_if(is.character,as.numeric)

# Preserve the order of factors
DataI$Bead=factor(DataI$Bead, levels=unique(DataI$Bead))

# Calculate the initial densities
DataI$IDens=round(DataI$Cells*DataI$Volu*DataI$Site*DataI$Dilu*DataI$Cove,0)
DataI=ddply(DataI, .(Strain,Conc,Bead), summarize, IDens=round(mean(IDens),0))

# Import the dataset
DataF=read.table("Data_FRBF.txt", h=T, dec=",")
names(DataF)
summary(DataF)

# Specify the variables as numeric numbers
DataF[,c(3,8)] %<>% mutate_if(is.numeric,as.character())
DataF[,c(4:7,9)] %<>% mutate_if(is.character,as.numeric)

# Preserve the order of factors
DataF$Bead=factor(DataF$Bead, levels=unique(DataF$Bead))

# Calculate the final densities
DataF$FDens=round(DataF$Cells*DataF$Volu*DataF$Site*DataF$Dilu*DataF$Cove,0)
DataF=ddply(DataF, .(Strain,Conc,Bead,Trial), summarize, FDens=round(mean(FDens),0))

# Create a complete dataset
Data=data.frame(Strain=DataF[,1], Conc=DataF[,2], Bead=DataF[,3], Trial=DataF[,4], IDens=rep(DataI[,4], each=3), FDens=DataF[,5])

# Convert density into ingestion rate
Data$Inges=(Data$IDens-Data$FDens)/(8*60*60)

# Calculate ingestion rate per rotifer
Data$Inges=Data$Inges/4

# Correct ingestion rate per volume
Data$Inges=Data$Inges/5

# Replace negative values by 0 values
Data$Inges=round(Data$Inges,4)
Data$Inges[Data$Inges<0]=0
Data[is.na(Data)]=0

# Calculate mean ingestion rates
Data2=setDT(Data)[, .(MeanInges=mean(Inges), IDens=mean(IDens)), by=list(Strain,Conc,Bead)]
Data2=as.data.frame(Data2)
Data2$IDens=Data2$IDens/10^5


##################################
### Functional response models ###
##################################

# Calculate bead densities
IBead=rep(NA,length(Data2[,1]))
for (i in 1:length(Data2[,1])) {
  if (Data2[i,3]=="C") Data2$IBead[i]=0.000*Data2[i,2]
  if (Data2[i,3]=="L") Data2$IBead[i]=0.063*Data2[i,2]
  if (Data2[i,3]=="M") Data2$IBead[i]=0.125*Data2[i,2]
  if (Data2[i,3]=="H") Data2$IBead[i]=0.250*Data2[i,2]
}

# Calculate total particle densities
Data2$IPart=Data2$IDens+Data2$IBead

# Split the dataset
SplitData2=split(Data2, list(Data2$Strain,Data2$Bead))

# Extract combinations of names
Names=unique(Data2[,c("Strain","Bead","Conc")])
Names=Names[order(Names$Bead,Names$Strain),]
Strain=Names$Strain; Bead=Names$Bead

# Holling I model
FuncHI=function(x, ModHI) {
  ModHI=nls(MeanInges ~ b + a * IDens, start=c(b=0.1, a=1.0), data=x)
  Param=coef(ModHI)
  b=Param[1]
  a=Param[2]

  Fit=ModHI$m$fitted()
  Res=ModHI$m$resid()
  MSS=sum((Fit-mean(Fit))^2)
  RSS=sum(Res^2)
  R2=MSS/(MSS+RSS)
  AIC=AIC(ModHI)
  
  Summary=data.frame(R2,AIC)
  Parameters=data.frame(b,a)
  Rates=data.frame(
    Dens = x$IDens,
    IngesP = b + a * x$IDens)
  
  Out=list(Model=ModHI, Summary=Summary, Parameters=Parameters, Rates=Rates)
  return(Out)
}
OutHI=lapply(SplitData2, FuncHI)

RateHI=bind_rows(lapply(OutHI, function (x) x[c("Rates")]))
RateHI=round(as.data.frame(do.call("rbind",RateHI)),4)
RateHI=cbind(Strain=Strain[c(1:144)],Bead=Bead[c(1:144)],Model="Holling I",RateHI)
rownames(RateHI)=c()
ParamHI=bind_rows(lapply(OutHI, function (x) x[c("Parameters")]))
ParamHI=round(as.data.frame(do.call("rbind",ParamHI)),4)
rownames(ParamHI)=c()
SummaHI=bind_rows(lapply(OutHI, function (x) x[c("Summary")]))
SummaHI=round(as.data.frame(do.call("rbind",SummaHI)),4)
rownames(SummaHI)=c()
ModelHI=unlist(lapply(OutHI, function (x) x[c("Model")]),recursive=F)


# Holling II model
FuncHII=function(x, ModHII) {
  ModHII=nls(MeanInges ~ (a * IDens) / (1 + a * h * IDens), start=c(a=1.0, h=0.1), data=x)
  Param=coef(ModHII)
  a=Param[1]
  h=Param[2]

  Fit=ModHII$m$fitted()
  Res=ModHII$m$resid()
  MSS=sum((Fit-mean(Fit))^2)
  RSS=sum(Res^2)
  R2=MSS/(MSS+RSS)
  AIC=AIC(ModHII)
  
  Summary=data.frame(R2,AIC)
  Parameters=data.frame(a,h)
  Rates=data.frame(
    Dens = x$IDens,
    IngesP = (a * x$IDens) / (1 + a * h * x$IDens))
  
  Out=list(Model=ModHII, Summary=Summary, Parameters=Parameters, Rates=Rates)
  return(Out)
}
OutHII=lapply(SplitData2, FuncHII)

RateHII=bind_rows(lapply(OutHII, function (x) x[c("Rates")]))
RateHII=round(as.data.frame(do.call("rbind",RateHII)),4)
RateHII=cbind(Strain=Strain[c(1:144)],Bead=Bead[c(1:144)],Model="Holling II",RateHII)
rownames(RateHII)=c()
ParamHII=bind_rows(lapply(OutHII, function (x) x[c("Parameters")]))
ParamHII=round(as.data.frame(do.call("rbind",ParamHII)),4)
rownames(ParamHII)=c()
SummaHII=bind_rows(lapply(OutHII, function (x) x[c("Summary")]))
SummaHII=round(as.data.frame(do.call("rbind",SummaHII)),4)
rownames(SummaHII)=c()
ModelHII=unlist(lapply(OutHII, function (x) x[c("Model")]),recursive=F)


# Holling III model
FuncHIII=function(x, ModHIII) {
  ModHIII=nls(MeanInges ~ (a * IDens^c) / (1 + a * h * IDens^c), start=c(a=1.0, h=0.1, c=1.0), data=x)
  Param=coef(ModHIII)
  a=Param[1]
  h=Param[2]
  c=Param[3]
  
  Fit=ModHIII$m$fitted()
  Res=ModHIII$m$resid()
  MSS=sum((Fit-mean(Fit))^2)
  RSS=sum(Res^2)
  R2=MSS/(MSS+RSS)
  AIC=AIC(ModHIII)
  
  Summary=data.frame(R2,AIC)
  Parameters=data.frame(a,h,c)
  Rates=data.frame(
    Dens = x$IDens,
    IngesP = (a * x$IDens^c) / (1 + a * h * x$IDens^c))
  
  Out=list(Model=ModHIII, Summary=Summary, Parameters=Parameters, Rates=Rates)
  return(Out)
}
OutHIII=lapply(SplitData2[-c(1,7,9,12,14,17:24)], FuncHIII)

RateHIII=bind_rows(lapply(OutHIII, function (x) x[c("Rates")]))
RateHIII=round(as.data.frame(do.call("rbind",RateHIII)),4)
RateHIII=cbind(Strain=Strain[-c(1:6,37:42,49:54,67:72,79:84,97:144)],Bead=Bead[-c(1:6,37:42,49:54,67:72,79:84,97:144)],Model="Holling III",RateHIII)
rownames(RateHIII)=c()
ParamHIII=bind_rows(lapply(OutHIII, function (x) x[c("Parameters")]))
ParamHIII=round(as.data.frame(do.call("rbind",ParamHIII)),4)
rownames(ParamHIII)=c()
SummaHIII=bind_rows(lapply(OutHIII, function (x) x[c("Summary")]))
SummaHIII=round(as.data.frame(do.call("rbind",SummaHIII)),4)
rownames(SummaHIII)=c()
ModelHIII=unlist(lapply(OutHIII, function (x) x[c("Model")]),recursive=F)


##########################################
### Fitting functional response models ###
##########################################

# Set the predicted dataset
Strain=rep(rep(unique(Data2$Strain), each=150),4)
Bead=rep(sort(unique(Data2$Bead)), each=150*6)
IDensP=as.numeric(rep(seq(0.1,15,by=0.1),6*4))
Data3=data.frame(Strain,Bead,IDensP)

# Holling I model
IngesP=list(); IngesPSD=list(); IngesPLCI=list(); IngesPUCI=list()
for (i in 1:length(ModelHI)) {IngesP[[i]]=predict(ModelHI[[i]], newdata=data.frame(IDens=unique(IDensP)))}
for (i in 1:length(ModelHI)) {IngesPSD[[i]]=predictNLS(ModelHI[[i]], newdata=data.frame(IDens=unique(IDensP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,3]}
for (i in 1:length(ModelHI)) {IngesPLCI[[i]]=predictNLS(ModelHI[[i]], newdata=data.frame(IDens=unique(IDensP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,5]}
for (i in 1:length(ModelHI)) {IngesPUCI[[i]]=predictNLS(ModelHI[[i]], newdata=data.frame(IDens=unique(IDensP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,6]}
DataHI=data.frame(Data3[c(1:3600),], IngesP=unlist(IngesP), IngesPLSD=unlist(IngesP)-unlist(IngesPSD), IngesPUSD=unlist(IngesP)+unlist(IngesPSD), IngesPLCI=unlist(IngesPLCI), IngesPUCI=unlist(IngesPUCI), Model="Holling I")

# Holling II model
IngesP=list(); IngesPSD=list(); IngesPLCI=list(); IngesPUCI=list()
for (i in 1:length(ModelHII)) {IngesP[[i]]=predict(ModelHII[[i]], newdata=data.frame(IDens=unique(IDensP)))}
for (i in 1:length(ModelHII)) {IngesPSD[[i]]=predictNLS(ModelHII[[i]], newdata=data.frame(IDens=unique(IDensP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,3]}
for (i in 1:length(ModelHII)) {IngesPLCI[[i]]=predictNLS(ModelHII[[i]], newdata=data.frame(IDens=unique(IDensP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,5]}
for (i in 1:length(ModelHII)) {IngesPUCI[[i]]=predictNLS(ModelHII[[i]], newdata=data.frame(IDens=unique(IDensP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,6]}
DataHII=data.frame(Data3[c(1:3600),], IngesP=unlist(IngesP), IngesPLSD=unlist(IngesP)-unlist(IngesPSD), IngesPUSD=unlist(IngesP)+unlist(IngesPSD), IngesPLCI=unlist(IngesPLCI), IngesPUCI=unlist(IngesPUCI), Model="Holling II")

# Holling III model
IngesP=list(); IngesPSD=list(); IngesPLCI=list(); IngesPUCI=list()
for (i in 1:length(ModelHIII)) {IngesP[[i]]=predict(ModelHIII[[i]], newdata=data.frame(IDens=unique(IDensP)))}
for (i in 1:length(ModelHIII)) {IngesPSD[[i]]=predictNLS(ModelHIII[[i]], newdata=data.frame(IDens=unique(IDensP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,3]}
for (i in 1:length(ModelHIII)) {IngesPLCI[[i]]=predictNLS(ModelHIII[[i]], newdata=data.frame(IDens=unique(IDensP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,5]}
for (i in 1:length(ModelHIII)) {IngesPUCI[[i]]=predictNLS(ModelHIII[[i]], newdata=data.frame(IDens=unique(IDensP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,6]}
DataHIII=data.frame(Data3[-c(1:150,901:1050,1201:1350,1651:1800,1951:2100,2401:3600),], IngesP=unlist(IngesP), IngesPLSD=unlist(IngesP)-unlist(IngesPSD), IngesPUSD=unlist(IngesP)+unlist(IngesPSD), IngesPLCI=unlist(IngesPLCI), IngesPUCI=unlist(IngesPUCI), Model="Holling III")


############################################
### Comparing functional response models ###
############################################

# Extract combinations of names
Strain=rep(unique(Data2[,c("Strain")]),4)
Bead=rep(unique(Data2[,c("Bead")]), each=6)
Strain1=Strain[-c(1,7,9,12,14,17:24)]
Bead1=Bead[-c(1,7,9,12,14,17:24)]

# Summaries of models
SummaHI$Strain=Strain; SummaHI$Bead=Bead; SummaHI$Model="Holling I"
SummaHII$Strain=Strain; SummaHII$Bead=Bead; SummaHII$Model="Holling II"
SummaHIII$Strain=Strain1; SummaHIII$Bead=Bead1; SummaHIII$Model="Holling III"

# Combine AIC coefficients
DataAIC=rbind(SummaHI,SummaHII,SummaHIII)
DataAIC=DataAIC[order(DataAIC$Strain),]

# Calculate likelihood ratios
LR1=list(); LR2=list()
for (i in 1:length(Strain)) {
  for (j in 1:length(Strain1)) {
    LR1[[i]]=lrtest(ModelHII[[i]],ModelHI[[i]])
    LR2[[j]]=lrtest(ModelHIII[[j]],ModelHI[-c(1,7,9,12,14,17:24)][[j]])
  }
}


###############################################
### Plotting predicted functional responses ###
###############################################

# Create a dataset
Data4=rbind(DataHI[,c(1:9)],DataHII[,c(1:9)],DataHIII[,c(1:9)])
Data4[,c(4:8)]=round(Data4[,c(4:8)],4)
Data4[,c(4:8)][Data4[,c(4:8)]<0]=0
Data4=Data4[order(Data4$Strain,Data4$Bead),]

# Include model selection
Data4$AIC=rep(DataAIC[,2], each=150)
Data5=as.data.frame(Data4 %>% group_by(Strain,Bead) %>% slice(which.min(AIC)))
Data4$Selection=ifelse(Data4[,1] %in% Data5[,1] & Data4[,2] %in% Data5[,2] & Data4[,10] %in% Data5[,10], "Accepted", "Rejected")

tiff('Functional Responses Beads Models.tiff', units="in", width=15, height=20, res=1000)
ggplot(Data4, aes(IDensP, IngesP)) +
  geom_line(data=subset(Data4, Selection=="Rejected"), aes(linetype=Model, size=Model), color="grey70") +
  geom_line(data=subset(Data4, Selection=="Accepted"), aes(color=Strain, linetype=Model, size=Model)) +
  geom_point(data=Data, aes(IDens/10^5, Inges, color=Strain), size=2, pch=16) +
  ylab(expression(italic('B. calyciflorus')~'ingestion rate'~'('*cells~sec^-1~rotifer^-1*')')) +
  xlab(expression(italic('C. reinhardtii')~'density'~'('*10^5~cells~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,0.4,by=0.1), limits=c(0,0.4)) +
  scale_x_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,15,by=3), limits=c(0,15)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
  scale_linetype_manual(values=c("Holling I"="dotted","Holling II"="solid","Holling III"="dotted","Ivlev II"="dashed")) +
  scale_size_manual(values=c("Holling I"=1,"Holling II"=1,"Holling III"=1.4,"Ivlev II"=1)) +
  theme(strip.text.x=element_blank()) +
  facet_wrap(Strain~Bead, scales="free", ncol=4, nrow=6) +
  theme(legend.position="none")
dev.off()
