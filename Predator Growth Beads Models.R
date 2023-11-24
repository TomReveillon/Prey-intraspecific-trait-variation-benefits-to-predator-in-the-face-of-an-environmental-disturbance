setwd("~/LIMNO 2019-2022/Experiments/Predator Growth Beads")

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

####################################################################################
####################################################################################
##### PREDATOR GROWTH RATE FOR PREDATOR-PREY SYSTEM WITH DIFFERENT CLONE TYPES #####
####################################################################################
####################################################################################

# Import the dataset
Data=read.table("Data_RGB.txt", h=T, dec=",")
names(Data)
summary(Data)

# Specify the variables as numeric or factor
Data[,c(3,5)] %<>% mutate_if(is.character,as.numeric)

# Preserve the order of factors
Data$Bead=factor(Data$Bead, levels=unique(Data$Bead))

# Sort the dataset by strain
Data=Data[order(Data$Strain),]

# Calculate mean densities of trials
Data2=setDT(Data)[, .(MeanDens=mean(Count), MeanEgg=mean(Eggs)), by=list(Strain,Bead,Day)]
Data2=as.data.frame(Data2)

# Calculate mean egg number per individual
Data2$AveEgg=Data2$MeanEgg/Data2$MeanDens


##########################
### Growth rate models ###
##########################

# Linear model
SetLi=function(mCall,LHS,data) {
  xy=sortedXyData(mCall[["x"]],LHS,data)
  b=min(xy$y); a=coef(lm(y~x,xy))[2] 
  Out=c(b,a); names(Out)=mCall[c("b","a")]
  return(Out)}
SSline=selfStart(as.formula("~b + a*x"), initial=SetLi, parameters=c("b","a"))

FuncLi=function(Model, Day) {
  Params=coef(Model)
  names(Params)=NULL
  b=Params[1]
  a=Params[2]
  
  Fit=Model$m$fitted()
  Res=Model$m$resid()
  MSS=sum((Fit-mean(Fit))^2)
  RSS=sum(Res^2)
  R2=MSS/(MSS+RSS)
  AIC=AIC(Model)
  
  Summary=data.frame(R2,AIC)
  Parameters=data.frame(b=b,a=a)
  Rates=data.frame(
    Time = Day,
    DensP = b + a * Day,
    AGR = a)
  
  Rates$RGR = Rates$AGR / Rates$DensP
  
  Out=list(Model=Model, Summary=Summary, Parameters=Parameters, Rates=Rates)
  return(Out)
}

# Exponential model
SetEx=function(mCall,LHS,data) {
  xy=sortedXyData(mCall[["x"]],LHS,data)
  b=min(xy$y); r=coef(lm(log(y)~x,xy))[2]
  Out=c(b,r); names(Out)=mCall[c("b","r")]
  return(Out)}
SSexpo=selfStart(as.formula("~b * exp(r * x)"), initial=SetEx, parameters=c("b","r"))

FuncEx=function(Model,Day) {
  Params=coef(Model)
  names(Params)=NULL
  b=Params[1]
  r=Params[2]
  
  Fit=Model$m$fitted()
  Res=Model$m$resid()
  MSS=sum((Fit-mean(Fit))^2)
  RSS=sum(Res^2)
  R2=MSS/(MSS+RSS)
  AIC=AIC(Model)
  
  Summary=data.frame(R2,AIC)
  Parameters=data.frame(b=b,r=r)
  Rates=data.frame(
    Time = Day,
    DensP = b * exp(r * Day),
    AGR = r * b * exp(r * Day))
  
  Rates$RGR = r
  
  Out=list(Model=Model, Summary=Summary, Parameters=Parameters, Rates=Rates)
  return(Out)
}

# Logistic model
SetLo=function(Coef){
  K=Coef[1]; r=1/(Coef[3]); b=K/(1 + exp(Coef[2]/Coef[3]))
  Out=c(b,K,r); names(Out)=c("b","K","r")
  return(Out)}

FuncLo=function(Model,Day) {
  Coef=coef(Model)
  Params=SetLo(Coef)
  b=Params[1] 
  K=Params[2] 
  r=Params[3]
  
  Fit=Model$m$fitted()
  Res=Model$m$resid()
  MSS=sum((Fit-mean(Fit))^2)
  RSS=sum(Res^2)
  R2=MSS/(MSS+RSS)
  AIC=AIC(Model)
  
  Summary=data.frame(R2,AIC)
  Parameters=data.frame(b=b,K=K,r=r)
  Rates=data.frame(
    Time = Day,
    DensP = K * (b * K) / (b + (K - b) * exp(-r * Day)),
    AGR = (r * b * K * (K - b) * exp(-r * Day)) / (b + (K - b)*exp(-r * Day))^2)

  Rates$RGR = Rates$AGR / Rates$DensP
  
  Out=list(Model=Model, Summary=Summary, Parameters=Parameters, Rates=Rates)
  return(Out)
}

# Extinction model
SetEc=function(mCall,LHS,data) {
  xy=sortedXyData(mCall[["x"]],LHS,data)
  b=5.0; r=3.0; K=5.0
  Out=c(b,r,K); names(Out)=mCall[c("b","r","K")]
  return(Out)}
SSexct=selfStart(as.formula("~b/(1 + exp(r * (log(x) - log(K))))"), initial=SetEc, parameters=c("b","r","K"))

FuncEc=function(Model,Day) {
  Params=coef(Model)
  b=Params[1]
  r=Params[2]
  K=Params[3]
  
  Fit=Model$m$fitted()
  Res=Model$m$resid()
  MSS=sum((Fit-mean(Fit))^2)
  RSS=sum(Res^2)
  R2=MSS/(MSS+RSS)
  AIC=AIC(Model)
  
  Summary=data.frame(R2,AIC)
  Parameters=data.frame(b=b,r=r,K=K)
  Rates=data.frame(
    Time = Day,
    DensP = b/(1 + exp(r * (log(Day) - log(K)))))
  
  Out=list(Model=Model, Summary=Summary, Parameters=Parameters, Rates=Rates)
  return(Out)
}


##################################
### Fitting growth rate models ###
##################################

# Split the dataset
SplitData2=split(Data2, list(Data2$Strain,Data2$Bead))

# Extract combinations of names
Names=unique(Data2[,c("Strain","Bead","Day")])
Names=Names[order(Names$Bead,Names$Strain),]
Strain=Names$Strain; Bead=Names$Bead; Time=Names$Day

# Linear model
ModLi=function(x) {
  FitLi=nls(MeanDens ~ SSline(Day, a, b), data=x)
  OutLi=FuncLi(FitLi, x$Day)} 
OutLi=lapply(SplitData2, ModLi) 

RateLi=bind_rows(lapply(OutLi, function (x) x[c("Rates")]))
RateLi=round(as.data.frame(do.call("rbind",RateLi)),4)
RateLi=cbind(Strain=Strain[c(1:144)],Bead=Bead[c(1:144)],Model="Linear",RateLi)
rownames(RateLi)=c()
ParamLi=bind_rows(lapply(OutLi, function (x) x[c("Parameters")]))
ParamLi=round(cbind(do.call("rbind",ParamLi)),4)
rownames(ParamLi)=c()
SummaLi=bind_rows(lapply(OutLi, function (x) x[c("Summary")]))
SummaLi=round(cbind(do.call("rbind",SummaLi)),4)
rownames(SummaLi)=c()
ModelLi=unlist(lapply(OutLi, function (x) x[c("Model")]),recursive=F)

# Exponential model
ModEx=function(x) {
  FitEx=nls(MeanDens ~ SSexpo(Day, r, b), data=x)
  OutEx=FuncEx(FitEx, x$Day)}
OutEx=lapply(SplitData2, ModEx)

RateEx=bind_rows(lapply(OutEx, function (x) x[c("Rates")]))
RateEx=round(as.data.frame(do.call("rbind",RateEx)),4)
RateEx=cbind(Strain=Strain[c(1:144)],Bead=Bead[c(1:144)],Model="Exponential",RateEx)
rownames(RateEx)=c()
ParamEx=bind_rows(lapply(OutEx, function (x) x[c("Parameters")]))
ParamEx=round(as.data.frame(do.call("rbind",ParamEx)),4)
rownames(ParamEx)=c()
SummaEx=bind_rows(lapply(OutEx, function (x) x[c("Summary")]))
SummaEx=round(as.data.frame(do.call("rbind",SummaEx)),4)
rownames(SummaEx)=c()
ModelEx=unlist(lapply(OutEx, function (x) x[c("Model")]),recursive=F)

# Logistic model
ModLo=function(x) {
  FitLo=nls(MeanDens ~ SSlogis(Day, b, K, r), data=x)
  OutLo=FuncLo(FitLo, x$Day)}
OutLo=lapply(SplitData2[c(1:7,11,13:14,19,21:22,24)], ModLo)

RateLo=bind_rows(lapply(OutLo, function (x) x[c("Rates")]))
RateLo=round(as.data.frame(do.call("rbind",RateLo)),4)
RateLo=cbind(Strain=Strain[c(1:42,61:66,73:84,109:114,121:132,139:144)],Bead=Bead[c(1:42,61:66,73:84,109:114,121:132,139:144)],Model="Logistic",RateLo)
rownames(RateLo)=c()
ParamLo=bind_rows(lapply(OutLo, function (x) x[c("Parameters")]))
ParamLo=round(cbind(do.call("rbind",ParamLo)),4)
rownames(ParamLo)=c()
SummaLo=bind_rows(lapply(OutLo, function (x) x[c("Summary")]))
SummaLo=round(cbind(do.call("rbind",SummaLo)),4)
rownames(SummaLo)=c()
ModelLo=unlist(lapply(OutLo, function (x) x[c("Model")]),recursive=F)

# Extinction model
ModEc=function(x) {
  FitEc=nls(MeanDens ~ b/(1 + exp(r * (log(Day) - log(K)))), start=c(b=5.0, r=3.0, K=5.0), data=x)
  OutEc=FuncEc(FitEc, x$Day)}
OutEc=lapply(SplitData2[c(19:24)], ModEc)

RateEc=bind_rows(lapply(OutEc, function (x) x[c("Rates")]))
RateEc=round(as.data.frame(do.call("rbind",RateEc)),4)
RateEc=cbind(Strain=Strain[c(109:144)],Bead=Bead[c(109:144)],Model="Extinction",RateEc)
rownames(RateEc)=c()
ParamEc=bind_rows(lapply(OutEc, function (x) x[c("Parameters")]))
ParamEc=round(cbind(do.call("rbind",ParamEc)),4)
rownames(ParamEc)=c()
SummaEc=bind_rows(lapply(OutEc, function (x) x[c("Summary")]))
SummaEc=round(cbind(do.call("rbind",SummaEc)),4)
rownames(SummaEc)=c()
ModelEc=unlist(lapply(OutEc, function (x) x[c("Model")]),recursive=F)

# Set the predicted dataset
Strain=rep(rep(unique(Data2$Strain), each=51),4)
Bead=rep(sort(unique(Data2$Bead)), each=51*6)
DayP=as.numeric(rep(seq(0,5,by=0.1),6*4))
Data3=data.frame(Strain,Bead,DayP)

# Linear model
DensP=list(); DensPSD=list(); DensPLCI=list(); DensPUCI=list()
for (i in 1:length(ModelLi)) {DensP[[i]]=predict(ModelLi[[i]], newdata=data.frame(Day=unique(DayP)))}
for (i in 1:length(ModelLi)) {DensPSD[[i]]=predictNLS(ModelLi[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,3]}
for (i in 1:length(ModelLi)) {DensPLCI[[i]]=predictNLS(ModelLi[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,5]}
for (i in 1:length(ModelLi)) {DensPUCI[[i]]=predictNLS(ModelLi[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,6]}
DataLi=data.frame(Data3[c(1:1224),], DensP=unlist(DensP), DensPLSD=unlist(DensP)-unlist(DensPSD), DensPUSD=unlist(DensP)+unlist(DensPSD), DensPLCI=unlist(DensPLCI), DensPUCI=unlist(DensPUCI), Model="Linear")

# Exponential model
DensP=list(); DensPSD=list(); DensPLCI=list(); DensPUCI=list()
for (i in 1:length(ModelEx)) {DensP[[i]]=predict(ModelEx[[i]], newdata=data.frame(Day=unique(DayP)))}
for (i in 1:length(ModelEx)) {DensPSD[[i]]=predictNLS(ModelEx[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,3]}
for (i in 1:length(ModelEx)) {DensPLCI[[i]]=predictNLS(ModelEx[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,5]}
for (i in 1:length(ModelEx)) {DensPUCI[[i]]=predictNLS(ModelEx[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,6]}
DataEx=data.frame(Data3[c(1:1224),], DensP=unlist(DensP), DensPLSD=unlist(DensP)-unlist(DensPSD), DensPUSD=unlist(DensP)+unlist(DensPSD), DensPLCI=unlist(DensPLCI), DensPUCI=unlist(DensPUCI), Model="Exponential")

# Logistic model
DensP=list(); DensPSD=list(); DensPLCI=list(); DensPUCI=list()
for (i in 1:length(ModelLo)) {DensP[[i]]=predict(ModelLo[[i]], newdata=data.frame(Day=unique(DayP)))}
for (i in 1:length(ModelLo)) {DensPSD[[i]]=predictNLS(ModelLo[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,3]}
for (i in 1:length(ModelLo)) {DensPLCI[[i]]=predictNLS(ModelLo[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,5]}
for (i in 1:length(ModelLo)) {DensPUCI[[i]]=predictNLS(ModelLo[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,6]}
DataLo=data.frame(Data3[c(1:357,511:561,613:714,919:969,1021:1122,1174:1224),], DensP=unlist(DensP), DensPLSD=unlist(DensP)-unlist(DensPSD), DensPUSD=unlist(DensP)+unlist(DensPSD), DensPLCI=unlist(DensPLCI), DensPUCI=unlist(DensPUCI), Model="Logistic")

# Extinction model
DensP=list(); DensPSD=list(); DensPLCI=list(); DensPUCI=list()
for (i in 1:length(ModelEc)) {DensP[[i]]=predict(ModelEc[[i]], newdata=data.frame(Day=unique(DayP)))}
for (i in 1:length(ModelEc)) {DensPSD[[i]]=predictNLS(ModelEc[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,3]}
for (i in 1:length(ModelEc)) {DensPLCI[[i]]=predictNLS(ModelEc[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,5]}
for (i in 1:length(ModelEc)) {DensPUCI[[i]]=predictNLS(ModelEc[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,6]}
DataEc=data.frame(Data3[c(919:1224),], DensP=unlist(DensP), DensPLSD=unlist(DensP)-unlist(DensPSD), DensPUSD=unlist(DensP)+unlist(DensPSD), DensPLCI=unlist(DensPLCI), DensPUCI=unlist(DensPUCI), Model="Extinction")


############################################
### Comparing functional response models ###
############################################

# Extract combinations of names
Strain=rep(unique(Data2[,c("Strain")]),4)
Bead=rep(unique(Data2[,c("Bead")]), each=6)
Strain1=Strain[-c(8:10,12,15:18,20,23)]
Bead1=Bead[-c(8:10,12,15:18,20,23)]
Strain2=Strain[-c(1:18)]
Bead2=Bead[-c(1:18)]

# Summaries of models
SummaLi$Strain=Strain; SummaLi$Bead=Bead; SummaLi$Model="Linear"
SummaEx$Strain=Strain; SummaEx$Bead=Bead; SummaEx$Model="Exponential"
SummaLo$Strain=Strain1; SummaLo$Bead=Bead1; SummaLo$Model="Logistic"
SummaEc$Strain=Strain2; SummaEc$Bead=Bead2; SummaEc$Model="Extinction"

# Combine AIC coefficients
DataAIC=rbind(SummaLi,SummaEx,SummaLo,SummaEc)
DataAIC=DataAIC[order(DataAIC$Strain),]

# Calculate likelihood ratios
LR1=list(); LR2=list(); LR3=list()
for (i in 1:length(Strain)) {
  for (j in 1:length(Strain1)) {
    for (k in 1:length(Strain2)) {
      LR1[[i]]=lrtest(ModelEx[[i]],ModelLi[[i]])
      LR2[[j]]=lrtest(ModelLo[[j]],ModelLi[-c(8:10,12,15:18,20,23)][[j]])
      LR3[[k]]=lrtest(ModelEc[[k]],ModelLi[-c(1:18)][[k]])
    }
  }
}


########################################
### Plotting predicted growth curves ###
########################################

# Create a dataset
Data4=rbind(DataLi[,c(1:9)],DataEx[,c(1:9)],DataLo[,c(1:9)],DataEc[,c(1:9)])
Data4[,c(4:8)]=round(Data4[,c(4:8)],4)
Data4[,c(4:8)][is.na(Data4[,c(4:8)])]=0
Data4=Data4[order(Data4$Strain,Data4$Bead),]

# Include model selection
Data4$AIC=rep(DataAIC[,2], each=51)
Data5=as.data.frame(Data4 %>% group_by(Strain,Bead) %>% slice(which.min(AIC)))
Data4$Selection=ifelse(Data4[,1] %in% Data5[,1] & Data4[,2] %in% Data5[,2] & Data4[,10] %in% Data5[,10], "Accepted", "Rejected")

tiff('Growth Curves Beads Beads Models.tiff', units="in", width=15, height=20, res=1000)
ggplot(Data4, aes(DayP, DensP)) +
  geom_line(data=subset(Data4, Selection=="Rejected"), aes(linetype=Model, size=Model), color="grey70") +
  geom_line(data=subset(Data4, Selection=="Accepted"), aes(color=Strain, linetype=Model, size=Model)) +
  geom_point(data=Data, aes(Day, Count, color=Strain), size=2, pch=16) +
  ylab(expression(italic('B. calyciflorus')~'density'~'('*rotifers~mL^-1*')')) + xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,25,by=5), limits=c(0,25)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,5,by=1), limits=c(0,5)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
  scale_linetype_manual(values=c("Linear"="dotted","Exponential"="dotted","Logistic"="solid","Gompertz"="dashed","Mortality"="longdash")) +
  scale_size_manual(values=c("Linear"=1,"Exponential"=1.4,"Logistic"=1,"Gompertz"=1,"Mortality"=1)) +
  theme(strip.text.x=element_blank()) +
  facet_wrap(Strain~Bead, scales="free", ncol=4, nrow=6) +
  theme(legend.position="none")
dev.off()
