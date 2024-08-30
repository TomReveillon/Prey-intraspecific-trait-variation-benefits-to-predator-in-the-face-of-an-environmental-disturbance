setwd("~/LIMNO 2019-2023/Experiments/Predator Growth Beads")

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


#########################################
### Plots of replicated growth curves ###
#########################################

ggplot(subset(Data, Strain=="CR1"), aes(Day, Count, group=interaction(Strain,Bead))) + 
  geom_smooth(method="loess", color="mediumpurple3", size=1, se=F) +
  ylab(expression(italic('B. calyciflorus')~'density'~'('*rotifers~mL^-1*')')) +
  xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,20,by=5), limits=c(0,20)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,5,by=1), limits=c(0,5)) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.position="none")

ggplot(subset(Data, Strain=="CR2"), aes(Day, Count, group=interaction(Strain,Bead))) +
  geom_smooth(method="loess", color="cornflowerblue", size=1, se=F) +
  ylab(expression(italic('B. calyciflorus')~'density'~'('*rotifers~mL^-1*')')) +
  xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,20,by=5), limits=c(0,20)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,5,by=1), limits=c(0,5)) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.position="none")

ggplot(subset(Data, Strain=="CR3"), aes(Day, Count, group=interaction(Strain,Bead))) +   
  geom_smooth(method="loess", color="chartreuse3", size=1, se=F) +
  ylab(expression(italic('B. calyciflorus')~'density'~'('*rotifers~mL^-1*')')) +
  xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,20,by=5), limits=c(0,20)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,5,by=1), limits=c(0,5)) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.position="none")

ggplot(subset(Data, Strain=="CR4"), aes(Day, Count, group=interaction(Strain,Bead))) +   
  geom_smooth(method="loess", color="gold2", size=1, se=F) +
  ylab(expression(italic('B. calyciflorus')~'density'~'('*rotifers~mL^-1*')')) +
  xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,20,by=5), limits=c(0,20)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,5,by=1), limits=c(0,5)) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.position="none")

ggplot(subset(Data, Strain=="CR5"), aes(Day, Count, group=interaction(Strain,Bead))) +   
  geom_smooth(method="loess", color="darkorange1", size=1, se=F) +
  ylab(expression(italic('B. calyciflorus')~'density'~'('*rotifers~mL^-1*')')) +
  xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,20,by=5), limits=c(0,20)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,5,by=1), limits=c(0,5)) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.position="none")

ggplot(subset(Data, Strain=="CR6"), aes(Day, Count, group=interaction(Strain,Bead))) +   
  geom_smooth(method="loess", color="tomato2", size=1, se=F) +
  ylab(expression(italic('B. calyciflorus')~'density'~'('*rotifers~mL^-1*')')) +
  xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,20,by=5), limits=c(0,20)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,5,by=1), limits=c(0,5)) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.position="none")


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
OutLi=lapply(SplitData2[c(7:9,13:18)], ModLi) 

RateLi=bind_rows(lapply(OutLi, function (x) x[c("Rates")]))
RateLi=round(as.data.frame(do.call("rbind",RateLi)),4)
RateLi=cbind(Strain=Strain[c(37:54,73:108)],Bead=Bead[c(37:54,73:108)],Model="Linear",RateLi)
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
OutEx=lapply(SplitData2[c(10:12)], ModEx)

RateEx=bind_rows(lapply(OutEx, function (x) x[c("Rates")]))
RateEx=round(as.data.frame(do.call("rbind",RateEx)),4)
RateEx=cbind(Strain=Strain[c(55:72)],Bead=Bead[c(55:72)],Model="Exponential",RateEx)
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
OutLo=lapply(SplitData2[c(1:6)], ModLo)

RateLo=bind_rows(lapply(OutLo, function (x) x[c("Rates")]))
RateLo=round(as.data.frame(do.call("rbind",RateLo)),4)
RateLo=cbind(Strain=Strain[c(1:30,31:36)],Bead=Bead[c(1:30,31:36)],Model="Logistic",RateLo)
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
DataLi=data.frame(Data3[c(307:459,613:918),], DensP=unlist(DensP), DensPLSD=unlist(DensP)-unlist(DensPSD), DensPUSD=unlist(DensP)+unlist(DensPSD), DensPLCI=unlist(DensPLCI), DensPUCI=unlist(DensPUCI), Model="Linear")

# Exponential model
DensP=list(); DensPSD=list(); DensPLCI=list(); DensPUCI=list()
for (i in 1:length(ModelEx)) {DensP[[i]]=predict(ModelEx[[i]], newdata=data.frame(Day=unique(DayP)))}
for (i in 1:length(ModelEx)) {DensPSD[[i]]=predictNLS(ModelEx[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,3]}
for (i in 1:length(ModelEx)) {DensPLCI[[i]]=predictNLS(ModelEx[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,5]}
for (i in 1:length(ModelEx)) {DensPUCI[[i]]=predictNLS(ModelEx[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,6]}
DataEx=data.frame(Data3[c(460:612),], DensP=unlist(DensP), DensPLSD=unlist(DensP)-unlist(DensPSD), DensPUSD=unlist(DensP)+unlist(DensPSD), DensPLCI=unlist(DensPLCI), DensPUCI=unlist(DensPUCI), Model="Exponential")

# Logistic model
DensP=list(); DensPSD=list(); DensPLCI=list(); DensPUCI=list()
for (i in 1:length(ModelLo)) {DensP[[i]]=predict(ModelLo[[i]], newdata=data.frame(Day=unique(DayP)))}
for (i in 1:length(ModelLo)) {DensPSD[[i]]=predictNLS(ModelLo[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,3]}
for (i in 1:length(ModelLo)) {DensPLCI[[i]]=predictNLS(ModelLo[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,5]}
for (i in 1:length(ModelLo)) {DensPUCI[[i]]=predictNLS(ModelLo[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,6]}
DataLo=data.frame(Data3[c(1:306),], DensP=unlist(DensP), DensPLSD=unlist(DensP)-unlist(DensPSD), DensPUSD=unlist(DensP)+unlist(DensPSD), DensPLCI=unlist(DensPLCI), DensPUCI=unlist(DensPUCI), Model="Logistic")

# Extinction model
DensP=list(); DensPSD=list(); DensPLCI=list(); DensPUCI=list()
for (i in 1:length(ModelEc)) {DensP[[i]]=predict(ModelEc[[i]], newdata=data.frame(Day=unique(DayP)))}
for (i in 1:length(ModelEc)) {DensPSD[[i]]=predictNLS(ModelEc[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,3]}
for (i in 1:length(ModelEc)) {DensPLCI[[i]]=predictNLS(ModelEc[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,5]}
for (i in 1:length(ModelEc)) {DensPUCI[[i]]=predictNLS(ModelEc[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,6]}
DataEc=data.frame(Data3[c(919:1224),], DensP=unlist(DensP), DensPLSD=unlist(DensP)-unlist(DensPSD), DensPUSD=unlist(DensP)+unlist(DensPSD), DensPLCI=unlist(DensPLCI), DensPUCI=unlist(DensPUCI), Model="Extinction")


########################################
### Plotting predicted growth curves ###
########################################

# Create a dataset
Data4=rbind(DataLi[,c(1:9)],DataEx[,c(1:9)],DataLo[,c(1:9)],DataEc[,c(1:9)])
Data4[,c(4:8)]=round(Data4[,c(4:8)],4)
Data4[,c(4:8)][is.na(Data4[,c(4:8)])]=0
Data4=Data4[order(Data4$Strain,Data4$Bead),]

# Correct standard errors and confidence intervals
Data4$DensPLSD=ifelse(Data4$DensPLSD==0, Data4$DensP-Data4$DensP*((lead(Data4$DensP)-lead(Data4$DensPLSD))/lead(Data4$DensP)), Data4$DensPLSD)
Data4$DensPUSD=ifelse(Data4$DensPUSD==0, Data4$DensP+Data4$DensP*((lead(Data4$DensPUSD)-lead(Data4$DensP))/lead(Data4$DensP)), Data4$DensPUSD)
Data4$DensPLCI=ifelse(Data4$DensPLCI==0, Data4$DensP-Data4$DensP*((lead(Data4$DensP)-lead(Data4$DensPLCI))/lead(Data4$DensP)), Data4$DensPLCI)
Data4$DensPUCI=ifelse(Data4$DensPUCI==0, Data4$DensP+Data4$DensP*((lead(Data4$DensPUCI)-lead(Data4$DensP))/lead(Data4$DensP)), Data4$DensPUCI)

# Correct standard errors and confidence intervals
Data4$DensPLSD=ifelse(Data4$DensPLSD < Data4$DensP*0.80, Data4$DensP*0.80, Data4$DensPLSD)
Data4$DensPUSD=ifelse(Data4$DensPUSD > Data4$DensP*1.20, Data4$DensP*1.20, Data4$DensPUSD)
Data4$DensPLCI=ifelse(Data4$DensPLCI < Data4$DensP*0.80, Data4$DensP*0.80, Data4$DensPLCI)
Data4$DensPUCI=ifelse(Data4$DensPUCI > Data4$DensP*1.20, Data4$DensP*1.20, Data4$DensPUCI)

# Export the dataset
Data4[,c(5:8)]=replace(Data4[,c(5:8)],Data4[,c(5:8)]<0,0)
write.table(Data4, file="Data_RGBP.txt", sep="\t", row.names=F)

# Create a dataset for labels
DataL=as.data.frame(Data4%>% group_by(Strain,Bead) %>% filter(DayP == max(DayP)))

tiff('Growth Curves Beads.tiff', units="in", width=15, height=8, res=1000)
ggplot(Data4, aes(DayP, DensP, group=interaction(Strain,Bead))) +
  geom_line(aes(color=Strain, linetype=Bead), size=1) +
  geom_point(data=Data, aes(Day, Count, color=Strain), size=2, pch=16, position=position_jitter(h=0, w=0.2)) +
  geom_text_repel(DataL, mapping=aes(label=Bead, color=Strain), size=5, nudge_y=0, nudge_x=0.5, direction="y", segment.color="transparent") +
  ylab(expression(italic('B. calyciflorus')~'density'~'('*rotifers~mL^-1*')')) + xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,24,by=6), limits=c(0,24)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,5,by=1), limits=c(0,5.5)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
  scale_linetype_manual(values=c("C"="solid","L"="solid","M"="solid","H"="solid")) +
  theme(strip.text.x=element_blank()) +
  facet_wrap(~Strain, scales="free", ncol=3, nrow=2) +
  theme(legend.position="none")
dev.off()

tiff('Growth Curves Beads Intervals.tiff', units="in", width=15, height=8, res=1000)
ggplot(Data4, aes(DayP, DensP, group=interaction(Strain,Bead))) +
  geom_line(aes(color=Strain, linetype=Bead), size=1) +
  geom_ribbon(data=subset(Data4, Bead=="C"), aes(ymin=DensPLSD, ymax=DensPUSD, fill=Strain, color=Strain), linetype="solid", size=0.5, alpha=0.3) +
  geom_point(data=Data, aes(Day, Count, color=Strain), size=2, pch=16, position=position_jitter(h=0, w=0.2)) +
  geom_text_repel(DataL, mapping=aes(label=Bead, color=Strain), size=5, nudge_y=0, nudge_x=0.5, direction="y", segment.color="transparent") +
  ylab(expression(italic('B. calyciflorus')~'density'~'('*rotifers~mL^-1*')')) + xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,24,by=6), limits=c(0,24)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,5,by=1), limits=c(0,5.5)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_fill_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
  scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
  scale_linetype_manual(values=c("C"="solid","L"="solid","M"="solid","H"="solid")) +
  theme(strip.text.x=element_blank()) +
  facet_wrap(~Strain, scales="free", ncol=3, nrow=2) +
  theme(legend.position="none")
dev.off()


########################################
### Plotting predicted growth curves ###
########################################

# Create a dataset
Data4=rbind(DataLi[,c(1:9)],DataEx[,c(1:9)],DataLo[,c(1:9)],DataEc[,c(1:9)])
Data4[,c(4:8)]=round(Data4[,c(4:8)],6)
Data4[,c(4:8)][is.na(Data4[,c(4:8)])]=0
Data4=Data4[order(Data4$Strain,Data4$Bead),]

# Correct standard errors and confidence intervals
Data4$DensPLSD=ifelse(Data4$DensPLSD==0, Data4$DensP-Data4$DensP*((lead(Data4$DensP)-lead(Data4$DensPLSD))/lead(Data4$DensP)), Data4$DensPLSD)
Data4$DensPUSD=ifelse(Data4$DensPUSD==0, Data4$DensP+Data4$DensP*((lead(Data4$DensPUSD)-lead(Data4$DensP))/lead(Data4$DensP)), Data4$DensPUSD)
Data4$DensPLCI=ifelse(Data4$DensPLCI==0, Data4$DensP-Data4$DensP*((lead(Data4$DensP)-lead(Data4$DensPLCI))/lead(Data4$DensP)), Data4$DensPLCI)
Data4$DensPUCI=ifelse(Data4$DensPUCI==0, Data4$DensP+Data4$DensP*((lead(Data4$DensPUCI)-lead(Data4$DensP))/lead(Data4$DensP)), Data4$DensPUCI)

# Correct standard errors and confidence intervals
Data4$DensPLSD=ifelse(Data4$DensPLSD < Data4$DensP*0.80, Data4$DensP*0.80, Data4$DensPLSD)
Data4$DensPUSD=ifelse(Data4$DensPUSD > Data4$DensP*1.20, Data4$DensP*1.20, Data4$DensPUSD)
Data4$DensPLCI=ifelse(Data4$DensPLCI < Data4$DensP*0.80, Data4$DensP*0.80, Data4$DensPLCI)
Data4$DensPUCI=ifelse(Data4$DensPUCI > Data4$DensP*1.20, Data4$DensP*1.20, Data4$DensPUCI)

# Arrange the datasets
Data=Data %>% arrange(factor(Bead, levels=unique(Data$Bead)))
Data4=Data4 %>% arrange(factor(Bead, levels=unique(Data4$Bead)))

# Include experimental points
Data4$DayE=c(Reduce(rbind,append(list(matrix(Data$Day,30)),rep(NA,21))))
Data4$DensE=c(Reduce(rbind,append(list(matrix(Data$Count,30)),rep(NA,21))))

# Include positions for labels
DataL=as.data.frame(Data4 %>% group_by(Bead,Strain) %>% filter(DayP == max(DayP)))
Data4$DayM=c(Reduce(rbind,append(list(matrix(DataL$DayP,1)),rep(NA,50))))
Data4$DensM=c(Reduce(rbind,append(list(matrix(DataL$DensP,1)),rep(NA,50))))

# Create names for labels
Data4$Labels=factor(Data4$Strain, levels=unique(Data4$Strain), labels=c(expression(C[R1]),expression(C[R2]),expression(C[R3]),expression(C[R4]),expression(C[R5]),expression(C[R6])))

# Split the dataset
SplitData4=split(Data4, list(Data4$Bead))

PlotFunc=function(x) {
  ggplot(x, aes(DayP, DensP, group=interaction(Strain,Bead))) +
    geom_line(aes(color=Strain, linetype=Bead), size=1) +
    geom_point(data=x, aes(DayE, DensE, color=Strain), size=2, pch=16, position=position_jitter(h=0, w=0.2)) +
    ylab(expression(italic('B. calyciflorus')~'density'~'('*rotifers~mL^-1*')')) + xlab(expression('Time (days)')) +
    theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
    theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
    theme(axis.title.y=element_blank()) +
    theme(axis.title.x=element_blank()) +
    scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,5,by=1), limits=c(0,5)) +
    theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
    scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
    scale_linetype_manual(values=c("C"="solid","L"="solid","M"="solid","H"="solid")) +
    theme(legend.position="none")
}

tiff('Growth Curves Beads.tiff', units="in", width=12, height=12, res=1000)
Panel=lapply(SplitData4, PlotFunc)
Panel[[1]]=Panel[[1]] + scale_y_continuous(labels=sprintf(seq(0,24,by=6), fmt="%.0f"), breaks=seq(0,24,by=6), limits=c(0,24))
Panel[[2]]=Panel[[2]] + scale_y_continuous(labels=sprintf(seq(0,20,by=5), fmt="%.0f"), breaks=seq(0,20,by=5), limits=c(0,20))
Panel[[3]]=Panel[[3]] + scale_y_continuous(labels=sprintf(seq(0,12,by=3), fmt="%.0f"), breaks=seq(0,12,by=3), limits=c(0,12))
Panel[[4]]=Panel[[4]] + scale_y_continuous(labels=sprintf(seq(0,8,by=2), fmt="%.0f"), breaks=seq(0,8,by=2), limits=c(0,8))
Panel[[1]]=Panel[[1]] + geom_text(mapping=aes(x=0, y=Inf, group=Bead), color="black", label="C", size=5, vjust=2.2, hjust=0)
Panel[[2]]=Panel[[2]] + geom_text(mapping=aes(x=0, y=Inf, group=Bead), color="black", label="L", size=5, vjust=2.2, hjust=0)
Panel[[3]]=Panel[[3]] + geom_text(mapping=aes(x=0, y=Inf, group=Bead), color="black", label="M", size=5, vjust=2.2, hjust=0)
Panel[[4]]=Panel[[4]] + geom_text(mapping=aes(x=0, y=Inf, group=Bead), color="black", label="H", size=5, vjust=2.2, hjust=0)
Panel[[1]]=Panel[[1]] + theme(plot.margin=unit(c(0.20,0.20,0.20,0.20),"cm"))
Panel[[2]]=Panel[[2]] + theme(plot.margin=unit(c(0.20,0.20,0.20,0.20),"cm"))
Panel[[3]]=Panel[[3]] + theme(plot.margin=unit(c(0.20,0.20,0.20,0.20),"cm"))
Panel[[4]]=Panel[[4]] + theme(plot.margin=unit(c(0.20,0.20,0.20,0.55),"cm"))
Yaxis=textGrob(expression(italic('B. calyciflorus')~'density'~'('*rotifers~mL^-1*')'), gp=gpar(fontface="bold", fontsize=18), rot=90)
Xaxis=textGrob(expression('Time (days)'), gp=gpar(fontface="bold", fontsize=18), rot=0)
grid.arrange(grobs=Panel, left=Yaxis, bottom=Xaxis, ncol=2, nrow=2)
dev.off()

PlotFunc=function(x) {
  ggplot(x, aes(DayP, DensP, group=interaction(Strain,Bead))) +
    geom_line(aes(color=Strain, linetype=Bead), size=1) +
    geom_ribbon(aes(ymin=DensPLSD, ymax=DensPUSD, fill=Strain, color=Strain), linetype="solid", size=0.5, alpha=0.3) +
    geom_point(data=x, aes(DayE, DensE, color=Strain), size=2, pch=16, position=position_jitter(h=0, w=0.2)) +
    geom_text(data=x, mapping=aes(x=-Inf, y=Inf, label=Bead), color="black", size=5, vjust=1.2, hjust=-1.5, parse=T) +
    ylab(expression(italic('B. calyciflorus')~'density'~'('*rotifers~mL^-1*')')) + xlab(expression('Time (days)')) +
    theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
    theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
    theme(axis.title.y=element_blank()) +
    theme(axis.title.x=element_blank()) +
    scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,5,by=1), limits=c(0,5)) +
    theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
    scale_fill_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
    scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
    scale_linetype_manual(values=c("C"="solid","L"="solid","M"="solid","H"="solid")) +
    theme(strip.text.x=element_blank()) + theme(panel.spacing.x=unit(0.5,"cm")) +
    facet_wrap(~Strain, scales="free", ncol=6, nrow=1) +
    theme(legend.position="none")
}

tiff('Growth Curves Beads Intervals.tiff', units="in", width=20, height=13, res=1000)
Panel=lapply(SplitData4, PlotFunc)
Panel[[1]]=Panel[[1]] + scale_y_continuous(labels=sprintf(seq(0,24,by=6), fmt="%.0f"), breaks=seq(0,24,by=6), limits=c(0,24))
Panel[[2]]=Panel[[2]] + scale_y_continuous(labels=sprintf(seq(0,20,by=5), fmt="%.0f"), breaks=seq(0,20,by=5), limits=c(0,20))
Panel[[3]]=Panel[[3]] + scale_y_continuous(labels=sprintf(seq(0,12,by=3), fmt="%.0f"), breaks=seq(0,12,by=3), limits=c(0,12))
Panel[[4]]=Panel[[4]] + scale_y_continuous(labels=sprintf(seq(0,8,by=2), fmt="%.0f"), breaks=seq(0,8,by=2), limits=c(0,8))
Panel[[1]]=Panel[[1]] + theme(plot.margin=unit(c(0.20,0.50,0.20,0.20),"cm"))
Panel[[2]]=Panel[[2]] + theme(plot.margin=unit(c(0.20,0.50,0.20,0.20),"cm"))
Panel[[3]]=Panel[[3]] + theme(plot.margin=unit(c(0.20,0.50,0.20,0.20),"cm"))
Panel[[4]]=Panel[[4]] + theme(plot.margin=unit(c(0.20,0.50,0.20,0.20),"cm"))
Yaxis=textGrob(expression(italic('B. calyciflorus')~'density'~'('*rotifers~mL^-1*')'), gp=gpar(fontface="bold", fontsize=18), rot=90)
Xaxis=textGrob(expression('Time (days)'), gp=gpar(fontface="bold", fontsize=18), rot=0)
grid.arrange(grobs=Panel, left=Yaxis, bottom=Xaxis, ncol=1, nrow=4)
dev.off()


################################
### Calculating growth rates ###
################################

# Split the dataset
SplitData4=split(Data4, list(Data4$Strain,Data4$Bead))

# Subset important columns
SplitData5=lapply(SplitData4, "[", c("DensP","DayP"))

# Consecutive per capita growth rates function
FuncGR=function(x) {
  FitGR=dynlm(formula=DensP~L(DayP), data=as.data.frame(x))
  OutGR=c(GrowP=summary(FitGR)$coef[2,1], GrowPSD=summary(FitGR)$coef[2,2])}

# Consecutive per capita growth rates model
ModelGR=function(y) {rollapply(y, width=3, FUN=FuncGR, by.column=F)}

# Calculate consecutive per capita growth rates
OutGR=lapply(SplitData5, ModelGR)
RateGR=round(as.data.frame(do.call("rbind",OutGR)),2)
RateGR$GrowPLSD=RateGR[,1]-RateGR[,2]
RateGR$GrowPUSD=RateGR[,1]+RateGR[,2]

# Select important columns
Strain=c(do.call("rbind",lapply(SplitData4, "[", c("Strain"))))[[1]]
Bead=c(do.call("rbind",lapply(SplitData4, "[", c("Bead"))))[[1]]
DayP=c(do.call("rbind",lapply(SplitData4, "[", c("DayP"))))[[1]]
Model=c(do.call("rbind",lapply(SplitData4, "[", c("Model"))))[[1]]

# Create a dataset
DataGR=subset(data.frame(Strain,Bead,DayP,Model), !DayP %in% c("4.8","4.9","5.0"))
DataGR=data.frame(DataGR,RateGR[,c(1,3,4)]); DataGR=DataGR[,c(1,2,3,5,6,7,4)]


#######################################
### Plotting predicted growth rates ###
#######################################

# Create a dataset
Data5=DataGR[,c(1:7)]
Data5[,c(4:6)]=round(Data5[,c(4:6)],2)
Data5=Data5[order(Data5$Strain),]

# Create a dataset for labels
DataL=as.data.frame(Data5 %>% group_by(Strain,Bead) %>% filter(DayP == max(DayP)))

tiff('Growth Rates Beads.tiff', units="in", width=15, height=8, res=1000)
ggplot(Data5, aes(DayP, GrowP, group=interaction(Strain,Bead))) + 
  geom_line(aes(color=Strain), linetype="solid", size=1) +
  geom_text_repel(DataL, mapping=aes(label=Bead, color=Strain), size=5, nudge_y=0, nudge_x=0.5, direction="y", segment.color="transparent") +
  ylab(expression(italic('B. calyciflorus')~'growth rate'~'('*day^-1*')')) + xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(-2,6,by=2), limits=c(-2,6)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,5,by=1), limits=c(0,5.5)) +
  theme(panel.background=element_rect(fill="white", colour="black", size=0.7, linetype="solid")) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
  theme(strip.text.x=element_blank()) +
  facet_wrap(~Strain, scales="free", ncol=3, nrow=2) +
  theme(legend.position="none")
dev.off()


##############################################
### Plotting predicted reproduction curves ###
##############################################

# Polynomial model
ModPo=function(x) {nls(AveEgg ~ a + b * Day + c * Day^2, start=c(a=0.5, b=0.5, c=-0.5), data=x)}
ModelPo=lapply(SplitData2, ModPo)

# Create a dataset
Strain=rep(rep(unique(Data2$Strain), each=51),4)
Bead=rep(sort(unique(Data2$Bead)), each=51*6)
DayP=as.numeric(rep(seq(0,5,by=0.1),6*4))
Data6=data.frame(Strain,Bead,DayP)

# Polynomial model
EggP=list(); EggPSD=list(); EggPLCI=list(); EggPUCI=list()
for (i in 1:length(ModelPo)) {EggP[[i]]=predict(ModelPo[[i]], newdata=data.frame(Day=unique(DayP)))}
for (i in 1:length(ModelPo)) {EggPSD[[i]]=predictNLS(ModelPo[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,3]}
for (i in 1:length(ModelPo)) {EggPLCI[[i]]=predictNLS(ModelPo[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,5]}
for (i in 1:length(ModelPo)) {EggPUCI[[i]]=predictNLS(ModelPo[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,6]}
DataPo=data.frame(Data6[c(1:1224),], EggP=unlist(EggP), EggPLSD=unlist(EggP)-unlist(EggPSD), EggPUSD=unlist(EggP)+unlist(EggPSD), EggPLCI=unlist(EggPLCI), EggPUCI=unlist(EggPUCI), Model="Polynomial")

# Calculate egg number per individual
Data$AveEgg=Data$Eggs/Data$Count
Data$AveEgg[is.na(Data$AveEgg)]=0

# Combine the datasets
Data7=rbind(DataPo[,c(1:9)])
Data7[,c(4:8)]=round(Data7[,c(4:8)],4)
Data7[,c(4:8)][Data7[,c(4:8)]<0]=0
Data7=Data7[order(Data7$Strain),]

# Create a dataset for labels
DataL=as.data.frame(Data7%>% group_by(Strain,Bead) %>% filter(DayP == max(DayP)))

tiff('Reproduction Curves Beads.tiff', units="in", width=15, height=8, res=1000)
ggplot(Data7, aes(DayP, EggP, group=interaction(Strain,Bead))) +
  geom_line(aes(color=Strain, linetype=Bead), size=1) +
  geom_point(data=Data, aes(Day, AveEgg, color=Strain), size=2, pch=16, position=position_jitter(h=0, w=0.2)) +
  geom_text_repel(DataL, mapping=aes(label=Bead, color=Strain), size=5, nudge_y=0, nudge_x=0.5, direction="y", segment.color="transparent") +
  ylab(expression(italic('B. calyciflorus')~'egg ratio'~'('*eggs~rotifer^-1*')')) + xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,1.2,by=0.3), limits=c(0,1.2)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,5,by=1), limits=c(0,5.5)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
  scale_linetype_manual(values=c("C"="solid","L"="solid","M"="solid","H"="solid")) +
  theme(strip.text.x=element_blank()) +
  facet_wrap(~Strain, scales="free", ncol=3, nrow=2) +
  theme(legend.position="none")
dev.off()


##############################################
### Plotting predicted reproduction curves ###
##############################################

# Arrange the datasets
Data=Data %>% arrange(factor(Bead, levels=unique(Data$Bead)))
Data7=Data7 %>% arrange(factor(Bead, levels=unique(Data7$Bead)))

# Include experimental points
Data7$DayE=c(Reduce(rbind,append(list(matrix(Data$Day,30)),rep(NA,21))))
Data7$EggE=c(Reduce(rbind,append(list(matrix(Data$AveEgg,30)),rep(NA,21))))

# Include positions for labels
DataL=as.data.frame(Data7 %>% group_by(Bead,Strain) %>% filter(DayP == max(DayP)))
Data7$DayM=c(Reduce(rbind,append(list(matrix(DataL$DayP,1)),rep(NA,50))))
Data7$EggM=c(Reduce(rbind,append(list(matrix(DataL$EggP,1)),rep(NA,50))))

# Create names for labels
Data7$Labels=factor(Data7$Strain, levels=unique(Data7$Strain), labels=c(expression(C[R1]),expression(C[R2]),expression(C[R3]),expression(C[R4]),expression(C[R6]),expression(C[R7])))

# Split the dataset
SplitData7=split(Data7, list(Data7$Bead))

PlotFunc=function(x) {
  ggplot(x, aes(DayP, EggP, group=interaction(Strain,Bead))) +
    geom_line(aes(color=Strain, linetype=Bead), size=1) +
    geom_point(data=x, aes(DayE, EggE, color=Strain), size=2, pch=16, position=position_jitter(h=0, w=0.2)) +
    ylab(expression(italic('B. calyciflorus')~'egg ratio'~'('*eggs~rotifer^-1*')')) + xlab(expression('Time (days)')) +
    theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
    theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
    theme(axis.title.y=element_blank()) +
    theme(axis.title.x=element_blank()) +
    scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,5,by=1), limits=c(0,5)) +
    theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
    scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
    scale_linetype_manual(values=c("C"="solid","L"="solid","M"="solid","H"="solid")) +
    theme(legend.position="none")
}

tiff('Reproduction Curves Beads.tiff', units="in", width=12, height=12, res=1000)
Panel=lapply(SplitData7, PlotFunc)
Panel[[1]]=Panel[[1]] + scale_y_continuous(labels=sprintf(seq(0,1.2,by=0.3), fmt="%.1f"), breaks=seq(0,1.2,by=0.3), limits=c(0,1.2))
Panel[[2]]=Panel[[2]] + scale_y_continuous(labels=sprintf(seq(0,1.2,by=0.3), fmt="%.1f"), breaks=seq(0,1.2,by=0.3), limits=c(0,1.2))
Panel[[3]]=Panel[[3]] + scale_y_continuous(labels=sprintf(seq(0,1.2,by=0.3), fmt="%.1f"), breaks=seq(0,1.2,by=0.3), limits=c(0,1.2))
Panel[[4]]=Panel[[4]] + scale_y_continuous(labels=sprintf(seq(0,1.2,by=0.3), fmt="%.1f"), breaks=seq(0,1.2,by=0.3), limits=c(0,1.2))
Panel[[1]]=Panel[[1]] + theme(plot.margin=unit(c(0.20,0.20,0.20,0.20),"cm"))
Panel[[2]]=Panel[[2]] + theme(plot.margin=unit(c(0.20,0.20,0.20,0.20),"cm"))
Panel[[3]]=Panel[[3]] + theme(plot.margin=unit(c(0.20,0.20,0.20,0.20),"cm"))
Panel[[4]]=Panel[[4]] + theme(plot.margin=unit(c(0.20,0.20,0.20,0.20),"cm"))
Yaxis=textGrob(expression(italic('B. calyciflorus')~'egg ratio'~'('*eggs~rotifer^-1*')'), gp=gpar(fontface="bold", fontsize=18), rot=90)
Xaxis=textGrob(expression('Time (days)'), gp=gpar(fontface="bold", fontsize=18), rot=0)
grid.arrange(grobs=Panel, left=Yaxis, bottom=Xaxis, ncol=2, nrow=2)
dev.off()

PlotFunc=function(x) {
  ggplot(x, aes(DayP, EggP, group=interaction(Strain,Bead))) +
    geom_line(aes(color=Strain, linetype=Bead), size=1) +
    geom_ribbon(aes(ymin=EggPLSD, ymax=EggPUSD, fill=Strain, color=Strain), linetype="solid", size=0.5, alpha=0.3) +
    geom_point(data=x, aes(DayE, EggE, color=Strain), size=2, pch=16, position=position_jitter(h=0, w=0.2)) +
    geom_text(data=x, mapping=aes(x=-Inf, y=Inf, label=Bead), color="black", size=5, vjust=1.2, hjust=-1.5, parse=T) +
    ylab(expression(italic('B. calyciflorus')~'egg ratio'~'('*eggs~rotifer^-1*')')) + xlab(expression('Time (days)')) +
    theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
    theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
    theme(axis.title.y=element_blank()) +
    theme(axis.title.x=element_blank()) +
    scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,5,by=1), limits=c(0,5)) +
    theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
    scale_fill_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
    scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
    scale_linetype_manual(values=c("C"="solid","L"="solid","M"="solid","H"="solid")) +
    theme(strip.text.x=element_blank()) + theme(panel.spacing.x=unit(0.5,"cm")) +
    facet_wrap(~Strain, scales="free", ncol=6, nrow=1) +
    theme(legend.position="none")
}

tiff('Reproduction Curves Beads Intervals.tiff', units="in", width=20, height=13, res=1000)
Panel=lapply(SplitData7, PlotFunc)
Panel[[1]]=Panel[[1]] + scale_y_continuous(labels=sprintf(seq(0,1.2,by=0.3), fmt="%.1f"), breaks=seq(0,1.2,by=0.3), limits=c(0,1.2))
Panel[[2]]=Panel[[2]] + scale_y_continuous(labels=sprintf(seq(0,1.2,by=0.3), fmt="%.1f"), breaks=seq(0,1.2,by=0.3), limits=c(0,1.2))
Panel[[3]]=Panel[[3]] + scale_y_continuous(labels=sprintf(seq(0,1.2,by=0.3), fmt="%.1f"), breaks=seq(0,1.2,by=0.3), limits=c(0,1.2))
Panel[[4]]=Panel[[4]] + scale_y_continuous(labels=sprintf(seq(0,1.2,by=0.3), fmt="%.1f"), breaks=seq(0,1.2,by=0.3), limits=c(0,1.2))
Panel[[1]]=Panel[[1]] + theme(plot.margin=unit(c(0.20,0.20,0.20,0.20),"cm"))
Panel[[2]]=Panel[[2]] + theme(plot.margin=unit(c(0.20,0.20,0.20,0.20),"cm"))
Panel[[3]]=Panel[[3]] + theme(plot.margin=unit(c(0.20,0.20,0.20,0.20),"cm"))
Panel[[4]]=Panel[[4]] + theme(plot.margin=unit(c(0.20,0.20,0.20,0.20),"cm"))
Yaxis=textGrob(expression(italic('B. calyciflorus')~'egg ratio'~'('*eggs~rotifer^-1*')'), gp=gpar(fontface="bold", fontsize=18), rot=90)
Xaxis=textGrob(expression('Time (days)'), gp=gpar(fontface="bold", fontsize=18), rot=0)
grid.arrange(grobs=Panel, left=Yaxis, bottom=Xaxis, ncol=1, nrow=4)
dev.off()
