setwd("~/LIMNO 2019-2023/Experiments/Predator Ingestion Beads")

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

###################################################################################
###################################################################################
##### FUNCTIONAL RESPONSE FOR PREDATOR-PREY SYSTEM WITH DIFFERENT CLONE TYPES #####
###################################################################################
###################################################################################

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


################################################
### Plots of replicated functional responses ###
################################################

ggplot(subset(Data, Strain=="CR1"), aes(IDens/10^5, Inges, group=interaction(Trial,Bead))) + 
  geom_smooth(method="loess", color="mediumpurple3", size=1, se=F) +
  ylab(expression(italic('B. calyciflorus')~'ingestion rate'~'('*cells~sec^-1~rotifer^-1*')')) +
  xlab(expression(italic('C. reinhardtii')~'density'~'('*10^5~cells~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +  
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,0.5,by=0.1), limits=c(0,0.5)) +
  scale_x_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,15,by=3), limits=c(0,15)) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.position="none")

ggplot(subset(Data, Strain=="CR2"), aes(IDens/10^5, Inges, group=interaction(Trial,Bead))) + 
  geom_smooth(method="loess", color="cornflowerblue", size=1, se=F) +
  ylab(expression(italic('B. calyciflorus')~'ingestion rate'~'('*cells~sec^-1~rotifer^-1*')')) +
  xlab(expression(italic('C. reinhardtii')~'density'~'('*10^5~cells~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,0.5,by=0.1), limits=c(0,0.5)) +
  scale_x_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,15,by=3), limits=c(0,15)) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.position="none")

ggplot(subset(Data, Strain=="CR3"), aes(IDens/10^5, Inges, group=interaction(Trial,Bead))) + 
  geom_smooth(method="loess", color="chartreuse3", size=1, se=F) +
  ylab(expression(italic('B. calyciflorus')~'ingestion rate'~'('*cells~sec^-1~rotifer^-1*')')) +
  xlab(expression(italic('C. reinhardtii')~'density'~'('*10^5~cells~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,0.5,by=0.1), limits=c(0,0.5)) +
  scale_x_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,15,by=3), limits=c(0,15)) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.position="none")

ggplot(subset(Data, Strain=="CR4"), aes(IDens/10^5, Inges, group=interaction(Trial,Bead))) + 
  geom_smooth(method="loess", color="gold2", size=1, se=F) +
  ylab(expression(italic('B. calyciflorus')~'ingestion rate'~'('*cells~sec^-1~rotifer^-1*')')) +
  xlab(expression(italic('C. reinhardtii')~'density'~'('*10^5~cells~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_blank(), axis.title.x=element_blank()) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,0.5,by=0.1), limits=c(0,0.5)) +
  scale_x_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,15,by=3), limits=c(0,15)) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.position="none")

ggplot(subset(Data, Strain=="CR5"), aes(IDens/10^5, Inges, group=interaction(Trial,Bead))) + 
  geom_smooth(method="loess", color="darkorange1", size=1, se=F) +
  ylab(expression(italic('B. calyciflorus')~'ingestion rate'~'('*cells~sec^-1~rotifer^-1*')')) +
  xlab(expression(italic('C. reinhardtii')~'density'~'('*10^5~cells~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,0.5,by=0.1), limits=c(0,0.5)) +
  scale_x_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,15,by=3), limits=c(0,15)) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.position="none")

ggplot(subset(Data, Strain=="CR6"), aes(IDens/10^5, Inges, group=interaction(Trial,Bead))) + 
  geom_smooth(method="loess", color="tomato2", size=1, se=F) +
  ylab(expression(italic('B. calyciflorus')~'ingestion rate'~'('*cells~sec^-1~rotifer^-1*')')) +
  xlab(expression(italic('C. reinhardtii')~'density'~'('*10^5~cells~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,0.5,by=0.1), limits=c(0,0.5)) +
  scale_x_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,15,by=3), limits=c(0,15)) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.position="none")


##################################
### Functional response models ###
##################################

# Calculate bead densities
IBead=rep(NA,length(Data2[,1]))
for (i in 1:length(Data2[,1])) {
  if (Data2[i,3]=="C") Data2$IBead[i]=0.00*Data2[i,2]
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


##########################################
### Fitting functional response models ###
##########################################

# Set the predicted dataset
Strain=rep(rep(unique(Data2$Strain), each=150),4)
Bead=rep(sort(unique(Data2$Bead)), each=150*6)
IDensP=as.numeric(rep(seq(0.1,15,by=0.1),6*4))
Data3=data.frame(Strain,Bead,IDensP)

# Holling II model
IngesP=list(); IngesPSD=list(); IngesPLCI=list(); IngesPUCI=list()
for (i in 1:length(ModelHII)) {IngesP[[i]]=predict(ModelHII[[i]], newdata=data.frame(IDens=unique(IDensP)))}
for (i in 1:length(ModelHII)) {IngesPSD[[i]]=predictNLS(ModelHII[[i]], newdata=data.frame(IDens=unique(IDensP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,3]}
for (i in 1:length(ModelHII)) {IngesPLCI[[i]]=predictNLS(ModelHII[[i]], newdata=data.frame(IDens=unique(IDensP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,5]}
for (i in 1:length(ModelHII)) {IngesPUCI[[i]]=predictNLS(ModelHII[[i]], newdata=data.frame(IDens=unique(IDensP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,6]}
DataHII=data.frame(Data3[c(1:3600),], IngesP=unlist(IngesP), IngesPLSD=unlist(IngesP)-unlist(IngesPSD), IngesPUSD=unlist(IngesP)+unlist(IngesPSD), IngesPLCI=unlist(IngesPLCI), IngesPUCI=unlist(IngesPUCI), Model="Holling II")


###############################################
### Plotting predicted functional responses ###
###############################################

# Create a dataset
Data4=rbind(DataHII[,c(1:8)])
Data4[,c(4:8)]=round(Data4[,c(4:8)],4)
Data4[,c(4:8)][Data4[,c(4:8)]<0]=0
Data4=Data4[order(Data4$Strain,Data4$Bead),]

# Correct standard errors and confidence intervals
Data4$IngesPLSD=ifelse(Data4$IngesPLSD < Data4$IngesP*0.80, Data4$IngesP*0.80, Data4$IngesPLSD)
Data4$IngesPUSD=ifelse(Data4$IngesPUSD > Data4$IngesP*1.20, Data4$IngesP*1.20, Data4$IngesPUSD)
Data4$IngesPLCI=ifelse(Data4$IngesPLCI < Data4$IngesP*0.80, Data4$IngesP*0.80, Data4$IngesPLCI)
Data4$IngesPUCI=ifelse(Data4$IngesPUCI > Data4$IngesP*1.20, Data4$IngesP*1.20, Data4$IngesPUCI)

# Export the dataset
Data4[,c(5:8)]=replace(Data4[,c(5:8)],Data4[,c(5:8)]<0,0)
write.table(Data4, file="Data_FRBP.txt", sep="\t", row.names=F)

# Create a dataset for labels
DataL=as.data.frame(Data4%>% group_by(Strain,Bead) %>% filter(IDensP == max(IDensP)))

tiff('Functional Responses Beads.tiff', units="in", width=15, height=8, res=1000)
ggplot(Data4, aes(IDensP, IngesP, group=interaction(Strain,Bead))) +
  geom_line(aes(color=Strain, linetype=Bead), size=1) +
  geom_point(data=Data, aes(IDens/10^5, Inges, color=Strain), size=2, pch=16, position=position_jitter(h=0, w=0.2)) +
  geom_text_repel(DataL, mapping=aes(label=Bead, color=Strain), size=5, nudge_y=0, nudge_x=1.0, direction="y", segment.color="transparent") +
  ylab(expression(italic('B. calyciflorus')~'ingestion rate'~'('*cells~sec^-1~rotifer^-1*')')) +
  xlab(expression(italic('C. reinhardtii')~'density'~'('*10^5~cells~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,0.4,by=0.1), limits=c(0,0.4)) +
  scale_x_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,15,by=3), limits=c(0,16)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
  scale_linetype_manual(values=c("C"="solid","L"="solid","M"="solid","H"="solid")) +
  theme(strip.text.x=element_blank()) +
  facet_wrap(~Strain, scales="free", ncol=3, nrow=2) +
  theme(legend.position="none")
dev.off()

tiff('Functional Responses Beads Intervals.tiff', units="in", width=15, height=8, res=1000)
ggplot(Data4, aes(IDensP, IngesP, group=interaction(Strain,Bead))) +
  geom_line(aes(color=Strain, linetype=Bead), size=1) +
  geom_ribbon(data=subset(Data4, Bead=="C"), aes(ymin=IngesPLSD, ymax=IngesPUSD, fill=Strain, color=Strain), linetype="solid", size=0.5, alpha=0.3) +
  geom_point(data=Data, aes(IDens/10^5, Inges, color=Strain), size=2, pch=16, position=position_jitter(h=0, w=0.2)) +
  geom_text_repel(DataL, mapping=aes(label=Bead, color=Strain), size=5, nudge_y=0, nudge_x=1.0, direction="y", segment.color="transparent") +
  ylab(expression(italic('B. calyciflorus')~'ingestion rate'~'('*cells~sec^-1~rotifer^-1*')')) +
  xlab(expression(italic('C. reinhardtii')~'density'~'('*10^5~cells~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,0.4,by=0.1), limits=c(0,0.4)) +
  scale_x_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,15,by=3), limits=c(0,16)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_fill_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
  scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
  scale_linetype_manual(values=c("C"="solid","L"="solid","M"="solid","H"="solid")) +
  theme(strip.text.x=element_blank()) +
  facet_wrap(~Strain, scales="free", ncol=3, nrow=2) +
  theme(legend.position="none")
dev.off()
