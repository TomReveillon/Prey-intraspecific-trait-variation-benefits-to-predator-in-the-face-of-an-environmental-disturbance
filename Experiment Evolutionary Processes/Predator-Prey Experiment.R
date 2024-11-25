setwd("~/LIMNO 2019-2023/Experiments/Population Dynamics")

rm(list=ls())

library(car)
library(cowplot)
library(data.table)
library(deSolve)
library(devtools)
library(directlabels)
library(dplyr)
library(dynlm)
library(emmeans)
library(foreach)
library(ggplot2)
library(ggrepel)
library(grid)
library(gridExtra)
library(gsubfn)
library(keras)
library(lme4)
library(lmtest)
library(magrittr)
library(MASS)
library(nlme)
library(plotly)
library(plyr)
library(pracma)
library(propagate)
library(reshape2)
library(reticulate)
library(secr)
library(tensorflow)
library(tidyverse)

####################################################################
####################################################################
##### INDIRECT EVOLUTIONARY PROCESSES FOR PREDATOR-PREY SYSTEM #####
####################################################################
####################################################################

# Import the dataset for predator
Data1=read.table("Data_RD2.txt", h=T, dec=",")
names(Data1)
summary(Data1)

# Import the dataset for preys
Data2=read.table("Data_AD2.txt", h=T, dec=",")
names(Data2)
summary(Data2)

# Specify the variables as numeric or factor
Data1[,c(4,6)] %<>% mutate_if(is.character,as.numeric)
Data2[,c(3:6)] %<>% mutate_if(is.character,as.numeric)

# Calculate the densities
Data1$Dens=Data1$Count/Data1$Dilu
Data2$Dens=Data2$Cells*Data2$Volu*Data2$Site*Data2$Dilu*Data2$Cove

# Preserve the order in the dataset
Data1$Strain=factor(Data1$Strain, levels=unique(Data1$Strain))
Data1$Bead=factor(Data1$Bead, levels=unique(Data1$Bead))


########################################################
########################################################
##### Phenotyping and calculating prey frequencies #####
########################################################
########################################################

# Load the training datasets
Folder1="Training 2 Strains/Day 1"
Folder2="Training 2 Strains/Day 10"
Folder3="Training 2 Strains/Day 20"
Folder4="Training 2 Strains/Day 30"
Folder5="Training 2 Strains/Day 40"
Folder6="Training 2 Strains/Day 50"
Folder=c(Folder1,Folder2,Folder3,Folder4,Folder5,Folder6)
Files=list.files(path=Folder, pattern="*.txt", full.names=T)
DataTrain=ldply(Files, read.table, sep="\t", fill=T, header=T, dec=",")

# Load the sample datasets
Folder1="Collection 2 Strains/Days 1-3"
Folder2="Collection 2 Strains/Days 4-6"
Folder3="Collection 2 Strains/Days 7-9"
Folder4="Collection 2 Strains/Days 10-12"
Folder5="Collection 2 Strains/Days 13-15"
Folder6="Collection 2 Strains/Days 16-18"
Folder7="Collection 2 Strains/Days 19-21"
Folder8="Collection 2 Strains/Days 22-24"
Folder9="Collection 2 Strains/Days 25-27"
Folder10="Collection 2 Strains/Days 28-30"
Folder11="Collection 2 Strains/Days 31-33"
Folder12="Collection 2 Strains/Days 34-36"
Folder13="Collection 2 Strains/Days 37-39"
Folder14="Collection 2 Strains/Days 40-42"
Folder15="Collection 2 Strains/Days 43-45"
Folder16="Collection 2 Strains/Days 46-48"
Folder17="Collection 2 Strains/Days 49-50"
Folder=c(Folder1,Folder2,Folder3,Folder4,Folder5,Folder6,Folder7,Folder8,Folder9,Folder10,Folder11,Folder12,Folder13,Folder14,Folder15,Folder16,Folder17)
Files=list.files(path=Folder, pattern="*.txt", full.names=T)
DataExpe=ldply(Files, read.table, sep="\t", fill=T, header=T, dec=",")

# Mutate dataset classes to numeric
DataTrain %<>% mutate_if(is.factor,as.numeric); DataTrain %<>% mutate_if(is.character,as.numeric); DataTrain %<>% mutate_if(is.integer,as.numeric)
DataExpe %<>% mutate_if(is.factor,as.numeric); DataExpe %<>% mutate_if(is.character,as.numeric); DataExpe %<>% mutate_if(is.integer,as.numeric)

# Change column names and select autofluoresence
colnames(DataTrain)=gsub("\\_", "\\.", colnames(DataTrain))
DataTrain=DataTrain %>% select_if(grepl("M01|M05|M06", names(.)))
colnames(DataTrain)=gsub("\\.M01", "\\.TL", colnames(DataTrain))
colnames(DataTrain)=gsub("\\.M05", "\\.RL", colnames(DataTrain))
colnames(DataTrain)=gsub("\\.M06", "\\.PL", colnames(DataTrain))
colnames(DataTrain)=gsub("\\.", " ", colnames(DataTrain))

colnames(DataExpe)=gsub("\\_", "\\.", colnames(DataExpe))
DataExpe=DataExpe %>% select_if(grepl("M01|M05|M06", names(.)))
colnames(DataExpe)=gsub("\\.M01", "\\.TL", colnames(DataExpe))
colnames(DataExpe)=gsub("\\.M05", "\\.RL", colnames(DataExpe))
colnames(DataExpe)=gsub("\\.M06", "\\.PL", colnames(DataExpe))
colnames(DataExpe)=gsub("\\.", " ", colnames(DataExpe))

# Create treatment identities
Strain=1; Bead=3; Day=50; Obs=300; Rep=3
DataTrain$Strain=rep(c(rep("CR2",1000*Rep),rep("CR4",1000*Rep)),6)
DataExpe$Strain=rep(c(rep("CR2-CR4",Obs)),Rep*Day*Bead)
DataExpe$Bead=rep(c(rep("C",Obs*Rep),rep("L",Obs*Rep),rep("H",Obs*Rep)),Day)
DataExpe$Trial=rep(c(rep(seq(1,3,by=1), each=Obs)),Day*Bead)
DataExpe$Day=rep(seq(1,50,by=1), each=Obs*Rep*Bead)

# Preserve the order in the dataset
DataTrain$Strain=factor(DataTrain$Strain, levels=unique(DataTrain$Strain))
DataExpe$Strain=factor(DataExpe$Strain, levels=unique(DataExpe$Strain))
DataExpe$Bead=factor(DataExpe$Bead, levels=unique(DataExpe$Bead))
DataExpe$Day=factor(DataExpe$Day, levels=unique(DataExpe$Day))


##############################################
### Preparing data for deep learning model ###
##############################################

# Create training and validation datasets
TrainCR2=subset(DataTrain, Strain=="CR2")[c(1:15000),]
ValidCR2=subset(DataTrain, Strain=="CR2")[c(15001:16500),]
TestCR2=subset(DataTrain, Strain=="CR2")[c(16501:18000),]

TrainCR4=subset(DataTrain, Strain=="CR4")[c(1:15000),]
ValidCR4=subset(DataTrain, Strain=="CR4")[c(15001:16500),]
TestCR4=subset(DataTrain, Strain=="CR4")[c(16501:18000),]

# Combine and randomize the datasets
DataTrain2=rbind(TrainCR2,TrainCR4)
DataTrain2$ID=c(rep(0,15000),rep(1,15000))
DataTrain2=DataTrain2[sample(nrow(DataTrain2)),]

DataValid2=rbind(ValidCR2,ValidCR4)
DataValid2$ID=c(rep(0,1500),rep(1,1500))
DataValid2=DataValid2[sample(nrow(DataValid2)),]

DataTest2=rbind(TestCR2,TestCR4)
DataTest2$ID=c(rep(0,1500),rep(1,1500))
DataTest2=DataTest2[sample(nrow(DataTest2)),]

# Normalize the datasets
DataTrain3=subset(DataTrain2, select=-c(Strain,ID))
DataTrain3=scale(DataTrain3, center=apply(DataTrain3,2,mean), scale=apply(DataTrain3,2,sd))
DataTrain3=as.data.frame(DataTrain3)
DataTrain3$ID=DataTrain2$ID

DataValid3=subset(DataValid2, select=-c(Strain,ID))
DataValid3=scale(DataValid3, center=apply(DataValid3,2,mean), scale=apply(DataValid3,2,sd))
DataValid3=as.data.frame(DataValid3)
DataValid3$ID=DataValid2$ID

DataTest3=subset(DataTest2, select=-c(Strain,ID))
DataTest3=scale(DataTest3, center=apply(DataTest3,2,mean), scale=apply(DataTest3,2,sd))
DataTest3=as.data.frame(DataTest3)
DataTest3$ID=DataTest2$ID

# Format the datasets
DataTrainX=as.matrix(DataTrain3[,-c(118)])
DataTrainY=to_categorical(DataTrain3$ID, 2)

DataValidX=as.matrix(DataValid3[,-c(118)])
DataValidY=to_categorical(DataValid3$ID, 2)

DataTestX=as.matrix(DataTest3[,-c(118)])
DataTestY=to_categorical(DataTest3$ID, 2)

# Deep learning model
Model=keras_model_sequential() 
Model %>% 
  layer_dense(units=100, activation="relu", input_shape=c(117)) %>%
  layer_dropout(rate=0.2) %>%  
  layer_dense(units=50, activation="relu", input_shape=c(117)) %>%
  layer_dropout(rate=0.2) %>%  
  layer_dense(units=10, activation="relu", input_shape=c(117)) %>%
  layer_dropout(rate=0.2) %>%  
  layer_dense(units=2, activation="softmax")

# Optimize model learning
Model %>% compile(
  loss="categorical_crossentropy",
  optimizer=optimizer_rmsprop(learning_rate=0.001),
  metrics="accuracy"
)

######################################################
### Training and validation of deep learning model ###
######################################################

# Plots of lost information and accuracy
Plot=Model %>% fit(DataTrainX, DataTrainY, epochs=50, batch_size=100, validation_split=0.2, validation_data=list(DataValidX,DataValidY))

# Evaluation of the model
Model %>% evaluate(DataTestX, DataTestY, verbose=0)

# Prediction of the model
DataPred=Model %>% predict(DataTestX)
colnames(DataPred)=c("CR2","CR4")
rownames(DataPred)=c(seq(1,3000,by=1))

# Summarize the predictions
DataPred2=data.frame(TrueID=DataTest3$ID+1, PredID=rep(NA,3000), Accuracy=rep(NA,3000), Confirmation=rep(NA,3000))
for (i in 1:length(DataPred[,1])){
  DataPred2$PredID[i]=which.max(DataPred[i,])
  DataPred2$Accuracy[i]=round(DataPred[i,which.max(DataPred[i,])],4)
  DataPred2$Confirmation[i]=DataPred2$PredID[i]==DataPred2$TrueID[i]
  if(DataPred2$TrueID[i]==1) {DataPred2$TrueID[i]="CR2"} 
  if(DataPred2$TrueID[i]==2) {DataPred2$TrueID[i]="CR4"} 
  if(DataPred2$PredID[i]==1) {DataPred2$PredID[i]="CR2"} 
  if(DataPred2$PredID[i]==2) {DataPred2$PredID[i]="CR4"} 
}


##########################################
### Application of deep learning model ###
##########################################

# Format the datasets
DataExpe2=subset(DataExpe, select=-c(Strain,Bead,Trial,Day))
DataExpe2=scale(DataExpe2, center=apply(DataExpe2,2,mean), scale=apply(DataExpe2,2,sd))
DataExpe2=as.data.frame(DataExpe2)
DataExpe2$Strain=DataExpe$Strain
DataExpe2$Bead=DataExpe$Bead
DataExpe2$Trial=DataExpe$Trial
DataExpe2$Day=DataExpe$Day

# Split the datasets
SplitExpe2=split(DataExpe2, list(DataExpe2$Strain,DataExpe2$Bead,DataExpe2$Trial,DataExpe2$Day))
SplitExpe2=SplitExpe2[grep("CR2-CR4", names(SplitExpe2))]  

# Normalize the datasets
SplitExpe2=lapply(SplitExpe2, "[", -c(118:121))
SplitExpe3=lapply(SplitExpe2, as.matrix)

# Predictions of the model
Func=function(x) {Out=Model %>% predict(x)}
ListPred=lapply(SplitExpe3, Func)

# Bind the datasets
DataPred=as.data.frame(do.call("rbind",ListPred))
colnames(DataPred)=c("CR2","CR4"); rownames(DataPred)=c()

# Summarize the predictions
DataPred2=data.frame(PredID=rep(NA,135000))
for (i in 1:length(DataPred[,1])){
  DataPred2$PredID[i]=which.max(DataPred[i,])
  if(DataPred2$PredID[i]==1) {DataPred2$PredID[i]="CR2"} 
  if(DataPred2$PredID[i]==2) {DataPred2$PredID[i]="CR4"} 
}

# Add column for treatments
DataPred2=cbind(DataPred2,subset(DataExpe2, Strain=="CR2-CR4")[,c(118:121)])
DataPred2=droplevels(DataPred2)

# Sum the predictions
DataPred3=data.frame(table(Strain=DataPred2[,2], Bead=DataPred2[,3], Trial=DataPred2[,4], Day=DataPred2[,5], Alga=DataPred2[,1]))
DataPred3=rbind(DataPred3)
colnames(DataPred3)[6]="Count"


#############################
### Preparing the dataset ###
#############################

# Merge the datasets
Data3=Data2[Data2$Day!=0,]
Data4=merge(Data3[,-c(3:6)], DataPred3, all=F)
Data4$Count[Data4$Count==300]=299
Data4$Count[Data4$Count==0]=1

# Create a dataset for monocultures
Data5=subset(Data2, Strain=="CR2"|Strain=="CR4")[,-c(3:6)]
Data5$Alga=Data5$Strain
Data5$Count=300

# Create a dataset for starting days
Data6=Data2[order(Data2$Strain),]
Data6=subset(Data6, !Strain %in% c("CR2","CR4") & Day==0)
Data6=Data6[rep(seq_len(nrow(Data6)), each=2),][,-c(3:6)]
Data6$Alga=c(rep(c("CR2","CR4"),9))
Data6$Count=150

# Create a complete dataset
Data7=rbind(Data4,Data5,Data6)

# Calculate prey frequencies
Data7$Freq=Data7$Count/300

# Calculate prey densities
Data7$DensA=Data7$Dens*Data7$Freq

# Calculate prey defense
Data7$Defense=as.numeric(as.character(revalue(Data7$Alga,c("CR2"=1.0, "CR4"=0.1))))

# Arrange the dataset
Data7=Data7 %>% arrange(factor(Strain, levels=c("CR2","CR4","CR2-CR4")), factor(Bead, levels=c("C","L","H")), factor(Trial, levels=c("1","2","3")), factor(Day, levels=seq(0,50,by=1)))

# Calculate mean prey defense
Data8=setDT(Data7)[, .(Defense=sum(Defense*Freq)/sum(Freq)), by=list(Strain,Bead,Day,Trial)]
Data8a=subset(Data8, c(Strain=="CR2"|Strain=="CR4"))
Data8b=subset(Data8, !c(Strain=="CR2"|Strain=="CR4"))

# Include mean prey defense
Data7$Defense=rbind(Data8a[,5], Data8b[rep(seq_len(nrow(Data8b)), each=2),][,5])
Data7=as.data.frame(Data7)

# Create a dataset for predator
Data9a=subset(Data1, c(Strain=="CR2"|Strain=="CR4"))
Data9b=subset(Data1, !c(Strain=="CR2"|Strain=="CR4"))
Data9b=Data9b[rep(seq_len(nrow(Data9b)), each=2),]
Data10=rbind(Data9a,Data9b); colnames(Data10)[7]="DensR"

# Arrange the dataset
Data10=Data10 %>% arrange(factor(Strain, levels=c("CR2","CR4","CR2-CR4")), factor(Bead, levels=c("C","L","H")), factor(Trial, levels=c("1","2","3")), factor(Day, levels=seq(0,50,by=1)))


###################################
### Analyze population dynamics ###
###################################

# Combine the datasets
Data11=data.frame(Data7[,c(1:4,7)], A=Data7[,10], R=Data10[,7])
Data11[,c(6:7)]=round(Data11[,c(6:7)],0)

# Melt the dataset
Data12=melt(Data11, id.vars=c("Strain","Bead","Day","Trial","Alga"), variable.name="Organism", value.name="Density")

# Rename the organisms
Organism=rep(NA,length(Data12[,1]))
for (i in 1:length(Data12[,1])) {
  if (Data12[i,5]=="CR2" & Data12[i,6]=="A")     
    Organism[i]="CR2"
  else if(Data12[i,5]=="CR4" & Data12[i,6]=="A")
    Organism[i]="CR4"
  else (Organism[i]="R")}
Data12$Organism=Organism

# Split the dataset
SplitData12=split(Data12, list(Data12$Strain,Data12$Bead,Data12$Trial,Data12$Organism))

# Find peaks and valleys for densities
Peaks=list(); Valleys=list()
FuncWei=function(x) {
  Peaks=findpeaks(x$Density, nups=2, ndowns=2)
  Valleys=findpeaks(-x$Density, nups=2, ndowns=2)
  Names=unique(x[1,c(1:2,4:6)])
  DataP=data.frame(Names[rep(seq_len(nrow(Names)), length(Peaks[,2])),], Peak=Peaks[,2])
  DataV=data.frame(Names[rep(seq_len(nrow(Names)), length(Valleys[,2])),], Valley=Valleys[,2])
  Data=list(Peaks=DataP, Valleys=DataV)
}
OutWei=lapply(SplitData12, FuncWei)

# Extract peaks and valleys for densities
DataP=bind_rows(lapply(OutWei, function (x) x[c("Peaks")]))
DataP=as.data.frame(do.call("rbind",DataP))
rownames(DataP)=c()

DataV=bind_rows(lapply(OutWei, function (x) x[c("Valleys")]))
DataV=as.data.frame(do.call("rbind",DataV))
rownames(DataV)=c()

# Include weights for origins or peaks and valleys
DataP=merge(Data12, DataP, all=T)
DataP=DataP[order(DataP$Strain,DataP$Bead,DataP$Trial,DataP$Alga,DataP$Organism,DataP$Day),]
DataP$Peak=ifelse(DataP[,6]==DataP[,8], 10, 1)
DataP$Peak=ifelse(DataP[,6]==0, 10, DataP[,8])
DataP$Peak=ifelse(is.na(DataP[,8]), 1, DataP[,8])
DataP=DataP[!duplicated(DataP[,c(1:6)]),]

DataV=merge(Data12, DataV, all=T)
DataV=DataV[order(DataV$Strain,DataV$Bead,DataV$Trial,DataV$Alga,DataV$Organism,DataV$Day),]
DataV$Valley=ifelse(DataV[,6]==DataV[,8], 10, 1)
DataV$Valley=ifelse(DataV[,6]==0, 10, DataV[,8])
DataV$Valley=ifelse(is.na(DataV[,8]), 1, DataV[,8])
DataV=DataV[!duplicated(DataV[,c(1:6)]),]

# Combine the datasets
Data12=merge(DataP, DataV, all=T)
Data12$Weight=ifelse(Data12[,8]==10|Data12[,9]==10, 10, 1)
Data12=Data12[,-c(8:9)]


#############################################
### Calculating coefficients of variation ###
#############################################

# Find extinction dynamics
Data13=subset(Data12, Day==50 & Organism=="R" & Density==0)
Data13a=subset(Data12, Strain %in% c(Data13[,1]) & Bead %in% c(Data13[,2]) & Trial %in% c(Data13[,3]))

# Remove extinction days
Data13a$Density=ifelse(Data13a$Density==0, NA, Data13a$Density)
Data13a=Data13a %>% group_by(Strain,Bead,Trial,Day) %>% filter(!any(is.na(Density)))
Data13a=as.data.frame(Data13a)

# Merge the dataset
Data13b=subset(Data12, !(Strain %in% c(Data13[,1]) & Bead %in% c(Data13[,2]) & Trial %in% c(Data13[,3])))
Data14=rbind(Data13a,Data13b)

# Prepare the dataset
Data14=Data14[!duplicated(Data14),]
Data14$Organism=ifelse(Data14$Organism=="CR2", "A", Data14$Organism)
Data14$Organism=ifelse(Data14$Organism=="CR4", "A", Data14$Organism)

# Calculate the coefficients
Data15=setDT(Data14)[, .(CV=round((sd(Density)/mean(Density))*100,2)), by=list(Strain,Bead,Trial,Organism)]
Data15=Data15[order(Data15$Strain,Data15$Bead,Data15$Trial),]
Data15$CV=ifelse(is.na(Data15$CV), 0, Data15$CV)
Data15=as.data.frame(Data15)


###############################################
### Plotting replicated population dynamics ###
###############################################

# Arrange the dataset
Data12=Data12 %>% arrange(factor(Strain, levels=c("CR2","CR4","CR2-CR4")), factor(Bead, levels=c("C","L","H")), factor(Trial, levels=c("1","2","3")), factor(Day, levels=seq(0,50,by=1)))
Data12$Strain=factor(Data12$Strain, levels=c("CR2","CR4","CR2-CR4"))
Data15$Strain=factor(Data15$Strain, levels=c("CR2","CR4","CR2-CR4"))
Data12$Bead=factor(Data12$Bead, levels=c("C","L","H"))
Data15$Bead=factor(Data15$Bead, levels=c("C","L","H"))
Data12$Trial=factor(Data12$Trial, levels=c("1","2","3"))
Data15$Trial=factor(Data15$Trial, levels=c("1","2","3"))

# Extract the dataset
write.table(Data12, file="~/Activité Professionnelle/LIMNO 2019-2023/Experiments/Population Dynamics/Data_EXP1.txt", sep="\t", row.names=F)

# Rescale densities
Plot12=Data12[,c(1:8)]
Plot12[Plot12$Organism=="CR2",]$Density=log((Plot12[Plot12$Organism=="CR2",]$Density)+1)-10
Plot12[Plot12$Organism=="CR4",]$Density=log((Plot12[Plot12$Organism=="CR4",]$Density)+1)-10
Plot12[Plot12$Organism=="R",]$Density=log((Plot12[Plot12$Organism=="R",]$Density)+1)

tiff('Population Dynamics 1.tiff', units="in", width=20, height=20, res=1000)
ggplot(subset(Plot12, Strain=="CR2"), aes(Day, Density)) +
  geom_point(aes(color=Organism, pch=Organism, size=Organism), alpha=0.6) +
  geom_smooth(aes(color=Organism, linetype=Organism, weight=Weight), method="loess", size=2, span=0.2, se=F) +
  geom_text(data=subset(Data15, Strain=="CR2" & Organism=="A"), mapping=aes(y=5.0, x=50-(50+0)*0.02, label=paste(format(CV,nsmall=2))), color="royalblue2", size=7, hjust=1) +
  geom_text(data=subset(Data15, Strain=="CR2" & Organism=="R"), mapping=aes(y=4.5, x=50-(50+0)*0.02, label=paste(format(CV,nsmall=2))), color="firebrick3", size=7, hjust=1) +
  ylab(expression(italic('B. calyciflorus')~'density'~'('*ln~individuals~mL^-1*')')) + xlab(expression('Time (days)')) +    
  theme(axis.text.y=element_text(face="plain", colour="black", size=24)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=24)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=24)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=24)) +
  theme(axis.title.y.left=element_text(face="plain", colour="black", angle=90)) +
  theme(axis.title.y.right=element_text(face="plain", colour="black", angle=90)) +
  scale_y_continuous(labels=function(y) sprintf("%.0f", y), breaks=seq(0,5,by=1), limits=c(0,5),
  sec.axis=sec_axis(~./5*5, expression(italic('C. reinhardtii')~'density'~'('*ln~cells~mL^-1*')'), labels=sprintf(seq(10,15,by=1), fmt="%.0f"), breaks=seq(0,5,by=1))) +
  scale_x_continuous(labels=function(y) sprintf("%.0f", y), breaks=seq(0,50,by=10), limits=c(0,50)) +
  theme(panel.background=element_blank(), panel.border=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.line=element_line(colour="black", linetype="solid", size=0.5)) +
  scale_color_manual(values=alpha(c("CR2"="royalblue2","CR4"="royalblue2","R"="firebrick3"),0.8), labels=c("CR2"=expression('Defended'~italic('C. reinhardtii')),"CR4"=expression('Undefended'~italic('C. reinhardtii')),"R"=expression(italic('B. calyciflorus')))) +
  scale_linetype_manual(values=c("CR2"="solid","CR4"="11","R"="solid"), labels=c("CR2"=expression('Defended'~italic('C. reinhardtii')),"CR4"=expression('Undefended'~italic('C. reinhardtii')),"R"=expression(italic('B. calyciflorus')))) +
  scale_shape_manual(values=c("CR2"=16,"CR4"=15,"R"=16), labels=c("CR2"=expression('Defended'~italic('C. reinhardtii')),"CR4"=expression('Undefended'~italic('C. reinhardtii')),"R"=expression(italic('B. calyciflorus')))) + 
  scale_size_manual(values=c("CR2"=2.0,"CR4"=1.8,"R"=2.0), labels=c("CR2"=expression('Defended'~italic('C. reinhardtii')),"CR4"=expression('Undefended'~italic('C. reinhardtii')),"R"=expression(italic('B. calyciflorus')))) +
  theme(panel.spacing.y=unit(1,'lines'), panel.spacing.x=unit(1,'lines')) + 
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color="black", size=0.9) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color="black", size=0.9) +
  facet_wrap(Bead~Trial, ncol=3, nrow=3) +
  theme(legend.position="none")
dev.off()

tiff('Population Dynamics 2.tiff', units="in", width=20, height=20, res=1000)
ggplot(subset(Plot12, Strain=="CR4"), aes(Day, Density)) +
  geom_point(aes(color=Organism, pch=Organism, size=Organism), alpha=0.6) +
  geom_smooth(aes(color=Organism, linetype=Organism, weight=Weight), method="loess", size=2, span=0.2, se=F) +
  geom_text(data=subset(Data15, Strain=="CR4" & Organism=="A"), mapping=aes(y=5.0, x=50-(50+0)*0.02, label=paste(format(CV,nsmall=2))), color="royalblue2", size=7, hjust=1) +
  geom_text(data=subset(Data15, Strain=="CR4" & Organism=="R"), mapping=aes(y=4.5, x=50-(50+0)*0.02, label=paste(format(CV,nsmall=2))), color="firebrick3", size=7, hjust=1) +
  ylab(expression(italic('B. calyciflorus')~'density'~'('*ln~individuals~mL^-1*')')) + xlab(expression('Time (days)')) +    
  theme(axis.text.y=element_text(face="plain", colour="black", size=24)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=24)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=24)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=24)) +
  theme(axis.title.y.left=element_text(face="plain", colour="black", angle=90)) +
  theme(axis.title.y.right=element_text(face="plain", colour="black", angle=90)) +
  scale_y_continuous(labels=function(y) sprintf("%.0f", y), breaks=seq(0,5,by=1), limits=c(0,5),
  sec.axis=sec_axis(~./5*5, expression(italic('C. reinhardtii')~'density'~'('*ln~cells~mL^-1*')'), labels=sprintf(seq(10,15,by=1), fmt="%.0f"), breaks=seq(0,5,by=1))) +
  scale_x_continuous(labels=function(y) sprintf("%.0f", y), breaks=seq(0,50,by=10), limits=c(0,50)) +
  theme(panel.background=element_blank(), panel.border=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.line=element_line(colour="black", linetype="solid", size=0.5)) +
  scale_color_manual(values=alpha(c("CR2"="royalblue2","CR4"="royalblue2","R"="firebrick3"),0.8), labels=c("CR2"=expression('Defended'~italic('C. reinhardtii')),"CR4"=expression('Undefended'~italic('C. reinhardtii')),"R"=expression(italic('B. calyciflorus')))) +
  scale_linetype_manual(values=c("CR2"="solid","CR4"="11","R"="solid"), labels=c("CR2"=expression('Defended'~italic('C. reinhardtii')),"CR4"=expression('Undefended'~italic('C. reinhardtii')),"R"=expression(italic('B. calyciflorus')))) +
  scale_shape_manual(values=c("CR2"=16,"CR4"=15,"R"=16), labels=c("CR2"=expression('Defended'~italic('C. reinhardtii')),"CR4"=expression('Undefended'~italic('C. reinhardtii')),"R"=expression(italic('B. calyciflorus')))) + 
  scale_size_manual(values=c("CR2"=2.0,"CR4"=1.8,"R"=2.0), labels=c("CR2"=expression('Defended'~italic('C. reinhardtii')),"CR4"=expression('Undefended'~italic('C. reinhardtii')),"R"=expression(italic('B. calyciflorus')))) +
  theme(panel.spacing.y=unit(1,'lines'), panel.spacing.x=unit(1,'lines')) + 
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color="black", size=0.9) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color="black", size=0.9) +
  facet_wrap(Bead~Trial, ncol=3, nrow=3) +
  theme(legend.position="none")
dev.off()

tiff('Population Dynamics 3.tiff', units="in", width=20, height=20, res=1000)
ggplot(subset(Plot12, Strain=="CR2-CR4"), aes(Day, Density)) +
  geom_point(aes(color=Organism, pch=Organism, size=Organism), alpha=0.6) +
  geom_smooth(aes(color=Organism, linetype=Organism, weight=Weight), method="loess", size=2, span=0.2, se=F) +
  geom_text(data=subset(Data15, Strain=="CR2-CR4" & Organism=="A"), mapping=aes(y=5.0, x=50-(50+0)*0.02, label=paste(format(CV,nsmall=2))), color="royalblue2", size=7, hjust=1) +
  geom_text(data=subset(Data15, Strain=="CR2-CR4" & Organism=="R"), mapping=aes(y=4.5, x=50-(50+0)*0.02, label=paste(format(CV,nsmall=2))), color="firebrick3", size=7, hjust=1) +
  ylab(expression(italic('B. calyciflorus')~'density'~'('*ln~individuals~mL^-1*')')) + xlab(expression('Time (days)')) +    
  theme(axis.text.y=element_text(face="plain", colour="black", size=24)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=24)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=24)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=24)) +
  theme(axis.title.y.left=element_text(face="plain", colour="black", angle=90)) +
  theme(axis.title.y.right=element_text(face="plain", colour="black", angle=90)) +
  scale_y_continuous(labels=function(y) sprintf("%.0f", y), breaks=seq(0,5,by=1), limits=c(0,5),
  sec.axis=sec_axis(~./5*5, expression(italic('C. reinhardtii')~'density'~'('*ln~cells~mL^-1*')'), labels=sprintf(seq(10,15,by=1), fmt="%.0f"), breaks=seq(0,5,by=1))) +
  scale_x_continuous(labels=function(y) sprintf("%.0f", y), breaks=seq(0,50,by=10), limits=c(0,50)) +
  theme(panel.background=element_blank(), panel.border=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.line=element_line(colour="black", linetype="solid", size=0.5)) +
  scale_color_manual(values=alpha(c("CR2"="royalblue2","CR4"="royalblue1","R"="firebrick3"),0.8), labels=c("CR2"=expression('Defended'~italic('C. reinhardtii')),"CR4"=expression('Undefended'~italic('C. reinhardtii')),"R"=expression(italic('B. calyciflorus')))) +
  scale_linetype_manual(values=c("CR2"="solid","CR4"="11","R"="solid"), labels=c("CR2"=expression('Defended'~italic('C. reinhardtii')),"CR4"=expression('Undefended'~italic('C. reinhardtii')),"R"=expression(italic('B. calyciflorus')))) +
  scale_shape_manual(values=c("CR2"=16,"CR4"=15,"R"=16), labels=c("CR2"=expression('Defended'~italic('C. reinhardtii')),"CR4"=expression('Undefended'~italic('C. reinhardtii')),"R"=expression(italic('B. calyciflorus')))) + 
  scale_size_manual(values=c("CR2"=2.0,"CR4"=1.8,"R"=2.0), labels=c("CR2"=expression('Defended'~italic('C. reinhardtii')),"CR4"=expression('Undefended'~italic('C. reinhardtii')),"R"=expression(italic('B. calyciflorus')))) +
  theme(panel.spacing.y=unit(1,'lines'), panel.spacing.x=unit(1,'lines')) + 
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color="black", size=0.9) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color="black", size=0.9) +
  facet_wrap(Bead~Trial, ncol=3, nrow=3) +
  theme(legend.position="none")
dev.off()


#####################################################
### Comparing ecology and evolution contributions ###
#####################################################

# Include prey mean defense and frequency
colnames(Data11)[6:7]=c("DensA","DensR")
Data11=cbind(Data11,Data7[,c(9,11)])
Data11[,c(8:9)]=round(Data11[,c(8:9)],2)

# Create a dataset
Data16=setDT(Data11)[, .(DensA=round(mean(DensA),0), DensR=round(mean(DensR),0)), by=list(Strain,Bead,Day,Trial,Alga,Defense,Freq)]
Data16=as.data.frame(Data16)

# Unmelt the dataset
Data17=subset(Data16, !Strain %in% c("CR2","CR4"))
Data17a=subset(Data17, Strain=="CR2-CR4" & Alga=="CR2")
Data17b=subset(Data17, Strain=="CR2-CR4" & Alga=="CR4")

# Rename column names
colnames(Data17a)[c(6:8)]=paste(colnames(Data17a)[c(6:8)], "DEF", sep="")
colnames(Data17b)[c(6:8)]=paste(colnames(Data17b)[c(6:8)], "UND", sep="")

# Create a dataset
Data18=cbind(Data17a[,c(1:4,6:8)],Data17b[,c(6:9)])

# Arrange the dataset
Data18=Data18 %>% arrange(factor(Strain, levels=c("CR2","CR4","CR2-CR4")), factor(Bead, levels=c("C","L","H")), factor(Trial, levels=c("1","2","3")), factor(Day, levels=seq(0,50,by=1)))

# Split the dataset
SplitData18=split(Data18, list(Data18$Strain,Data18$Trial,Data18$Bead))
SplitData18=SplitData18[sapply(SplitData18, function(x) dim(x)[1]) > 0]

# Subset important columns
SplitData18=lapply(SplitData18, "[", c("Day","DensR","DefenseDEF"))

# Instantaneous predator growth rate function
FuncGR=function(x) {
  FitGR=dynlm(formula=log(DensR+1)~L(Day), data=as.data.frame(x))
  OutGR=c(RateR=summary(FitGR)$coef[2,1])}

# Calculate instantaneous predator growth rates
ModelGR=function(y) {rollapply(y, width=2, FUN=FuncGR, by.column=F)}

# Calculate consecutive predator growth rates
OutGR=round(as.data.frame(do.call("cbind",lapply(SplitData18, ModelGR))),4)
GrowR=data.frame(OutGR=unlist(OutGR,use.names=F))[,1]

# Include predator growth rates
Data18$GrowR=c(t(cbind(rep(NA,9), matrix(GrowR, ncol=50, byrow=T))))
Data18$CombA=Data18$FreqDEF*(Data18$DensADEF+Data18$DensAUND)

# Split the dataset
SplitData18=split(Data18, list(Data18$Strain,Data18$Trial,Data18$Bead))
SplitData18=SplitData18[sapply(SplitData18, function(x) dim(x)[1]) > 0]

# Eco-evolutionary models
ModelEE1=function(x) {lm(x$GrowR~x$CombA)}
OutEE1=lapply(SplitData18, ModelEE1)

# Extract coefficients
CoefEE=lapply(OutEE1,function (x) {c(Coeff1=coef(summary(x))[1,2], Coeff2=coef(summary(x))[2,2])})
CoefEE=as.data.frame(do.call("rbind",CoefEE))
Data18$Coeff1=rep(CoefEE[,1], each=51)
Data18$Coeff2=rep(CoefEE[,2], each=51)

# Split the dataset
SplitData18=split(Data18, list(Data18$Strain,Data18$Trial,Data18$Bead))
SplitData18=SplitData18[sapply(SplitData18, function(x) dim(x)[1]) > 0]

# Extract combinations of names
Names=unique(Data18[,c("Strain","Bead","Trial","Day")])
Names=Names[order(Names$Bead,Names$Strain),]
Strain=Names$Strain; Bead=Names$Bead; Trial=Names$Trial

# Calculate the contributions
ModelEE2=function(x) {
  a=x[,13] + x[,14] * x[,6] * (x[,7] + x[,10])
  b=x[,13] + x[,14] * lead(x[,6]) * (x[,7] + x[,10])
  c=x[,13] + x[,14] * x[,6] * (lead(x[,7]) + lead(x[,10]))
  d=x[,13] + x[,14] * lead(x[,6]) * (lead(x[,7]) + lead(x[,10]))
  Eco=((b-a)+(d-c))/2
  Evo=((c-a)+(d-b))/2
  Coef=data.frame(Eco=Eco, Evo=Evo, EcoEvo=Eco/Evo)
}
OutEE2=lapply(SplitData18, ModelEE2)

RateEE=round(as.data.frame(do.call("rbind",OutEE2)),4)
RateEE=cbind(Strain=Data18[,1],Bead=Data18[,2],Day=Data18[,3],Trial=Data18[,4],RateEE)
rownames(RateEE)=c()

# Include the contributions
Data18$Eco=RateEE[,5]
Data18$Evo=RateEE[,6]
Data18$EcoEvo=RateEE[,7]

# Correct the contributions
Data18$EcoEvo=log(abs(Data18$EcoEvo))

# Replace infinite values
Data18$Eco=ifelse(is.infinite(Data18$Eco), 0, Data18$Eco)
Data18$Evo=ifelse(is.infinite(Data18$Evo), 0, Data18$Evo)
Data18$EcoEvo=ifelse(is.infinite(Data18$EcoEvo), 0, Data18$EcoEvo)

# Calculate the mean contributions
Data18=setDT(na.omit(Data18))[, MeanEcoEvo := mean(EcoEvo), by=list(Strain,Bead,Day)]
Data18=setDT(na.omit(Data18))[, HarmEcoEvo := mean(EcoEvo), by=list(Strain,Bead)]
Data18=as.data.frame(Data18)


#############################################
### Calculating coefficients of variation ###
#############################################

# Scale the dataset
Data19=na.omit(Data18)
Data19$ScaledEcoEvo=Data19$MeanEcoEvo+abs(min(Data19$MeanEcoEvo))

# Calculate the coefficients
Data20=setDT(Data19)[, .(CV=round((sd(ScaledEcoEvo)/mean(ScaledEcoEvo))*100,2), EcoEvo=round(HarmEcoEvo,2)), by=list(Strain,Bead)]
Data20=Data20[order(Data20$Strain,Data20$Bead),]
Data20=as.data.frame(Data20)


####################################################
### Plotting ecology and evolution contributions ###
####################################################

# Arrange the dataset
Data18=Data18 %>% arrange(factor(Strain, levels=c("CR2-CR4")), factor(Bead, levels=c("C","L","H")), factor(Trial, levels=c("1","2","3")), factor(Day, levels=seq(0,50,by=1)))
Data18$Strain=factor(Data18$Strain, levels=c("CR2-CR4"))
Data20$Strain=factor(Data20$Strain, levels=c("CR2-CR4"))
Data18$Bead=factor(Data18$Bead, levels=c("C","L","H"))
Data20$Bead=factor(Data20$Bead, levels=c("C","L","H"))
Data18$Trial=factor(Data18$Trial, levels=c("1","2","3"))

# Extract the dataset
write.table(Data18, file="~/Activité Professionnelle/LIMNO 2019-2023/Experiments/Population Dynamics/Data_EXP2.txt", sep="\t", row.names=F)

tiff('Contribution Ecology Evolution.tiff', units="in", width=14, height=5, res=1000)
ggplot(Data18, aes(Day, EcoEvo)) + coord_cartesian(clip="off") +
  geom_hline(aes(yintercept=0), color="black", linetype="11", size=1) + 
  geom_point(aes(Day, MeanEcoEvo), color="navyblue", size=2, pch=16) +
  geom_smooth(aes(Day, MeanEcoEvo), method="loess", color="navyblue", linetype="solid", size=1.5, span=0.2, se=F) +
  geom_text(data=Data20, mapping=aes(y=3.0, x=0+(50+0)*0.02, label=paste(format(EcoEvo,nsmall=2))), color="navyblue", size=5, hjust=0) +
  ylab(expression('Ratio of contributions')) + xlab(expression('Time (days)')) +    
  theme(axis.text.y.right=element_blank()) + 
  theme(axis.text.y.left=element_text(face="plain", colour="black", size=18)) +
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y.right=element_text(face="plain", colour="navyblue", size=14, angle=90)) +
  theme(axis.title.y.left=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(y) sprintf("%.1f", y), breaks=seq(-3.0,3.0,by=1.5), limits=c(-3.0,3.0), 
  sec.axis=sec_axis(~./1*1, expression('Evolution > Ecology'~~~~~~~~~'Ecology > Evolution'))) +
  scale_x_continuous(labels=function(y) sprintf("%.0f", y), breaks=seq(0,50,by=10), limits=c(0,50)) +
  theme(panel.background=element_blank(), panel.border=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(axis.line.y.right=element_blank()) +
  theme(axis.line.y.left=element_line(colour="black", linetype="solid", size=0.5)) +
  theme(axis.line.x=element_line(colour="black", linetype="solid", size=0.5)) +
  theme(axis.ticks.y.right=element_blank()) +
  theme(axis.ticks.y.left=element_line(colour="black", linetype="solid", size=0.5)) +
  scale_color_manual(values=alpha(c("1"="navyblue","2"="navyblue","3"="navyblue"),0.3)) +
  theme(panel.spacing.y=unit(1,'lines'), panel.spacing.x=unit(1,'lines')) + 
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=0.5) +
  facet_wrap(~Bead, ncol=3, nrow=1) +
  theme(legend.position="none")
dev.off()


#########################
### Statistical tests ###
#########################

### Rotifer densities ###

# Rename column names
Data21a=subset(Data11, Strain %in% c("CR2","CR4","CR2-CR4") & Alga=="CR2")
Data21b=subset(Data11, Strain %in% c("CR2","CR4","CR2-CR4") & Alga=="CR4")
colnames(Data21a)[c(6,8)]=paste(colnames(Data21a)[c(6,8)], "DEF", sep="")
colnames(Data21b)[c(6,8)]=paste(colnames(Data21b)[c(6,8)], "UND", sep="")

# Unmelt the dataset
Data21=merge(Data21a, Data21b, all=T)
Data18$Day=as.factor(Data18$Day)
Data21$Day=as.factor(Data21$Day)
Data21[is.na(Data21)]=0

# Subset the dataset
Data22=subset(Data21, Strain=="CR2-CR4" & Alga=="CR4")

# Normality of values
shapiro.test(Data21$DensR)

# Generalized linear mixed effect regression
Model1=glmer(log(DensR+1)~Strain*Bead + Day + (1|Trial), family=gaussian(link="identity"), data=Data21)
hist(residuals(Model1))
Anova(Model1)

# Posthoc pairwise regression
emmeans(Model1, pairwise~Strain)
emmeans(Model1, pairwise~Bead)


### Alga frequencies ###

# Normality of values
shapiro.test(Data22$FreqUND)

# Generalized linear mixed effect regression
Model2=glmer(invlogit(FreqUND)~Bead*Day + (1|Trial), family=Gamma(link="logit"), data=Data22)
hist(residuals(Model2))
Anova(Model2)

# Posthoc pairwise regression
emmeans(Model2, pairwise~Bead)


### Relative contributions ###

# Normality of values
shapiro.test(Data18$EcoEvo)

# Generalized linear mixed effect regression
Model3=glmer(EcoEvo~Bead*Day + (1|Trial), family=gaussian(link="identity"), data=Data18)
hist(residuals(Model3))
Anova(Model3)


### Coefficients of variation ###

# Normality of values
shapiro.test(subset(Data15, Organism=="A")$CV)
shapiro.test(subset(Data15, Organism=="R")$CV)

# Generalized linear regression
Model4=lm(CV~Strain*Bead, data=subset(Data15, Organism=="A"))
shapiro.test(residuals(Model4))
Anova(Model4)

# Generalized linear regression
Model5=lm(CV~Strain*Bead, data=subset(Data15, Organism=="R"))
shapiro.test(residuals(Model5))
Anova(Model5)

# Calculate the mean coefficients
Data23=setDT(Data15)[, .(CV=round(mean(CV),2), CVLSD=round(mean(CV)-sd(CV),2), CVUSD=round(mean(CV)+sd(CV),2)), by=list(Strain,Bead,Organism)]
Data23=as.data.frame(Data23)
