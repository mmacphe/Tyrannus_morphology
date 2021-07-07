rm(list=ls())

### Set source directory to the folder this file came from within RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

morpho<-read.csv('./Output Files/Tyrannus morphology data.csv')

#What we're trying to do is to see if, within each taxon:
#1) morphologies differ between the sexes, and 
#2) morphologies differ between age classes

### Standardize notation for Sex column
names(table(morpho$Sex))
morpho<-morpho[!morpho$Sex=="Unknown",] #removes unknown sex individuals
names(table(morpho$Sex)) #should now only have "Female" and "Male"

### Standardize notation for Age column
names(table(morpho$Age)) #to see which assignments are in the Age column
morpho$Age[morpho$Age=="SecondYear"]<-"Adult"
morpho$Age[morpho$Age=="Fledgling"]<-"Juvenile"
morpho$Age[morpho$Age=="HatchYear"]<-"Juvenile"
morpho$Age[morpho$Age=="Nestling"]<-"Juvenile"
morpho<-morpho[!morpho$Age=="Unknown",] #removes unknown age individuals
names(table(morpho$Age)) #should now only have "Adult" and "Juvenile"

### For Use in later steps: Make separate male and female voucher tables
morpho_Females<-morpho[!morpho$Sex=="Male",]
write.csv(morpho_Females, file="Voucher_Table_Females.csv")

morpho_Males<-morpho[!morpho$Sex=="Female",]
write.csv(morpho_Males, file="Voucher_Table_Males.csv")

###############################################
#####/\/\/\/\/\/\/TRY MCMCGLM/\/\/\/\/\/\/\####
###############################################
library(MCMCglmm)
library(caper) #for comparative.data() but doesn't allow duplicate 'row names'(i.e., individuals from the same taxon)

## Installing the package mulTree
if(!require(devtools)) install.packages("devtools")
library(devtools)
install_github("TGuillerme/mulTree", ref = "release")
library(mulTree)

phy<-read.tree('./Output Files/Tyrannus_phylogeny.tre')

options(na.action = "na.fail")
morpho1<-subset(morpho, select=c(BL.Average,BW.Average,BD.Average,Tarsus.Average,Kipp.s.Distance,WC.Average,Tail,Age,Sex,Kipp.s.Index,tip.label)) #remove the Mass column for the first comparisons bc it is unnecessary for the next step and has NAs that will affect the outcome
colnames(morpho1) #look at the column names to pick out the right ones for response_var
morpho1<-na.omit(morpho1) #we have to remove all the NAs to run the dredge()
#response_var<-colnames(morpho1)[c(1,2,3,4,5,6,7,10)]

comp_data<-as.mulTree(data=morpho1, tree=phy, taxa="tip.label")

my_formula_Sex_BL<-BL.Average~Sex
my_formula_Sex_BW<-BW.Average~Sex
my_formula_Sex_BD<-BD.Average~Sex
my_formula_Sex_Tarsus<-Tarsus.Average~Sex
my_formula_Sex_KippsD<-Kipp.s.Distance~Sex
my_formula_Sex_WC<-WC.Average~Sex
my_formula_Sex_Tail<-Tail~Sex
my_formula_Sex_KI<-Kipp.s.Index~Sex

my_parameters<-c(100000,10,1000) #The MCMC parameters (generations, sampling, burnin). Used default settings from https://github.com/TGuillerme/mulTree/blob/master/doc/mulTree-manual.pdf

my_priors<-list(R=list(V=1/2,nu=0.002),G=list(G1=list(V=1/2,nu=0.002)))

## Run the MCMCglmm on the comp_data object
LM_Sex_BL<-mulTree(mulTree.data=comp_data, formula=my_formula_Sex_BL,priors=my_priors,parameters=my_parameters,output="Sex_BL",ESS=50,chains=1) #chains=1 bc only 1 tree and not multiple trees
LM_Sex_BW<-mulTree(mulTree.data=comp_data, formula=my_formula_Sex_BW,priors=my_priors,parameters=my_parameters,output="Sex_BW",ESS=50,chains=1) #chains=1 bc only 1 tree and not multiple trees
LM_Sex_BD<-mulTree(mulTree.data=comp_data, formula=my_formula_Sex_BD,priors=my_priors,parameters=my_parameters,output="Sex_BD",ESS=50,chains=1) #chains=1 bc only 1 tree and not multiple trees
LM_Sex_Tarsus<-mulTree(mulTree.data=comp_data, formula=my_formula_Sex_Tarsus,priors=my_priors,parameters=my_parameters,output="Sex_Tarsus",ESS=50,chains=1) #chains=1 bc only 1 tree and not multiple trees
LM_Sex_KippsD<-mulTree(mulTree.data=comp_data, formula=my_formula_Sex_KippsD,priors=my_priors,parameters=my_parameters,output="Sex_KippsD",ESS=50,chains=1) #chains=1 bc only 1 tree and not multiple trees
LM_Sex_WC<-mulTree(mulTree.data=comp_data, formula=my_formula_Sex_WC,priors=my_priors,parameters=my_parameters,output="Sex_WC",ESS=50,chains=1) #chains=1 bc only 1 tree and not multiple trees
LM_Sex_Tail<-mulTree(mulTree.data=comp_data, formula=my_formula_Sex_Tail,priors=my_priors,parameters=my_parameters,output="Sex_Tail",ESS=50,chains=1) #chains=1 bc only 1 tree and not multiple trees
LM_Sex_KI<-mulTree(mulTree.data=comp_data, formula=my_formula_Sex_KI,priors=my_priors,parameters=my_parameters,output="Sex_KI",ESS=50,chains=1) #chains=1 bc only 1 tree and not multiple trees

#Read each single, specific model
Sex_BL<-read.mulTree("Sex_BL",model=TRUE)
#summarise results
summary_Sex_BL<-summary(Sex_BL)

Sex_BW<-read.mulTree("Sex_BW",model=TRUE)
#summarise results
summary_Sex_BW<-summary(Sex_BW)

Sex_BD<-read.mulTree("Sex_BD",model=TRUE)
#summarise results
summary_Sex_BD<-summary(Sex_BD)

Sex_Tarsus<-read.mulTree("Sex_Tarsus",model=TRUE)
#summarise results
summary_Sex_Tarsus<-summary(Sex_Tarsus)

Sex_KippsD<-read.mulTree("Sex_KippsD",model=TRUE)
#summarise results
summary_Sex_KippsD<-summary(Sex_KippsD)

Sex_WC<-read.mulTree("Sex_WC",model=TRUE)
#summarise results
summary_Sex_WC<-summary(Sex_WC)

Sex_Tail<-read.mulTree("Sex_Tail",model=TRUE)
#summarise results
summary_Sex_Tail<-summary(Sex_Tail)

Sex_KI<-read.mulTree("Sex_KI",model=TRUE)
#summarise results
summary_Sex_KI<-summary(Sex_KI)

LM1_output<-list(summary_Sex_BL, summary_Sex_BW, summary_Sex_BD, summary_Sex_Tarsus,
                 summary_Sex_KippsD, summary_Sex_WC, summary_Sex_Tail, summary_Sex_KI) #Affect of Sex Class

### Write out output ###
sink("MCMCglmm_LMOutput_Sex.txt")
for(i in 1:length(LM1_output)){
  print(LM1_output[[i]])
  cat("########################################################")
}
sink()

### Repeat for Age classes
my_formula_Age_BL<-BL.Average~Age
my_formula_Age_BW<-BW.Average~Age
my_formula_Age_BD<-BD.Average~Age
my_formula_Age_Tarsus<-Tarsus.Average~Age
my_formula_Age_KippsD<-Kipp.s.Distance~Age
my_formula_Age_WC<-WC.Average~Age
my_formula_Age_Tail<-Tail~Age
my_formula_Age_KI<-Kipp.s.Index~Age

LM_Age_BL<-mulTree(mulTree.data=comp_data, formula=my_formula_Age_BL,priors=my_priors,parameters=my_parameters,output="Age_BL",ESS=50,chains=1) #chains=1 bc only 1 tree and not multiple trees
LM_Age_BW<-mulTree(mulTree.data=comp_data, formula=my_formula_Age_BW,priors=my_priors,parameters=my_parameters,output="Age_BW",ESS=50,chains=1) #chains=1 bc only 1 tree and not multiple trees
LM_Age_BD<-mulTree(mulTree.data=comp_data, formula=my_formula_Age_BD,priors=my_priors,parameters=my_parameters,output="Age_BD",ESS=50,chains=1) #chains=1 bc only 1 tree and not multiple trees
LM_Age_Tarsus<-mulTree(mulTree.data=comp_data, formula=my_formula_Age_Tarsus,priors=my_priors,parameters=my_parameters,output="Age_Tarsus",ESS=50,chains=1) #chains=1 bc only 1 tree and not multiple trees
LM_Age_KippsD<-mulTree(mulTree.data=comp_data, formula=my_formula_Age_KippsD,priors=my_priors,parameters=my_parameters,output="Age_KippsD",ESS=50,chains=1) #chains=1 bc only 1 tree and not multiple trees
LM_Age_WC<-mulTree(mulTree.data=comp_data, formula=my_formula_Age_WC,priors=my_priors,parameters=my_parameters,output="Age_WC",ESS=50,chains=1) #chains=1 bc only 1 tree and not multiple trees
LM_Age_Tail<-mulTree(mulTree.data=comp_data, formula=my_formula_Age_Tail,priors=my_priors,parameters=my_parameters,output="Age_Tail",ESS=50,chains=1) #chains=1 bc only 1 tree and not multiple trees
LM_Age_KI<-mulTree(mulTree.data=comp_data, formula=my_formula_Age_KI,priors=my_priors,parameters=my_parameters,output="Age_KI",ESS=50,chains=1) #chains=1 bc only 1 tree and not multiple trees

Age_BL<-read.mulTree("Age_BL",model=TRUE)
#summarise results
summary_Age_BL<-summary(Age_BL)

Age_BW<-read.mulTree("Age_BW",model=TRUE)
#summarise results
summary_Age_BW<-summary(Age_BW)

Age_BD<-read.mulTree("Age_BD",model=TRUE)
#summarise results
summary_Age_BD<-summary(Age_BD)

Age_Tarsus<-read.mulTree("Age_Tarsus",model=TRUE)
#summarise results
summary_Age_Tarsus<-summary(Age_Tarsus)

Age_KippsD<-read.mulTree("Age_KippsD",model=TRUE)
#summarise results
summary_Age_KippsD<-summary(Age_KippsD)

Age_WC<-read.mulTree("Age_WC",model=TRUE)
#summarise results
summary_Age_WC<-summary(Age_WC)

Age_Tail<-read.mulTree("Age_Tail",model=TRUE)
#summarise results
summary_Age_Tail<-summary(Age_Tail)

Age_KI<-read.mulTree("Age_KI",model=TRUE)
#summarise results
summary_Age_KI<-summary(Age_KI)

LM2_output<-list(summary_Age_BL, summary_Age_BW, summary_Age_BD, summary_Age_Tarsus,
                 summary_Age_KippsD, summary_Age_WC, summary_Age_Tail, summary_Age_KI) #Affect of Age Class

### Write out output ###
sink("MCMCglmm_LMOutput_Age.txt")
for(i in 1:length(LM2_output)){
  print(LM2_output[[i]])
  cat("########################################################")
}
sink()

#################################################
#####/\/\/\/\/\/\/TRIED MCMCglmm/\/\/\/\/\/\/\####
#################################################

##Create empty list to store output of each for loop iteration
#LM1_output<-list() #Affect of Sex Class
#LM2_output<-list() #Affect of Age Class

#options(na.action = "na.fail")
#morpho1<-subset(morpho, select=c(BL.Average,BW.Average,BD.Average,Tarsus.Average,Kipp.s.Distance,WC.Average,Tail,Age,Sex,Kipp.s.Index,tip.label)) #remove the Mass column for the first comparisons bc it is unnecessary for the next step and has NAs that will affect the outcome
#colnames(morpho1) #look at the column names to pick out the right ones for response_var
#morpho1<-na.omit(morpho1) #we have to remove all the NAs to run the dredge()
#response_var<-colnames(morpho1)[c(1,2,3,4,5,6,7)]

#for(i in 1:length(response_var)){
#  print(paste0("Character",i," -- ",response_var[i]))#Report character
#  character_of_interest<-morpho1[,response_var[i]] #Extract character of interest (coi)
#  names(character_of_interest)<-morpho1$tip.label #Add names to vector of coi 
#  birds.fullmodel_Sex<-lm(morpho1[,response_var[i]]~morpho1$Sex, data = morpho1)
#  birds.fullmodel_Age<-lm(morpho1[,response_var[i]]~morpho1$Age, data = morpho1)
#  LM1_output[[i]]<-dredge(birds.fullmodel_Sex, rank = "AIC")
#  LM2_output[[i]]<-dredge(birds.fullmodel_Age, rank = "AIC")
#}

#names(LM1_output)<-response_var 
#names(LM2_output)<-response_var

#### Write out output ###
#sink("Dredge_LMOutput_Sex.txt")
#for(i in 1:length(LM1_output)){
#  cat(paste0("Character ",i," -- ",response_var[i]))
#  cat("\n\n")
#  print(LM1_output[[i]])
#  cat("#############")
#  cat("\n\n")
#}
#sink()


#sink("Dredge_LMOutput_Age.txt")
#for(i in 1:length(LM2_output)){
#  cat(paste0("Character ",i," -- ",response_var[i]))
#  cat("\n\n")
#  print(LM2_output[[i]])
#  cat("#############")
#  cat("\n\n")
#}
#sink()


### Make morphology datasets that do not include juveniles
### Remove juvenile individuals
morpho_Adults<-morpho[!morpho$Age=="Juvenile",]

### Split the morphology data sets up into categories by tip.label
morpho_Adults$tip.label
names(table(morpho_Adults$tip.label))

morpho_split<-split(morpho_Adults,morpho_Adults$tip.label)
#names(morpho_split)<-gsub(" ","_",names(morpho_split))
names(morpho_split) 

### Read in phylogeny 
phy<-read.tree('./Output Files/Tyrannus_phylogeny.tre')

phy$tip.label[!phy$tip.label %in% names(morpho_split)] #Check that this returns "character(0)"
names(morpho_split)[! names(morpho_split) %in% phy$tip.label] #Check that this also returns "character(0)"


### Create summary output df for data frame
morpho_trim<-lapply(morpho_split,function(x) x[colnames(x) %in% c("BL.Average","BW.Average","BD.Average", "Kipp.s.Distance", "Kipp.s.Index", "WC.Average", "Tail", "Tarsus.Average","Mass")])

### Get average values for each morphometric
otu_avg<-list()
otu_sd<-list()
otu_cv<-list()
otu_n<-list()

for(i in 1:length(morpho_trim)){
  otu_avg[[i]]<-apply(morpho_trim[[i]],2,function(x) mean(na.omit(as.numeric(x))))
  otu_sd[[i]]<-apply(morpho_trim[[i]],2,function(x) sd(na.omit(as.numeric(x))))
  otu_cv[[i]]<-otu_sd[[i]]/otu_avg[[i]]
  otu_n[[i]]<-apply(morpho_trim[[i]],2,function(x) length(which(!is.na(x))))
}

otu_avg_df<-do.call(rbind,otu_avg)
otu_sd_df<-do.call(rbind,otu_sd)
otu_cv_df<-do.call(rbind,otu_cv)
otu_n_df<-do.call(rbind,otu_n)

rownames(otu_avg_df)<-names(morpho_split)

### Write out .csv file with mean and SD for each OTU for table in manuscript
morphology_summary<-data.frame(matrix(paste(round(otu_avg_df,2)," ± ",round(otu_sd_df,2),"\n(",otu_n_df,")",sep=""),ncol=ncol(otu_avg_df)))
colnames(morphology_summary)<-colnames(otu_avg_df)
rownames(morphology_summary)<-names(morpho_split)
write.csv(morphology_summary,file="Adults_morphology_summary_table.csv")

### Write out .csv file just with the mean for each OTU to be used in further analyses
Tyrannus_OTU_averages<-data.frame(matrix(paste(round(otu_avg_df,2)),ncol=ncol(otu_avg_df)))
colnames(Tyrannus_OTU_averages)<-colnames(otu_avg_df)
rownames(Tyrannus_OTU_averages)<-names(morpho_split)

### Add migratory strategy as a column in the Tyrannus OTU averages .csv ###
migratory_data<-read.csv("Tyrannus_subspecies_MigrationStrategies.csv",row.names=1)
migratory_data[rownames(Tyrannus_OTU_averages),]
Tyrannus_data<-merge(Tyrannus_OTU_averages, migratory_data, by=0)
write.csv(Tyrannus_data, file="Tyrannus_Adults_data.csv")

### Write out .csv file with coefficient of variation for each OTU to be used in further analyses
cv_summary<-otu_cv_df
rownames(cv_summary)<-names(morpho_split)
cv_summary<-as.data.frame(cv_summary)
migratory_data[rownames(cv_summary),]
CV_data<-merge(cv_summary, migratory_data, by=0)
write.csv(CV_data, file = "cv_summary_Adults_table.csv")

### What we're trying to do now is to test whether tarsus length is the best approximation of body mass.
morpho2<-subset(morpho, select=c(Mass,BL.Average,BW.Average,BD.Average,Tarsus.Average,Kipp.s.Distance,WC.Average,Tail,Kipp.s.Index,tip.label)) #remove the Mass column for the first comparisons bc it is unnecessary for the next step and has NAs that will affect the outcome
colnames(morpho2) #check the column names to make sure they're all there
morpho2<-na.omit(morpho2) #we have to remove all the NAs to run the dredge(); leaves 190 observations
#birds.fullmodel_Tarsus<-lm(Mass~(BL.Average+BW.Average+BD.Average+Tarsus.Average+Kipp.s.Distance+WC.Average+Tail+Kipp.s.Index)*tip.label, data = morpho2)

#LM3_output<-dredge(birds.fullmodel_Tarsus, rank = "AIC")
#head(LM3_output,25) #the top 25 models
#nrow(LM3_output[LM3_output$delta<7, ]) #models with delta AIC < 7 (returns [1] 49)
#w<-Weights(LM3_output)
#w

##Refit best linear model
#bestmodel<-get.models(LM3_output, 1)[[1]]
#z<-lm(bestmodel,data=morpho2)
#summary(z) #the best model has tarsus length as the most important factor (p=0.000448)

#Using mulTree
comp_data<-as.mulTree(data=morpho2, tree=phy, taxa="tip.label")

my_formula<-Mass~BL.Average+BW.Average+BD.Average+Tarsus.Average+Kipp.s.Distance+WC.Average+Tail

my_parameters<-c(100000,10,1000) #The MCMC parameters (generations, sampling, burnin). Used default settings from https://github.com/TGuillerme/mulTree/blob/master/doc/mulTree-manual.pdf

#my_priors<-list(R=list(V=1/2,nu=0.002),G=list(G1=list(V=1/2,nu=0.002))) #tutorial values from mulTree
my_priors<-list(R=list(V=1,nu=0),G=list(G1=list(V=1,nu=0))) #defaults from MCMCglmm that are not automatically specified in mulTree. Defaults listed at https://cran.r-project.org/web/packages/MCMCglmm/MCMCglmm.pdf

## Run the MCMCglmm on the comp_data object
LM<-mulTree(mulTree.data=comp_data, formula=my_formula, priors=my_priors, parameters=my_parameters,output="Tarsus_Mass_MCMC",chains=1) #chains=1 bc only 1 tree and not multiple trees. ESS = effective sample size; default is 1000.

Tarsus_Mass<-read.mulTree("Tarsus_Mass_MCMC",model=TRUE)
#summarise results
summary_Tarsus_Mass<-summary(Tarsus_Mass)
summary_Tarsus_Mass

### Write out output ###
sink("MCMCglmm_LMOutput_Tarsus_Mass.txt")
print(summary_Tarsus_Mass)
sink()
