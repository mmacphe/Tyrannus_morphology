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

### For Use in later steps: Make separate male and female voucher tables
morpho_Females<-morpho[!morpho$Sex=="Male",]
write.csv(morpho_Females, file="Voucher_Table_Females.csv")

morpho_Males<-morpho[!morpho$Sex=="Female",]
write.csv(morpho_Males, file="Voucher_Table_Males.csv")

#response_var<-colnames(morpho)[c(18,19,20,21,22,23,24)]
response_var<-colnames(morpho)[c(19,20,21,22,23,24,25)]

#Create empty list to store output of each for loop iteration
LM1_output<-list()
LM2_output<-list()

for(i in 1:length(response_var)){
  print(paste0("Character",i," -- ",response_var[i]))#Report character
  character_of_interest<-morpho[,response_var[i]] #Extract character of interest (coi)
  names(character_of_interest)<-morpho$tip.label #Add names to vector of coi
  LM1_output[[i]]<-AIC(lm(morpho[,response_var[i]]~morpho$tip.label, data=morpho)) #Calculate AIC scores
  LM2_output[[i]]<-AIC(lm(morpho[,response_var[i]]~morpho$tip.label*morpho$Sex, data=morpho)) #Calculate AIC scores while accounting for sex class
}

names(LM1_output)<-response_var
names(LM2_output)<-response_var 

AIC_Diff<-list()

for(i in 1:length(response_var)){
  print(paste0("Character",i," -- ",response_var[[i]]))
  AIC_Diff[[i]]<-LM1_output[[i]]-LM2_output[[i]]
}

names(AIC_Diff)<-response_var

### Write out output ###
sink("Sex_Influence.txt")
for(i in 1:length(AIC_Diff)){
  cat(paste0("Character ",i," -- ",response_var[i]))
  cat("\n\n")
  print(AIC_Diff[[i]])
  cat("#############")
  cat("\n\n")
}
sink()
#Note: For all values >2 AIC units, one model is better than the other.
#Here, we subtracted LM2 (accounted for sex class) from LM1 (did not account for sex class).
#All negative values indicate that the model was improved by the accounting of sex class, and positive values that it was not.

###/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\###
###############################################################
###/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\###

### Empty global environment to assess the role of Age
#Note: we do this because some individuals may have unknown sexes (which were previously removed, above) but known ages.
rm(list=ls())

morpho<-read.csv('./Output Files/Tyrannus morphology data.csv')

### Remove any individuals with unresolved taxonomic identities
#morpho<-morpho[-grep("unknown",morpho$subspecies),]

### Classify each taxon using trinomial naming
#morpho$subspecies<-trimws(morpho$subspecies)
#morpho$trinomial<-paste(morpho$Species,morpho$subspecies)
#morpho$trinomial<-trimws(morpho$trinomial)
#names(table(morpho$trinomial)) #see that unknown subspecies were removed

### Standardize notation for Age column
names(table(morpho$Age)) #to see which assignments are in the Age column
morpho$Age[morpho$Age=="SecondYear"]<-"Adult"
morpho$Age[morpho$Age=="Fledgling"]<-"Juvenile"
morpho$Age[morpho$Age=="HatchYear"]<-"Juvenile"
morpho$Age[morpho$Age=="Nestling"]<-"Juvenile"
morpho<-morpho[!morpho$Age=="Unknown",] #removes unknown age individuals
names(table(morpho$Age)) #should now only have "Adult" and "Juvenile"

response_var<-colnames(morpho)[c(19,20,21,22,23,24,25)]

#Create empty list to store output of each for loop iteration
LM1_output<-list()
LM2_output<-list()

for(i in 1:length(response_var)){
  print(paste0("Character",i," -- ",response_var[i]))#Report character
  character_of_interest<-morpho[,response_var[i]] #Extract character of interest (coi)
  names(character_of_interest)<-morpho$tip.label #Add names to vector of coi
  LM1_output[[i]]<-AIC(lm(morpho[,response_var[i]]~morpho$tip.label, data=morpho)) #Calculate AIC scores
  LM2_output[[i]]<-AIC(lm(morpho[,response_var[i]]~morpho$tip.label*morpho$Age, data=morpho))
}

names(LM1_output)<-response_var
names(LM2_output)<-response_var 

AIC_Diff<-list()

for(i in 1:length(response_var)){
  print(paste0("Character",i," -- ",response_var[[i]]))
  AIC_Diff[[i]]<-LM1_output[[i]]-LM2_output[[i]]
}

names(AIC_Diff)<-response_var

### Write out output ###
sink("Age_Influence.txt")
for(i in 1:length(AIC_Diff)){
  cat(paste0("Character ",i," -- ",response_var[i]))
  cat("\n\n")
  print(AIC_Diff[[i]])
  cat("#############")
  cat("\n\n")
}
sink()
#Note: For all values >2 AIC units, one model is better than the other.
#Here, we subtracted LM2 (accounted for age class) from LM1 (did not account for age class).
#All negative values indicate that the model was improved by the accounting of age class, and positive values that it was not.