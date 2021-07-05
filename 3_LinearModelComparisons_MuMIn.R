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

#Create empty list to store output of each for loop iteration
LM1_output<-list() #Affect of Sex Class
LM2_output<-list() #Affect of Age Class

options(na.action = "na.fail")
morpho1<-na.omit(morpho) #we have to remove all the NAs to run the dredge()
response_var<-colnames(morpho1)[c(19,20,21,22,23,24,25)]

for(i in 1:length(response_var)){
  print(paste0("Character",i," -- ",response_var[i]))#Report character
  character_of_interest<-morpho1[,response_var[i]] #Extract character of interest (coi)
  names(character_of_interest)<-morpho1$tip.label #Add names to vector of coi 
  birds.fullmodel_Sex<-lm(morpho1[,response_var[i]]~morpho1$Sex, data = morpho1)
  birds.fullmodel_Age<-lm(morpho1[,response_var[i]]~morpho1$Age, data = morpho1)
  LM1_output[[i]]<-dredge(birds.fullmodel_Sex, rank = "AIC")
  LM2_output[[i]]<-dredge(birds.fullmodel_Age, rank = "AIC")
}

names(LM1_output)<-response_var 
names(LM2_output)<-response_var

### Write out output ###
sink("Dredge_LMOutput_Sex.txt")
for(i in 1:length(LM1_output)){
  cat(paste0("Character ",i," -- ",response_var[i]))
  cat("\n\n")
  print(LM1_output[[i]])
  cat("#############")
  cat("\n\n")
}
sink()


sink("Dredge_LMOutput_Age.txt")
for(i in 1:length(LM2_output)){
  cat(paste0("Character ",i," -- ",response_var[i]))
  cat("\n\n")
  print(LM2_output[[i]])
  cat("#############")
  cat("\n\n")
}
sink()

#What we're trying to do is to test whether tarsus length is the best approximation of body mass.

