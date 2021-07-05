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
morpho1<-subset(morpho, select=c(BL.Average,BW.Average,BD.Average,Tarsus.Average,Kipp.s.Distance,WC.Average,Tail,Age,Sex,Kipp.s.Index,tip.label)) #remove the Mass column for the first comparisons bc it is unnecessary for the next step and has NAs that will affect the outcome
colnames(morpho1) #look at the column names to pick out the right ones for response_var
morpho1<-na.omit(morpho1) #we have to remove all the NAs to run the dredge()
response_var<-colnames(morpho1)[c(1,2,3,4,5,6,7)]

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

# What we're trying to do is to test whether tarsus length is the best approximation of body mass.
morpho2<-subset(morpho, select=c(Mass,BL.Average,BW.Average,BD.Average,Tarsus.Average,Kipp.s.Distance,WC.Average,Tail,Kipp.s.Index,tip.label)) #remove the Mass column for the first comparisons bc it is unnecessary for the next step and has NAs that will affect the outcome
colnames(morpho2) #check the column names to make sure they're all there
morpho2<-na.omit(morpho2) #we have to remove all the NAs to run the dredge(); leaves 190 observations
birds.fullmodel_Tarsus<-lm(Mass~(BL.Average+BW.Average+BD.Average+Tarsus.Average+Kipp.s.Distance+WC.Average+Tail+Kipp.s.Index)*tip.label, data = morpho2)

LM3_output<-dredge(birds.fullmodel_Tarsus, rank = "AIC")
head(LM3_output,25) #the top 25 models
nrow(LM3_output[LM3_output$delta<7, ]) #models with delta AIC < 7 (returns [1] 49)
w<-Weights(LM3_output)
w

#Refit best linear model
bestmodel<-get.models(LM3_output, 1)[[1]]
z<-lm(bestmodel,data=morpho2)
summary(z) #the best model has tarsus length as the most important factor (p=0.000448)