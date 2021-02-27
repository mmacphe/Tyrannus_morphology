rm(list=ls())

### Set source directory to the folder this file came from within RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

morpho<-read.csv("Tyrannus_voucher_table.csv")

#What we're trying to do is to see if, within each taxon:
#1) morphologies differ between the sexes, and 
#2) morphologies differ between age classes

### Remove any individuals with unresolved taxonomic identities
morpho<-morpho[-grep("unknown",morpho$subspecies),]

### Classify each taxon using trinomial naming
morpho$subspecies<-trimws(morpho$subspecies)
morpho$trinomial<-paste(morpho$Species,morpho$subspecies)
morpho$trinomial<-trimws(morpho$trinomial)
names(table(morpho$trinomial)) #see that unknown subspecies were removed

### Standardize notation for Sex column
names(table(morpho$Sex))
morpho<-morpho[!morpho$Sex=="Unknown",] #removes unknown sex individuals
names(table(morpho$Sex)) #should now only have "Female" and "Male"

### For Use in later steps: Make separate male and female voucher tables
morpho_Females<-morpho[!morpho$Sex=="Male",]
write.csv(morpho_Females, file="Voucher_Table_Females.csv")

morpho_Males<-morpho[!morpho$Sex=="Female",]
write.csv(morpho_Males, file="Voucher_Table_Males.csv")

response_var<-colnames(morpho)[c(18,19,20,21,22,23,24)]

#Create empty list to store output of each for loop iteration
LM1_output<-list()
LM2_output<-list()

#Conduct the phylogenetic ANOVA he following code works without assuming y is the same order as tree$tip.label
for(i in 1:length(response_var)){
  print(paste0("Character",i," -- ",response_var[i]))#Report character
  character_of_interest<-morpho[,response_var[i]] #Extract character of interest (coi)
  names(character_of_interest)<-morpho$trinomial #Add names to vector of coi
  LM1_output[[i]]<-AIC(lm(morpho[,response_var[i]]~morpho$trinomial, data=morpho)) #Perform phylogenetic residuals
  LM2_output[[i]]<-AIC(lm(morpho[,response_var[i]]~morpho$trinomial*morpho$Sex, data=morpho))
}

names(LM1_output)<-response_var
names(LM2_output)<-response_var 

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
#Here, we subtracted LM2 (accounted for sex) from LM1 (did not account for sex).
#All negative values indicate that the model was improved by the accounting of sex, and positive values that it was not.


#/\/\/\/\/#Q2) Do the morphologies differ between age classes?
### Standardize notation for Age column
names(table(morpho$Age)) #to see which assignments are in the Age column
morpho$Age[morpho$Age=="SecondYear"]<-"Adult"
morpho$Age[morpho$Age=="Fledgling"]<-"Juvenile"
morpho$Age[morpho$Age=="HatchYear"]<-"Juvenile"
morpho$Age[morpho$Age=="Nestling"]<-"Juvenile"
morpho<-morpho[!morpho$Age=="Unknown",] #removes unknown age individuals
names(table(morpho$Age)) #should now only have "Adult" and "Juvenile"

#There is only one age category (Adult) so LM1 should always be better than LM2

lm1<-lm(Morph_data$`BL Average` ~ Morph_data$trinomial, data=Morph_data)
lm2<-lm(Morph_data$`BL Average`~ Morph_data$trinomial*Morph_data$Age)

AIC(lm1, lm2)
#    df      AIC
#lm1 29 4072.592
#lm2 57 9657.571

lm1<-lm(Morph_data$`BW Average`~ Morph_data$trinomial, data=Morph_data)
lm2<-lm(Morph_data$`BW Average`~ Morph_data$trinomial*Morph_data$Age)

AIC(lm1, lm2)
#    df      AIC
#lm1 29 2409.433
#lm2 57 9657.571

lm1<-lm(Morph_data$`BD Average`~ Morph_data$trinomial, data=Morph_data)
lm2<-lm(Morph_data$`BD Average`~ Morph_data$trinomial*Morph_data$Age)

AIC(lm1, lm2)
#    df      AIC
#lm1 29 1483.731
#lm2 57 9657.571

lm1<-lm(Morph_data$`Kipp's Average`~ Morph_data$trinomial, data=Morph_data)
lm2<-lm(Morph_data$`Kipp's Average`~ Morph_data$trinomial*Morph_data$Age)

AIC(lm1, lm2)
#    df      AIC
#lm1 29 5903.397
#lm2 57 9657.571

lm1<-lm(Morph_data$`WC Average` ~ Morph_data$trinomial, data=Morph_data)
lm2<-lm(Morph_data$`WC Average`~ Morph_data$trinomial*Morph_data$Age)

AIC(lm1, lm2)
#    df      AIC
#lm1 29 7943.603
#lm2 57 9657.571

lm1<-lm(Morph_data$Tail ~ Morph_data$trinomial, data=Morph_data)
lm2<-lm(Morph_data$Tail ~ Morph_data$trinomial*Morph_data$Age)

AIC(lm1, lm2)
#    df       AIC
#lm1 29 10195.159
#lm2 57  9657.571

#/\/\/\/\/#Q3) Do the morphologies differ with tarsus length?
lm1<-lm(Morph_data$`BL Average` ~ Morph_data$trinomial, data=Morph_data)
lm2<-lm(Morph_data$`BL Average`~ Morph_data$trinomial*Morph_data$`Tarsus Average`)

AIC(lm1, lm2)
#    df      AIC
#lm1 29 4072.592
#lm2 57 3839.143

lm1<-lm(Morph_data$`BW Average` ~ Morph_data$trinomial, data=Morph_data)
lm2<-lm(Morph_data$`BW Average`~ Morph_data$trinomial*Morph_data$`Tarsus Average`)

AIC(lm1, lm2)
#    df      AIC
#lm1 29 2409.433
#lm2 57 2218.043

lm1<-lm(Morph_data$`BD Average` ~ Morph_data$trinomial, data=Morph_data)
lm2<-lm(Morph_data$`BD Average`~ Morph_data$trinomial*Morph_data$`Tarsus Average`)

AIC(lm1, lm2)
#    df      AIC
#lm1 29 1483.731
#lm2 57 1337.364

lm1<-lm(Morph_data$`Kipp's Average` ~ Morph_data$trinomial, data=Morph_data)
lm2<-lm(Morph_data$`Kipp's Average`~ Morph_data$trinomial*Morph_data$`Tarsus Average`)

AIC(lm1, lm2)
#    df      AIC
#lm1 29 5903.397
#lm2 57 5896.600

lm1<-lm(Morph_data$`WC Average` ~ Morph_data$trinomial, data=Morph_data)
lm2<-lm(Morph_data$`WC Average`~ Morph_data$trinomial*Morph_data$`Tarsus Average`)

AIC(lm1, lm2)
#    df      AIC
#lm1 29 7943.603
#lm2 57 7366.929

lm1<-lm(Morph_data$Tail ~ Morph_data$trinomial, data=Morph_data)
lm2<-lm(Morph_data$Tail ~ Morph_data$trinomial*Morph_data$`Tarsus Average`)

AIC(lm1, lm2)
#    df       AIC
#lm1 29 10195.159
#lm2 57  9557.213