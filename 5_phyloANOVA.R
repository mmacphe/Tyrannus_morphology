rm(list=ls())

require(phytools)

### Set source directory to the folder this file came from within RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

phy<-read.tree('./Output Files/Tyrannus_phylogeny.tre')
#phy<-read.tree('./Output Files/Tyrannus_phylogeny_0.2.tre')
#phy<-read.tree('./Output Files/Tyrannus_phylogeny_0.93.tre')

### Bring in morphology dataset with average morphology values for each OTU and PC scores
otu_avg<-read.csv('./Output Files/Tyrannus morphology + PCA avg_Adults.csv', row.names = 1)

###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

### Note: If assessing sexes separately, use one of the following lines of code instead of the above line, then proceed.
#otu_avg<-read.csv('./Output Files/Tyrannus morphology + PCA avg_Females_Adults.csv', row.names = 1)
### OR
#otu_avg<-read.csv('./Output Files/Tyrannus morphology + PCA avg_Males_Adults.csv', row.names = 1)

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###

species<-rownames(otu_avg)
otu_avg<-cbind(species,otu_avg)
View(otu_avg)

### Check that tree$tip.label is the same as otu_avg
phy$tip.label[!phy$tip.label %in% otu_avg$species]

### phylANOVA on OTU morphometric means while controlling for body size ###
response_var<-colnames(otu_avg)[c(4,5,6,8,9,10)]

### Set up predictor variable for phylANOVA ###
migration<-otu_avg$Strategy
names(migration)<-rownames(otu_avg)

### Set up tarsus for size correction #
otu_avg$Tarsus.Average <-as.numeric(otu_avg$Tarsus.Average)
tarsus.length<-otu_avg$Tarsus.Average
names(tarsus.length)<-rownames(otu_avg)

### Create empty list to store output of each for loop iteration
phyl.resid_output<-list()
phylANOVA_output<-list() 

### Conduct the phylogenetic ANOVA he following code works without assuming y is the same order as tree$tip.label
for(i in 1:length(response_var)){
	print(paste0("Character ",i," -- ",response_var[i]))#Report character
	character_of_interest<-otu_avg[,response_var[i]] #Extract character of interest (coi)
	names(character_of_interest)<-rownames(otu_avg) #Add names to vector of coi
	phyl.resid_output[[i]]<-phyl.resid(phy,tarsus.length, character_of_interest,method="lambda") #Perform phylogenetic residuals
	phylANOVA_output[[i]]<-phylANOVA(phy, x=migration, y=phyl.resid_output[[i]]$resid[,1], p.adj="none")
}

names(phyl.resid_output)<-response_var
names(phylANOVA_output)<-response_var 

### Write out output ###
sink("phylANOVA_output.txt")
#sink("phylANOVA_output_0.2.txt")
#sink("phylANOVA_output_0.93.txt")

###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

### Note: If assessing sexes separately, use one of the following lines of code instead of the above line, then proceed.
#sink("phylANOVA_output_Females.txt")
### OR
#sink("phylANOVA_output_Males.txt")

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###

for(i in 1:length(phylANOVA_output)){
  cat(paste0("Character ",i," -- ",response_var[i]))
  cat("\n\n")
  print(phylANOVA_output[[i]])
  cat("#############")
  cat("\n\n")
}
sink()

### Conduct the phylogenetic ANVOA with p.adj="bonferroni" to account for multiple hypothesis testing
## Create empty list to store output of each for loop iteration
phyl.resid_output<-list()
phylANOVA_output<-list() 

for(i in 1:length(response_var)){
  print(paste0("Character ",i," -- ",response_var[i]))#Report character
  character_of_interest<-otu_avg[,response_var[i]] #Extract character of interest (coi)
  names(character_of_interest)<-rownames(otu_avg) #Add names to vector of coi
  phyl.resid_output[[i]]<-phyl.resid(phy,tarsus.length, character_of_interest,method="lambda") #Perform phylogenetic residuals
  phylANOVA_output[[i]]<-phylANOVA(phy, x=migration, y=phyl.resid_output[[i]]$resid[,1], p.adj="bonferroni")
}
names(phyl.resid_output)<-response_var
names(phylANOVA_output)<-response_var 

### Write out output ###
sink("phylANOVA_bonferonni_output.txt")

###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

### Note: If assessing sexes separately, use one of the following lines of code instead of the above line, then proceed.
#sink("phylANOVA_bonferonni_Females.txt")
### OR
#sink("phylANOVA_bonferonni_Males.txt")

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###

for(i in 1:length(phylANOVA_output)){
  cat(paste0("Character ",i," -- ",response_var[i]))
  cat("\n\n")
  print(phylANOVA_output[[i]])
  cat("#############")
  cat("\n\n")
}
sink()

###################################################################################
### Conduct the phylogenetic ANOVA for metrics that do not need to be size-corrected ###
phy$tip.label[!phy$tip.label %in% rownames(otu_avg)]
response_var2<-colnames(otu_avg)[c(12,13)]

phylANOVA_PCoutput<-list() 

for(i in 1:length(response_var2)){
  print(paste0("Character ",i," -- ",response_var2[i]))#Report character
  character_of_interest<-otu_avg[,response_var2[i]] #Extract character of interest (coi)
  names(character_of_interest)<-rownames(otu_avg) #Add names to vector of coi
  phylANOVA_PCoutput[[i]]<-phylANOVA(tree=phy, x=migration, y=character_of_interest, p.adj="none") 
}

names(phylANOVA_PCoutput)<-response_var2 

### Write out output ###
sink("phylANOVA_output_nosizecorrectionneeded.txt")
#sink("phylANOVA_output_nosizecorrectionneeded_0.2.txt")
#sink("phylANOVA_output_nosizecorrectionneeded_0.93.txt")

###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

### Note: If assessing sexes separately, use one of the following lines of code instead of the above line, then proceed.
#sink("phylANOVA_output_PCscores_Females.txt")
### OR
#sink("phylANOVA_output_PCscores_Males.txt")

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###

for(i in 1:length(phylANOVA_PCoutput)){
	cat(paste0("Character ",i," -- ",response_var2[i]))
	cat("\n\n")
	print(phylANOVA_PCoutput[[i]])
	cat("#############")
	cat("\n\n")
}
sink()

## Conduct the phylogenetic ANVOA on metrics that do not need to be size-corrected with p.adj="bonferroni" to account for multiple hypothesis testing
phylANOVA_PCoutput<-list() 

for(i in 1:length(response_var2)){
  print(paste0("Character ",i," -- ",response_var2[i]))#Report character
  character_of_interest<-otu_avg[,response_var2[i]] #Extract character of interest (coi)
  names(character_of_interest)<-rownames(otu_avg) #Add names to vector of coi
  phylANOVA_PCoutput[[i]]<-phylANOVA(tree=phy, x=migration, y=character_of_interest, p.adj="bonferroni") 
}

names(phylANOVA_PCoutput)<-response_var2 

### Write out output ###
sink("phylANOVA_output_bonferroni_nosizecorrectionneeded.txt")

###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

### Note: If assessing sexes separately, use one of the following lines of code instead of the above line, then proceed.
#sink("phylANOVA_output_bonferroni_PCscores_Females.txt")
### OR
#sink("phylANOVA_output_bonferroni_PCscores_Males.txt")

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###

for(i in 1:length(phylANOVA_PCoutput)){
  cat(paste0("Character ",i," -- ",response_var2[i]))
  cat("\n\n")
  print(phylANOVA_PCoutput[[i]])
  cat("#############")
  cat("\n\n")
}
sink()

########################################################
### Repeat for CV data
otu_cv_avg<-read.csv('./Output Files/cv_summary_Adults.csv', row.names = 1) 

###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

### Note: If assessing sexes separately, use one of the following lines of code instead of the above line, then proceed.
#otu_cv_avg<-read.csv('./Output Files/cv_summary_Females_Adults.csv', row.names = 1) 
### OR
#otu_cv_avg<-read.csv('./Output Files/cv_summary_Males_Adults.csv', row.names = 1) 

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###

species<-rownames(otu_cv_avg)
otu_cv_avg<-cbind(species,otu_cv_avg)
View(otu_cv_avg)

### Check that tree$tip.label is the same as otu_avg
phy$tip.label[!phy$tip.label %in% otu_cv_avg$species]

response_var3<-colnames(otu_cv_avg)[c(4,5,6,7,8,9,10,12,13)]

###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

### Note: If assessing females, use one of the following lines of code instead of the above line, then proceed.
#response_var3<-colnames(otu_cv_avg)[c(4,5,6,7,8,9,10,11,13,14)] #run anova for tails separately (see below) for females bc T. s sanctamartae has NA

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###

### Set up predictor variable for phylANOVA ###
migration<-otu_cv_avg$Strategy
names(migration)<-rownames(otu_cv_avg)

#Create empty list to store output of each for loop iteration
phylANOVA_CVoutput<-list() 

### Conduct the phylogenetic ANOVA on coefficients of variation ###
for(i in 1:length(response_var3)){
  print(paste0("Character ",i," -- ",response_var3[i]))#Report character
  character_of_interest<-otu_cv_avg[,response_var3[i]] #Extract character of interest (coi)
  names(character_of_interest)<-rownames(otu_cv_avg) #Add names to vector of coi
  phylANOVA_CVoutput[[i]]<-phylANOVA(phy, x=migration, y=character_of_interest, p.adj="none")
}

names(phylANOVA_CVoutput)<-response_var3

### Write out output ###
sink("phylANOVA_cv_output.txt")
#sink("phylANOVA_cv_output_0.2.txt")
#sink("phylANOVA_cv_output_0.93.txt")

###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

#Note: If assessing sexes separately, use one of the following lines of code instead of the above line, then proceed.
#sink("phylANOVA_cv_output_Females.txt")
### OR
#sink("phylANOVA_cv_output_Males.txt")

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###

for(i in 1:length(phylANOVA_CVoutput)){
  cat(paste0("Character ",i," -- ",response_var3[i]))
  cat("\n\n")
  print(phylANOVA_CVoutput[[i]])
  cat("#############")
  cat("\n\n")
}
sink()

#Do phylogenetic ANOVA again on CV values with bonferonni correction for multiple hypothesis testing
phylANOVA_CVoutput<-list() 

### Conduct the phylogenetic ANOVA on coefficients of variation ###
for(i in 1:length(response_var3)){
  print(paste0("Character ",i," -- ",response_var3[i]))#Report character
  character_of_interest<-otu_cv_avg[,response_var3[i]] #Extract character of interest (coi)
  names(character_of_interest)<-rownames(otu_cv_avg) #Add names to vector of coi
  phylANOVA_CVoutput[[i]]<-phylANOVA(phy, x=migration, y=character_of_interest, p.adj="bonferroni")
}

names(phylANOVA_CVoutput)<-response_var3

### Write out output ###
sink("phylANOVA_cv_bonferroni_output.txt")

###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

#Note: If assessing sexes separately, use one of the following lines of code instead of the above line, then proceed.
#sink("phylANOVA_cv_bonferroni_output_Females.txt")
### OR
#sink("phylANOVA_cv_bonferroni_output_Males.txt")

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###

for(i in 1:length(phylANOVA_CVoutput)){
  cat(paste0("Character ",i," -- ",response_var3[i]))
  cat("\n\n")
  print(phylANOVA_CVoutput[[i]])
  cat("#############")
  cat("\n\n")
}
sink()


###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

#Note: Run the following code in addition to conduct phylANOVA on tail length in females.
### Tail for females 
#otu_cv_avg<-read.csv('./Output Files/cv_summary_Females_Adults.csv', row.names = 1) 
#otu_cv_avg_females<-otu_cv_avg[c(1,2,3,4,5,6,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,24,25,26,27,28),]
#View(otu_cv_avg_females)

#response_var4<-colnames(otu_cv_avg_females)[9]

### Set up predictor variable for phylANOVA ###
#migration<-otu_cv_avg_females$Strategy
#names(migration)<-rownames(otu_cv_avg_females)

### Remove species with no variation in tail
#phy_trim<-drop.tip(phy, c("Tyrannus_savana_sanctaemartae","Tyrannus_caudifasciatus_jamaicensis"))
#plot(phy_trim)

### Check that tree$tip.label is the same as otu_avg
#phy_trim$tip.label[!phy_trim$tip.label %in% otu_cv_avg_females$species]

### Create empty list to store output of each for loop iteration
#phylANOVA_CVoutput<-list() 

### Conduct the phylogenetic ANOVA on coefficients of variation ###
#for(i in 1:length(response_var4)){
#  print(paste0("Character ",i," -- ",response_var4[i]))#Report character
#  character_of_interest<-otu_cv_avg_females[,response_var4[i]] #Extract character of interest (coi)
#  names(character_of_interest)<-rownames(otu_cv_avg_females) #Add names to vector of coi
#  phylANOVA_CVoutput[[i]]<-phylANOVA(phy_trim, x=migration, y=character_of_interest, p.adj="none")
#}

#names(phylANOVA_CVoutput)<-response_var4

### Write out output ###
#sink("phylANOVA_cv_output_Females_tail.txt")

#for(i in 1:length(phylANOVA_CVoutput)){
#  cat(paste0("Character ",i," -- ",response_var4[i]))
#  cat("\n\n")
#  print(phylANOVA_CVoutput[[i]])
#  cat("#############")
#  cat("\n\n")
#}
#sink()

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###

### Build residuals dataframe for a figure that is size and phylogeny corrected
phyl_resid<-as.data.frame(phyl.resid_output$BL.Average,row.names = phy$tip.label)
phyl_resid$BL<-phyl_resid$resid
phyl_resid2<-as.data.frame(phyl.resid_output$BW.Average, row.names=phy$tip.label)
phyl_resid$BW<-phyl_resid2$resid
phyl_resid3<-as.data.frame(phyl.resid_output$BD.Average, row.names=phy$tip.label)
phyl_resid$BD<-phyl_resid3$resid
phyl_resid4<-as.data.frame(phyl.resid_output$Kipp.s.Distance, row.names=phy$tip.label)
phyl_resid$KD<-phyl_resid4$resid
phyl_resid5<-as.data.frame(phyl.resid_output$WC.Average, row.names=phy$tip.label)
phyl_resid$WC<-phyl_resid5$resid
phyl_resid6<-as.data.frame(phyl.resid_output$Tail, row.names=phy$tip.label)
phyl_resid$TL<-phyl_resid6$resid

residuals<-merge(phyl_resid, migration, by=0, row.names=TRUE) #the taxonomic order can be checked visually by comparin residuals to phyl.resid_output$ coi
residuals$migration<-residuals$y

write.csv(residuals, "phylANOVA_tarsus-corrected_residuals.csv")

###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

### Note: If assessing sexes separately, use one of the following lines of code instead of the above line, then proceed.
#write.csv(residuals, "phylANOVA_tarsus-corrected_residuals_Females.csv")
### OR
#write.csv(residuals, "phylANOVA_tarsus-corrected_residuals_Males.csv")

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###
