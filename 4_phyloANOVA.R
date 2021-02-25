rm(list=ls())

require(phytools)

phy<-read.tree('./Output Files/Tyrannus_phylogeny.tre')

### Bring in morphology dataset with average morphology values for each OTU and PC scores
otu_avg<-read.csv('./Output Files/Tyrannus morphology + PCA avg.csv', row.names = 1)
species<-rownames(otu_avg)
otu_avg<-cbind(species,otu_avg)
View(otu_avg)

### Check that tree$tip.label is the same as otu_avg
phy$tip.label[!phy$tip.label %in% otu_avg$species]

### phylANOVA on OTU morphometric means while controlling for body size ###
response_var<-colnames(otu_avg)[c(3,4,5,7,8,9)]

### Set up predictor variable for phylANOVA ###
migration<-otu_avg$Strategy
names(migration)<-rownames(otu_avg)

### Set up tarsus for size correction #
otu_avg$Tarsus.Average <-as.numeric(otu_avg$Tarsus.Average)
tarsus.length<-otu_avg$Tarsus.Average
names(tarsus.length)<-rownames(otu_avg)

#Create empty list to store output of each for loop iteration
phyl.resid_output<-list()
phylANOVA_output<-list() 

#Conduct the phylogenetic ANOVA he following code works without assuming y is the same order as tree$tip.label
for(i in 1:length(response_var)){
	print(paste0("Character ",i," -- ",response_var[i]))#Report character
	character_of_interest<-otu_avg[,response_var[i]] #Extract character of interest (coi)
	names(character_of_interest)<-rownames(otu_avg) #Add names to vector of coi
	phyl.resid_output[[i]]<-phyl.resid(phy,character_of_interest,tarsus.length,method="lambda") #Perform phylogenetic residuals
	phylANOVA_output[[i]]<-phylANOVA(phy, migration, nsim=100000, phyl.resid_output[[i]]$resid[,1])
}

names(phyl.resid_output)<-response_var
names(phylANOVA_output)<-response_var 

### Write out output ###
sink("phylANOVA_output.txt")
for(i in 1:length(phylANOVA_output)){
  cat(paste0("Character ",i," -- ",response_var[i]))
  cat("\n\n")
  print(phylANOVA_output[[i]])
  cat("#############")
  cat("\n\n")
}
sink()

### Conduct the phylogenetic ANOVA for PC scores that are not size-corrected ###
phy$tip.label[!phy$tip.label %in% rownames(otu_avg)]
response_var2<-colnames(otu_avg)[c(11,12,13,14)]
#migration<-otu_avg$Strategy
#names(migration)<-rownames(otu_avg)

phylANOVA_PCoutput<-list() 
#Question: How to label y? (the tip.labels so that it doesn't assume same order as tree$tip.label)
for(i in 1:length(response_var2)){
  print(paste0("Character ",i," -- ",response_var2[i]))#Report character
  character_of_interest<-otu_avg[,response_var2[i]] #Extract character of interest (coi)
  names(character_of_interest)<-rownames(otu_avg) #Add names to vector of coi
  phylANOVA_PCoutput[[i]]<-phylANOVA(tree=phy, x=migration, y=character_of_interest, nsim=100000) 
}

names(phylANOVA_PCoutput)<-response_var2 

### Write out output ###
sink("phylANOVA_output_PCscores.txt")

for(i in 1:length(phylANOVA_PCoutput)){
	cat(paste0("Character ",i," -- ",response_var2[i]))
	cat("\n\n")
	print(phylANOVA_PCoutput[[i]])
	cat("#############")
	cat("\n\n")
}
sink()

### Repeat for CV data: Bring in cv_summary_table_v4
otu_cv_avg<-read.csv('./Output Files/cv_summary_v5.csv', row.names = 1) 
species<-rownames(otu_cv_avg)
otu_cv_avg<-cbind(species,otu_cv_avg)
View(otu_cv_avg)

### Check that tree$tip.label is the same as otu_avg
phy$tip.label[!phy$tip.label %in% otu_cv_avg$species]

response_var3<-colnames(otu_cv_avg)[c(3,4,5,6,7,8,9,11,12,13,14)]

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
  phylANOVA_CVoutput[[i]]<-phylANOVA(phy, migration, y=character_of_interest, nsim=100000)
}

names(phylANOVA_CVoutput)<-response_var3

### Write out output ###
sink("phylANOVA_cv_output.txt")

for(i in 1:length(phylANOVA_CVoutput)){
  cat(paste0("Character ",i," -- ",response_var3[i]))
  cat("\n\n")
  print(phylANOVA_CVoutput[[i]])
  cat("#############")
  cat("\n\n")
}
sink()
