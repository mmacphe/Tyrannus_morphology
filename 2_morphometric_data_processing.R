rm(list=ls())

require(phytools)

### Set source directory to the folder this file came from within RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### Read in morphological data set
morpho <- read.csv("Tyrannus_voucher_table.csv")

###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

#Note: If assessing sexes separately, use one of the following lines of code instead of the above line, then proceed.
#morpho<-read.csv('./Output Files/Voucher_Table_Females.csv')
#OR
#morpho<-read.csv('./Output Files/Voucher_Table_Males.csv')

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###

### Remove individuals with unknown subspecies ID 
### Note: This does not remove species that lack subspecies, only individuals that were not identified to the lowest taxonomic division
names(table(morpho$subspecies)) #run following line of code if there are "unknown" subspecies
morpho<-morpho[-grep("unknown",morpho$subspecies),] #do not run this for sex-specific analyses

### Read in phylogeny 
phy<-read.tree('./Output Files/Tyrannus_phylogeny.tre')

### Match morphometric data to names in phylogeny 
morpho<-morpho[!morpho$Species=="",]

morpho$subspecies[is.na(morpho$subspecies)]<-""
morpho$subspecies<-trimws(morpho$subspecies) #removes leading/trailing whitespace
morpho$Species<-trimws(morpho$Species)

morpho$Species<-gsub(" ","_",paste(morpho$Species))
morpho$tip.label<-gsub(" ","_",paste(morpho$Species,morpho$subspecies,sep="_"))
morpho$tip.label[morpho$subspecies==""]<-morpho$Species[morpho$subspecies==""]

names(table(morpho$tip.label))

phy$tip.label[!phy$tip.label %in% morpho$tip.label] #Check that returns "character(0)"

### Save morpho as a .csv ###
write.csv(morpho, file="Tyrannus morphology data.csv")

###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

#Note: If assessing sexes separately, use one of the following lines to write out the data file, then proceed.
#write.csv(morpho, file="Tyrannus morphology data_Females.csv")
#OR
#write.csv(morpho, file="Tyrannus morphology data_Males.csv")

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###

### Split the morphology data sets up into categories by tip.label
morpho$tip.label
names(table(morpho$tip.label))

morpho_split<-split(morpho,morpho$tip.label)
names(morpho_split)<-gsub(" ","_",names(morpho_split))
names(morpho_split) 

phy$tip.label[!phy$tip.label %in% names(morpho_split)] #Check that this returns "character(0)"
names(morpho_split)[! names(morpho_split) %in% phy$tip.label] #Check that this also returns "character(0)"

### Create summary output df for data frame
morpho_trim<-lapply(morpho_split,function(x) x[colnames(x) %in% c("BL.Average","BW.Average","BD.Average", "Kipp.s.Average", "WC.Average", "Tail", "Tarsus.Average")])

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
write.csv(morphology_summary,file="morphology_summary_table.csv")

###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

#Note: If assessing sexes separately, use one of the following lines to write out the data file, then proceed.
#write.csv(morphology_summary, file="morphology_summary_table_Females.csv")
#OR
#write.csv(morpho, file="morphology_summary_table_Males.csv")

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###

### Write out .csv file just with the mean for each OTU to be used in further analyses
Tyrannus_OTU_averages<-data.frame(matrix(paste(round(otu_avg_df,2)),ncol=ncol(otu_avg_df)))
colnames(Tyrannus_OTU_averages)<-colnames(otu_avg_df)
rownames(Tyrannus_OTU_averages)<-names(morpho_split)

### Add migratory strategy as a column in the Tyrannus OTU averages .csv ###
migratory_data<-read.csv("Tyrannus_subspecies_MigrationStrategies.csv",row.names=1)
migratory_data[rownames(Tyrannus_OTU_averages),]
Tyrannus_data<-merge(Tyrannus_OTU_averages, migratory_data, by=0)
write.csv(Tyrannus_data, file="Tyrannus_data.csv")

###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

#Note: If assessing sexes separately, use one of the following lines to write out the data file, then proceed.
#write.csv(Tyrannus_data, file="Tyrannus_data_Females.csv")
#OR
#write.csv(Tyrannus_data, file="Tyrannus_data_Males.csv")

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###

### Write out .csv file with coefficient of variation for each OTU to be used in further analyses
cv_summary<-otu_cv_df
rownames(cv_summary)<-names(morpho_split)
cv_summary<-as.data.frame(cv_summary)
migratory_data[rownames(cv_summary),]
CV_data<-merge(cv_summary, migratory_data, by=0)
write.csv(CV_data, file = "cv_summary_table.csv")

###NOTE/\/\/\/\/\/\/\/\/\/\/\###
###########################
###/\/\/\/\/\/\/\/\/\/\/\###

#Note: If assessing sexes separately, use one of the following lines to write out the data file, then proceed.
#write.csv(CV_data, file="cv_summary_table_Females.csv")
#OR
#write.csv(morpho, file="cv_summary_table_Males.csv")

###/\/\/\/\/\/\/\/\/\/\/\###
###########################
###END/\/\/\/\/\/\/\/\/\/\/\###