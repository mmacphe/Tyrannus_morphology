rm(list=ls())

require("phytools")

### Set source directory to the folder this file came from within RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### Bring in base phylogeny ###
tyrannus_base<-read.tree("MacPherson_Tyrannus_base.tre")

### make tree ultrametric ###
is.ultrametric(tyrannus_base)
tyrannus_um<-force.ultrametric(tyrannus_base,method="nnls")

#Rooted; includes branch lengths. 
plot(tyrannus_um)

### Find edge lengths 
require(ape)
node.depth.edgelength(tyrannus_um)

### read in list of missing taxa sensu clements 2019 ###
taxa2add<-readLines(file("Tyrannus_taxa2add.txt"))
taxa2add<-gsub(" ","_",taxa2add)

### Add in T. savana monachus tips based on Gomez-Bahamon phylogeny ###
taxa2add<- taxa2add[-10] #remove T. savana monachus as a single OTU because we are separating Central American and South American populations
taxa2add
tyrannus_mod <-bind.tip(tyrannus_um,tip.label="Tyrannus_savana_monachus_CA",where=which(tyrannus_um $tip.label=="Tyrannus_savana_savana"),position=1.4)
tyrannus_mod <-bind.tip(tyrannus_mod,tip.label="Tyrannus_savana_monachus_SA",where=which(tyrannus_mod $tip.label=="Tyrannus_savana_monachus_CA"),position=(3682/1e6))
plot(tyrannus_mod)

### Add missing tips to the phylogeny assuming a MRCA of 0.5 mya ### Note: change the 'position' argument below to add missing tips at different branch lengths
for(i in 1:length(taxa2add)){
  taxa2add[i]
  binomial<-paste(strsplit(taxa2add[i],"_")[[1]][1:2],collapse="_")
  recipient<-tyrannus_mod $tip.label[grep(binomial, tyrannus_mod $tip.label)]
  
  if(length(recipient)>1){
    recipient<-findMRCA(tyrannus_mod,recipient)	
    tyrannus_mod <-bind.tip(tyrannus_mod,taxa2add[i],where=recipient)
  }else{
    tyrannus_mod <-bind.tip(tyrannus_mod,taxa2add[i],where=which(tyrannus_mod $tip.label== recipient),position=0.5) #the position argument is what determines the terminal branch length for the added taxon in millions of years
  }
}

plot(tyrannus_mod)

### Replace any 0-length branches with a very small value ###
tyrannus_mod$edge.length[tyrannus_mod$edge.length==0]<-0.001
plot(tyrannus_mod)

### Write out phylogeny ###
write.tree(tyrannus_mod,file="Tyrannus_phylogeny_0.93.tre")
