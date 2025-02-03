renv::restore()

library(ape)
library(phytools)
library(splits)
library(tidyverse)

rm(list = ls())

cyanistes <- read.nexus("C:/Users/daank/OneDrive - University of Twente/Documents/Rijksuniversiteit Groningen/Year 1 - 24-25/Island Biology/Practicals/Practical 1 - Cyanistes/Cyanistes.tre")
cyanistes


####################################################
#Visualize and handle trees
cyanistes$tip.label

plot(cyanistes)
plot(cyanistes, cex = 0.4)

par(mar = c(0,0,0,0))
plot(ladderize(cyanistes), cex = 0.4, no.margin = TRUE)

is.ultrametric(cyanistes)

plot(ladderize(cyanistes), cex = 0.4, no.margin = TRUE)
nodelabels(cex = 0.3, frame = "circle", bg = "grey")

plot(ladderize(cyanistes), type = "fan", cex = 0.3, no.margin = TRUE)

subtreeplot(ladderize(cyanistes), cex = 0.4)

cyanistes_ingroup <- drop.tip(cyanistes, c("Parus_major_major_KP759174.1",
                                        "Parus_major_DQ792786.1", 
                                        "Parus_major_DQ792787.1",
                                        "Parus_major_EU167009.1", 
                                        "Parus_major_KJ456375.1"))

par(mfrow = c(1,1))
plot(ladderize(cyanistes_ingroup), cex = 0.4, no.margin = TRUE)


####################################################
#Ancestral area reconstruction
island_data<-read.csv("C:/Users/daank/OneDrive - University of Twente/Documents/Rijksuniversiteit Groningen/Year 1 - 24-25/Island Biology/Practicals/Practical 1 - Cyanistes/Cyanistes_distribution.csv",header=T)

View(island_data)

island_d<-as.data.frame(island_data$Distribution)
taxa<-as.data.frame(island_data$Species)
islands<-as.data.frame(island_d[match(cyanistes$tip.label,taxa[,1]),])
islands<-t(islands)
islands<-as.character(islands)
names(islands)<-cyanistes$tip.label

set.seed(3)
cyanistes_simmap<-make.simmap(cyanistes,islands,model="ER",nsim=1000)

pd<-summary(cyanistes_simmap,plot=FALSE)

par(oma=c(0,0,0,0))
cols<-setNames(palette()[1:length(unique(islands))],sort(unique(islands)))
plot(cyanistes_simmap[[1]],cols,fsize=0.8)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=25,fsize=0.8)

nodelabels(pie=pd$ace,piecol=cols,cex=0.3)

tiplabels(pie=to.matrix(islands,sort(unique(islands))),piecol=cols,cex=0.1)

#Questions
#Based on this tree, can you infer how many colonisation events of the Canary Islands there have been?
  
#Combining with Figtree to check for the node ages, can you estimate when these colonisations took place?
  
#Do you see any evidence for back-colonisation in the cytochrome-B data? (this would be a species or clade of mainland individuals, nested within an island clade).


####################################################
#Species delimitation
cyanistes_gmyc <- gmyc(cyanistes_ingroup)

summary(cyanistes_gmyc)

spec.list(cyanistes_gmyc)

plot(cyanistes_gmyc)


#Questions
#How many “species” does the model identify?
##9 species -> GMYC_spec indicates 9
  
#Do the clusters match species based on the tip names?
##No, sometimes the tip names are spread over multiple clades.
  
#Would you recommend a taxonomic revision of this group?
##Yes, this taxonomy is based on just one model and should be confirmed by other methods.
  
#How many endemic species of Cyanistes are there in the Canary Islands?
##4
  
#How many endemic subspecies?
##0
  
#What other sources of evidence would be useful to define species of blue tits in the Canarian archipelago?
##Song, morphology

#What species would you give priority to conservation (if any)?
##


