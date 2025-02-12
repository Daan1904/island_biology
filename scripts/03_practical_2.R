renv::restore()

install.packages('DAISIE',dependencies = TRUE)
install.packages('DAISIEprep',dependencies = TRUE)
install.packages('ape',dependencies = TRUE)

library(DAISIE)
library(ape)
library(DAISIEprep)
library(tidyverse)

rm(list = ls())

###Questions to answer in the report
##How many times has Insula colonised the Caribbean islands?
# Only once. Look at the tree on GitHub and see in which clades there is mainland included. That is the case nowhere, only the outgroup is mainland. That makes that there is just 1 colonization.

##How many times has Insula colonised Jamaica?
#There have been 5 colonizations. Count in the tree made in R.

##How many radiations have occurred on Jamaica?
#There have been 2 radiations. There are 2 clades with recent cladogenesis.

##What are the rates of colonisation, speciation and extinction for Insula in Jamaica?
#The outcomes have been taken from model 3 as that was found as the best model via AIC
#Colonization = gamma = 0.02762404
#Speciation = lambda_c = 6.76418
#Speciation = lambda_a = 0
#Extinction = mu = 8.775854

##Is there evidence for diversity-dependence in the Insula species of Jamaica?
#No, there is no evidence for diversity dependence. Diversity dependence is indicated by the carrying capacity K and that variable is infinite in the best model (M3).

##Is there evidence for equilibrium dynamics on the island of Jamaica?
#Yes, there is evidence since the number of species over time is stable in the simulation. A plateau is shown in the graph for both the endemic and non-endemic species. The dynamics consist of speciation and extinction and if the the number of species is stable, the dynamics are in equilibrium.

##How has the diversity of Insula on Jamaica changed through time (according to DAISIE simulations)?
#The diversity of Insula on Jamaica showed a very quick increase in the first 1 million years, after which the endemic diversity stabilizes around 10 species. The diversity of the non-endemic species is stable around 3 species.

##Modify one or two parameters of your choice and run simulations again under those parameters. How does this affect the simulations? Describe.
#Option 1: low colonization / regular extinction -> equilibrium
#Option 2: high colonization / low extinction -> increase in diversity
#Option 3: low colonization / high extinction -> decrease in diversity
#Base the answer on the outcomes of these 3 simulations shown at the bottom of this script.


######################################################
####Loading the data
###load xlsx data
#install.packages("readxl")
library(readxl)

insula_checklist <- read_excel("C:/Users/daank/OneDrive - University of Twente/Documents/Rijksuniversiteit Groningen/Year 1 - 24-25/Island Biology/Practicals/Practical 2 - DAISIE/Insula_checklist.xlsx") |>
  as_tibble()
insula_checklist

###Load tree data
insula_tree <- read.nexus("C:/Users/daank/OneDrive - University of Twente/Documents/Rijksuniversiteit Groningen/Year 1 - 24-25/Island Biology/Practicals/Practical 2 - DAISIE/Insula.tre")
insula_tree


######################################################
####DAISIEprep
###Create the data / phylogeny
#Define tip labels
insula_tree$tip.label

#Convert phylogeny to phylo4 format, to easily work with tip labels
phylo <- phylobase::phylo4(insula_tree)
phylobase::plot(phylo)

#Add distribution information to the tips of the tree
#h is nonendemic, f and e are endemic, all others are not found on the island
island_species <- data.frame(
  tip_labels = c("Spec_9",
                 "Spec_19",
                 "Spec_25",
                 "Spec_26",
                 "Spec_29",
                 "Spec_33",
                 "Spec_38",
                 "Spec_39",
                 "Spec_40",
                 "Spec_41",
                 "Spec_42",
                 "Spec_43",
                 "Spec_47",
                 "Spec_48",
                 "Spec_24"),
  tip_endemicity_status = c("endemic",
                            "endemic",
                            "endemic",
                            "endemic",
                            "endemic",
                            "endemic",
                            "endemic",
                            "endemic",
                            "endemic",
                            "endemic",
                            "endemic",
                            "endemic",
                            "endemic",
                            "endemic",
                            "nonendemic"))

#Add endemicity status to the tips of the tree
endemicity_status <- create_endemicity_status(
  phylo = phylo,
  island_species = island_species)

phylod <- phylobase::phylo4d(phylo, endemicity_status)

#Answer question B
plot_phylod(phylod = phylod)

###Extract data using DAISIEprep
#Create object to store island colonization events
island_tbl <- island_tbl()
island_tbl

#Fill the data frame with data
island_tbl_min <- extract_island_species(
  phylod = phylod,
  extraction_method = "min")

island_tbl_min

#Find probability of the ancestors of the island species being present on the island to determine the time of colonisation
phylod <- add_asr_node_states(phylod = phylod, asr_method = "mk")

#Plotting the phylogeny with indicated which ancestors were present on the island
plot_phylod(phylod = phylod)
#Optional probability plotting
plot_phylod(phylod = phylod, node_pies = TRUE)

#Extract island colonisation and diversification times from the phylogeny with asr using the reconstructed ancestral states of island presence/absence
island_tbl_asr <- extract_island_species(
  phylod = phylod,
  extraction_method = "asr")

island_tbl_asr

#Compare min and asr tables
all.equal(island_tbl_min,island_tbl_asr)

#Equal outcome, we will here use asr, so rename logically
island_tbl<-island_tbl_asr

#Visualize the table
View(island_tbl)

#For example, the code below shows the names of the species included in each lineage
island_tbl@island_tbl$species

###Adding missing species
##Adding an entire lineage when a phylogeny is not available for the lineage
#Adding lineage with 1 species
island_tbl <- add_island_colonist(
  island_tbl = island_tbl,
  clade_name = "Spec_51",
  status = "endemic",
  # clade with just 1 species, missing_species = 0
  # because adding the lineage already counts as 1
  missing_species = 0,
  col_time = NA_real_,
  col_max_age = FALSE,
  branching_times = NA_real_,
  min_age = NA_real_,
  clade_type = 1,
  species = "Spec_51")

##Adding missing species to an island clade that has been sampled in the phylogeny
island_tbl <- add_missing_species(
  island_tbl = island_tbl,
  # num_missing_species equals total species missing
  num_missing_species = 1,
  # name of a sampled species you want to "add" the missing to
  # it can be any in the clade
  species_to_add_to = "Spec_42")

##Prepare object for analyses in DAISIE
#Convert to data object to be used in DAISIE
data_list <- create_daisie_data(
  data = island_tbl,
  island_age = 5,
  num_mainland_species = 1000,
  precise_col_time = TRUE)
data_list[1]
data_list[2]

################################################
####DAISIE
DAISIE_plot_island(data_list)
title(ylab = "Time (Mya)", xlab = "Clades", main = "Species through time")

#Plot age versus diversity
DAISIE_plot_age_diversity(data_list)

##Fitting DAISIE models
#Fit a full DAISIE model (M1)
M1_results <- DAISIE_ML(
  datalist = data_list,
  initparsopt = c(5.8,7.7,20,0.02,3.2),
  ddmodel = 11,
  idparsopt = 1:5,
  parsfix = NULL,
  idparsfix = NULL)

M1_results
#  lambda_c       mu       K      gamma lambda_a    loglik df conv
#  5.90395 7.821783 2947851 0.02523945  3.19606 -37.10743  5    0

#Fit model with no carrying-capacity (M2) (since very high in M1, make it a free parameter in this model)
M2_results <- DAISIE_ML(
  datalist = data_list,
  initparsopt = c(5.8,7.7,0.02,3.2),
  idparsopt = c(1,2,4,5),
  parsfix = Inf,
  idparsfix = 3,
  ddmodel = 0)

M2_results
#  lambda_c       mu   K      gamma lambda_a    loglik df conv
#  5.904042 7.822047 Inf 0.02522024  3.19712 -37.10743  4    0

#Fit model with no carrying capacity AND no anagenesis (M3) (since anagenesis (labda_a) very low in M2)
M3_results <- DAISIE_ML(
  datalist = data_list,
  initparsopt = c(5.8,7.7,0.02),
  idparsopt = c(1,2,4),
  parsfix = c(Inf,0),
  idparsfix = c(3,5),
  ddmodel = 0)

M3_results
#  lambda_c       mu   K      gamma lambda_a    loglik df conv
#  6.76418 8.775854 Inf 0.02762404        0 -37.30016  3    0

#Select the best model using AIC
#Create new function to compute AIC values base on the likelihood values and number of parameters for each model
AIC_compare <- function(LogLik,k){
  aic <- (2 * k) - (2 * LogLik)
  return(aic)
}
#Fill in the values from the data
AICs <- AIC_compare(LogLik = c(M1_results$loglik,M2_results$loglik,M3_results$loglik),
                    k = c(M1_results$df,M2_results$df,M3_results$df))
names(AICs) <- c('M1','M2','M3')
AICs
#>    M1    M2    M3 
#> 84.21 82.21 80.60
#Lowest AIC is preferred -> M3 is the best model

##Simulate islands
#Simulate islands with the parameters estimated from the best model for the Galápagos bird data (takes a while to run, you can reduce number of replicates if you want)
jamaica_sims <- DAISIE_sim(
  time = 5,
  M = 1000,
  pars = c(6.76418, 8.775854, Inf, 0.02762404,0),
  replicates = 100,
  plot_sims = FALSE)

#Plot the species-through-time plots resulting from the simulations
DAISIE_plot_sims(jamaica_sims)

#Now lets try this for model 1, which includes diversity dependence
jamaica_sims_1 <- DAISIE_sim(
  time = 5,
  M = 1000,
  pars = c(5.90395,7.821783,2947851,0.02523945,3.19606),
  replicates = 100,
  plot_sims = FALSE)

#Plot the species-through-time plots resulting from the simulations
DAISIE_plot_sims(jamaica_sims_1)

#Try the model with high colonization and low extinction
jamaica_sims_growth <- DAISIE_sim(
  time = 5,
  M = 1000,
  pars = c(1, 0.02, Inf, 0.05,0),
  replicates = 100,
  plot_sims = FALSE)

#Plot the species-through-time plots resulting from the simulations
DAISIE_plot_sims(jamaica_sims_growth)

#Try the model with low colonization and high extinction
jamaica_sims_loss <- DAISIE_sim(
  time = 5,
  M = 1000,
  pars = c(1, 5, Inf, 0.02,0),
  replicates = 100,
  plot_sims = FALSE)

#Plot the species-through-time plots resulting from the simulations
DAISIE_plot_sims(jamaica_sims_loss)







