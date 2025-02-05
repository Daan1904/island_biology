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
#

##What are the rates of colonisation, speciation and extinction for Insula in Jamaica?
#

##Is there evidence for diversity-dependence in the Insula species of Jamaica?
#

##Is there evidence for equilibrium dynamics on the island of Jamaica?
#

##How has the diversity of Insula on Jamaica changed through time (according to DAISIE simulations)?
#

##Modify one or two parameters of your choice and run simulations again under those parameters. How does this affect the simulations? Describe.
#


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
####DAISIE explanation
###Fitting DAISIE models using the Galápagos birds as example
##Load the required packages
library(DAISIE)
library(ape)

rm(list = ls())

##Load and visualise Galapágos bird data
#load(file = "C:/Users/daank/OneDrive - University of Twente/Documents/Rijksuniversiteit Groningen/Year 1 - 24-25/Island Biology/Practicals/Practical 2 - DAISIE/galapagos_datalist.Rdata")

#View data list
galapagos_datalist
View(galapagos_datalist)
?Galapagos_datalist 

galapagos_datalist[[1]]
#> $island_age
#> [1] 4
#> 
#> $not_present
#> [1] 992

#To view just the Mimus colonisation
galapagos_datalist[[4]]
#> $colonist_name
#> [1] "Mimus"
#> 
#> $branching_times
#> [1] 4.00000 3.99999 3.68000 2.93000 0.29000
#> 
#> $stac
#> [1] 6
#> 
#> $missing_species
#> [1] 0
#> 
#> $type1or2
#> [1] 1

#Visualise Galápagos data
#Reset the margins of the plot window, this gives the default options
dev.off()

DAISIE_plot_island(galapagos_datalist)

#Plot age versus diversity
DAISIE_plot_age_diversity(galapagos_datalist)

##Fitting DAISIE models
#Fit a full DAISIE model (M1)
M1_results <- DAISIE_ML(
  datalist = galapagos_datalist,
  initparsopt = c(1.5,1.1,20,0.009,1.1),
  ddmodel = 11,
  idparsopt = 1:5,
  parsfix = NULL,
  idparsfix = NULL)

M1_results
#>   lambda_c       mu       K       gamma     lambda_a    loglik df conv
#> 1 1.258226 1.136924 9968506 0.004957378 1.251793e-06 -84.78145  5    0

#Fit model with no carrying-capacity (M2) (since very high in M1, make it a free parameter in this model)
M2_results <- DAISIE_ML(
  datalist = galapagos_datalist,
  initparsopt = c(1.5,1.1,0.009,1.1),
  idparsopt = c(1,2,4,5),
  parsfix = Inf,
  idparsfix = 3,
  ddmodel = 0)

M2_results
#>   lambda_c       mu   K      gamma     lambda_a    loglik df conv
#> 1 1.264389 1.149378 Inf 0.00505558 1.662578e-05 -84.78088  4    0

#Fit model with no carrying capacity AND no anagenesis (M3) (since anagenesis (labda_a) very low in M2)
M3_results <- DAISIE_ML(
  datalist = galapagos_datalist,
  initparsopt = c(1.5,1.1,0.009),
  idparsopt = c(1,2,4),
  parsfix = c(Inf,0),
  idparsfix = c(3,5),
  ddmodel = 0)

M3_results
#>   lambda_c       mu   K       gamma lambda_a    loglik df conv
#> 1 1.263034 1.146225 Inf 0.005040353        0 -84.78082  3    0

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
#>       M1       M2       M3 
#> 179.5629 177.5618 175.5616
#Lowest AIC is preferred -> M3 is the best model

##Simulate islands
#Simulate islands with the parameters estimated from the best model for the Galápagos bird data (takes a while to run, you can reduce number of replicates if you want)
Galapagos_sims <- DAISIE_sim(
  time = 4,
  M = 1000,
  pars = c(1.26, 1.146, Inf, 0.005,0),
  replicates = 100,
  plot_sims = FALSE)

#Plot the species-through-time plots resulting from the simulations
DAISIE_plot_sims(Galapagos_sims)







