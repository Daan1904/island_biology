renv::restore()

library(DAISIEprep)
library(ape)
################################################
####DAISIEprep explanation
###Create the data / phylogeny
#Make a phylogeny
set.seed(
  4,
  kind = "Mersenne-Twister",
  normal.kind = "Inversion",
  sample.kind = "Rejection")
phylo <- ape::rcoal(10)

#Define tip labels
phylo$tip.label <- c("Plant_a", "Plant_b", "Plant_c", "Plant_d", "Plant_e",
                     "Plant_f", "Plant_g", "Plant_h", "Plant_i", "Plant_j")

#Convert phylogeny to phylo4 format, to easily work with tip labels
phylo <- phylobase::phylo4(phylo)
phylobase::plot(phylo)

#Add distribution information to the tips of the tree
#h is nonendemic, f and e are endemic, all others are not found on the island
island_species <- data.frame(
  tip_labels = c("Plant_h",
                 "Plant_f",
                 "Plant_e"),
  tip_endemicity_status = c("nonendemic","endemic","endemic"))

#Add endemicity status to the tips of the tree
endemicity_status <- create_endemicity_status(
  phylo = phylo,
  island_species = island_species)

phylod <- phylobase::phylo4d(phylo, endemicity_status)

plot_phylod(phylod = phylod)

###Extract data using DAISIEprep
#Create object to store island colonization events
island_tbl <- island_tbl()
island_tbl
#> Class:  Island_tbl 
#> [1] clade_name      status          missing_species col_time       
#> [5] col_max_age     branching_times min_age         species        
#> [9] clade_type     
#> <0 rows> (or 0-length row.names)

#Fill the data frame with data
island_tbl_min <- extract_island_species(
  phylod = phylod,
  extraction_method = "min")

island_tbl_min
#> Class:  Island_tbl 
#>   clade_name     status missing_species  col_time col_max_age branching_times
#> 1    Plant_e    endemic               0 0.1925924       FALSE    0.154371....
#> 2    Plant_h nonendemic               0 0.1233674       FALSE              NA
#>   min_age      species clade_type
#> 1      NA Plant_e,....          1
#> 2      NA      Plant_h          1

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
#> Class:  Island_tbl 
#>   clade_name     status missing_species  col_time col_max_age branching_times
#> 1    Plant_e    endemic               0 0.1925924       FALSE    0.154371....
#> 2    Plant_h nonendemic               0 0.1233674       FALSE              NA
#>   min_age      species clade_type
#> 1      NA Plant_e,....          1
#> 2      NA      Plant_h          1

#Compare min and asr tables
all.equal(island_tbl_min,island_tbl_asr)
#> [1] TRUE

#Equal outcome, we will here use asr, so rename logically
island_tbl<-island_tbl_asr

#Visualize the table
View(island_tbl)

#For example, the code below shows the names of the species included in each lineage
island_tbl@island_tbl$species
#> [[1]]
#> [1] "Plant_e" "Plant_f"
#> 
#> [[2]]
#> [1] "Plant_h"

###Adding missing species
##Adding missing species to an island clade that has been sampled in the phylogeny
#Find out which species are already sampled in each lineage using
island_tbl@island_tbl$species
#> [[1]]
#> [1] "Plant_e" "Plant_f"
#> 
#> [[2]]
#> [1] "Plant_h"

island_tbl <- add_missing_species(
  island_tbl = island_tbl,
  # num_missing_species equals total species missing
  num_missing_species = 2,
  # name of a sampled species you want to "add" the missing to
  # it can be any in the clade
  species_to_add_to = "Plant_e")

#The new island table now has missing species added to the lineage we wanted to
island_tbl@island_tbl$missing_species
#> [1] 2 0

#With the new missing species added to the island_tbl we can repeat the conversion steps above using create_daisie_data() to produce data accepted by the DAISIE model
data_list <- create_daisie_data(
  data = island_tbl,
  island_age = 12,
  num_mainland_species = 100,
  precise_col_time = TRUE)

##Adding an entire lineage when a phylogeny is not available for the lineage
#Example for adding lineage with 1 species
island_tbl <- add_island_colonist(
  island_tbl = island_tbl,
  clade_name = "Plant_y",
  status = "endemic",
  # clade with just 1 species, missing_species = 0
  # because adding the lineage already counts as 1
  missing_species = 0,
  col_time = NA_real_,
  col_max_age = FALSE,
  branching_times = NA_real_,
  min_age = NA_real_,
  clade_type = 1,
  species = "Plant_a")

#Example for adding lineage with 5 species (“Plant_a”, “Plant_b”, “Plant_c”, “Plant_d”, “Plant_e”) to a lineage called “Plant_radiation”
island_tbl <- add_island_colonist(
  island_tbl = island_tbl,
  clade_name = "Plant_radiation",
  status = "endemic",
  # the total species is 5 and all are missing
  # but we add missing_species = 4 because
  # adding the lineage already counts as 1
  missing_species = 4,
  col_time = NA_real_,
  col_max_age = FALSE,
  branching_times = NA_real_,
  min_age = NA_real_,
  clade_type = 1,
  species = c("Plant_a", "Plant_b", "Plant_c",
              "Plant_d", "Plant_e"))

##Prepare object for analyses in DAISIE
#Convert to data object to be used in DAISIE
data_list <- create_daisie_data(
  data = island_tbl,
  island_age = 12,
  num_mainland_species = 100,
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

#Comparison with bird data from the Azores islands
data(Macaronesia_datalist)
Azores <- Macaronesia_datalist[[1]]

#Visualise Azores data
DAISIE_plot_island(Azores)

#Three of the species are extinct and only known from fossils
#Simulate Azores with pre-identified ML parameters
Azores_sims <- DAISIE_sim(
  time = 6.3,
  M = 300,
  pars = c(0,1.053151832,Inf,0.052148979,0.512939011),
  replicates = 100,
  plot_sims = FALSE)

#Plot the species-through-time plot for the Azores from the resulting from the simulations
DAISIE_plot_sims(Azores_sims)










