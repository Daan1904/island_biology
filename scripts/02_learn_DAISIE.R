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

##Load and visualise Galapágos bird data
load(file = "C:/Users/daank/OneDrive - University of Twente/Documents/Rijksuniversiteit Groningen/Year 1 - 24-25/Island Biology/Practicals/Practical 2 - DAISIE/galapagos_datalist.Rdata")

#View data list
galapagos_datalist
View(galapagos_datalist)
?Galapagos_datalist 













