## For questions, contact Kyle Dahlin, kydahlin@gmail.com
## Originated January 2023
##
## Title: Analysis code for "Pace of life and Mosquito-borne pathogens" ########
##
## Project: Global zoonoses - spillover of mosquito-borne pathogens
##
## Purpose: Build data frames for: 1) vector traits, 2) all traits and model
##          outputs, and 3) thermal characteristics of transmission
##
## Contents: 0) Load in necessary packages, functions, settings
##           1) Build vector trait data frame
##           2) Build all traits and model outputs data frame
##           3) Build thermal characteristics data frame
##
##
## Settings:  Options for reducing the resolution of variables for memory
##            allocation and figure plotting purposes
##
##
## Outputs: (in ./results)
##          1) VectorTraits.csv - vector trait data frame
##          2) AllOutputs.csv - all traits and model outputs data frame
##          3) ThermalCharacteristics.csv - thermal characteristics data frame
##

# 0) Load in necessary packages, functions, settings ###########################
# Packages
library(tidyverse)
library(foreach)
library(doParallel)
library(progress)

# Functions for computing transmission measures
source("code/output-functions.R")

# Function: Slice data to optimize usage of parallel processing and memory
slice_func<-function(x,n) {
  N<-length(x);
  lapply(seq(1,N,n),function(i) x[i:min(i+n-1,N)])
}

###* Settings-------------------------------------------------------------------

# Resolution settings: (set these to Inf to keep original resolution)
# Temperature: Factor by which to reduce resolution
Temp_vec_length <- 301

# Pace of life history: 
pace_vec_length <- 100

pace_seq <- seq(0, 1, length.out = pace_vec_length)

# Biting tolerance values: values to consider for biting tolerance
sigmaH_value <- 100
sigmaH_vec <- c(10, 100, Inf)

# 1) Build vector trait data frame #############################################

## Assign names to mosquitoes and pathogens
# Define mosquito species names and codes:
mosquito_names <- c(
  "Aedes aegypti", "Aedes albopictus",
  "Culex quinquefasciatus", "Anopheles spp."
)
mosquito_code <- c('AE','AL','CQ', 'AN')
names(mosquito_names) <- mosquito_code

###* Build table of mosquito parameters ----

## NB: Carrying capacity for larval mosquitoes
# In the absence of good estimates for each species or temperature-dependence of
# this trait, we assume that this parameter is constant. It can be used as a 
# scaling parameter for overall mosquito abundance 
# (it could alternately be used to fix the maximum adult mosquito density across species)
# This value ensures that vector-host ratio is bounded between approx. 10 and 100, which is a commonly-assumed range
larval_mosquito_carrying_capacity <- 10

# Get mosquito, pathogen trait thermal performance curve (TPC) data from:
#   "results/MosquitoThermalResponse.csv"
VectorTraits_df <- read_rds("data/clean/parameter_TPCs.rds") %>%
  # Lower the resolution of temperature to length Temp_vec_length
  filter(Temperature %in% res_reduce(Temperature, Temp_vec_length)) %>%
  # Larval mosquito carrying capacity
  mutate(KL = larval_mosquito_carrying_capacity) %>%
  # Average adult lifespan
  mutate(muV = 1 / lf) %>%
  # Compute equilibrium mosquito abundance (V0) from other traits
  mutate(V0 = compute.V0(.)) %>% 
  ungroup() 

### * Save data frame to use for trait TPCs figure ----
write_rds(VectorTraits_df, "results/VectorTraits.rds", compress = 'gz')

# 2) Build pathogen trait data frame ###########################################

# Define pathogen types
Pathogen_transmissibility <- c("High","Low")
Pathogen_type <- c("Chronic","Acute")
Immune_variation = c("High", "Low", "None")
PathogenTraits_df <- expand.grid(Pathogen_transmissibility = Pathogen_transmissibility,
                                 Pathogen_type = Pathogen_type,
                                 Immune_variation = Immune_variation)

# 3) Build host trait data frame ###############################################
# Pace of life history (POLH) settings
family_in = "rodent"
source("code/trait-functions.R")

HostTraits_df <- expand_grid(getHostParms(pace_seq), 
                             sigmaH = sigmaH_vec) %>% 
  mutate(Model = if_else(is.finite(sigmaH), "CCH", "RM"),
         K.c = K.c,
         mu.c = mu.c,
         Family = family_in) %>% 
  filter(sigmaH %in% c(10, Inf))

family_in = "primate"
source("code/trait-functions.R")
HostTraits_df <- rbind(HostTraits_df, 
                       expand_grid(getHostParms(pace_seq), 
                                   sigmaH = sigmaH_vec) %>% 
                         mutate(Model = if_else(is.finite(sigmaH), "CCH", "RM"),
                                K.c = K.c,
                                mu.c = mu.c,
                                Family = family_in) %>%
                         filter(sigmaH %in% c(100, Inf)))

# 4) Build all traits and model outputs data frame #############################

###* Build table with all parameter combinations ----

# Combine into table
HostPathogen_df <- expand_grid(HostTraits_df, 
                               PathogenTraits_df) %>%
  # compute susceptibility as a function of Pathogen_type and PoLH
  mutate(betaH = get.Host.beta(.)) %>%
  # compute infectious period as a function of Pathogen_type and PoLH
  mutate(alphaH = get.Host.alpha(.)) %>%
  # compute recovery rate from infectious period  and mortality rate
  mutate(gammaH = muH/alphaH - muH) 

# Initialize data frame
AllOutputs_df <- tibble(mosquito_species = c(), pathogen = c(), 
                        Temperature = c(), Pace = c(), Model = c(),
                        Pathogen_transmissibility = c(), Pathogen_type = c(),
                        Immune_variation = c(), 
                        sigmaH = c(), lambda_R0 = c(), V0 = c(), R0 = c(),
                        Family = c())

VectorTraits_df <- ungroup(VectorTraits_df)
HostPathogen_df <- ungroup(HostPathogen_df)

# Start new cluster for doParallel
cluster_size <- parallel::detectCores()-1
my.cluster <- parallel::makeCluster(cluster_size, type = "PSOCK"
)
# Register cluster for doParallel
doSNOW::registerDoSNOW(cl = my.cluster)

# Set up iteration grid
Pace_slices = slice_func(unique(HostPathogen_df$Pace), 10)
Temp_slices = slice_func(unique(VectorTraits_df$Temperature), 10)

iter_grid = expand_grid(Pace = Pace_slices, Temperature = Temp_slices)

# Set up progress bar
iterations <- dim(iter_grid)[1]

# Collect R0, V0, lambda_R0 across systems and host trait values
file_index = 1
for (system_name in unique(VectorTraits_df$system_ID)) {
  print(paste0("(",file_index, "/", length(unique(VectorTraits_df$system_ID)),") R0/V0/lambdaR0: ", system_name))
  pb <- progress_bar$new(
    format = ":spin progress = :percent [:bar] elapsed: :elapsed | eta: :eta",
    total = iterations,
    width = 100)                                                                                                                         
  progress <- function(n){pb$tick()}
  opts <- list(progress = progress)            
  
  
  
  AllOutputs_df <- foreach(index =  1:dim(iter_grid)[1],
                           .packages = "tidyverse",
                           .combine = 'rbind',
                           .export = c("mu.z0"),
                           .options.snow = opts) %dopar% {
                             index_Pace = iter_grid$Pace[index]
                             index_Temp = iter_grid$Temperature[index]
                             
                             VectorTraits_df %>% 
                               dplyr::filter(system_ID == system_name) %>% 
                               dplyr::filter(Temperature %in% index_Temp[[1]]) %>% 
                               expand_grid(dplyr::filter(HostPathogen_df, Pace %in% index_Pace[[1]])) %>% 
                               ## Compute outputs ##
                               # Basic reproduction number
                               mutate(R0 = sqrt(compute.RH(.) * compute.RV(.))) %>%
                               # # Total mosquito abundance
                               mutate(V0 = compute.V0(.)) %>%
                               # monotonicity of R0 wrt POLH
                               mutate(lambda_R0 = compute.lambdaR0(.)) %>% 
                               pivot_longer(cols = c(R0, V0, lambda_R0), 
                                            names_to = "variable", 
                                            values_to = "value") %>% 
                               select(system_ID, Temperature, Pace, Model, 
                                      Pathogen_transmissibility, 
                                      Pathogen_type, Immune_variation, 
                                      sigmaH, Family, sample_num, variable, value) %>% 
                               ungroup() %>% 
                               group_by(system_ID, Temperature, Pace, Model, 
                                        sigmaH, Pathogen_transmissibility, 
                                        Pathogen_type, Immune_variation, Family, 
                                        variable) %>% 
                               # partition(cluster) %>%
                               summarise(
                                 lowHCI = quantile(value, 0.055),
                                 highHCI = quantile(value, 0.945),
                                 mean = mean(value),
                                 median = median(value)
                               ) %>% 
                               # collect() %>%
                               # separate system_ID into mosquito species and pathogen labels
                               separate(system_ID, c("mosquito_species", "pathogen"), sep = " / ") %>% 
                               ungroup() %>% 
                               distinct()
                           }
  pb$terminate()
  
  write_rds(AllOutputs_df, paste0("results/AllOutputs_", file_index,".rds"))
  rm(AllOutputs_df)
  gc()
  AllOutputs_df <- tibble(mosquito_species = c(), pathogen = c(), 
                          Temperature = c(), Pace = c(), Model = c(),
                          Pathogen_transmissibility = c(), Pathogen_type = c(),
                          Immune_variation = c(), Family = c(),
                          sigmaH = c(), lambda_R0 = c(), V0 = c(), R0 = c())
  
  file_index <- file_index + 1
}
pb$terminate()

# Combine data frames for the species
for (file_index in 1:5) {
  file_name <- paste0("results/AllOutputs_", file_index,".rds")
  
  AllOutputs_df <- rbind(AllOutputs_df, 
                         read_rds(file_name))
  file.remove(file_name)
}

###* Save data frame ----

write_rds(ungroup(AllOutputs_df) %>% 
            filter(Immune_variation == "None"), "results/AllOutputs_none.rds", compress = 'gz')

write_rds(ungroup(AllOutputs_df) %>% 
            filter(Immune_variation == "Low"), "results/AllOutputs_low.rds", compress = 'gz')

write_rds(ungroup(AllOutputs_df) %>% 
            filter(Immune_variation == "High"), "results/AllOutputs_high.rds", compress = 'gz')
