## Kyle Dahlin, University of Georgia, kydahlin@gmail.com
## Started Oct 2020
##
## Title: Trait functions
##
## Project: Global zoonoses, pace-of-life and mosquito-borne pathogens
##
## Purpose: Functions for evaluating trait values for hosts, mosquitoes, and 
##          pathogens
##
## Contents: 1) Set-up,load packages, get data, etc.
##           2) Build parameter sets for host, mosquito, and pathogen
##              a) Define parameter-functions (of appropriate variable)
##              b) Build parameter table
##
## Inputs:  PaceHyperParameters:
##            for each host group considered, a set of hyperparameters which 
##            describe the range of values taken for each life history trait 
##            (recruitment rate, mortality rate, carrying capacity, recovery 
##            period, susceptibility to infection)
##
##          MosquitoThermalParameters:
##            for the three mosquito species, defines the functional form and 
##            relevant parameters describing how traits vary with temperature
##
## Outputs: getHostParms(p):
##            get recruitment (lambda_H), mortality (mu_H), carrying capacity (K_H), 
##            recovery period (Gamma), and susceptbility (beta_H) of hosts as a 
##            function of a pace-of-life (p) between 0 and 1
##
##          getVecParms(SpeciesName,T):
##             get maximum biting frequency (sigma_V), fecundity (f), probability of egg
##             survival to adulthood (delta_L), immature mosquito development rate (rho_L), and
##             adult mosquito mortality rate (mu_V) as a function of SpeciesName (AE = Aedes 
##             aegypti), (AL = Aedes albopictus), (CQ = Culex quinquefasciatus) and temperature
##             (T in Celsius)
##          getPathParms(N):
##             get in-mosquito incubation rate (eta_V), host recovery rate (gamma_H), mosquito-
##             to-host transmission probability (beta_H) and host-to-mosquito transmission
##             probability (beta_V) as a function of sample number (N), where the total number 
##             of samples is given by the value of Path.NumberSamples
##______________________________________________________________________________

#_____________________________________________
# 1) Set-up,load packages, get data, etc. ----
#_____________________________________________
# Load Libraries ----
library(tidyverse)
library(reshape2)
library(cowplot)

#_____________________________
# 2) Build parameter sets ----
#_____________________________

###* Host Trait functions ----

# Function representing a trait varying exponentially with pace
## z0 is the minimum value (when s=0), and c is the logarithmic slope
trait.exp.function <- function(s,z0,c) {
  z0*exp(c*s)
}

# Get hyper-parameters for trait functions of pace with supporting data

HyperParms <- read_csv("data/clean/trait-hyperparameters.csv") %>%
  filter(Family == family_in)

lambda.z0 <- filter(HyperParms, variable == "lambda")$z0
lambda.c <- filter(HyperParms, variable == "lambda")$c
mu.z0 <- filter(HyperParms, variable == "mu")$z0
mu.c <- filter(HyperParms, variable == "mu")$c
K.z0 <- filter(HyperParms, variable == "K")$z0
K.c <- filter(HyperParms, variable == "K")$c

# Host parameter functions
## recruitment rate (lambda)
get.Host.lambda <- function(pace) {trait.exp.function(pace,lambda.z0,lambda.c)}

## mortality rate (mu)
get.Host.mu <- function(pace) {trait.exp.function(pace,mu.z0,mu.c)}

## carrying capacity (K)
get.Host.K <- function(pace) {trait.exp.function(pace,K.z0,K.c)}

## infectious period as proportion of lifespan (alpha) [pathogen dependent]
# basepoint
get.Host.alpha.z0 <- function(input) {
  with(input, {
    alpha.z0 <- case_when(
      (Pathogen_type == "Chronic") ~ 0.01, # 5% of host lifespan
      (Pathogen_type == "Acute") ~ 2*mu.z0) # 2 days
  }
  )
}
# exponential rate of increase
get.Host.alpha.c <- function(input) {
  with(input, {
    alpha.c <- case_when(
      (Pathogen_type == "Chronic" & Immune_variation == "High") ~ log(0.05/0.01), # range of 1% - 5% of host lifespan
      (Pathogen_type == "Chronic" & Immune_variation == "Low") ~ log(0.03/0.01), # range of 1% - 3% of host lifespan
      (Pathogen_type == "Chronic" & Immune_variation == "None") ~ 0, # range of 5% - 5% of host lifespan
      (Pathogen_type == "Acute" & Immune_variation == "High") ~ log((12*mu.z0*exp(mu.c))/(2*mu.z0)), # range of 2 days to 2 weeks
      (Pathogen_type == "Acute" & Immune_variation == "Low") ~ log((7*mu.z0*exp(mu.c))/(2*mu.z0)), # range of 2 days to 1 weeks
      (Pathogen_type == "Acute" & Immune_variation == "None") ~ 0 # range of 2 days to 2 days
    )
  }
  )
}
# function
get.Host.alpha <- function(input) {
  with(input, {
    alpha.z0 <- get.Host.alpha.z0(input)
    alpha.c <- get.Host.alpha.c(input)
    
    alpha.z0*exp(alpha.c*Pace)
  }
  )
}

## recovery rate (gamma = mu/alpha - mu)
get.Host.gamma <- function(pace) {get.Host.mu(pace)/get.Host.alpha(pace) - get.Host.mu(pace)}

## susceptibility (beta) [pathogen dependent]
# basepoints
get.Host.beta.z0 <- function(input) {
  with(input, {
    beta.z0 <- case_when(
      (Pathogen_transmissibility == "High") ~ 0.6,
      (Pathogen_transmissibility == "Low") ~ 0.05)
  }
  )
}
# exponential rates of increase
get.Host.beta.c <- function(input) {
  with(input, {
    beta.c <-  case_when(
      # these are chosen to span the same breadth of probabilities (40%) but at 
      # different magnitudes
      (Pathogen_transmissibility == "High" & Immune_variation == "High") ~ log(1/0.6), 
      (Pathogen_transmissibility == "High" & Immune_variation == "Low") ~ log(0.8/0.6),  
      (Pathogen_transmissibility == "High" & Immune_variation == "None") ~ 0, 
      (Pathogen_transmissibility == "Low" & Immune_variation == "High") ~ log(0.45/0.05), 
      (Pathogen_transmissibility == "Low" & Immune_variation == "Low") ~ log(0.25/0.05),
      (Pathogen_transmissibility == "Low" & Immune_variation == "None") ~ 0
    )
  }
  )
}
# function
get.Host.beta <- function(input) {
  with(input, {
    beta.z0 <- get.Host.beta.z0(input)
    beta.c <- get.Host.beta.c(input)
    
    beta.z0*exp(beta.c*Pace)
  }
  )
}

###* Produce table of all host parameters as a function of pace-of-life ----
# outputs dataframe of Host parameter values for the given values of 'pace'
getHostParms <- function(pace) {
  tibble(Pace = pace,lambdaH = get.Host.lambda(pace), muH = get.Host.mu(pace), 
         KH = get.Host.K(pace))
}

###* Mosquito traits ----
# See './code/get-mosquito-traits.R
