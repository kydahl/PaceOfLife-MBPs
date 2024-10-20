## Kyle Dahlin, University of Georgia, kyle.dahlin@uga.edu
## Started Oct 2020
##
## Title: Trait data analysis
##
## Project: Global zoonoses: pace-of-life and mosquito-borne pathogens
##
## Purpose: Load host trait databases and derive "pace hyper-parameters"
##
## Contents: 1) Set-up,load packages and data
##           2) Estimate pace hyper-parameters from trait data
##           3) Output hyper-parameters to be used by 'trait-functions.R'
##
## Inputs:
##          Host trait database: "data/clean/primate_imputed.csv"
##
##
## Outputs:
##          Trait hyper-parameters: z0 and c for lambda, mu, and K
##
##
## ______________________________________________________________________________

# _____________________________________________
# 1) Set-up,load packages, get data, etc. ----
# _____________________________________________
###* Load Libraries ----
library(tidyverse)
library(reshape2)
library(cowplot)
library(ggpmisc)
library(broom)
library(ggfortify)

### Functions for fitting parameters to pace ----

# Function: get coefficients (z0, c) for the best fit curve of the model
#     log(parameter) = z0 + c * pace
extract_pace_fit <- function(input, var_name, measure_ID) {
  with(as.list(input), {
    df = filter(input, var == var_name, Measure == measure_ID) %>% 
      arrange(value)
    pace.vec <- seq(0, 1, length.out = dim(df)[1])
    # Assume that log(parameter) = z0 + c * pace
    var.fit = lm(log(value) ~ pace, data = df)
    return(as.vector(c(var.fit$coefficients,summary(var.fit)$r.squared)))
  })
}

# Function: create a curve illustrating the parameter as a function of 
# pace of life history
create_pace_curve <- function(input, var_name, measure_ID) {
  with(as.list(input), {
    coeffs = extract_pace_fit(input, var_name, measure_ID)
    curve = exp(coeffs[1] + coeffs[2] * pace)
  })
}

### Function: plot fit of each parameter to the pace of life measure ----
fitplot_func <- function(in_df, Measure_ID) {
  with(as.list(in_df), {
    # specify the pace of life history measure
    df <- filter(in_df, Measure == Measure_ID)
    
    # make separate dataframe with R2 values
    R2_df <- df %>% 
      select(Measure, var, R2) %>% 
      mutate(R2 = sprintf("R^2 ==%.2f", R2)) %>% 
      distinct()
    
    df %>% 
      ggplot(mapping = aes(x = pace, y = (value)), group = c(var, Family)) +
      # plot empirical values as points
      geom_point() +
      # plot fitted curve in blue
      geom_line(aes(x = pace, y = (curve)), lwd = 2, color = "blue") +
      # label R^2 value
      geom_label(data = filter(R2_df, Measure == Measure_ID), aes(label = R2),
                 x = -Inf, y = Inf, hjust = 0, vjust = 1.5, parse = TRUE) +
      xlab("Pace of life history") +
      facet_wrap(Family ~ var, scales = "free", strip.position = "left",
                  labeller = as_labeller(c(
                    K = "log(Population density)",
                    lambda = "log(Recruitment rate)",
                    mu = "log(Mortality rate)"
                  ))) +
      scale_y_continuous(trans = "log10") +
      ylab(NULL) +
      theme_half_open(12) +
      theme(
        strip.background = element_blank(),
        strip.placement = "outside"
      )
  })
}


# _______________________________________________________
# 2) Estimate pace hyper-parameters from trait data ----
# _______________________________________________________

## Notes:
### There are no hard and fast rules for how host traits should vary with the "pace of life history" as,
### at the end of the day, this is just a "made-up" quantity. That said, for our analysis, a linear or
### exponential relationship would allow for the greatest potential for mathematical tractibility.
### Therefore, we start by looking at the trait data directly to evaluate whether these functional forms
###  are somewhat (read qualitatively) close to the actual way that host traits vary across species.
### To do this, we first define the "pace" metric, then order species by this metric. We will fit each
### function type to these data to verify (in an ad-hoc manner) that they reasonable represent the
### variation.


### a) Read in primate traits ----
## The following traits should increase with pace of life history:
## recruitment rate (lambda) and mortality rate (mu, inverse of lifespan)
## Carrying capacity (K) might also increase with pace of life history
## Immunological traits will also vary with pace, see './code/trait-functions.R'

HyperParms <- tibble(variable = double(), z0 = double(), c = double(),
                     Family = character(), Measure = character())

FitsTable <- tibble(Spp = double(), Measure = c(), pace = c(), var = c(), 
                    value = c(), coeff1 = c(), coeff2 = c(), R2 = c(), curve = c(),
                    Family = c())


for (family in c("rodent", "primate")) {
  
  if (family == "rodent") {
    base_data <- read_csv("data/clean/rodent_taxa_imputed.csv", col_types = cols())
  } else if (family == "primate") {
    base_data <- read_csv("data/clean/primate_taxa_imputed.csv", col_types = cols())
  }
  
  # Load in imputed primate trait data
  Traits <- base_data %>% 
    # Remove the only human primate
    mutate(Spp = MSW05_Binomial) %>%
    # Rename column headers for easier reference
    mutate(AdultMass = `X5.1_AdultBodyMass_g`) %>%
    mutate(LitterSize = `X15.1_LitterSize`) %>%
    mutate(LittersPerDay = `X16.1_LittersPerYear` / 365.25) %>% # convert from per year to per day
    mutate(Longevity = `X17.1_MaxLongevity_m` * 365.25 / 12) %>% # convert from months to days
    mutate(PopulationDensity = `X21.1_PopulationDensity_n.km2` / 100) %>% # convert from square kilometers to hectares
    mutate(InterbirthInterval = `X14.1_InterbirthInterval_d`) %>%
    mutate(SexualMaturityAge = `X23.1_SexualMaturityAge_d`) %>%
    mutate(AgeatFirstBirth = `X3.1_AgeatFirstBirth_d`) %>% 
    # Select only relevant traits
    select(
      Spp, AdultMass, LitterSize, LittersPerDay, Longevity,
      PopulationDensity, InterbirthInterval, SexualMaturityAge, AgeatFirstBirth
    )
  
  ### b) Translate traits to life history parameters used in modeling ----
  Parameters <- Traits %>%
    # mortality rate (per day)
    mutate(mu = 1 / Longevity) %>%
    # population density (indviduals per hectare)
    mutate(K = PopulationDensity) %>%
    # recruitment rate (per day)
    mutate(lambda = LitterSize * LittersPerDay) %>%
    # Generation time taken to be the "doubling time" for a growth rate defined as the difference between recruitment and mortality
    mutate(GenerationTime = log(2) / (lambda - mu)) %>% 
    mutate(SexualMaturityPeriod = Longevity - SexualMaturityAge) %>%
    # maximum number of reproductive events in a single lifetime
    mutate(MaxNumReps = SexualMaturityPeriod / InterbirthInterval) %>%
    # Maximum Lifetime Reproductive Output (MLRO)
    mutate(MLRO = LitterSize * MaxNumReps) %>%
    mutate(MLROperMass = MLRO / AdultMass) %>%
    # Normalized recruitment rate
    mutate(norm.lambda = percent_rank(lambda)) %>%
    # Normalized mortality rate
    mutate(norm.mu = percent_rank(mu)) %>%
    # Normalized population density
    mutate(norm.K = percent_rank(K))
  
  ### c) Compute pace of life history measures----
  # See below for complete definition of measures
  Measures <- c("1","2","3","4")
  
  # correct 'direction' of PCA used for pace of life history
  # we choose for pace of life history to be the direction of increasing reproduction and decreasing longevity
  family_mult <- ifelse(family == "rodent", -1, 1)
  
  # PCA measure for pace of life history
  PCA_df <- base_data %>% 
    # Remove constant columns
    select(-c(`...1`,`X2.1_AgeatEyeOpening_d`, `X13.3_WeaningHeadBodyLen_mm`)) %>%
    select(where(is.numeric)) %>% 
    prcomp(., scale.=TRUE) %>%
    augment(Traits) %>% 
    mutate(Measure = 5) %>% 
    mutate(pace = percent_rank(family_mult * .fittedPC1)) %>% 
    right_join(Parameters, by = "Spp") %>% 
    select(Spp, Measure, pace, lambda, mu, K)
  
  # Compute pace of life history measures for each species
  PaceTable <- expand_grid(Parameters, Measure = Measures) %>%
    mutate(pace = case_when(
      # Measure 1 = generation time
      Measure == 1 ~ percent_rank(-GenerationTime),
      # Measure 2 = maximum lifetime reproductive output
      Measure == 2 ~ percent_rank(MLRO),
      # Measure 3 = maximum lifetime reproductive output per unit mass
      Measure == 3 ~ percent_rank(MLROperMass),
      # Measure 4 = sum of all life history traits considered in Van de Walle et al 2023, Proceedings B
      Measure == 4 ~ percent_rank((percent_rank(AgeatFirstBirth^(-1)) + percent_rank(Longevity) + percent_rank(MaxNumReps / Longevity) + percent_rank(LittersPerDay)))
    )
    ) %>% 
    select(Spp, Measure, pace, lambda, mu, K) %>% 
    rbind(PCA_df)
  
  
  # Pace measure 0: generation time
  # average age at first birth
  
  # Pace measure 1: maximum lifetime reproductive output per unit mass
  # this definition has an empirical basis
  
  # Pace measure 2: 2-norm of normalized reproductive rate and mortality rate
  # this definition comes from assuming that pace of life history is simply a
  # description of the relative rates of recruitment and mortality
  # "live fast and die hard"
  
  # Pace measure 3: 2-norm of normalized reproductive, mortality rates and pop. density
  # this definition is similar to Pace measure 2, with the additional assumption
  # that population density increases with pace of life history
  
  #________________________________________________________
  # 3) Fit parameters to pace of life history measures ----
  #________________________________________________________
  
  ### Fit parameters as functions of pace of life history ----
  tempFits <- PaceTable %>%
    pivot_longer(cols = c("lambda", "mu", "K"), names_to = "var") %>%
    group_by(var, Measure) %>%
    # c value
    mutate(coeff1 = extract_pace_fit(., var, Measure)[1]) %>% 
    # z0 value
    mutate(coeff2 = extract_pace_fit(., var, Measure)[2]) %>% 
    # R2 of the fit
    mutate(R2 = extract_pace_fit(., var, Measure)[3]) %>%
    # value of the parameter curve for the specified value of pace
    mutate(curve = exp(coeff1 + coeff2 * pace)) %>% 
    mutate(Family = family)
  
  FitsTable <- rbind(FitsTable, tempFits)
  #__________________________________________________________________
  # 4) Output hyper-parameters to be used by 'trait-functions.R' ----
  #__________________________________________________________________
  
  # Output coefficients as named variables (to be added to the appropriate table later)
  
  temp_df <- tempFits %>% 
    # Choose the measure we think gives the "best fit"
    # See './code/compare-pace-measures.R' for comparison analysis
    filter(Measure == 1) %>% 
    rename(variable = var) %>%
    select(variable, coeff1, coeff2) %>% 
    distinct() %>%
    mutate(z0 = exp(coeff1), 
           c = coeff2,
           Family = family,
           .keep = 'unused')
  
  HyperParms <- rbind(HyperParms, temp_df)
  
}
# Save fits and hyperparameters
write_csv(FitsTable, "data/clean/fits-table.csv")
write_csv(HyperParms, "data/clean/trait-hyperparameters.csv")