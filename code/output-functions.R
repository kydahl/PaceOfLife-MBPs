## Kyle Dahlin, University of Georgia, kydahlin@gmail.com
## Started Oct 2020
##
## Title: Functions for getting model outputs (R0, pV, IV, and E0) #############
##
## Project: Global zoonoses - spillover of mosquito-borne pathogens
##
## Purpose: Provide functions for evaluating several model properties and equilibrium values
##
## Contents: 1) Set-up,load packages, get data, etc.
##           2) Build parameter sets
##           3) Create functions for getting parameters
##
##
## Inputs:
##          A dataframe containing all the variables needed to compute these quantities including:
##            All host parameters
##            Mosquito species considered
##            Temperature
##              Mosquito parameters associated with these
##            Model type ("RM" = Ross-MacDonald or "CCH" = Chitnis dynamic contact rate model)
##
## Outputs: A collection of functions with the above inputs. Functions labeled "compute" have continuous
##          outputs and those with "assign" have discrete outputs. Remaining functions are more general
##
##          compute.RH: Evaluate the host type reproduction number
##
##          compute.RV: Evaluate the mosquito type reproduction number
##
##          compute.R0: Evaluate the basic reproduction number (actually just sqrt(RH*RV))
##
##          compute.pH: Evaluate equilibrium prevalence in hosts
##
##          compute.pV: Evaluate equilibrium prevalence in mosquitoes
##
##          compute.V0: Evaluate equilibrium total abundance of mosquitoes

# 0) Helper functions ----
###* Reduce resolution of a vector by sub-sampling to length new_length
res_reduce <- function(df, new_length) {
  old_length <- length(unique(df))
  if (old_length < new_length) {
    warning("new_length is larger than old length. Vector will be unchanged")
    new_length <- old_length
  }
  ret <- unique(df)[seq(1, length(unique(df)), length.out = new_length)]
}

# 1) Functions for computing model outputs ----

###* Equilibrium total vector population ---------------------------------------
compute.V0 <- function(input) {
  with(input, {
    eps <- .Machine$double.eps
    temp <- sigmaV_f * deltaL
    # check if net reproduction rate exceeds mortality rate
    temp_bool <- temp > (1 / lf)
    ifelse(temp_bool,
           V0 <- KL * rhoL * lf * (1 - 1 / (lf * temp + eps)) ,
           V0 <- 0 # if mortality exceeds reproduction, set to zero
    )
  })
}

### Host-to-mosquito reproduction number----------------------------------------
compute.RH <- function(input) {
  with(input, {
    eps <- .Machine$double.eps
    # Derive intermediate quantities
    bV <- ifelse(is.infinite(sigmaH),
                 sigmaV, # Ross-Macdonald model
                 sigmaV * sigmaH * KH / (sigmaH * KH + sigmaV * V0 + eps)
    )
    RH <- ifelse(V0 == 0, 
                 0, 
                 bV * betaV * V0 * exp(-1 / (lf * etaV)) / (KH * (gammaH + muH) + eps))
  })
}

### Mosquito-to-host reproduction number----------------------------------------
compute.RV <- function(input) {
  with(input, {
    eps <- .Machine$double.eps
    RV <- ifelse(is.infinite(sigmaH),
                 sigmaV * betaH / (1 / (lf + eps)), # Ross-Macdonald
                 sigmaH * sigmaV * betaH * KH / ((1 / (lf + eps)) * (sigmaH * KH + sigmaV * V0))
    )
  })
}

### Basic reproduction number---------------------------------------------------
# Computed via the type reproduction numbers (see above)
compute.R0 <- function(input) {
  RH <- compute.RH(input)
  RV <- compute.RV(input)
  # Return R0
  return(sqrt(RH * RV))
}

###* LambdaR0 (derivative of R0 wrt pace) --------------------------------------
compute.lambdaR0 <- function(input) {
  with(as.list(input), {
    # Bite supply
    Bs <- sigmaH * KH
    # Bite demand
    Bd <- sigmaV * V0
    # Bite supply-demand ratio
    B_ratio <- Bd / Bs
    #
    beta.c <- get.Host.beta.c(input)
    #
    alpha.c <- get.Host.alpha.c(input)
    #
    r.p <- beta.c + alpha.c - mu.c
    #
    lambdaR0 <- ifelse(Model == "CCH",
      ((B_ratio - 1) / (B_ratio + 1)) * K.c - mu.c + alpha.c + beta.c,
      -K.c + r.p
    )
  })
}

###* Reactivity of the disease-free equilibrium---------------------------------
compute.E0 <- function(input) {
  RH <- compute.RH(input)
  RV <- compute.RV(input)
  with(input, {
    # Return E0
    return(0.5 * (muV * RV + (gammaH + muH) * RH) / sqrt((gammaH + muH) * muV))
  })
}

###* Equilibrium Host prevalence--------------------------------------------------
compute.pH <- function(input) {
  # V0 = compute.V0(input)
  R0 <- compute.R0(input)
  RV <- compute.RV(input)
  with(input, {
    # Derive intermediate quantities
    M <- muV * V0 / (muH * KH)
    # Return pV
    pH <- (muH / (gammaH + muH)) * (M * RV / (M * RV + 1)) * (1 - (1 / R0^2))
    return(pH)
  })
}

###* Equilibrium Vector prevalence------------------------------------------------
compute.pV <- function(input) {
  R0 <- compute.R0(input)
  # V0 <- compute.V0(input)
  RH <- compute.RH(input)
  with(input, {
    # Derive intermediate quantities
    M <- muV * V0 / (muH * KH)
    pV <- ((RH / M) / ((RH / M) + 1)) * (1 - (1 / R0^2))
  })
}


###* Equilibrium Infected vector abundance----------------------------------------
compute.IV <- function(input) {
  V0 <- compute.V0(input)
  pV <- compute.pV(input)
  # return IV
  pV * V0
}

###* Full Jacobian Matrix at the DFE----------------------------------------------
compute.J0 <- function(input) {
  with(input, {
    # Derive intermediate quantities
    V0 <- KL * rhoL * (1 - muV / (sigmaV * fecundity * deltaL)) / muV
    deltaV <- etaV / (etaV + muV)
    zeta <- sigmaV * sigmaH / (sigmaV * V0 + sigmaH * KH)
    CH <- zeta * betaV * V0
    CV <- zeta * betaH * KH

    matrix(c(
      -muH, 0, CV,
      CH, -(etaV + muV), 0,
      0, etaV, -muV
    ),
    3, 3,
    byrow = TRUE
    )
    # return(J0)
  })
}