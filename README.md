# Fast-lived Vertebrate Hosts Exhibit Higher Potential for Mosquito-borne Parasite Transmission

## Background:
### Purpose: investigate the pace of life hypothesis in the context of mosquito-borne pathogens
The resurgence and emergence of vector-borne zoonoses is a pressing concern across several regions. Predicting hotspots (areas with an increased risk of the emergence or resurgence of such zoonoses) will be an essential preparation and prevention tool for countries and communities likely to be affected by these diseases. A spillover event is one in which a pathogen which is usually transmitted within one population enters a novel host population. Many studies have sought to estimate the risk of spillover of pathogens from various wildlife populations into human populations in order to identify areas to focus surveillance or to implement pre-emptive control policies.

Mosquito-borne pathogens of wildlife (in particular, viruses) have been a recent focus because they have led to several outbreaks and epidemics (e.g. Dengue, Zika, Chikungunya). Estimating spillover risk is a challenge in vector-borne disease systems because the transmission cycle involves the ecology of hosts, vectors, and pathogens. Many studies have estimated risk entirely based on properties of the mosquito vector. While mosquitoes would ultimately the be the culprit of any spillover event to humans, the ability of the wildlife host population to effectively maintain the transmission of the pathogen (a reservoir or maintenance host) or to amplify transmission risk (an amplification host) is also vital. Integrating information on wildlife hosts has mostly been accomplished through incorporating data on presence and susceptibility. However, collecting these data is costly and time-intensive and thus data are not available for most potential host species.

We investigate the generality of the pace of life hypothesis to mosquito-borne pathogen transmission systems. This hypothesis suggests "that species that live fast have the highest equilibrium prevalence" (Han et al., Ecology Letters 2020). We explore this hypothesis by representing host species according to a pace of life parameter (which determines the species birth rate, longevity, and carrying capacity density) and evaluating invariants of an ODE model of disease transmission relevant to spillover risk or reservoir status (e.g. equilibrium vector prevalence or basic reproduction number). We investigate these values across four medically-important mosquito species and along a broad parameter space of pathogen traits.

Thus we obtain a relationship between pace of life and reservoir status (R_0 > 1) and determine when it is that R_0 increases with pace of life. The extent to which host pace of life and these epidemiologically important predictions are correlated can inform managers of target species for surveillance of emerging vector-borne diseases in wildlife populations. We also further the investigation into paradigm hypotheses in disease ecology and compare predictions in mosquito-borne transmission systems to those of directly- and environmentally-transmitted pathogens.

## Getting started:
* Run `manuscript-figures.Rmd` to generate the figures and supplementary figures from the manuscript
* If you wish to replicate the full analysis, run the following scripts in the `code` folder
  * `trait-imputation.R` to create the imputed data sets for primate and rodent traits from PANTHERIA
  * `trait-hyperparameters.R` to obtain the model parameters as functions of pace of life
  * `get-analysis.R` for obtain all model outputs (R_0, lambda_R0, etc.) as functions of pace of life and the other variables

## Contents:
* `code`: contains all code necessary to obtain all study results
* `data`: includes original (`data/raw`) and modified (`data/clean`) data used in the study
* `results`: contains final results of the study, to be used to generate figures with `manuscript-figures.Rmd`
* `doc`: contains outputs from  `manuscript-figures.Rmd`
## Authors:
 Kyle Dahlin ([kydahlin@gmail.com](kydahlin@gmail.com), Virginia Tech), Suzanne O'Regan (Paperpile), JP Schmidt (University of Georgia), Barbara Han (Cary Institute of Technology), John Drake (University of Georgia)

## License:
This project is licensed under the MIT License - see the LICENSE.md file for details.