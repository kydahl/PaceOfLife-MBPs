# Code for generating results and figures
To reproduce the study results, run the following code in this order:
  * `trait-imputation.R` creates the imputed data sets for primate and rodent traits from PANTHERIA
  * `trait-hyperparameters.R` provides the model parameters as functions of pace of life
  * `get-analysis.R` for obtain all model outputs (R_0, lambda_R0, etc.) as functions of pace of life and the other variables

These scripts contain functions used for analysis:
  * `output-functions.R`: functions for calculating model outputs such as the basic reproduction number R_0
  * `trait-functions.R`: functions describing the relationships between pace of life and immunological parameters