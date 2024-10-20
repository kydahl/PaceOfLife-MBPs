# Results of study analyses
* `AllOutputs_high.rds`: all model outputs when immune trait variability is wide (`..._low.Rds` corresponds to narrow variability and `..._none.Rds` to no variability)
  * NB: These files are kept separate to ensure the file size is small enough to post to GitHub
* `criticalB_table.csv`: values of the vector-host ratio which must be exceed in order for R_0 to increase with pace of life
* `flow_diagram.png`: a figure illustrating the ODE model as a compartmental diagram
* `VectorTraits.rds`: data frame containing just the mosquito thermal performance curves used in this analysis (reduced from `data/clean/parameter_TPCs.rds`)