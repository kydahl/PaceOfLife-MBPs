# Raw and formatted data files

The `raw` folder contains data directly obtained from the PANTHERIA database. 

The `clean` folder contains modified data:
* `primate_imputed.csv` and `primate_taxa_imputed.csv` are primate trait data sets obtained by imputing missing data from the original PANTHERIA dataset (for ropdents, these data are in `rodent_imputed.csv` and `rodent_taxa_imputed.csv`)
* `trait-hyperparameters.csv` contains the values for hyperparameters used to determine the functional relationship between pace of life and the life history parameters of the model
* `parameter_TPCs.rds` is a dataset of fitted thermal performance curves for the mosquito traits used in the model