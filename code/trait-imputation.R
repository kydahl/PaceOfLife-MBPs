# Contributors: Kyle Dahlin, JP Schmidt

library(tidyverse)
library(doParallel)
library(caret)
library(missForest)

cl <- makeCluster(parallel::detectCores() - 1)
registerDoParallel(cl)

# primate data from PanTheria
primate <- read_tsv("./data/raw/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  # Filter to primates
  filter(MSW05_Order == "Primates") %>% 
  # Remove human primates
  filter(MSW05_Binomial != "Homo sapiens") %>% 
  # Remove references column
  select(-c("References")) %>% 
  # Replace NA tokens (-999) with NAs
  mutate(across(where(is.numeric), ~ na_if(., -999)))

# primate = data.frame(apply(primate, 2, as.factor))

# rodent data from PanTheria
rodent <- read_tsv("./data/raw/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  # Filter to primates
  filter(MSW05_Order == "Rodentia") %>% 
  # Remove references column
  select(-c("References")) %>% 
  # Replace NA tokens (-999) with NAs
  mutate(across(where(is.numeric), ~ na_if(., -999)))

# get taxonomic columns
primate_taxa <- primate[, 3:4]
primate_test <- rep(1, nrow(primate_taxa))
primate_taxa <- cbind(primate_taxa, primate_test)

rodent_taxa <- rodent[, 3:4]
rodent_test <- rep(1, nrow(rodent_taxa))
rodent_taxa <- cbind(rodent_taxa, rodent_test)

# turn taxonomic groupings into dummy binary variables
primate_dummies <- dummyVars(primate_test ~ ., data = primate_taxa)
primate_taxa <- predict(primate_dummies, newdata = primate_taxa)

rodent_dummies <- dummyVars(rodent_test ~ ., data = rodent_taxa)
rodent_taxa <- predict(rodent_dummies, newdata = rodent_taxa)

# impute using missForest
set.seed(82)
primate_imp <- missForest(data.frame(select(primate, c(6:35))), 
                          maxiter = 100, ntree = 1000, 
                          verbose = TRUE, variablewise = TRUE, 
                          parallelize = "variables") 

primate_imp_w_taxa <- missForest(data.frame(cbind(primate[, 6:35], primate_taxa)), 
                                 maxiter = 100, ntree = 1000, 
                                 verbose = TRUE, variablewise = TRUE,
                                 parallelize = "variables")

set.seed(82)
rodent_imp <- missForest(data.frame(select(rodent, c(6:35))), 
                         maxiter = 100, ntree = 1000, 
                         verbose = TRUE, variablewise = TRUE,
                         parallelize = "variables") 

rodent_imp_w_taxa <- missForest(data.frame(cbind(rodent[, 6:35], rodent_taxa)), 
                                maxiter = 100, ntree = 1000, 
                                verbose = TRUE, variablewise = TRUE,
                                parallelize = "variables")

# bring back taxonomic columns
primate2 <- cbind(primate[, 1:5], primate_imp$ximp[, 1:29])
primate2_taxa <- cbind(primate[, 1:5], primate_imp_w_taxa$ximp[, 1:29])


rodent2 <- cbind(rodent[, 1:5], rodent_imp$ximp[, 1:29])
rodent2_taxa <- cbind(rodent[, 1:5], rodent_imp_w_taxa$ximp[, 1:29])

# write out imputed data
write.csv(primate2, "data/clean/primate_imputed.csv")
write.csv(primate2_taxa, "data/clean/primate_taxa_imputed.csv")


write.csv(rodent2, "data/clean/rodent_imputed.csv")
write.csv(rodent2_taxa, "data/clean/rodent_taxa_imputed.csv")