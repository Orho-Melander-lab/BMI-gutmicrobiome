#!/ludc/Tools/Software/R/4.2.2/bin/Rscript

cat('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
cat('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
cat(paste0('Sim started on ', Sys.time(), '\n'))


### nice reference for transforming compositional data before applying penalized regression: https://upcommons.upc.edu/bitstream/handle/2117/340287/lqaa029.pdf?sequence=1
set.seed(1)

library(rio)
library(BiocParallel) 
library(glmnet)
library(caret)
library(tidyverse)
library(ggplot2)
library(readxl)
library(here)


### nice reference for transforming compositional data before applying penalized regression: https://upcommons.upc.edu/bitstream/handle/2117/340287/lqaa029.pdf?sequence=1
data <- readxl::read_xlsx(here::here("data_nds_filt_CLRtrans_9993_updated250208.xlsx"))


data <- data %>% 
  mutate(male = sex == 1)


rownames(data) <- data$Data_sample_code
BMI <- data$BMI


## Get the MGSs
# identify species columns
start <- which(colnames(data) == 'shannon_nds') + 1
species_cols <- colnames(data)[start:(ncol(data)-1)]

# species matrix
tax <- data[, species_cols]

# covariates
covariates <- data.frame(
  age = data$age,
  male = data$male
)

# phenotype
pheno <- data$BMI


var.fun.species <- function(pheno, species_vec, covariates) {
  
  data <- data.frame(
    pheno = pheno,
    species = species_vec,
    covariates
  )
  
  cv_folds <- createFolds(seq_len(nrow(data)), k = 10)
  
  r2_full_vec <- numeric(length(cv_folds))
  r2_null_vec <- numeric(length(cv_folds))
  r2_diff_vec <- numeric(length(cv_folds))
  
  for (j in seq_along(cv_folds)) {
    
    idx <- cv_folds[[j]]
    
    train <- data[-idx, ]
    test  <- data[idx, ]
    
    fit_full <- glm(species ~ pheno + age + male, data = train)
    fit_null <- glm(species ~ age + male, data = train)
    
    pred_full <- predict(fit_full, test)
    pred_null <- predict(fit_null, test)
    
    r2_full <- cor(pred_full, test$species)^2
    r2_null <- cor(pred_null, test$species)^2
    
    r2_full_vec[j] <- r2_full
    r2_null_vec[j] <- r2_null
    r2_diff_vec[j] <- r2_full - r2_null
  }
  
  return(c(
    r2_full = mean(r2_full_vec),
    r2_null = mean(r2_null_vec),
    r2_diff = mean(r2_diff_vec)
  ))
}

var_exp_est <- data.frame(
  r2_full = numeric(length(species_cols)),
  r2_null = numeric(length(species_cols)),
  r2_diff = numeric(length(species_cols)),
  species = species_cols
)

for (i in seq_along(species_cols)) {
  
  sp <- species_cols[i]
  
  res <- var.fun.species(
    pheno,
    tax[[sp]],
    covariates
  )
  
  var_exp_est[i, "r2_full"] <- res["r2_full"]
  var_exp_est[i, "r2_null"] <- res["r2_null"]
  var_exp_est[i, "r2_diff"] <- res["r2_diff"]
}


bootstrap_var_exp <- function(pheno, tax, covariates, species_cols,
                              n_boot = 500,
                              save_every = 50,
                              file_prefix = "bootstrap_results") {
  
  n <- length(pheno)
  n_species <- length(species_cols)
  
  res_boot <- data.frame(
    boot = integer(n_boot * n_species),
    species = character(n_boot * n_species),
    r2_diff = numeric(n_boot * n_species),
    stringsAsFactors = FALSE
  )
  
  row_id <- 1
  
  for (b in seq_len(n_boot)) {
    
    boot_ix <- sample.int(n, replace = TRUE)
    
    pheno_b <- pheno[boot_ix]
    tax_b <- tax[boot_ix, , drop = FALSE]
    cov_b <- covariates[boot_ix, , drop = FALSE]
    
    for (i in seq_along(species_cols)) {
      
      sp <- species_cols[i]
      
      res <- var.fun.species(
        pheno_b,
        tax_b[[sp]],
        cov_b
      )
      
      res_boot[row_id, "boot"] <- b
      res_boot[row_id, "species"] <- sp
      res_boot[row_id, "r2_diff"] <- res["r2_diff"]
      
      row_id <- row_id + 1
    }
    
    print(paste("Finished bootstrap", b))
    
    # SAVE CHECKPOINT
    if (b %% save_every == 0) {
      
      # keep only filled rows
      current_res <- res_boot[1:(row_id - 1), ]
      
      file_name <- paste0(file_prefix, "_up_to_boot_", b, ".csv")
      
      write.csv(current_res, 
                file.path("results","results_raclr", file_name), 
                row.names = FALSE)
      
      print(paste("Saved checkpoint:", file_name))
    }
  }
  
  return(res_boot[1:(row_id - 1), ])
}


res_boot_df <- bootstrap_var_exp(
  pheno = pheno,
  tax = tax,
  covariates = covariates,
  species_cols = species_cols,
  n_boot = 500,          
  save_every = 50,      # save after every 3 bootstraps
  file_prefix = "bootstrap_results"
) 

boot_ci <- res_boot_df %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(
    lower_CI = quantile(r2_diff, 0.025), 
    upper_CI = quantile(r2_diff, 0.975),
    mean_boot = mean(r2_diff),
    .groups = "drop"
  )

var_exp_all <- merge(var_exp_est, boot_ci, by = "species")

export(var_exp_all, file = here::here("results/results_raclr/VE_by_BMI_in_each_speciesRACLR_with_CI.xlsx"))


cat(paste0('Sim finished on ', Sys.time(), '\n'))
cat('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
cat('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
