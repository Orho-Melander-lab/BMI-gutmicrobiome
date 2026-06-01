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
library(vegan)

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
Y <- as.matrix(tax)

# covariates
covariates <- data.frame(
  age = data$age,
  male = data$male
)

# phenotype
pheno <- data$BMI

X <- data.frame(
  BMI = pheno,
  age = covariates$age,
  male = covariates$male
)

rda_bmi <- rda(Y ~ BMI + Condition(age + male), data = X)

# significance test, BMI effect only
#anova_res <- anova(rda_bmi, by = "term", permutations = 999) #0.001 *** BMI was associated with gut microbial composition (permutation test, p = 0.001), explaining a small but statistically robust fraction of variance (adjusted R² = 0.46%)


##### bootstrapping
rda_boot_fun <- function(Y, X, idx) {
  
  Yb <- Y[idx, , drop = FALSE]
  Xb <- X[idx, , drop = FALSE]
  
  model <- rda(Yb ~ BMI + Condition(age + male), data = Xb)
  
  RsquareAdj(model)$adj.r.squared
}


n_boot <- 500
n <- nrow(Y)

boot_r2 <- numeric(n_boot)

# ensure directory exists
if (!dir.exists("results")) dir.create("results")

for (b in seq_len(n_boot)) {
  
  idx <- sample.int(n, replace = TRUE)
  
  boot_r2[b] <- rda_boot_fun(Y, X, idx)
  
  # progress print
  if (b %% 10 == 0) print(paste("Bootstrap", b))
  
  # save every 50 iterations
  if (b %% 50 == 0) {
    
    current <- data.frame(
      iteration = seq_len(b),
      adj_r2 = boot_r2[1:b]
    )
    
    file_name <- paste0("rda_bootstrap_up_to_", b, ".csv")
    
    write.csv(
      current,
      file.path("results","results_rda", file_name),
      row.names = FALSE
    )
    
    print(paste("Saved checkpoint at", b))
  }
}


boot_r2_clean <- boot_r2[!is.na(boot_r2)]

ci <- quantile(boot_r2_clean, probs = c(0.025, 0.975))
mean_r2 <- mean(boot_r2_clean)

rda_summary_df <- data.frame(
  term = "BMI",
  adj_r2 = RsquareAdj(rda_bmi)$adj.r.squared,
  mean_r2_boot = mean_r2,
  ci_low = ci[1],
  ci_high = ci[2]
)


export(rda_summary_df, file = here::here("results/results_rda/VE_by_BMI_in_all_species_RDA_with_CI.xlsx"))


cat(paste0('Sim finished on ', Sys.time(), '\n'))
cat('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
cat('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
