# Marlena Maziarz and Mariam Miari
# #MM = edited by Mariam Miari
# 250205
# This script takes a data frame of est and p-values from linear, logistic and linear present only analyses
# and finds families that are over-represented in the stat sig results
# exports results for each, for all combined and for the most significant only

# the metric to look at is most promising (TRUE if stat sig based on Fisher or chisquared and overrepresented among the stat sig results for individual mgs's)
# as well as teh ratio_sig_vs_all (proportion in the stat sig vs proportion in all)



rm(list = ls())

library(rio)
library(tidyverse)
library(here)
library(epitools)

d_org <- import(here::here('data.xlsx'))

str(d_org)
summary(d_org)

d <- d_org %>%
    transmute(mgs = maintaxa,
              family = family,
              linear_all_est = est_linear_all,
              linear_all_p_adj = p_value_adj_linear_all,
              logistic_all_est = est_logistic_all,
              logistic_all_p_adj = p_value_adj_logistic_all,
              linear_present_only_est = est_presence,
              linear_present_only_p_adj = p_value_adj_presence)



# input
# x_table = output of table()
get_or <- function(x){
    or <- (x[1, 1]*x[2, 2])/(x[1, 2]*x[2, 1])
    return(or)
}



# model = 'linear_all', 'logistic_all', 'linear_present_only'
# subset = 'all', 'positive', 'negative'

get_stat_sig_families <- function(d, 
                                  model = 'linear_all', 
                                  est_subset = 'all_est', # subset of the estimates
                                  p_val_thr = 0.05,
                                  perc_digits = 2,
                                  or_digits = 2,
                                  p_val_digits = 4){
  
  if(model == 'linear_all'){
    d <- d %>%
      mutate(est = linear_all_est,
             p_adj = linear_all_p_adj)
    est_thr <- 0
    
  }else if(model == 'logistic_all'){
    d <- d %>%
      mutate(est = logistic_all_est,
             p_adj = logistic_all_p_adj)
    est_thr <- 0 # this should be changed to 1 if est is an OR
    
  }else if(model == 'linear_present_only'){
    d <- d %>%
      mutate(est = linear_present_only_est,
             p_adj = linear_present_only_p_adj)
    est_thr <- 0
    
  }else{
    print('No such option for model, exiting...')
    return(NULL)
  }
  
  # store n_family from all_est
  d_all_est <- d  #MM
  n_family_all <- d_all_est %>%
    group_by(family) %>%
    summarise(n_family = n(), .groups = 'drop') %>%
    rename(n_family_all_est = n_family)
  
  # Filter based on est_subset
#  if(est_subset == 'positive_est'){
 #   d <- d %>% filter(est >= est_thr)
 # }else if(est_subset == 'negative_est'){
  #  d <- d %>% filter(est < est_thr)
  #}else if(est_subset != 'all_est'){
   # print('No such option for subset, exiting')
    #return(NULL)
 # }
  
  families <- sort(unique(d$family))
  n_families <- length(families) 
  
  results <- data.frame(model = rep(model, n_families),
                        est_subset = rep(est_subset, n_families),
                        family = families,
                        most_promising = rep(NA, n_families), # 1 if ratio_sig_vs_all > 1 and p_val < 0.05
                        ratio_sig_vs_all = rep(NA, n_families), # perc_sig_in_family/perc_in_family (calc before rounding)
                        p_val = rep(NA, n_families), 
                        n_sig_family = rep(NA, n_families), 
                        n_sig_total = rep(sum(d$p_adj < p_val_thr), n_families), 
                        n_family = rep(NA, n_families),
                        n_total = rep(nrow(d), n_families),
                        perc_sig_in_family = rep(NA, n_families), # n_sig_family/n_sig_total*100
                        perc_in_family = rep(NA, n_families),  # n_family/n_total*100
                        or = rep(NA, n_families), 
                        test_used = rep(NA, n_families))
  
  if(n_families == 0){
    print('THIS SHOULD NEVER HAPPEN')
    return(results)
  }
  
  for(i in 1:n_families){
    if(est_subset == 'all_est'){
    	x <- table(d$family == families[i], d$p_adj < p_val_thr)
    }else if(est_subset == 'positive_est'){
	x <- table(d$family == families[i], d$p_adj < p_val_thr & d$est >=est_thr)
    }else if (est_subset == 'negative_est'){
	x <- table(d$family == families[i], d$p_adj < p_val_thr & d$est < est_thr)
    }  
    results[i, 'n_sig_family'] <- x[2,2]
    
    # n_family used from all_est
    results[i, 'n_family'] <- n_family_all$n_family_all_est[n_family_all$family == families[i]] #MM
    
    results[i, 'or'] <- get_or(x)
    
    # are any of the cells < 5? 
    ## run Fisher's exact test
    if(any(x) <= 5){
      p_val <- fisher.test(x)$p.value
      test <- 'fisher'
      ## run Chi-squared
    } else {
      p_val <- chisq.test(x)$p.value
      test <- 'chi-squared'
    }
    
    results[i, 'p_val'] <- p_val
    results[i, 'test_used'] <- test
  }
  
  # recalculate percentages using fixed n_family
  results$perc_sig_in_family <- round(results$n_sig_family / results$n_sig_total * 100, perc_digits)
  results$perc_in_family <- round(results$n_family / results$n_total * 100, perc_digits)
  results$ratio_sig_vs_all <- round((results$n_sig_family / results$n_sig_total) / (results$n_family / results$n_total), or_digits)
  results$or <- round(results$or, or_digits)
  results$p_val <- signif(results$p_val, p_val_digits)
  results$most_promising <- results$ratio_sig_vs_all > 1 & results$p_val < 0.05
  
  return(results)
}

results_list <- NULL

for(curr_model in c('linear_all', 'logistic_all', 'linear_present_only')){
  for(curr_est_subset in c('all_est', 'positive_est', 'negative_est')){
    results_list <- c(results_list, list(get_stat_sig_families(d, model = curr_model, est_subset = curr_est_subset)))
  }
}




results_list <- NULL

for(curr_model in c('linear_all', 'logistic_all', 'linear_present_only')){
    for(curr_est_subset in c('all_est', 'positive_est', 'negative_est')){
        results_list <- c(results_list, list(get_stat_sig_families(d, model = curr_model, est_subset = curr_est_subset)))
    }
}


results_all <- NULL

counter <- 1
for(curr_model in c('linear_all', 'logistic_all', 'linear_present_only')){
    for(curr_est_subset in c('all_est', 'positive_est', 'negative_est')){
        results_all <- rbind(results_all, results_list[[counter]])
        export(results_list[[counter]], file = paste0('results/overrepresented_families_in_sig_results_', curr_model, '_', curr_est_subset, '_250205.xlsx'))
        counter <- counter + 1
        
    }
}


export(results_all, file = 'results/overrepresented_families_in_sig_results_all_models_and_subsets_250205.xlsx')

export(results_all %>% filter(most_promising), file = 'results/_overrepresented_families_in_sig_results_MOST_PROMISING_FROM_ALL_RESULTS_250205.xlsx')


