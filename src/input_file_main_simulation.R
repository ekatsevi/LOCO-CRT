######### SIMULATION PARAMETERS ###################
# problem parameters
n = 1000                          # number of observations
p = 500                           # number of variables
sigma = 1                         # noise standard deviation
reps = 25                         # number of repetitions
rho_vals = c(0.2, 0.5, 0.8)       # AR(1) parameter for covariance of X
k_vals = c(25,50,100)             # number of nonzero coefficients
A_vals_1 = seq(0.5, 5, length.out = 10) # specification of signal amplitudes
A_vals_2 = seq(0.7, 7, length.out = 10)
A_vals_3 = seq(1,10, length.out = 10)
A_vals = unique(c(A_vals_1, A_vals_2, A_vals_3))
A_vals_ranks = vector("list", 3)
names(A_vals_ranks) = rho_vals
indices = c(1,10,2,9,3,8,4,7,5,6)
A_vals_ranks[["0.2"]] = A_vals_1[indices]
A_vals_ranks[["0.5"]] = A_vals_2[indices]
A_vals_ranks[["0.8"]] = A_vals_3[indices]

# multiple testing parameters
global_error_rates = c("FWER", "FDR") # global error rates to control
alpha_FDR    = 0.1                    # FDR control level
alpha_FWER   = 0.05                   # FWER control level

# simulation parameters
parameters = cross_df(list(k = k_vals, A = A_vals, rep = 1:reps, rho = rho_vals)) %>%
  filter((rho == rho_vals[1]) & (A %in% A_vals_1) | 
         (rho == rho_vals[2]) & (A %in% A_vals_2) |
         (rho == rho_vals[3]) & (A %in% A_vals_3)) %>%
  group_by(rep, rho) %>%
  mutate(A_rank = match(A, A_vals_ranks[[as.character(unique(rho))]])) %>%
  arrange(A_rank, .by_group = TRUE) %>%
  ungroup() %>%
  select(-A_rank) 

######### METHODS TO COMPARE ###################

loco_crt_methods = tibble(inner_method = "quadratic", stopping_rule = "1se", 
                          fit_type = "in_sample", method_name = "LOCO CRT")
marginal_crt_methods = tibble(inner_method = "quadratic", method_name = "Marginal CRT") 

hrt_methods = tibble(learning_method = "relaxed_lasso", inner_method = "quadratic_fixed", 
                     stopping_rule = "1se", holdout_prop = 0.2, method_name = "HRT") 

oracle_hrt_methods = tibble(inner_method = "quadratic_fixed", holdout_prop = 0.2, 
                            method_name = "Oracle HRT")

knockoff_methods = tibble(test_stat = "glmnet_coefdiff", knockoff_realization = 1, 
                          method_name = "Knockoffs") 

full_crt_methods = tibble()

# list of all methods
methods = c(loco_crt_methods$method_name,  
            marginal_crt_methods$method_name,
            hrt_methods$method_name, 
            oracle_hrt_methods$method_name, 
            knockoff_methods$method_name)