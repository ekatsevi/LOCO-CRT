######### SIMULATION PARAMETERS ###################
# problem parameters
n = 1000                          # number of observations
p = 500                           # number of variables
sigma = 1                         # noise standard deviation
reps = 25                         # number of repetitions
k_vals = c(25,50,100)             # number of nonzero coefficients
rho_vals = c(0.2, 0.5, 0.8)       # AR(1) parameter for covariance of X
A_vals_1 = seq(2, 4.5, by = 0.5)  # specification of signal amplitudes
A_vals_2 = seq(2.5, 6.25, by = 0.75)
A_vals_3 = seq(4,9, by = 1)
A_vals = unique(c(A_vals_1, A_vals_2, A_vals_3))
A_vals_ranks = vector("list", 3)
names(A_vals_ranks) = rho_vals
indices = c(1,6,2,5,3,4)
A_vals_ranks[["0.2"]] = A_vals_1[indices]
A_vals_ranks[["0.5"]] = A_vals_2[indices]
A_vals_ranks[["0.8"]] = A_vals_3[indices]


# multiple testing parameters
global_error_rates = "FDR"        # global error rates to control
alpha_FDR    = 0.1                # FDR control level

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

loco_crt_methods = tibble()
marginal_crt_methods = tibble()
hrt_methods = tibble()
oracle_hrt_methods = tibble()
full_crt_methods = tibble()

knockoff_methods = cross_df(list(test_stat = "glmnet_coefdiff", knockoff_realization = 1:50)) %>%
  mutate(method_name = sprintf("knockoffs_%s_%d", test_stat, knockoff_realization))  

# list of all methods
methods = knockoff_methods$method_name