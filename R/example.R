# Generating simulated dataset
set.seed(12)
data.sim = simulate_group_inter(N = 50, n_mvnorm = 3, grouped = TRUE,
                                timepoints = 3:5, nonpara_inter = TRUE,
                                sample_from = seq(0,52,13), cst_ni = FALSE, 
                                cos = FALSE, A_vec = c(1, 1.5))

sim1 = data.sim$sim

# lambdas <- round(exp(seq(log(0.1), log(1 * 0.0001),
#                              length.out = 10
# )), digits = 4)
# 
# gammas <- c(
#   0.000001, 0.0000001, 0.00000001, 0.000000001
# )

lambdas = c(0.0046, 0.0001)
gammas <- 0.00000001

x = as.matrix(sim1[,-1:-3])
y = sim1$y
series = sim1$series
t = sim1$position

bases = create_bases(t, keep = NULL)


lasso_sim1 = plmm_lasso(x, y, series, t, name_group_var = "group", bases$bases,
                        gamma = gammas[1], lambda = lambdas[1], timexgroup = TRUE,
                        criterion = "BIC")


# Run model on a grid of hyperparameter and retrieve the model with the best BIC
tuned_plmm <- tune_plmm(x, y, series, t, name_group_var = "group", bases$bases,
  gamma_vec = gammas, lambda_vec = lambdas, timexgroup = TRUE,
  criterion = "BIC")

tuned_plmm$hyperparameters
tuned_plmm$lasso_output$theta

# Visualize overall fit and nonlinear functions 
plot_fit(x, y, series, t,  name_group_var = "group", 
  plmm_output, predicted = FALSE)

# Get debiased fixed-effects and pvalues
posi = debias_plmm(x, y, series, tuned_plmm, a = 1, Z = NULL)
posi


# Inference for the nonlinear functions
test_nonlinear_functions = test_f(x, y, series, t,  name_group_var = "group",
                                  tuned_plmm, n_boot = 10, predicted = F)

# Test of equality of the functions H0: f1 = f2
test_nonlinear_functions[[1]]

head(test_nonlinear_functions[[2]])

# Continuous joint confidence bands for the difference between f1 and f2
test_nonlinear_functions[[3]]

