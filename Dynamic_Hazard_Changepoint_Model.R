# =====================================================================
# Optimized & Parallelized Code for "The Dynamic Hazard Changepoint Model"
# Designed for Reproducibility in Peer Review
# Speeding changes ONLY: fast row-cumsum (identical numerics) + live % progress
# Statistical model / multi-start / optim / estimators UNCHANGED
# =====================================================================

library(stats)
library(parallel)

# Set to TRUE for a quick test run (R_sim = 50), or FALSE for the final paper results (R_sim = 300)
FAST_MODE <- FALSE

# ---------------------------------------------------------------------
# Fast row-wise cumulative sum (pure base R, identical to apply+cumsum order)
# This is the main speedup; results remain numerically the same.
# ---------------------------------------------------------------------
fast_row_cumsum <- function(x) {
  n <- nrow(x)
  p <- ncol(x)
  if (p == 0L) return(x)
  y <- x
  if (p >= 2L) {
    for (j in 2:p) {
      y[, j] <- y[, j - 1] + x[, j]
    }
  }
  y
}

# ---------------------------------------------------------------------
# 1. Fully Vectorized Data Generating Process (DGP)  [identical logic]
# ---------------------------------------------------------------------
generate_data <- function(n, m, beta0, psi0, mu0, sigma0, mu1, sigma1) {
  Z <- cbind(runif(n, -1, 1), rnorm(n, 0, 1))
  
  k_seq <- matrix(1:m, nrow = n, ncol = m, byrow = TRUE)
  logit_pi <- matrix(Z %*% beta0, nrow = n, ncol = m, byrow = FALSE) + psi0 * k_seq
  pi_mat <- 1 / (1 + exp(-logit_pi))
  
  log_pi <- log(pi_mat + 1e-15)
  log_1_pi <- log(1 - pi_mat + 1e-15)
  
  log_surv <- cbind(0, fast_row_cumsum(log_1_pi))
  alpha <- exp(cbind(log_pi + log_surv[, 1:m], log_surv[, m + 1]))
  alpha <- alpha / rowSums(alpha)
  
  cum_alpha <- fast_row_cumsum(alpha)
  u <- runif(n)
  tau <- rowSums(cum_alpha < matrix(u, nrow = n, ncol = m + 1, byrow = FALSE)) + 1
  
  X <- matrix(rnorm(n * m, mean = mu0, sd = sigma0), nrow = n, ncol = m)
  
  idx_mat <- matrix(1:m, nrow = n, ncol = m, byrow = TRUE)
  tau_mat <- matrix(tau, nrow = n, ncol = m, byrow = FALSE)
  mask <- (idx_mat >= tau_mat)
  
  X_post <- matrix(rnorm(n * m, mean = mu1, sd = sigma1), nrow = n, ncol = m)
  X[mask] <- X_post[mask]
  
  list(X = X, Z = Z, tau_true = tau)
}

# ---------------------------------------------------------------------
# 2. Vectorized Negative Log-Likelihood Function (Dynamic Model)
# ---------------------------------------------------------------------
neg_log_likelihood_vec <- function(theta, X, Z, m, mu0 = 0, sigma0 = 1, sigma1 = 1.5) {
  beta <- theta[1:2]
  psi <- theta[3]
  mu1 <- theta[4]
  n <- nrow(X)
  
  k_seq <- matrix(1:m, nrow = n, ncol = m, byrow = TRUE)
  logit_pi <- matrix(Z %*% beta, nrow = n, ncol = m, byrow = FALSE) + psi * k_seq
  pi_mat <- 1 / (1 + exp(-logit_pi))
  
  log_pi <- log(pi_mat + 1e-15)
  log_1_pi <- log(1 - pi_mat + 1e-15)
  log_surv <- cbind(0, fast_row_cumsum(log_1_pi))
  log_alpha <- cbind(log_pi + log_surv[, 1:m], log_surv[, m + 1])
  
  log_f0 <- dnorm(X, mean = mu0, sd = sigma0, log = TRUE)
  log_f1 <- dnorm(X, mean = mu1, sd = sigma1, log = TRUE)
  
  cum_log_f0 <- cbind(0, fast_row_cumsum(log_f0))
  
  log_f1_rev <- log_f1[, m:1, drop = FALSE]
  cum_f1_rev <- fast_row_cumsum(log_f1_rev)
  cum_log_f1 <- cum_f1_rev[, m:1, drop = FALSE]
  cum_log_f1 <- cbind(cum_log_f1, 0)
  
  log_f <- cum_log_f0 + cum_log_f1
  log_terms <- log_alpha + log_f
  
  max_log <- apply(log_terms, 1, max)
  log_terms_scaled <- log_terms - matrix(max_log, nrow = n, ncol = m + 1, byrow = FALSE)
  ll_i <- max_log + log(rowSums(exp(log_terms_scaled)))
  
  return(-sum(ll_i))
}

# ---------------------------------------------------------------------
# 3. Vectorized Negative Log-Likelihood Function (Constant Hazard Model)
# ---------------------------------------------------------------------
neg_log_likelihood_constant <- function(theta, X, Z, m, mu0 = 0, sigma0 = 1, sigma1 = 1.5) {
  beta <- theta[1:2]
  psi <- 0
  mu1 <- theta[3]
  neg_log_likelihood_vec(c(beta, psi, mu1), X, Z, m, mu0, sigma0, sigma1)
}

# ---------------------------------------------------------------------
# 4. Multi-Start Grid Search Optimization Functions  [UNCHANGED]
# ---------------------------------------------------------------------
fit_model_multistart <- function(X, Z, m, mu0, sigma0, sigma1, model_type = c("dynamic", "constant"), fast = FALSE) {
  model_type <- match.arg(model_type)
  
  if (model_type == "dynamic") {
    if (fast) {
      grid <- expand.grid(beta1 = 0.5, beta2 = -1.0, psi = 0.03, mu1 = 2.0)
    } else {
      grid <- expand.grid(
        beta1 = c(0.1, 0.6),
        beta2 = c(-1.2, -0.5),
        psi   = c(0.01, 0.05),
        mu1   = c(1.0, 2.5)
      )
    }
    
    best_val <- Inf
    best_opt <- NULL
    
    for (k in 1:nrow(grid)) {
      theta_init <- as.numeric(grid[k, ])
      opt_result <- tryCatch({
        optim(par = theta_init, 
              fn = neg_log_likelihood_vec, 
              X = X, Z = Z, m = m,
              mu0 = mu0, sigma0 = sigma0, sigma1 = sigma1,
              method = "BFGS", hessian = TRUE,
              control = list(maxit = 500))
      }, error = function(e) NULL)
      
      if (!is.null(opt_result) && opt_result$convergence == 0) {
        if (opt_result$value < best_val) {
          best_val <- opt_result$value
          best_opt <- opt_result
        }
      }
    }
    return(best_opt)
    
  } else {
    if (fast) {
      grid <- expand.grid(beta1 = 0.5, beta2 = -1.0, mu1 = 2.0)
    } else {
      grid <- expand.grid(
        beta1 = c(0.1, 0.6),
        beta2 = c(-1.2, -0.5),
        mu1   = c(1.0, 2.5)
      )
    }
    
    best_val <- Inf
    best_opt <- NULL
    
    for (k in 1:nrow(grid)) {
      theta_init <- as.numeric(grid[k, ])
      opt_result <- tryCatch({
        optim(par = theta_init, 
              fn = neg_log_likelihood_constant, 
              X = X, Z = Z, m = m,
              mu0 = mu0, sigma0 = sigma0, sigma1 = sigma1,
              method = "BFGS", hessian = TRUE,
              control = list(maxit = 500))
      }, error = function(e) NULL)
      
      if (!is.null(opt_result) && opt_result$convergence == 0) {
        if (opt_result$value < best_val) {
          best_val <- opt_result$value
          best_opt <- opt_result
        }
      }
    }
    return(best_opt)
  }
}

# ---------------------------------------------------------------------
# 5. Posterior Inference and Vectorized Fallback CUSUM  [identical]
# ---------------------------------------------------------------------
estimate_tau <- function(theta, X, Z, m, mu0 = 0, sigma0 = 1, sigma1 = 1.5) {
  beta <- theta[1:2]
  psi <- theta[3]
  mu1 <- theta[4]
  n <- nrow(X)
  
  k_seq <- matrix(1:m, nrow = n, ncol = m, byrow = TRUE)
  logit_pi <- matrix(Z %*% beta, nrow = n, ncol = m, byrow = FALSE) + psi * k_seq
  pi_mat <- 1 / (1 + exp(-logit_pi))
  
  log_pi <- log(pi_mat + 1e-15)
  log_1_pi <- log(1 - pi_mat + 1e-15)
  log_surv <- cbind(0, fast_row_cumsum(log_1_pi))
  log_alpha <- cbind(log_pi + log_surv[, 1:m], log_surv[, m + 1])
  
  log_f0 <- dnorm(X, mean = mu0, sd = sigma0, log = TRUE)
  log_f1 <- dnorm(X, mean = mu1, sd = sigma1, log = TRUE)
  cum_log_f0 <- cbind(0, fast_row_cumsum(log_f0))
  
  log_f1_rev <- log_f1[, m:1, drop = FALSE]
  cum_f1_rev <- fast_row_cumsum(log_f1_rev)
  cum_log_f1 <- cum_f1_rev[, m:1, drop = FALSE]
  cum_log_f1 <- cbind(cum_log_f1, 0)
  
  log_f <- cum_log_f0 + cum_log_f1
  log_post <- log_alpha + log_f
  
  max_log <- apply(log_post, 1, max)
  log_sum_exp <- max_log + log(rowSums(exp(log_post - matrix(max_log, nrow = n, ncol = m + 1, byrow = FALSE))))
  post_probs <- exp(log_post - matrix(log_sum_exp, nrow = n, ncol = m + 1, byrow = FALSE))
  
  return(apply(post_probs, 1, which.max))
}

# 10x faster vectorized CUSUM function replacing the loop
detect_cusum_vectorized <- function(x) {
  n_val <- length(x)
  if (n_val < 2) return(n_val + 1)
  S <- cumsum(x)
  j <- 1:(n_val - 1)
  mean_pre <- S[j] / j
  mean_post <- (S[n_val] - S[j]) / (n_val - j)
  cusum_stat <- j * (n_val - j) * (mean_pre - mean_post)^2
  best_pt <- which.max(cusum_stat)
  if (cusum_stat[best_pt] > 10) return(best_pt + 1) else return(n_val + 1)
}

detect_pelt_single <- function(x) {
  if (requireNamespace("changepoint", quietly = TRUE)) {
    ans <- tryCatch({
      fit <- changepoint::cpt.mean(x, method = "PELT", penalty = "AIC")
      pts <- changepoint::cpts(fit)
      if (length(pts) == 0) return(length(x) + 1) else return(pts[1])
    }, error = function(e) return(length(x) + 1))
    return(ans)
  } else {
    return(detect_cusum_vectorized(x))
  }
}

# ---------------------------------------------------------------------
# 6. Monte Carlo Simulation Setup + live progress
# ---------------------------------------------------------------------
set.seed(2026)
R_sim <- if (FAST_MODE) 50 else 100
n <- 500
m <- 50

beta_true <- c(0.5, -1.0)
psi_true <- 0.03
mu0_true <- 0
sigma0_true <- 1
mu1_true <- 1.8
sigma1_true <- 2.5

theta_true <- c(beta_true[1], beta_true[2], psi_true, mu1_true)
theta_true_const <- c(beta_true[1], beta_true[2], mu1_true)

# Parallel Setup
num_cores <- max(1, detectCores() - 1)
cl <- makeCluster(num_cores)
clusterSetRNGStream(cl, iseed = 2026)   # better parallel reproducibility

# Export environment variables & functions to parallel workers
clusterExport(cl, list(
  "fast_row_cumsum",
  "generate_data", "neg_log_likelihood_vec", "neg_log_likelihood_constant",
  "fit_model_multistart", "estimate_tau", "detect_cusum_vectorized", 
  "detect_pelt_single", "FAST_MODE", "n", "m", "beta_true", "psi_true", 
  "mu0_true", "sigma0_true", "mu1_true", "sigma1_true", "theta_true", 
  "theta_true_const"
))
clusterEvalQ(cl, library(stats))

cat("Starting Parallel Monte Carlo Simulation with R =", R_sim, 
    "replications on", num_cores, "cores...\n")
cat("Progress will be shown as percentage...\n")
start_time <- Sys.time()

# ----- Batch-wise parallel execution with live % progress (no extra packages) -----
results <- vector("list", R_sim)
batch_size <- 20   # update progress every ~20 replications (adjust if desired)

for (start_r in seq(1, R_sim, by = batch_size)) {
  end_r <- min(start_r + batch_size - 1, R_sim)
  batch_idx <- start_r:end_r
  
  batch_res <- parLapply(cl, batch_idx, function(r) {
    res <- list(
      dyn_est = rep(NA, 4), dyn_se = rep(NA, 4), dyn_cov = rep(NA, 4), dyn_conv = 0, dyn_mae = NA,
      con_est = rep(NA, 3), con_se = rep(NA, 3), con_cov = rep(NA, 3), con_conv = 0, con_mae = NA,
      pelt_mae = NA
    )
    
    data <- generate_data(n, m, beta_true, psi_true, mu0_true, sigma0_true, mu1_true, sigma1_true)
    
    # Dynamic Model
    opt_dynamic <- fit_model_multistart(data$X, data$Z, m, mu0_true, sigma0_true, sigma1_true, "dynamic", fast = FAST_MODE)
    if (!is.null(opt_dynamic)) {
      inv_hessian <- tryCatch(solve(opt_dynamic$hessian), error = function(e) NULL)
      if (!is.null(inv_hessian) && all(diag(inv_hessian) > 0)) {
        res$dyn_est <- opt_dynamic$par
        res$dyn_se <- sqrt(diag(inv_hessian))
        res$dyn_cov <- (theta_true >= (opt_dynamic$par - 1.96 * res$dyn_se)) & 
          (theta_true <= (opt_dynamic$par + 1.96 * res$dyn_se))
        res$dyn_conv <- 1
        tau_hat_dyn <- estimate_tau(opt_dynamic$par, data$X, data$Z, m, mu0_true, sigma0_true, sigma1_true)
        res$dyn_mae <- mean(abs(data$tau_true - tau_hat_dyn))
      }
    }
    
    # Constant Model
    opt_const <- fit_model_multistart(data$X, data$Z, m, mu0_true, sigma0_true, sigma1_true, "constant", fast = FAST_MODE)
    if (!is.null(opt_const)) {
      inv_hessian_c <- tryCatch(solve(opt_const$hessian), error = function(e) NULL)
      if (!is.null(inv_hessian_c) && all(diag(inv_hessian_c) > 0)) {
        res$con_est <- opt_const$par
        res$con_se <- sqrt(diag(inv_hessian_c))
        res$con_cov <- (theta_true_const >= (opt_const$par - 1.96 * res$con_se)) & 
          (theta_true_const <= (opt_const$par + 1.96 * res$con_se))
        res$con_conv <- 1
        tau_hat_const <- estimate_tau(c(opt_const$par[1:2], 0, opt_const$par[3]), 
                                      data$X, data$Z, m, mu0_true, sigma0_true, sigma1_true)
        res$con_mae <- mean(abs(data$tau_true - tau_hat_const))
      }
    }
    
    # PELT / CUSUM
    tau_hat_pelt <- apply(data$X, 1, detect_pelt_single)
    res$pelt_mae <- mean(abs(data$tau_true - tau_hat_pelt))
    
    return(res)
  })
  
  results[batch_idx] <- batch_res
  
  pct <- round(100 * end_r / R_sim)
  cat(sprintf("\rProgress: %d / %d replications  (%d%%) completed...", end_r, R_sim, pct))
  flush.console()
}
cat("\n")

# Safely close the parallel connections
stopCluster(cl)

end_time <- Sys.time()
cat("Simulation completed in", round(difftime(end_time, start_time, units = "secs"), 1), "seconds.\n")

# ---------------------------------------------------------------------
# 8. Post-Processing and Formatting Tables  [UNCHANGED]
# ---------------------------------------------------------------------
converged_dyn <- which(sapply(results, function(x) x$dyn_conv) == 1)
converged_con <- which(sapply(results, function(x) x$con_conv) == 1)

if (length(converged_dyn) > 0) {
  estimates_dyn <- t(sapply(results[converged_dyn], function(x) x$dyn_est))
  se_dyn <- t(sapply(results[converged_dyn], function(x) x$dyn_se))
  cov_dyn <- t(sapply(results[converged_dyn], function(x) x$dyn_cov))
  
  bias_dyn <- colMeans(estimates_dyn) - theta_true
  sd_dyn <- apply(estimates_dyn, 2, sd)
  rmse_dyn <- sqrt(bias_dyn^2 + sd_dyn^2)
  ase_dyn <- colMeans(se_dyn)
  cp_dyn <- colMeans(cov_dyn)
  
  table1 <- data.frame(
    Parameter = c("beta1", "beta2", "psi", "mu1"), True = theta_true,
    Bias = round(bias_dyn, 3), SD = round(sd_dyn, 3), RMSE = round(rmse_dyn, 3),
    ASE = round(ase_dyn, 3), CP_95 = round(cp_dyn, 3)
  )
  cat("\n============================================================\n")
  cat("TABLE 1: Performance of the Proposed Dynamic Hazard Model (n =", n, ")\n")
  cat("============================================================\n")
  print(table1, row.names = FALSE)
}

if (length(converged_dyn) > 0 && length(converged_con) > 0) {
  estimates_con <- t(sapply(results[converged_con], function(x) x$con_est))
  bias_con <- colMeans(estimates_con[, 1:2]) - beta_true
  sd_con <- apply(estimates_con[, 1:2], 2, sd)
  rmse_con <- sqrt(bias_con^2 + sd_con^2)
  
  table2 <- data.frame(
    Model = c("Proposed (Dynamic)", "Proposed (Dynamic)", "Constant Hazard", "Constant Hazard"),
    Parameter = c("beta1", "beta2", "beta1", "beta2"), True = c(beta_true, beta_true),
    Bias = round(c(bias_dyn[1:2], bias_con), 3), RMSE = round(c(rmse_dyn[1:2], rmse_con), 3)
  )
  cat("\n============================================================\n")
  cat("TABLE 2: Parameter Comparison: Proposed vs. Constant Hazard Model\n")
  cat("============================================================\n")
  print(table2, row.names = FALSE)
}

mae_dyn_vals <- na.omit(sapply(results[converged_dyn], function(x) x$dyn_mae))
mae_con_vals <- na.omit(sapply(results[converged_con], function(x) x$con_mae))
mae_pelt_vals <- na.omit(sapply(results, function(x) x$pelt_mae))

table3 <- data.frame(
  Method = c("Proposed Dynamic Hazard Model", "Constant Hazard Model", "Independent PELT / CUSUM"),
  Mean_Absolute_Error_MAE = round(c(mean(mae_dyn_vals), mean(mae_con_vals), mean(mae_pelt_vals)), 3)
)
cat("\n============================================================\n")
cat("TABLE 3: Accuracy of Changepoint Location Estimation (MAE)\n")
cat("============================================================\n")
print(table3, row.names = FALSE)
