library(pbmcapply)

# check contraction constant
make_matrix <- function(a, N) {
  v <- a^(0:(N - 1))
  toeplitz(v)
}


# --- 2. Efficient T builder with precomputed powers ---
build_T_fast <- function(a, N) {
  M <- make_matrix(a, N)
  R <- matrix(0, N, N)
  coef <- (a^2)^(1:N)  # correct power order
  for (k in N:1) {
    R[1:k, 1:k] <- R[1:k, 1:k] + coef[N - k + 1] * M[1:k, 1:k]
  }
  R
}



# --- 3. Parameter grid ---
T_final <- 2
tau_vector <- c(0.0001, 0.05, 0.1, 0.5)
N_vector <- T_final / tau_vector
kappa_vector <- c(0.5, 1, 10)
beta_vector <- c(0.5, 0.7, 0.9)
mu_vector <- c(0.1, 1, 10)

results <- expand.grid(
  tau = tau_vector,
  kappa = kappa_vector,
  beta = beta_vector,
  mu = mu_vector
)

# --- 4. Parallel computation with caching (A depends on ω and N only) ---
cache <- new.env(hash = TRUE)

results$L_c <- unlist(pbmclapply(
  1:nrow(results),
  function(i) {
    tau <- results$tau[i]
    kappa <- results$kappa[i]
    beta <- results$beta[i]
    mu <- results$mu[i]
    
    omega <- 1 / (1 + tau * kappa^(2 * beta))
    N <- as.integer(T_final / tau)
    key <- paste0(round(omega, 10), "_", N)
    
    if (exists(key, envir = cache)) {
      lambda_max <- get(key, envir = cache)
    } else {
      A <- build_T_fast(omega, N)
      lambda_max <- max(eigen(A, symmetric = TRUE, only.values = TRUE)$values)
      assign(key, lambda_max, envir = cache)
    }
    
    (tau^2 * lambda_max) / mu
  },
  mc.cores = parallel::detectCores()
))

save(results, file = here::here("data_files/contraction_constant_results.RData"))