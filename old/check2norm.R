source("functionality.R")

T_final <- 2
tau <- 0.02
h <- 0.02
kappa <- 4
alpha <- 1.8
m = 2
beta <- alpha/2
mu <- 0.1

scale.factor <- kappa^2
scaling <- scale.factor^beta
roots <- my.get.roots(m, beta)
poles_rs_k <- compute.partial.fraction.param(roots$factor, roots$pr_roots, roots$pl_roots, tau, scaling)

poles <- poles_rs_k$p
residues <- poles_rs_k$r

omega1 <- 1/(1+tau*kappa^(2*beta))
omega1
omega <- sum(residues/(1-poles))
omega


for (m in 1:8) {
  roots <- my.get.roots(m, beta)
  poles_rs_k <- compute.partial.fraction.param(roots$factor, roots$pr_roots, roots$pl_roots, tau, scaling)
  
  poles <- poles_rs_k$p
  residues <- poles_rs_k$r
  
  omega_m <- sum(residues/(1-poles))
  print(omega_m)
}

N <- T_final/tau

assemble_matrix <- function(N, omega) {
  i <- row(matrix(0, N, N))
  j <- col(matrix(0, N, N))
  A <- ifelse(j >= i, (j - i + 1) * omega^(j - i + 2), 0)
  return(A)
}

N <- 10
omega <- 1
A <- assemble_matrix(N, omega)

library(RSpectra)
eig <- eigs_sym(t(A) %*% A, 1)   
norm2 <- sqrt(eig$values)







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


build_T_fast(2,5)

T_final <- 2
tau_vector <- c(0.0001, 0.05, 0.1, 0.5)
N_vector <- T_final/tau_vector
kappa_vector <- c(0.5, 1, 10)
beta_vector <- c(0.5, 0.7, 0.9)
mu_vector <- c(0.1, 1, 10)

results <- expand.grid(tau = tau_vector, kappa = kappa_vector, beta = beta_vector, mu = mu_vector)

library(parallel)
results$L_c <- mclapply(1:nrow(results), function(i) {
  tau <- results$tau[i]
  kappa <- results$kappa[i]
  beta <- results$beta[i]
  mu <- results$mu[i]
  
  omega <- 1 / (1 + tau * kappa^(2 * beta))
  N <- T_final / tau
  print(i)
  A <- build_T_fast(omega, N)
  norm2 <- max(eigen(A, symmetric = TRUE, only.values = TRUE)$values)
  (tau^2 * norm2) / mu
}, mc.cores = parallel::detectCores())

save(results, file = here::here("data_files/contraction_constant_results.RData"))


load(here::here("data_files/contraction_constant_results.RData"))
all_results <- results %>% 
  mutate(upperbound = 1/(mu*kappa^(4*beta))) %>%
  mutate(limit = upperbound * (1 - exp(-T_final * kappa^(2*beta))*(1+T_final*kappa^(2*beta)))) %>%
  mutate(limit = limit/2) %>%
  mutate(T_final = T_final) %>%
  mutate(N = T_final/tau) %>%
  mutate(L_cless1 = L_c < 1) %>%
  mutate(Tlessmukappa2beta = T_final < (mu * kappa^(2*beta))) %>%
  mutate(onelessmukappa4beta  = upperbound < 1) 


library(ggplot2)
library(scales)

p <- ggplot(all_results, aes(x = N, y = L_c,
                             color = factor(kappa),
                             shape = factor(beta),
                             linetype = factor(mu))) +
  geom_point(size = 3) +
  geom_line(aes(group = interaction(kappa, beta, mu)), alpha = 1) +
  # scale_x_log10() +
  # scale_y_log10(
  #   breaks = trans_breaks("log10", function(x) 10^x),
  #   labels = label_number(accuracy = 0.0001, trim = TRUE)
  # ) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  labs(x = "N",
       y = expression(L[c]),
       color = expression(kappa),
       shape = expression(beta),
       linetype = expression(mu)) +
  theme_minimal(base_size = 14)

q <- ggplot(all_results, aes(x = N, y = limit, # upperbound
                             color = factor(kappa),
                             shape = factor(beta),
                             linetype = factor(mu))) +
  geom_point(size = 3) +
  geom_line(aes(group = interaction(kappa, beta, mu)), alpha = 1) +
  # scale_x_log10() +
  # scale_y_log10(
  #   breaks = trans_breaks("log10", function(x) 10^x),
  #   labels = label_number(accuracy = 0.0001, trim = TRUE)
  # ) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  labs(x = "N",
       y = expression(frac(1, mu * kappa^{4 * beta})),
       color = expression(kappa),
       shape = expression(beta),
       linetype = expression(mu)) +
  theme_minimal(base_size = 14)

library(patchwork)

combined <- (p | q) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
combined

ggsave(here::here("fixedpointconvergence_sharp.png"), width = 7, height = 7, plot = p, dpi = 300)

ggsave(here::here("fixedpointconvergence_approx.png"), width = 7, height = 7, plot = q, dpi = 300)

ggsave(
  here::here("fixedpointconvergence_combined.png"),
  width = 12, height = 6, plot = combined, dpi = 300
)
