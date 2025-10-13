library(pbmcapply)
library(ggplot2)
library(scales)
library(dplyr)
library(patchwork)
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
tau_vector <- c(0.001, 0.01, 0.1)
N_vector <- T_final / tau_vector
kappa_vector <- c(4)
beta_vector <- c(0.5, 0.7, 0.9)
mu_vector <- c(1)

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

save(results, file = here::here("data_files/contraction_constant_results1.RData"))





T_final <- 2

load(here::here("data_files/contraction_constant_results1.RData"))

all_results <- results %>% 
  mutate(upperbound = 1/(mu*kappa^(4*beta))) %>%
  mutate(T_final = T_final) %>%
  mutate(N = T_final/tau) %>%
  mutate(L_cless1 = L_c < 1) %>%
  mutate(Tlessmukappa2beta = T_final < (mu * kappa^(2*beta))) %>%
  mutate(onelessmukappa4beta  = upperbound < 1) 


# Find combined y-limits
ymin <- min(all_results$L_c, all_results$upperbound, na.rm = TRUE)
ymax <- max(all_results$L_c, all_results$upperbound, na.rm = TRUE)


p <- ggplot(all_results, aes(x = N, y = L_c,
                             color = factor(beta),
                             shape = factor(kappa),
                             linetype = factor(mu))) +
  geom_point(size = 3) +
  geom_line(aes(group = interaction(beta, kappa, mu)), alpha = 1) +
  scale_x_continuous() +
  scale_y_continuous(
    limits = c(ymin, ymax),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = label_number(accuracy = 0.0001, trim = TRUE)
  ) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  labs(x = "N",
       y = expression(gamma),
       color = expression(beta),
       shape = expression(kappa),
       linetype = expression(mu)) +
  theme_minimal(base_size = 14, base_family = "Palatino")

q <- ggplot(all_results, aes(x = N, y = upperbound,
                             color = factor(beta),
                             shape = factor(kappa),
                             linetype = factor(mu))) +
  geom_point(size = 3) +
  geom_line(aes(group = interaction(beta, kappa, mu)), alpha = 1) +
  scale_x_continuous() +
  scale_y_continuous(
    limits = c(ymin, ymax),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = label_number(accuracy = 0.0001, trim = TRUE),
    position = "right"
  ) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  labs(x = "N",
       y = expression(1 / mu * kappa^{4 * beta}),
       color = expression(beta),
       shape = expression(kappa),
       linetype = expression(mu)) +
  theme_minimal(base_size = 14, base_family = "Palatino")


combined1 <- (p | q) +
  plot_annotation(
    title = expression("Contraction constant " * gamma * " and its upper bound " * 1 / mu * kappa^{4 * beta}),
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, 
                                            family = "Palatino"))
  ) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

combined1

ggsave(
  here::here("data_files/fixedpointconvergence_combined1.png"),
  width = 12, height = 6, plot = combined1, dpi = 300
)



# --- 3. Parameter grid ---
T_final <- 2
tau_vector <- c(0.001, 0.01, 0.1)
N_vector <- T_final / tau_vector
kappa_vector <- c(1, 4, 16)
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

save(results, file = here::here("data_files/contraction_constant_results2.RData"))



T_final <- 2

load(here::here("data_files/contraction_constant_results2.RData"))

all_results <- results %>% 
  mutate(upperbound = 1/(mu*kappa^(4*beta))) %>%
  mutate(T_final = T_final) %>%
  mutate(N = T_final/tau) %>%
  mutate(L_cless1 = L_c < 1) %>%
  mutate(Tlessmukappa2beta = T_final < (mu * kappa^(2*beta))) %>%
  mutate(onelessmukappa4beta  = upperbound < 1) 


# Find combined y-limits
ymin <- min(all_results$L_c, all_results$upperbound, na.rm = TRUE)
ymax <- max(all_results$L_c, all_results$upperbound, na.rm = TRUE)

p <- ggplot(all_results, aes(x = N, y = L_c,
                             color = factor(kappa),
                             shape = factor(beta),
                             linetype = factor(mu))) +
  geom_point(size = 3) +
  geom_line(aes(group = interaction(kappa, beta, mu)), alpha = 1) +
  scale_x_log10() +
  scale_y_log10(
    limits = c(ymin, ymax),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = label_scientific()
  ) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  labs(x = "N",
       y = expression(gamma),
       color = expression(kappa),
       shape = expression(beta),
       linetype = expression(mu)) +
  theme_minimal(base_size = 14, base_family = "Palatino")

q <- ggplot(all_results, aes(x = N, y = upperbound,
                             color = factor(kappa),
                             shape = factor(beta),
                             linetype = factor(mu))) +
  geom_point(size = 3) +
  geom_line(aes(group = interaction(kappa, beta, mu)), alpha = 1) +
  scale_x_log10() +
  scale_y_log10(
    limits = c(ymin, ymax),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = label_scientific(),
    position = "right"
  ) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  labs(x = "N",
       y = expression(1 / mu * kappa^{4 * beta}),
       color = expression(kappa),
       shape = expression(beta),
       linetype = expression(mu)) +
  theme_minimal(base_size = 14, base_family = "Palatino")


combined2 <- (p | q) + 
  plot_annotation(
    title = expression("Contraction constant " * gamma * " and its upper bound " * 1 / mu * kappa^{4 * beta}),
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, 
                                            family = "Palatino"))
  ) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
combined2


ggsave(
  here::here("data_files/fixedpointconvergence_combined2.png"),
  width = 12, height = 7, plot = combined2, dpi = 300
)















