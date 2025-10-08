library(dplyr)
make_matrix <- function(a, N) {
  M <- outer(1:N, 1:N, function(i, j) a^abs(i - j))
  return(M)
}

build_T <- function(a, N){
  M <- make_matrix(a, N) 
  R <- M*0
  for (k in N:1) {
    aux <- diag(c(rep(1, k), rep(0, N - k)))
    temp <- a^(2*(N - k + 1)) * aux %*% M %*% aux
    R <- R + temp
  }
  return(R)
}
build_block_matrix <- function(egs, Nh) {
  N <- length(egs)
  # total size
  m <- N * Nh
  AA <- matrix(0, nrow = m, ncol = m)
  
  for (j in 1:N) {
    rows <- ((j - 1) * Nh + 1):(j * Nh)
    cols <- ((j - 1) * Nh + 1):(j * Nh)
    AA[rows, cols] <- build_T(egs[j], Nh)
  }
  return(AA)
}
build_perm_matrix_general <- function(N, Nh) {
  m <- N * Nh
  # target indices
  perm <- as.vector(sapply(1:Nh, function(i) i + Nh*(0:(N-1))))
  P <- diag(m)[perm, ]
  return(P)
}

build_perm_matrix <- function(N, Nh) {
  # total size
  m <- N * Nh
  
  # initialize zero matrix
  P <- matrix(0, nrow = m, ncol = m)
  
  # fill according to permutation rule
  for (j in 1:N) {
    for (i in 1:Nh) {
      k <- (j - 1) * Nh + i     # row index
      l <- (i - 1) * N + j      # column index
      P[k, l] <- 1
    }
  }
  return(P)
}

Nh <- 3
egs <- c(1, 3)
P <- build_perm_matrix_general(length(egs), Nh)

AA <- build_block_matrix(egs, Nh)
AA
P %*% AA %*% t(P)
P
build_perm_matrix(length(egs), Nh)
# P <- matrix(c(
#   1,0,0,0,0,0,0,0,0,
#   0,0,0,1,0,0,0,0,0,
#   0,0,0,0,0,0,1,0,0,
#   0,1,0,0,0,0,0,0,0,
#   0,0,0,0,1,0,0,0,0,
#   0,0,0,0,0,0,0,1,0,
#   0,0,1,0,0,0,0,0,0,
#   0,0,0,0,0,1,0,0,0,
#   0,0,0,0,0,0,0,0,1
# ), nrow = 9, byrow = TRUE)
# 
# P

Nh <- 3
eg1 <- 1
eg2 <- 3
eg3 <- 5
AA <- cbind(rbind(build_T(eg1,Nh), build_T(eg1,Nh)*0, build_T(eg1,Nh)*0),
            rbind(build_T(eg2,Nh)*0, build_T(eg2,Nh), build_T(eg2,Nh)*0),
            rbind(build_T(eg3,Nh)*0, build_T(eg3,Nh)*0, build_T(eg3,Nh)))
AA

P %*% AA %*% P



T_final <- 2


tau_vector   <- c(0.005, 0.01, 0.5)
kappa_vector <- c(1, 10, 20)
beta_vector  <- c(0.5, 0.7, 0.9)
mu_vector    <- c(1, 5, 10)

results <- expand.grid(tau = tau_vector,
                       kappa = kappa_vector,
                       beta = beta_vector,
                       mu = mu_vector) %>%
  mutate(N = T_final/tau)

for (i in 1:nrow(results)) {
  N <- results$N[i]
  tau <- results$tau[i]
  kappa <- results$kappa[i]
  beta <- results$beta[i]
  mu <- results$mu[i]
  
  omega <- 1/(1+tau * kappa^(2*beta))
  mat <- build_T(omega, N)
  results$L_c[i] <- (tau^2 * norm(mat, type = "2")) / mu
}




library(ggplot2)
library(scales)

p <- ggplot(results, aes(x = N, y = L_c,
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

p

ggsave(here::here("data_files/fixedpointconvergence_sharp.png"), width = 7, height = 7, plot = p, dpi = 300)
