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
    print(norm(aux, type = "2"))
    temp <- a^(2*(N - k + 1)) * aux %*% M %*% aux
    #print(temp)
    R <- R + temp
  }
  return(R)
}
build_T(2,4)


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

Nh <- 2
egs <- c(3, 5)
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


# check norm identity
source(here::here("functionality.R"))
T_final <- 2
time_step <- 0.001 
h <- 1
kappa <- 15
alpha <- 0.5 
m = 1
beta <- alpha/2

graph <- gets.graph.tadpole(h = h)
graph$compute_fem()
G <- graph$mesh$G
C <- graph$mesh$C
L <- kappa^2*C + G
I <- Matrix::Diagonal(nrow(C))

# Numerical solution
obj <- my.fractional.operators.frac(L, beta, C, scale.factor = kappa^2, m = m, time_step)
partial_fraction_terms <- obj$partial_fraction_terms
residues <- obj$residues
output <- I*0
for (i in 1:(m+1)) {output <- output + residues[i] * solve(partial_fraction_terms[[i]], I)}
R <- output

C_sqrt <- expm::sqrtm(C)       # matrix square root
Omega <- C_sqrt %*% R %*% C_sqrt
n <- nrow(Omega)
Omega2 <- Omega %*% Omega
Omega3 <- Omega2 %*% Omega
Omega4 <- Omega2 %*% Omega2
Omega5 <- Omega3 %*% Omega2
Omega6 <- Omega3 %*% Omega3
B11 <- matrix(0, nrow = n, ncol = n)
B12 <- Omega2 + Omega4 + Omega6
B13 <- Omega3 + Omega5
B14 <- Omega4

B21 <- matrix(0, nrow = n, ncol = n)
B22 <- Omega3 + Omega5
B23 <- Omega2 + Omega4
B24 <- Omega3

B31 <- matrix(0, nrow = n, ncol = n)
B32 <- Omega4
B33 <- Omega3
B34 <- Omega2

B41 <- matrix(0, nrow = n, ncol = n)
B42 <- matrix(0, nrow = n, ncol = n)
B43 <- matrix(0, nrow = n, ncol = n)
B44 <- matrix(0, nrow = n, ncol = n)
library(Matrix)  # optional, for bdiag if needed

big_matrix <- rbind(
  cbind(B11, B12, B13, B14),
  cbind(B21, B22, B23, B24),
  cbind(B31, B32, B33, B34),
  cbind(B41, B42, B43, B44)
)

omega <- 1/(1+time_step * kappa^(2*beta))

TT <- build_T(omega, 3)

norm(TT, type = "2")
norm(big_matrix, type = "2")

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



# CHECK NESTEDNESS

V_d_aux <- C %*% u_d
f_aux <- elliptic_u %*% t(psi_prime) + elliptic_f %*% t(psi) - z_bar
F_proj_aux <- C %*% f_aux

WW <- V_d - V_d_aux
VV <- F_proj - F_proj_aux
sqrt(as.double(time_step * sum(WW * (C %*% WW))))
sqrt(as.double(time_step * sum(VV * (C %*% VV))))

g1 <- C %*% elliptic_f
g2 <- R %*% overkill_elliptic_f
SS <- g1 - g2
sqrt(as.double(1 * sum(SS * (C %*% SS))))


C_aux <- R %*% Psi
norm(C - C_aux, type = "2")


ai <- as.vector(Psi %*% EIGENFUN[,3])
bi <- as.vector(overkill_EIGENFUN[,3])
SS <- ai - bi
sqrt(as.double(1 * sum(SS * (Cstar %*% SS))))
library(ggplot2)

# Example: assume ai and bi are numeric vectors of the same length
df <- data.frame(
  x = seq_along(ai),
  ai = ai,
  bi = bi
)

p <- ggplot(df, aes(x = x)) +
  geom_line(aes(y = ai, color = "ai")) +
  geom_line(aes(y = bi, color = "bi")) +
  scale_color_manual(values = c("ai" = "red", "bi" = "blue")) +
  theme_minimal() +
  labs(x = "Index", y = "Value", color = "Series")
plotly::ggplotly(p)
