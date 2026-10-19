## ------------------------------------------------------------------------------------------------------------------------------------------------
# remotes::install_github("davidbolin/rspde", ref = "devel")
# remotes::install_github("davidbolin/metricgraph", ref = "devel")
library(rSPDE)
library(MetricGraph)
library(grateful)

library(ggplot2)
library(reshape2)
library(plotly)


## ------------------------------------------------------------------------------------------------------------------------------------------------
# Functions to compute the exact solution to the fractional diffusion equation
g_linear <- function(r, A, lambda_j_alpha_half) {
  return(A * exp(-lambda_j_alpha_half * r))
  }
G_linear <- function(t, A) {
  return(A * t)
  }
g_exp <- function(r, A, mu) {
  return(A * exp(mu * r))
  }
G_exp <- function(t, A, lambda_j_alpha_half, mu) {
  exponent <- lambda_j_alpha_half + mu
  return(A * (exp(exponent * t) - 1) / exponent)
  }
g_poly <- function(r, A, n) {
  return(A * r^n)
}
G_poly <- function(t, A, lambda_j_alpha_half, n) {
  t <- as.vector(t)
  k_vals <- 0:n
  sum_term <- sapply(t, function(tt) {
    sum(((-lambda_j_alpha_half * tt)^k_vals) / factorial(k_vals))
  })
  coeff <- ((-1)^(n + 1)) * factorial(n) / (lambda_j_alpha_half^(n + 1))
  return(A * coeff * (1 - exp(lambda_j_alpha_half * t) * sum_term))
}
g_sin <- function(r, A, omega) {
  return(A * sin(omega * r))
}
G_sin <- function(t, A, lambda_j_alpha_half, omega) {
  denom <- lambda_j_alpha_half^2 + omega^2
  numerator <- exp(lambda_j_alpha_half * t) * (lambda_j_alpha_half * sin(omega * t) - omega * cos(omega * t)) + omega
  return(A * numerator / denom)
}
g_cos <- function(r, A, theta) {
  return(A * cos(theta * r)) 
}
G_cos <- function(t, A, lambda_j_alpha_half, theta) {
  denom <- lambda_j_alpha_half^2 + theta^2
  numerator <- exp(lambda_j_alpha_half * t) * (lambda_j_alpha_half * cos(theta * t) + theta * sin(theta * t)) - lambda_j_alpha_half
  return(A * numerator / denom)
}

