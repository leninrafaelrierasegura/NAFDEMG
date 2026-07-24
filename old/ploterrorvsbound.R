library(slackr)
source(here::here("keys.R"))
slackr_setup(token = token) # token comes from keys.R

time_step_vector <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1)
h_vector <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1)
m_vector <- c(1, 2, 3, 4)
alpha_vector <- c(0.5, 0.75, 1, 1.25, 1.5, 1.75, 2)

df <- expand.grid(
  time_step = time_step_vector,
  h = h_vector,
  m = m_vector,
  alpha = alpha_vector
)

df$error <- NA_real_
df$boundsharp <- NA_real_
df$boundloose <- NA_real_


for (i in seq_len(nrow(df))) {
  
  time_step <- df$time_step[i]
  h <- df$h[i]
  m <- df$m[i]
  alpha <- df$alpha[i]
  

  T_final <- 2
  kappa <- 15
  beta <- alpha/2
  N_finite = 4 # choose even
  adjusted_N_finite <- N_finite + N_finite/2 + 1
  # Coefficients for u_0 and f
  coeff_U_0 <- 50*(1:adjusted_N_finite)^-1
  coeff_U_0[-5] <- 0
  coeff_FF <- rep(0, adjusted_N_finite)
  coeff_FF[7] <- 10
  
  
  time_seq <- seq(0, T_final, by = time_step)
  graph <- gets.graph.tadpole(h = h)
  graph$compute_fem()
  L <- graph$mesh$G
  C <- graph$mesh$C
  L <- kappa^2*C + L
  
  HH <- gets.eigen.params(N_finite = N_finite, kappa = kappa, alpha = alpha, graph = graph)
  EIGENVAL_ALPHA <- HH$EIGENVAL_ALPHA
  HH <- HH$EIGENFUN
  
  AAA = 1
  OMEGA = pi
  
  U_0 <- HH %*% coeff_U_0
  
  U_true <- HH %*% 
    outer(1:length(coeff_U_0), 
          1:length(time_seq), 
          function(i, j) (coeff_U_0[i] + coeff_FF[i] * G_sin(t = time_seq[j], A = AAA, lambda_j_alpha_half = EIGENVAL_ALPHA[i], omega = OMEGA) ) * exp(-EIGENVAL_ALPHA[i] * time_seq[j]))
  
  overkill_graph <- gets.graph.tadpole(h = 0.001)
  overkill_graph$compute_fem()
  AUX <- gets.eigen.params(N_finite = N_finite, kappa = kappa, alpha = alpha, graph = overkill_graph)$EIGENFUN
  
  
  AUX <- t(graph$fem_basis(overkill_graph$get_mesh_locations())) %*%
    overkill_graph$mesh$C %*%
    AUX %*%
    outer(1:length(coeff_FF), 
          1:length(time_seq), 
          function(i, j) coeff_FF[i] * g_sin(r = time_seq[j], A = AAA, omega = OMEGA))
  
  # Numerical solution
  my_op_frac <- my.fractional.operators.frac(L, beta, C, scale.factor = kappa^2, m = m, time_step)
  
  
  AUX <- solve_fractional_evolution(my_op_frac, time_step, time_seq, val_at_0 = U_0, RHST = AUX)
  AUX <- U_true - AUX
  
  error <- sqrt(as.double(time_step * sum(AUX * (C %*% AUX))))
  df$error[i] <- error
  df$boundsharp[i] <- time_step + h^alpha + h^(alpha - 2) * exp(- 2 * pi * sqrt((1 - alpha / 2) * m))
  df$boundloose[i] <- time_step + h^alpha * abs(log(time_step)) + time_step^(-2) * h^(alpha - 2) * exp(- 2 * pi * sqrt((1 - alpha / 2) * m))
  slackr_msg(text = paste0("i = ", i, ", m = ", m, ", alpha = ", alpha, ", h = ", h, ", time_step = ", time_step), channel = "#research")
  print(paste0("i = ", i, ", m = ", m, ", alpha = ", alpha, ", h = ", h, ", time_step = ", time_step))
}


save(df, file = here::here("old/edf_error_vs_bounds.RData"))

load(here::here("old/edf_error_vs_bounds.RData"))

df$index <- seq_len(nrow(df))


library(ggplot2)


ggplot(df, aes(x = index)) +
  geom_line(aes(y = error, color = "Error"), linewidth = 0.8) +
  geom_line(aes(y = boundsharp, color = "Sharp bound"), linewidth = 0.8) +
  geom_line(aes(y = boundloose, color = "Loose bound"), linewidth = 0.8) +
  geom_point(aes(y = error, color = "Error"), size = 1) +
  geom_point(aes(y = boundsharp, color = "Sharp bound"), size = 1) +
  geom_point(aes(y = boundloose, color = "Loose bound"), size = 1) +
  scale_y_log10() +
  labs(
    x = "Index",
    y = "Value",
    color = NULL
  ) +
  theme_bw()


library(dplyr)
library(tidyr)
library(ggplot2)

df_long <- df %>%
  pivot_longer(
    cols = c(error, boundsharp, boundloose),
    names_to = "curve",
    values_to = "value"
  )

for (a in unique(df$alpha)) {
  
  p <- df_long %>%
    filter(alpha == a) %>%
    ggplot(aes(x = h,
               y = value,
               color = curve,
               linetype = curve)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    scale_x_log10() +
    scale_y_log10() +
    facet_grid(
      m ~ time_step,
      labeller = label_bquote(
        rows = m == .(m),
        cols = tau == .(time_step)
      )
    ) +
    labs(
      title = bquote(alpha == .(a)),
      x = "h",
      y = "Value",
      color = "",
      linetype = ""
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey95"),
      panel.grid.minor = element_blank()
    )
  
  print(p)
  
  ggsave(
    filename = paste0("alpha_", a, ".pdf"),
    plot = p,
    width = 18,
    height = 10
  )
}



library(dplyr)
library(tidyr)
library(ggplot2)

df_long <- df %>%
  pivot_longer(
    cols = c(error, boundsharp, boundloose),
    names_to = "curve",
    values_to = "value"
  )

for (a in unique(df$alpha)) {
  
  p <- df_long %>%
    filter(alpha == a) %>%
    ggplot(aes(x = time_step,
               y = value,
               color = curve,
               linetype = curve)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    scale_x_log10() +
    scale_y_log10() +
    facet_grid(
      m ~ h,
      labeller = label_bquote(
        rows = m == .(m),
        cols = h == .(h)
      )
    ) +
    labs(
      title = bquote(alpha == .(a)),
      x = expression(tau),
      y = "Value",
      color = "",
      linetype = ""
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey95"),
      panel.grid.minor = element_blank()
    )
  
  print(p)
  
  ggsave(
    filename = paste0("alpha_", a, "_vs_tau.pdf"),
    plot = p,
    width = 18,
    height = 10
  )
}


library(dplyr)
library(tidyr)
library(ggplot2)

df_long <- df %>%
  pivot_longer(
    cols = c(error, boundsharp, boundloose),
    names_to = "curve",
    values_to = "value"
  )

for (a in unique(df$alpha)) {
  
  p <- df_long %>%
    filter(alpha == a) %>%
    ggplot(aes(x = m,
               y = value,
               color = curve,
               linetype = curve,
               group = curve)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    scale_y_log10() +
    facet_grid(
      time_step ~ h,
      labeller = label_bquote(
        rows = tau == .(time_step),
        cols = h == .(h)
      )
    ) +
    labs(
      title = bquote(alpha == .(a)),
      x = "m",
      y = "Value",
      color = "",
      linetype = ""
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey95"),
      panel.grid.minor = element_blank()
    )
  
  print(p)
  
  ggsave(
    filename = paste0("alpha_", a, "_vs_m.pdf"),
    plot = p,
    width = 18,
    height = 10
  )
}


files <- c(
  paste0("alpha_", alpha_vector, "_vs_m.pdf")
)

for (f in files) {
  slackr_upload(
    file = f,
    channels = "#research"
  )
}


