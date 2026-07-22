
## ----eval = TRUE---------------------------------------------------------------------------------------------------------------------------------------------------------
# remotes::install_github("davidbolin/rspde", ref = "devel")
# remotes::install_github("davidbolin/metricgraph", ref = "devel")
library(rSPDE)
library(MetricGraph)
library(grateful)

library(ggplot2)
library(reshape2)
library(plotly)
library(slackr)
source(here::here("keys.R"))
slackr_setup(token = token) # token comes from keys.R

dividerh_vector <- c(12)
dividertau_vector <- c(21, 24, 30)

for (iii in 1:length(dividerh_vector)) {
  for (jjj in 1:length(dividertau_vector)) {
    dividerh <- dividerh_vector[iii]
    dividertau <- dividertau_vector[jjj]

## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Parameters
T_final <- 2
kappa <- 4
N_finite = 4 # choose even
adjusted_N_finite <- N_finite + N_finite/2 + 1
# Coefficients for u_0 and f
coeff_U_0 <- 50*(1:adjusted_N_finite)^-1
coeff_U_0[-5] <- 0
coeff_FF <- rep(0, adjusted_N_finite)
coeff_FF[7] <- 10

AAA = 1
OMEGA = pi

# Time step and mesh size
m_vector <- c(1, 2, 3, 4)
alpha_vector <- seq(1, 1.8, by = 0.2)
epsilon <- 1e-10


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Overkill parameters
overkill_time_step <- 5e-6 #0.1 * 2^-14
overkill_h <- 9e-4#(0.1 * 2^-14)^(1/2)
# Finest time and space mesh
overkill_time_seq <- seq(0, T_final, length.out = ((T_final - 0) / overkill_time_step + 1))
overkill_graph <- gets.graph.tadpole(h = overkill_h)
# Compute the weights on the finest mesh
overkill_graph$compute_fem() # This is needed to compute the weights
overkill_C <- overkill_graph$mesh$C


# ----eval = FALSE, class.source = "fold-show"----------------------------------------------------------------------------------------------------------------------------
# Create a matrix to store the errors
errors_projected <- matrix(NA, nrow = length(m_vector), ncol = length(alpha_vector))
for (j in 1:length(alpha_vector)) {
  alpha <- alpha_vector[j] # from 0.5 to 2
  beta <- alpha / 2

  # Compute the eigenvalues and eigenfunctions on the finest mesh
  auxX <- gets.eigen.params(N_finite = N_finite, kappa = kappa, alpha = alpha, graph = overkill_graph)
  EIGENVAL_ALPHA <- auxX$EIGENVAL_ALPHA # Eigenvalues (they are independent of the meshes)
  auxY <- auxX$EIGENFUN # Eigenfunctions on the finest mesh

  # Compute the true solution on the finest mesh
  auxX <- auxY %*%
    outer(1:length(coeff_U_0),
          1:length(overkill_time_seq),
          function(i, j) (coeff_U_0[i] + coeff_FF[i] * G_sin(t = overkill_time_seq[j], A = AAA, lambda_j_alpha_half = EIGENVAL_ALPHA[i], omega = OMEGA)) * exp(-EIGENVAL_ALPHA[i] * overkill_time_seq[j]))

  for (i in 1:length(m_vector)) {
    m <- m_vector[i]
    sol <- solve_htaunew(alpha, m, epsilon, tol = 1e-16)
    h <- sol$h/dividerh
    time_step <- sol$tau/dividertau
    graph <- gets.graph.tadpole(h = h)
    graph$compute_fem()
    L <- graph$mesh$G
    C <- graph$mesh$C
    L <- kappa^2*C + L
    U_0 <- gets.eigen.params(N_finite = N_finite, kappa = kappa, alpha = alpha, graph = graph)$EIGENFUN
    U_0 <- U_0 %*% coeff_U_0 # Compute U_0 on the current mesh
    Psi <- graph$fem_basis(overkill_graph$get_mesh_locations())

    time_seq <- seq(0, T_final, length.out = ((T_final - 0) / time_step + 1))
    my_op_frac <- my.fractional.operators.frac(L, beta, C, scale.factor = kappa^2, m = m, time_step)
    # Compute matrix F with columns F^k
    U_approx <- t(Psi) %*% overkill_C %*% auxY %*%
      outer(1:length(coeff_FF),
            1:length(time_seq),
        function(i, j) coeff_FF[i] * g_sin(r = time_seq[j], A = AAA, omega = OMEGA))

    U_approx <- Psi %*% solve_fractional_evolution(my_op_frac, time_step, time_seq, val_at_0 = U_0, RHST = U_approx)
    U_approx <- auxX - construct_piecewise_projection(U_approx, time_seq, overkill_time_seq)
    errors_projected[i,j] <- sqrt(as.double(overkill_time_step * sum(U_approx * (overkill_C %*% U_approx)))) / (time_step^(-2) * h^(alpha - 2))
    slackr_msg(text = paste0("m =", m, ", alpha =", alpha, ", h =", h, ", time_step =", time_step), channel = "#research")
  }
}

## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
observed_rates <- numeric(length(alpha_vector))
for (u in 1:length(alpha_vector)) {
  observed_rates[u] <- coef(lm(log(errors_projected[, u]) ~ sqrt(m_vector)))[2]
}

theoretical_rates <- - pi * 2 * sqrt((1 - alpha_vector / 2))

p_m_tadpole_with_logtau_and_tauminustwo <- error.convergence.plotter(x_axis_vector = m_vector, 
                                 alpha_vector, 
                                 errors_projected, 
                                 theoretical_rates, 
                                 observed_rates,
                                 line_equation_fun = exp_line_equation,
                                 fig_title = expression("Convergence in " * italic(m)),
                                 x_axis_label = expression(italic(m)),
                                 apply_sqrt = TRUE)

ggsave(here::here("old/conv_rates_m_tadpole_with_logtau_and_tauminustwo.png"), width = 4, height = 5, plot = p_m_tadpole_with_logtau_and_tauminustwo, dpi = 300)


## ----eval = TRUE, echo = FALSE-------------------------------------------------------------------------------------------------------------------------------------------
initial_comment <- paste0("Here’s the latest plot update for the convergence in m")# with kappa = ", kappa, " and divider = ", divider)
# option 1
slackr_upload(
  filename = here::here("old/conv_rates_m_tadpole_with_logtau_and_tauminustwo.png"),        # path to your image
  initial_comment = initial_comment,
  channels = "#research"
)
get_numeric_scalars <- function() {
  vars <- ls(envir = .GlobalEnv)
  scalars <- sapply(vars, function(v) {
    val <- get(v, envir = .GlobalEnv)
    is.numeric(val) && length(val) == 1
  })
  mget(vars[scalars], envir = .GlobalEnv)
}

# Collect the scalars
num_scalars <- get_numeric_scalars()

# Format scalars
scalar_msg <- paste(
  sprintf("%s = %s", names(num_scalars), unlist(num_scalars)),
  collapse = "\n"
)

# Send one Slack message
slackr::slackr_msg(
  text = paste(
    "Numeric scalars in this run:\n",
    scalar_msg
  ),
  channel = "#research"
)

  }
}
