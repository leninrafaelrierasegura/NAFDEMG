## ---------------------------------------------------------------------------
convergence_settings <- function(study, profile = "pilot") {
  study <- match.arg(study, c("h", "tau", "m"))
  profile <- match.arg(profile, c("pilot", "full"))

  # Shared scientific experiment: change these here, not in individual notebooks.
  T_final <- 0.5
  diffusion <- 1e-4
  kappa <- 4
  N_finite <- 4L
  active_modes <- c(5L, 7L)
  exact_initial <- c(10, 0)
  exact_forcing <- c(0, 10)
  AAA <- 1
  OMEGA <- pi
  alpha_vector <- c(0.5, 0.9, 1.3)
  range_window <- c(1e-4, 1e-2)
  range_points <- if (profile == "pilot") 6L else 8L
  study_range <- 10^seq(log10(range_window[1]), log10(range_window[2]),
                        length.out = range_points)

  # Shared numerical accuracy and execution settings.
  solver_backend <- "uniform-modal"
  root_interpolation <- "cubic"
  quadrature_order <- 4L
  step_sum_order <- if (profile == "pilot") 64L else 128L
  run_refinement_checks <- TRUE
  refinement_tolerance <- 0.05
  workers <- if (profile == "pilot") 4L else 8L
  temporal_cache_mib <- 512
  memory_budget_gib <- NULL

  # Auxiliary grids isolate the variable studied; they need not be identical.
  m_fixed <- 4L
  h_study_tau <- if (profile == "pilot") 1e-10 else 5e-11
  tau_study_h <- 1e-5
  m_study_h <- if (profile == "pilot") 1e-6 else 5e-7
  m_study_tau <- if (profile == "pilot") 1e-11 else 5e-12
  m_vector <- 1:4

  h_vector <- uniform_tadpole_h(study_range)
  tau_targets <- study_range
  time_grids <- lapply(tau_targets, function(dt) experiment_grid(T_final, dt))
  tau_vector <- vapply(time_grids, function(g) g$time_step, numeric(1))
  stopifnot(length(unique(tau_vector)) == range_points,
            max(tau_vector)/min(tau_vector) >= 100*(1-1e-12),
            max(h_vector)/min(h_vector) >= 100*(1-1e-12))
  h_fixed <- uniform_tadpole_h(if (study == "tau") tau_study_h else m_study_h)
  time_grid <- experiment_grid(T_final, if (study == "h") h_study_tau else m_study_tau)
  axis_values <- switch(study, h = h_vector, tau = tau_vector, m = m_vector)
  anchor_h <- if (study == "h") min(h_vector) else h_fixed
  anchor_grid <- if (study == "tau") time_grids[[which.min(tau_vector)]] else time_grid
  anchor_m <- if (study == "m") max(m_vector) else m_fixed

  shared_parameters <- list(T_final = T_final, diffusion = diffusion, kappa = kappa,
    modes = active_modes, N_finite = N_finite, alpha = alpha_vector,
    initial = exact_initial, forcing = exact_forcing, A = AAA, omega = OMEGA,
    range_window = range_window, range_points = range_points,
    solver_backend = solver_backend, root_interpolation = root_interpolation,
    quadrature_order = quadrature_order, step_sum_order = step_sum_order,
    refinement_tolerance = refinement_tolerance, run_refinement_checks = run_refinement_checks,
    auxiliary = list(h_study_tau = h_study_tau, tau_study_h = tau_study_h,
                     m_study_h = m_study_h, m_study_tau = m_study_tau, m_fixed = m_fixed),
    error_method = "uniform-modal-projection-v2")
  suite_id <- digest::digest(shared_parameters, algo = "sha256")
  fixed <- switch(study,
    h = list(m = m_fixed, time_step = time_grid$time_step),
    tau = list(h = h_fixed, m = m_fixed),
    m = list(h = h_fixed, time_step = time_grid$time_step))
  experiment_parameters <- list(study = study, profile = profile,
    diffusion = diffusion, shared = shared_parameters, suite_id = suite_id,
    requested_range = study_range, axis = axis_values, fixed = fixed)
  execution_settings <- list(workers = workers, temporal_cache_mib = temporal_cache_mib,
                             memory_budget_gib = memory_budget_gib)
  list(experiment_profile = profile, T_final = T_final, diffusion = diffusion, kappa = kappa, N_finite = N_finite,
    active_modes = active_modes, exact_initial = exact_initial, exact_forcing = exact_forcing,
    AAA = AAA, OMEGA = OMEGA, alpha_vector = alpha_vector, study_range = study_range,
    h_vector = h_vector, tau_targets = tau_targets, tau_vector = tau_vector,
    m_vector = m_vector, h_fixed = h_fixed, m_fixed = m_fixed, time_grid = time_grid,
    time_grids = time_grids, axis_values = axis_values, anchor_h = anchor_h,
    anchor_grid = anchor_grid, anchor_m = anchor_m,
    solver_backend = solver_backend, root_interpolation = root_interpolation,
    quadrature_order = quadrature_order, step_sum_order = step_sum_order,
    run_refinement_checks = run_refinement_checks, refinement_tolerance = refinement_tolerance,
    workers = workers, temporal_cache_mib = temporal_cache_mib,
    memory_budget_gib = memory_budget_gib, shared_parameters = shared_parameters,
    experiment_parameters = experiment_parameters, execution_settings = execution_settings)
}

