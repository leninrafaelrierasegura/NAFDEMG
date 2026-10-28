## -----------------------------------------------------------------------------
isolated_discrete_spectrum <- function(C, L, max_dimension = 1200L,
                                       max_dense_mib = 1024,
                                       residual_tolerance = 1e-8,
                                       orthogonality_tolerance = 1e-10) {
  n <- nrow(C)
  stopifnot(n >= 1L, ncol(C) == n, identical(dim(L), dim(C)),
            isSymmetric(C, tol = 1e-12), isSymmetric(L, tol = 1e-12),
            length(max_dimension) == 1L, is.finite(max_dimension), max_dimension >= 1,
            length(max_dense_mib) == 1L, is.finite(max_dense_mib), max_dense_mib > 0,
            is.finite(residual_tolerance), residual_tolerance > 0,
            is.finite(orthogonality_tolerance), orthogonality_tolerance > 0)
  estimated_mib <- 20 * 8 * as.double(n)^2 / 1024^2
  if (n > max_dimension || estimated_mib > max_dense_mib)
    stop("Complete eigendecomposition refused before dense allocation: dimension ", n,
         ", estimated temporary storage ", round(estimated_mib, 1),
         " MiB. Choose a moderate isolated-study mesh or explicitly raise the guards.")
  started <- proc.time()[["elapsed"]]
  dense_C <- as.matrix(C)
  dense_L <- as.matrix(L)
  stopifnot(all(is.finite(dense_C)), all(is.finite(dense_L)))
  R <- chol(dense_C)
  Rinv <- backsolve(R, diag(n))
  symmetric <- crossprod(Rinv, dense_L %*% Rinv)
  ee <- eigen((symmetric + t(symmetric)) / 2, symmetric = TRUE)
  stopifnot(all(is.finite(ee$values)), all(ee$values > 0))
  V <- backsolve(R, ee$vectors)
  CV <- C %*% V
  LV <- L %*% V
  residual <- as.matrix(LV - sweep(as.matrix(CV), 2L, ee$values, `*`))
  denominator <- sqrt(colSums(as.matrix(LV)^2)) +
    abs(ee$values) * sqrt(colSums(as.matrix(CV)^2))
  residual_error <- max(sqrt(colSums(residual^2)) / pmax(denominator, .Machine$double.xmin))
  orthogonality_error <- max(abs(crossprod(V, CV) - diag(n)))
  if (residual_error > residual_tolerance || orthogonality_error > orthogonality_tolerance)
    stop("Numerical eigensystem failed validation: residual = ", residual_error,
         ", C-orthogonality = ", orthogonality_error)
  list(values = ee$values, vectors = V,
       diagnostics = data.frame(dimension = n, eigenvalue_min = min(ee$values),
         eigenvalue_max = max(ee$values), relative_residual = residual_error,
         C_orthogonality_error = orthogonality_error, estimated_dense_mib = estimated_mib,
         elapsed_seconds = proc.time()[["elapsed"]] - started))
}

isolated_exact_be_reference <- function(C, spectrum, initial, spatial_load, time_seq,
                                        alpha, diffusion = 1, A = 1, omega = pi,
                                        max_history_mib = 256) {
  n <- nrow(C)
  stopifnot(length(initial) == n, length(spatial_load) == n,
            all(is.finite(initial)), all(is.finite(spatial_load)),
            length(time_seq) >= 2L, all(is.finite(time_seq)), all(diff(time_seq) > 0),
            length(alpha) == 1L, is.finite(alpha), alpha > 0, alpha <= 2,
            length(diffusion) == 1L, is.finite(diffusion), diffusion > 0,
            length(A) == 1L, is.finite(A), length(omega) == 1L, is.finite(omega),
            is.finite(max_history_mib), max_history_mib > 0,
            identical(dim(spectrum$vectors), c(n, n)), length(spectrum$values) == n)
  dt <- diff(time_seq)[1L]
  stopifnot(max(abs(diff(time_seq) - dt)) <=
              128 * .Machine$double.eps * max(1, abs(tail(time_seq, 1L))))
  history_mib <- 2 * 8 * as.double(n) * length(time_seq) / 1024^2
  if (history_mib > max_history_mib)
    stop("Exact reference history exceeds max_history_mib; use a moderate time grid.")
  V <- spectrum$vectors
  x <- as.numeric(crossprod(V, C %*% initial))
  f <- as.numeric(crossprod(V, spatial_load))
  rates <- diffusion * spectrum$values^(alpha / 2)
  mu <- 1 / (1 + dt * rates)
  g <- average_sine_forcing(time_seq, A = A, omega = omega)
  states <- coefficients <- matrix(0, n, length(time_seq))
  states[, 1L] <- initial
  coefficients[, 1L] <- x
  y <- x
  for (k in seq_len(length(time_seq) - 1L)) {
    y <- mu * (y + dt * f * g[k + 1L])
    coefficients[, k + 1L] <- y
    states[, k + 1L] <- as.numeric(V %*% y)
  }
  list(states = states, coefficients = coefficients, initial_coefficients = x,
       forcing_coefficients = f, exact_rates = rates, g_average = g,
       dt = dt, time_seq = time_seq, alpha = alpha, diffusion = diffusion,
       history_mib = history_mib)
}


## -----------------------------------------------------------------------------
isolated_tadpole_data <- function(h, kappa = 4, initial_amplitude = 10,
                                 forcing_amplitude = 10) {
  stopifnot(length(h) == 1L, is.finite(h), h > 0,
            length(kappa) == 1L, is.finite(kappa), kappa > 0,
            is.finite(initial_amplitude), is.finite(forcing_amplitude))
  graph <- gets.graph.tadpole(h)
  graph$compute_fem()
  C <- graph$mesh$C
  L <- kappa^2 * C + graph$mesh$G
  PtE <- graph$mesh$PtE
  lengths <- graph$get_edge_lengths()
  stopifnot(length(lengths) == 2L, max(abs(lengths - c(1, 2))) < 1e-12,
            nrow(C) == 2L + nrow(PtE))
  load <- numeric(nrow(C))
  element_lengths <- numeric()
  for (edge in 1:2) {
    idx <- which(PtE[, 1L] == edge)
    idx <- idx[order(PtE[idx, 2L])]
    coordinates <- c(0, lengths[edge] * PtE[idx, 2L], lengths[edge])
    nodes <- c(if (edge == 1L) 1L else 2L, 2L + idx, 2L)
    element_lengths <- c(element_lengths, diff(coordinates))
    for (i in seq_len(length(coordinates) - 1L)) {
      pair <- trig_hat_integrals(coordinates[i], coordinates[i + 1L], 2 * pi,
                                 A = forcing_amplitude * sqrt(2 / 3))
      load[nodes[i]] <- load[nodes[i]] + pair[1L]
      load[nodes[i + 1L]] <- load[nodes[i + 1L]] + pair[2L]
    }
  }
  list(C = C, L = L, initial = project_u0_tadpole(graph, initial_amplitude),
       spatial_load = load, kappa = kappa, h_requested = h,
       h_min = min(element_lengths), h_max = max(element_lengths))
}


## -----------------------------------------------------------------------------
isolated_default_coefficient_file <- function() {
  here::here("matlab", "computed_fractional_roots.csv")
}

isolated_coefficient_source <- function(parameters) {
  value <- parameters$coefficient_source
  if (identical(value, "computed")) return("computed")
  if (value %in% c("package", "original rSPDE package tables with linear interpolation"))
    return("package")
  stop("Unknown isolated-study coefficient source: ", value)
}

isolated_validate_roots <- function(roots, m) {
  values <- c(roots$pr_roots, roots$pl_roots)
  if (length(roots$factor) != 1L || !is.numeric(roots$factor) ||
      !is.finite(roots$factor) || roots$factor <= 0 ||
      length(roots$pr_roots) != m || length(roots$pl_roots) != m + 1L ||
      !is.numeric(values) || is.complex(values) || any(!is.finite(values)) || any(values >= 0))
    stop("Fractional coefficients need one positive factor and finite real negative roots of the expected counts.")
  numerator <- sort(roots$pr_roots)
  denominator <- sort(roots$pl_roots)
  if (any(denominator[-length(denominator)] >= numerator) ||
      any(numerator >= denominator[-1L]))
    stop("Fractional numerator and denominator roots must strictly interlace.")
  roots
}

isolated_read_computed_coefficients <- function(file = isolated_default_coefficient_file()) {
  if (length(file) != 1L || is.na(file) || !file.exists(file))
    stop("Missing computed coefficient file: ", file)
  tab <- read.csv(file, stringsAsFactors = FALSE)
  if (!identical(names(tab), c("m", "beta", "kind", "index", "value")) ||
      !nrow(tab) ||
      !all(vapply(tab[c("m", "beta", "index", "value")], is.numeric, logical(1))) ||
      any(!is.finite(as.matrix(tab[c("m", "beta", "index", "value")])) ) ||
      any(!tab$m %in% 1:8) || any(tab$beta < .25 | tab$beta >= 1) ||
      any(!tab$kind %in% c("factor", "rc", "rb")) ||
      any(tab$index != floor(tab$index)) || anyDuplicated(tab[c("m", "beta", "kind", "index")]))
    stop("Malformed computed coefficient artifact.")
  tab
}

isolated_fractional_roots <- function(m, beta, coefficient_source = "package",
                                      coefficient_file = isolated_default_coefficient_file(),
                                      coefficient_table = NULL) {
  coefficient_source <- match.arg(coefficient_source, c("package", "computed"))
  orders <- if (coefficient_source == "computed") 1:8 else 1:4
  if (length(m) != 1L || !is.finite(m) || !(m %in% orders))
    stop("Coefficient source '", coefficient_source, "' supports m = ", paste(orders, collapse = ", "))
  if (length(beta) != 1L || !is.finite(beta) || beta < .25 || beta >= 1)
    stop("Fractional coefficient beta must satisfy 0.25 <= beta < 1.")
  if (coefficient_source == "package") return(isolated_validate_roots(my.get.roots(m, beta), m))
  tab <- if (is.null(coefficient_table)) isolated_read_computed_coefficients(coefficient_file) else coefficient_table
  rows <- tab[tab$m == m & abs(tab$beta - beta) <=
    32 * .Machine$double.eps * max(1, abs(beta)), ]
  if (!nrow(rows)) stop("No computed coefficients for m = ", m, ", beta = ", beta,
    ". Available beta: ", paste(sort(unique(tab$beta)), collapse = ", "),
    ". Computed coefficients are never interpolated or replaced by package coefficients.")
  values <- function(kind, indices) {
    part <- rows[rows$kind == kind, ]
    if (nrow(part) != length(indices) || !setequal(part$index, indices))
      stop("Incomplete computed coefficient block for m = ", m, ", beta = ", beta)
    part$value[match(indices, part$index)]
  }
  isolated_validate_roots(list(factor = values("factor", 0),
    pr_roots = values("rc", seq_len(m)), pl_roots = values("rb", seq_len(m + 1L))), m)
}


## -----------------------------------------------------------------------------
isolated_rational_rates <- function(eigenvalues, alpha, m, kappa, diffusion = 1,
                                     coefficient_source = "package",
                                     coefficient_file = isolated_default_coefficient_file(), roots = NULL) {
  coefficient_source <- match.arg(coefficient_source, c("package", "computed"))
  orders <- if (coefficient_source == "computed") 1:8 else 1:4
  stopifnot(all(is.finite(eigenvalues)), all(eigenvalues > 0),
            length(alpha) == 1L, alpha >= .5, alpha <= 2,
            length(m) == 1L, m %in% orders, is.finite(diffusion), diffusion > 0,
            is.finite(kappa), kappa > 0)
  if (alpha == 2) return(diffusion * eigenvalues)
  if (is.null(roots)) roots <- isolated_fractional_roots(m, alpha / 2, coefficient_source, coefficient_file)
  stopifnot(is.finite(roots$factor), roots$factor > 0,
            all(is.finite(c(roots$pr_roots, roots$pl_roots))),
            all(c(roots$pr_roots, roots$pl_roots) < 0))
  z <- eigenvalues / kappa^2
  log_inverse <- log(roots$factor) + vapply(z, function(x)
    sum(log1p(-roots$pr_roots * x)) - sum(log1p(-roots$pl_roots * x)), numeric(1))
  rates <- diffusion * kappa^alpha * exp(-log_inverse)
  stopifnot(all(is.finite(rates)), all(rates > 0))
  rates
}

isolated_beta_bracket <- function(m, beta) {
  table <- get(paste0("m", m, "table"), envir = asNamespace("rSPDE"))
  tolerance <- 32 * .Machine$double.eps * max(1, abs(beta))
  exact <- any(abs(table$beta - beta) <= tolerance)
  list(exact = exact, lower = max(table$beta[table$beta <= beta + tolerance]),
       upper = min(table$beta[table$beta >= beta - tolerance]))
}

isolated_rational_error <- function(data, spectrum, reference, m,
                                    excitation_tolerance = 1e-8, coefficient_source = "package",
                                    coefficient_file = isolated_default_coefficient_file(),
                                    coefficient_table = NULL) {
  coefficient_source <- match.arg(coefficient_source, c("package", "computed"))
  started <- proc.time()[["elapsed"]]
  stopifnot(nrow(reference$states) == nrow(data$C),
            ncol(reference$states) == length(reference$time_seq))
  roots <- if (reference$alpha == 2) NULL else isolated_fractional_roots(
    m, reference$alpha / 2, coefficient_source, coefficient_file, coefficient_table)
  op <- my.fractional.operators.frac(
    L = data$L, beta = reference$alpha / 2, C = data$C,
    scale.factor = data$kappa^2, m = m, time_step = reference$dt,
    diffusion = reference$diffusion,
    fractional_roots = if (coefficient_source == "computed") roots else NULL)
  rates <- isolated_rational_rates(spectrum$values, reference$alpha, m,
                                   data$kappa, reference$diffusion, coefficient_source, coefficient_file, roots)
  mu <- 1 / (1 + reference$dt * rates)
  y <- reference$initial_coefficients
  u <- as.numeric(data$initial)
  total <- spectral_total <- compensation <- spectral_compensation <- 0
  for (k in seq_len(length(reference$time_seq) - 1L)) {
    forcing <- reference$g_average[k + 1L]
    rhs <- data$C %*% u + reference$dt * data$spatial_load * forcing
    u <- as.numeric(my.solver.frac(op, rhs))
    delta <- reference$states[, k + 1L] - u
    value <- reference$dt * sum(delta * (data$C %*% delta))
    adjusted <- value - compensation
    updated <- total + adjusted
    compensation <- (updated - total) - adjusted
    total <- updated
    y <- mu * (y + reference$dt * reference$forcing_coefficients * forcing)
    dy <- reference$coefficients[, k + 1L] - y
    value <- reference$dt * sum(dy^2)
    adjusted <- value - spectral_compensation
    updated <- spectral_total + adjusted
    spectral_compensation <- (updated - spectral_total) - adjusted
    spectral_total <- updated
  }
  sparse_error <- sqrt(max(0, total))
  product_error <- sqrt(max(0, spectral_total))
  gap <- u - as.numeric(spectrum$vectors %*% y)
  final_gap <- sqrt(max(0, sum(gap * (data$C %*% gap))))
  excited <- abs(reference$initial_coefficients) > excitation_tolerance *
    max(abs(reference$initial_coefficients)) |
    abs(reference$forcing_coefficients) > excitation_tolerance *
    max(abs(reference$forcing_coefficients))
  relative_rate_error <- rates / reference$exact_rates - 1
  bracket <- if (reference$alpha == 2) list(exact = TRUE, lower = 1, upper = 1) else
    if (coefficient_source == "computed") list(exact = TRUE, lower = reference$alpha / 2,
                                               upper = reference$alpha / 2) else
      isolated_beta_bracket(m, reference$alpha / 2)
  data.frame(alpha = reference$alpha, diffusion = reference$diffusion, m = m,
    error = sparse_error, product_error = product_error,
    absolute_error_disagreement = abs(sparse_error - product_error),
    relative_error_disagreement = abs(sparse_error - product_error) /
      max(product_error, .Machine$double.xmin), final_state_gap_C = final_gap,
    max_relative_rate_error_full = max(abs(relative_rate_error)),
    max_relative_rate_error_excited = if (any(excited)) max(abs(relative_rate_error[excited])) else 0,
    excited_numerical_modes = sum(excited), coefficient_source = coefficient_source,
    coefficient_beta_exact = bracket$exact,
    beta_exact_in_package_table = if (coefficient_source == "package") bracket$exact else NA,
    beta_table_lower = bracket$lower, beta_table_upper = bracket$upper,
    elapsed_seconds = proc.time()[["elapsed"]] - started)
}


## -----------------------------------------------------------------------------
isolated_m_parameters <- function(coefficient_source = "package",
                                   coefficient_file = isolated_default_coefficient_file()) {
  coefficient_source <- match.arg(coefficient_source, c("package", "computed"))
  parameters <- list(experiment = "isolated_discrete_fractional_m", T_final = 1e-3,
       diffusion = 1, kappa = 4, initial_amplitude = 10, forcing_amplitude = 10,
       A = 1, omega = pi, alpha = c(.5, .9, 1.3), m = 1:4,
       scenarios = data.frame(scenario = c("base", "half_h", "half_tau", "half_h_half_tau"),
         h = c(.01, .005, .01, .005), target_step = c(1e-5, 1e-5, 5e-6, 5e-6)),
       coefficient_source = "original rSPDE package tables with linear interpolation",
       load_assembly = "analytic P1 loads for 10e5 and 10e7",
       thresholds = list(error_agreement_relative = 1e-3, error_agreement_absolute = 1e-11,
         zero_control_absolute = 1e-9, error_refinement_relative = .02,
         slope_refinement_absolute = .05, resolution_ratio = 100),
       package_versions = vapply(c("Matrix", "MetricGraph", "rSPDE"),
         function(p) as.character(utils::packageVersion(p)), character(1)))
  if (coefficient_source == "computed") {
    isolated_read_computed_coefficients(coefficient_file)
    parameters$m <- 1:8
    parameters$coefficient_source <- "computed"
    parameters$coefficient_file <- normalizePath(coefficient_file, mustWork = TRUE)
    parameters$coefficient_sha256 <- digest::digest(file = coefficient_file, algo = "sha256")
    parameters$coefficient_generator <- if (identical(parameters$coefficient_sha256,
      "d299b6583a35d4613fa25d14318e7d0065842efb2b6efb9d4be2857b772b4596"))
      list(repository = "NAFDEMG", branch = "diffusion-parameter",
        script = "tests/audit_rational_generation.m",
        git_blob = "5d0f78d134c12ea7a3534f2142d05c2a9a5ec504",
        production_orders = 1:8) else list(description = "Custom artifact; generator not verified")
  }
  parameters
}

isolated_fit_slopes <- function(errors, parameters) {
  rows <- list()
  for (scenario in parameters$scenarios$scenario) for (alpha in parameters$alpha) {
    values <- errors[errors$scenario == scenario & errors$alpha == alpha, ]
    values <- values[order(values$m), ]
    valid <- nrow(values) >= 2L && all(is.finite(values$error)) && all(values$error > 0)
    fit <- if (valid) lm(log(error) ~ sqrt(m), data = values) else NULL
    resolved <- if ("resolved_above_diagnostic_noise" %in% names(values))
      values[values$resolved_above_diagnostic_noise & values$error > 0, ] else values
    resolved_fit <- if (nrow(resolved) >= 2L) lm(log(error) ~ sqrt(m), data = resolved) else NULL
    rows[[length(rows) + 1L]] <- data.frame(scenario = scenario, alpha = alpha,
      observed_slope = if (valid) unname(coef(fit)[2L]) else NA_real_,
      bound_benchmark = -2 * pi * sqrt(1 - alpha / 2),
      resolved_order_count = nrow(resolved),
      resolved_only_slope = if (!is.null(resolved_fit)) unname(coef(resolved_fit)[2L]) else NA_real_,
      r_squared = if (valid) summary(fit)$r.squared else NA_real_,
      max_log_fit_residual = if (valid) max(abs(residuals(fit))) else NA_real_,
      strictly_decreasing = all(diff(values$error) < 0),
      increasing_pairs = sum(diff(values$error) > 0),
      last_to_previous_error = if (nrow(values) >= 2L)
        tail(values$error, 1L) / tail(values$error, 2L)[1L] else NA_real_)
  }
  do.call(rbind, rows)
}

run_isolated_m_study <- function(parameters = isolated_m_parameters(), workers = 4L,
                                 max_dimension = 1200L, max_dense_mib = 1024,
                                 max_history_mib = 256, verbose = TRUE) {
  source <- isolated_coefficient_source(parameters)
  orders <- if (source == "computed") 1:8 else 1:4
  coefficient_file <- if (source == "computed") parameters$coefficient_file else isolated_default_coefficient_file()
  coefficient_table <- NULL
  if (source == "computed") {
    coefficient_table <- isolated_read_computed_coefficients(coefficient_file)
    if (is.null(parameters$coefficient_sha256) || !identical(parameters$coefficient_sha256,
        digest::digest(file = coefficient_file, algo = "sha256")))
      stop("Computed artifact changed after parameters were prepared; create new parameters and a fresh run.")
  }
  stopifnot(length(parameters$T_final) == 1L, is.finite(parameters$T_final), parameters$T_final > 0,
            length(parameters$diffusion) == 1L, is.finite(parameters$diffusion), parameters$diffusion > 0,
            length(parameters$alpha) >= 1L, all(is.finite(parameters$alpha)),
            all(parameters$alpha >= .5 & parameters$alpha < 2), !anyDuplicated(parameters$alpha),
            length(parameters$m) >= 2L, all(parameters$m %in% orders), !anyDuplicated(parameters$m),
            nrow(parameters$scenarios) >= 1L, !anyDuplicated(parameters$scenarios$scenario),
            parameters$scenarios$scenario[1L] == "base",
            all(is.finite(parameters$scenarios$h)), all(parameters$scenarios$h > 0),
            all(is.finite(parameters$scenarios$target_step)), all(parameters$scenarios$target_step > 0))
  for (alpha in parameters$alpha) for (m in parameters$m)
    isolated_fractional_roots(m, alpha / 2, source, coefficient_file, coefficient_table)
  started <- proc.time()[["elapsed"]]
  meshes <- unique(parameters$scenarios$h)
  spatial <- lapply(meshes, function(h) {
    data <- isolated_tadpole_data(h, parameters$kappa, parameters$initial_amplitude,
                                  parameters$forcing_amplitude)
    spectrum <- isolated_discrete_spectrum(data$C, data$L, max_dimension, max_dense_mib)
    if (verbose) message(sprintf("Complete FEM eigensystem: h target %.5g, dimension %d, %.2f s",
      h, nrow(data$C), spectrum$diagnostics$elapsed_seconds))
    list(data = data, spectrum = spectrum)
  })
  spectral_diagnostics <- do.call(rbind, lapply(seq_along(meshes), function(i)
    data.frame(h_requested = meshes[i], h_min = spatial[[i]]$data$h_min,
      h_max = spatial[[i]]$data$h_max, spatial[[i]]$spectrum$diagnostics)))
  errors <- controls <- reference_setup <- list()
  for (i in seq_len(nrow(parameters$scenarios))) {
    scenario <- parameters$scenarios[i, ]
    cache <- spatial[[match(scenario$h, meshes)]]
    data <- cache$data
    spectrum <- cache$spectrum
    grid <- make_time_grid(parameters$T_final, scenario$target_step)
    references <- lapply(parameters$alpha, function(alpha) {
      begin <- proc.time()[["elapsed"]]
      result <- isolated_exact_be_reference(data$C, spectrum, data$initial,
        data$spatial_load, grid$time_seq, alpha, parameters$diffusion,
        parameters$A, parameters$omega, max_history_mib)
      reference_setup[[length(reference_setup) + 1L]] <<- data.frame(
        scenario = scenario$scenario, alpha = alpha, history_mib = result$history_mib,
        elapsed_seconds = proc.time()[["elapsed"]] - begin)
      result
    })
    jobs <- expand.grid(alpha_index = seq_along(parameters$alpha), m_index = seq_along(parameters$m))
    if (verbose) message(sprintf("Isolated cases: %s, %d steps, %d jobs, up to %d workers",
      scenario$scenario, grid$N, nrow(jobs), min(workers, nrow(jobs))))
    rows <- run_reference_cases(seq_len(nrow(jobs)), function(job) {
      isolated_rational_error(data, spectrum, references[[jobs$alpha_index[job]]],
                              parameters$m[jobs$m_index[job]], coefficient_source = source,
                              coefficient_file = coefficient_file, coefficient_table = coefficient_table)
    }, workers = workers)
    values <- do.call(rbind, rows)
    errors[[i]] <- data.frame(scenario = scenario$scenario, h_requested = scenario$h,
      h_min = data$h_min, h_max = data$h_max, time_step = grid$time_step,
      n_steps = grid$N, dimension = nrow(data$C), values)
    zero_reference <- isolated_exact_be_reference(data$C, spectrum, data$initial,
      data$spatial_load, grid$time_seq, 2, parameters$diffusion, parameters$A,
      parameters$omega, max_history_mib)
    controls[[i]] <- data.frame(scenario = scenario$scenario,
      isolated_rational_error(data, spectrum, zero_reference, 1L,
        coefficient_source = source, coefficient_file = coefficient_file,
        coefficient_table = coefficient_table))
    if (verbose) message(sprintf("Completed %s; alpha=2 zero-rational-error control: %.3g",
      scenario$scenario, controls[[i]]$error))
  }
  errors <- do.call(rbind, errors)
  controls <- do.call(rbind, controls)
  noise <- controls$error[match(errors$scenario, controls$scenario)]
  errors$resolution_ratio <- errors$product_error /
    pmax(errors$absolute_error_disagreement, noise, .Machine$double.xmin)
  threshold <- parameters$thresholds
  errors$implementation_agreement_passed <-
    errors$absolute_error_disagreement <= threshold$error_agreement_absolute |
    errors$relative_error_disagreement <= threshold$error_agreement_relative
  errors$resolved_above_diagnostic_noise <- errors$resolution_ratio >= threshold$resolution_ratio
  slopes <- isolated_fit_slopes(errors, parameters)
  base <- errors[errors$scenario == "base", c("alpha", "m", "error")]
  refinements <- lapply(setdiff(parameters$scenarios$scenario, "base"), function(scenario) {
    refined <- errors[errors$scenario == scenario, c("alpha", "m", "error")]
    joined <- merge(base, refined, by = c("alpha", "m"), suffixes = c("_base", "_refined"))
    joined$scenario <- scenario
    joined$relative_change <- abs(joined$error_refined / joined$error_base - 1)
    joined$passed <- joined$relative_change <= threshold$error_refinement_relative
    joined
  })
  refinement_errors <- if (length(refinements)) do.call(rbind, refinements) else data.frame()
  base_slopes <- slopes[slopes$scenario == "base", c("alpha", "observed_slope")]
  slope_checks <- lapply(setdiff(parameters$scenarios$scenario, "base"), function(scenario) {
    refined <- slopes[slopes$scenario == scenario, c("alpha", "observed_slope")]
    joined <- merge(base_slopes, refined, by = "alpha", suffixes = c("_base", "_refined"))
    joined$scenario <- scenario
    joined$absolute_change <- abs(joined$observed_slope_refined - joined$observed_slope_base)
    joined$passed <- joined$absolute_change <= threshold$slope_refinement_absolute
    joined
  })
  refinement_slopes <- if (length(slope_checks)) do.call(rbind, slope_checks) else data.frame()
  status <- list(numerical_checks_passed = all(errors$implementation_agreement_passed) &&
      all(controls$error <= threshold$zero_control_absolute),
    reference_error_resolution_passed = all(errors$resolved_above_diagnostic_noise),
    mesh_time_sensitivity_passed = all(refinement_errors$passed) && all(refinement_slopes$passed),
    every_curve_strictly_decreasing = all(slopes$strictly_decreasing))
  list(parameters = parameters, errors = errors, slopes = slopes,
       refinement_errors = refinement_errors, refinement_slopes = refinement_slopes,
       zero_controls = controls, eigensystem_diagnostics = spectral_diagnostics,
       reference_setup = do.call(rbind, reference_setup), status = status,
       elapsed_seconds = proc.time()[["elapsed"]] - started)
}

isolated_m_diagnosis <- function(results) {
  base <- results$errors[results$errors$scenario == "base", ]
  base_slopes <- results$slopes[results$slopes$scenario == "base", ]
  words <- c(
    "The measured quantity is the difference between exact-discrete and rational FEM/backward Euler solutions on the same grid.",
    "Every primary rational trajectory uses the shifted sparse solver. The diagnostic root-product trajectory uses the complete numerical FEM eigenbasis.",
    sprintf("Numerical implementation checks passed: %s. All measured errors resolved above the diagnostic noise indicator: %s.",
      results$status$numerical_checks_passed, results$status$reference_error_resolution_passed),
    sprintf("Maximum normalized generalized-eigen residual: %.4g; maximum C-orthogonality error: %.4g.",
      max(results$eigensystem_diagnostics$relative_residual),
      max(results$eigensystem_diagnostics$C_orthogonality_error)),
    sprintf("Largest sparse/product error disagreement: %.4g absolute, %.4g relative. Largest alpha=2 control error: %.4g.",
      max(results$errors$absolute_error_disagreement), max(results$errors$relative_error_disagreement),
      max(results$zero_controls$error)),
    sprintf("Mesh/time sensitivity checks passed: %s. Largest relative error change: %.4g; largest absolute fitted-slope change: %.4g.",
      results$status$mesh_time_sensitivity_passed,
      if (nrow(results$refinement_errors)) max(results$refinement_errors$relative_change) else 0,
      if (nrow(results$refinement_slopes)) max(results$refinement_slopes$absolute_change) else 0))
  for (alpha in results$parameters$alpha) {
    rows <- base[base$alpha == alpha, ]
    rows <- rows[order(rows$m), ]
    fit <- base_slopes[base_slopes$alpha == alpha, ]
    words <- c(words, sprintf(
      "alpha=%g: observed slope %.6f, bound benchmark %.6f, R-squared %.5f; strictly decreasing: %s; final-order/previous-order error ratio %.6g.",
      alpha, fit$observed_slope, fit$bound_benchmark, fit$r_squared,
      fit$strictly_decreasing, fit$last_to_previous_error))
    if (any(diff(rows$error) > 0)) {
      pairs <- which(diff(rows$error) > 0)
      for (j in pairs) words <- c(words, sprintf(
        "For alpha=%g the error increases from m=%g to m=%g. The corresponding maximum excited-mode relative rate errors are %.6g and %.6g; this checks whether the scalar approximation itself shows the same deterioration.",
        alpha, rows$m[j], rows$m[j + 1L], rows$max_relative_rate_error_excited[j],
        rows$max_relative_rate_error_excited[j + 1L]))
    }
  }
  if (isolated_coefficient_source(results$parameters) == "package" &&
      any(!base$beta_exact_in_package_table)) words <- c(words,
    "Some beta values are interpolated between package rows. The beta brackets are recorded. Interpolation may contribute to a finite-order plateau; these runs do not identify it as the sole cause.")
  words <- c(words, paste0("Coefficient source: ", isolated_coefficient_source(results$parameters), "."))
  if (!results$status$reference_error_resolution_passed) words <- c(words,
    "Some orders are below the diagnostic numerical-resolution threshold. Their errors remain plotted, but the all-order fit includes unresolved points. Consult resolved_order_count and resolved_only_slope; the restricted fit is descriptive, not an asymptotic claim. Where unresolved points are involved, large relative refinement changes are inconclusive for mesh/time sensitivity.")
  c(words,
    "A fit across the tested orders is a finite-order measurement. The upper-bound benchmark does not require equality of slopes, linearity on the square-root scale, or strict monotonicity of these computed approximants.",
    "The numerical controls are sensitivity and consistency checks, not certified roundoff bounds or a uniform fractional-operator approximation theorem.")
}


## -----------------------------------------------------------------------------
save_isolated_m_results <- function(results, run_dir) {
  destination <- file.path(run_dir, "isolated_m_results.rds")
  if (file.exists(destination)) stop("Refusing to overwrite isolated results: ", destination)
  stopifnot(identical(readRDS(file.path(run_dir, "parameters.rds")), results$parameters))
  saveRDS(results, destination)
  for (name in c("errors", "slopes", "refinement_errors", "refinement_slopes",
                 "zero_controls", "eigensystem_diagnostics", "reference_setup"))
    write.csv(results[[name]], file.path(run_dir, paste0(name, ".csv")), row.names = FALSE)
  paragraphs <- isolated_m_diagnosis(results)
  writeLines(c("# Isolated m-study diagnosis", "", as.vector(rbind(paragraphs, ""))),
             file.path(run_dir, "diagnosis.md"))
  invisible(destination)
}

isolated_m_error_plot <- function(results) {
  errors <- results$errors
  errors$alpha_label <- factor(errors$alpha)
  errors$scenario <- factor(errors$scenario, levels = results$parameters$scenarios$scenario)
  guides <- do.call(rbind, lapply(split(errors, list(errors$scenario, errors$alpha), drop = TRUE),
    function(rows) {
      rows <- rows[order(rows$m), ]
      rows$error <- rows$error[1L] * exp(-2 * pi * sqrt(1 - rows$alpha[1L] / 2) *
        (sqrt(rows$m) - sqrt(rows$m[1L])))
      rows
    }))
  ggplot2::ggplot(errors, ggplot2::aes(sqrt(m), error, color = alpha_label, group = alpha_label)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(ggplot2::aes(shape = resolved_above_diagnostic_noise), size = 2) +
    ggplot2::scale_shape_manual(name = "Resolution indicator", values = c(`FALSE` = 1, `TRUE` = 16),
      labels = c(`FALSE` = "below indicator", `TRUE` = "resolved")) +
    ggplot2::geom_line(data = guides, linetype = "dashed", alpha = .7) +
    ggplot2::facet_wrap(~scenario, ncol = 2) +
    ggplot2::scale_x_continuous(breaks = sqrt(results$parameters$m), labels = results$parameters$m) +
    ggplot2::scale_y_continuous(trans = scales::log_trans(base = exp(1)),
      breaks = 10^seq(floor(log10(min(errors$error))), ceiling(log10(max(errors$error)))),
      labels = scales::label_scientific(digits = 2)) +
    ggplot2::labs(x = "Rational order m (square-root positions)",
      y = "Isolated space-time rational error (natural-log scale)", color = "alpha",
      subtitle = "Solid: sparse-solver error. Dashed: bound exponent, anchored at first order.\nOpen circles: below the numerical resolution indicator.") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.position = "bottom")
}

render_isolated_m_notebook <- function(parameters = NULL, run_id = NULL,
                                       execute = TRUE, workers = 4L,
                                       input = here::here("conv_in_m_isolated.Rmd"),
                                       max_dimension = 1200L, max_dense_mib = 1024,
                                       max_history_mib = 256, encoding = "UTF-8",
                                       coefficient_source = NULL, coefficient_file = NULL) {
  defaults <- rmarkdown::yaml_front_matter(input)$params
  if (is.null(parameters)) {
    if (is.null(coefficient_source)) coefficient_source <- defaults$coefficient_source
    if (is.null(coefficient_source)) coefficient_source <- "package"
    if (is.null(coefficient_file)) coefficient_file <- defaults$coefficient_file
    if (is.null(coefficient_file)) coefficient_file <- isolated_default_coefficient_file()
    parameters <- isolated_m_parameters(coefficient_source, coefficient_file)
  }
  stopifnot(length(execute) == 1L, !is.na(execute), is.logical(execute))
  run_dir <- reference_experiment_run("m", parameters, run_id)
  file <- file.path(run_dir, "isolated_m_results.rds")
  if (execute && file.exists(file)) stop("Run already has results; set execute = FALSE or use a new run ID.")
  if (!execute && !file.exists(file)) stop("No isolated results to reload in this run.")
  intermediates <- file.path(run_dir, "render")
  dir.create(intermediates, showWarnings = FALSE)
  rmarkdown::render(input, output_dir = run_dir, intermediates_dir = intermediates,
    knit_root_dir = dirname(normalizePath(input)), encoding = encoding,
    params = list(settings = parameters, run_id = basename(run_dir), execute = execute,
      coefficient_source = isolated_coefficient_source(parameters),
      coefficient_file = if (is.null(parameters$coefficient_file)) isolated_default_coefficient_file() else parameters$coefficient_file,
      workers = workers, max_dimension = max_dimension, max_dense_mib = max_dense_mib,
      max_history_mib = max_history_mib), envir = new.env(parent = globalenv()), quiet = TRUE)
}
