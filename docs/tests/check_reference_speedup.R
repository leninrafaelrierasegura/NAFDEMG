# Run from the repository root: Rscript --vanilla tests/check_reference_speedup.R
source("tests/reference_test_helpers.R")
functions <- load_reference_functions()
checks <- 0L
check <- function(condition, label) {
  if (!isTRUE(condition)) stop(label)
  checks <<- checks + 1L
}
close <- function(x, y, tolerance = 2e-11) {
  max(abs(x - y)) <= tolerance * max(1, max(abs(x)), max(abs(y)))
}
fails <- function(expression) inherits(try(force(expression), silent = TRUE), "try-error")

# This fixture has a nonuniform mesh, nonidentity mass, nonnested reference
# partition, and two exact modes whose FEM projections need not be eigenvectors.
grids <- list(sort(unique(c(seq(0, .2, length.out = 14), .04, .08, .12, .16))),
              c(0, .001, .153, .2), seq(0, .2, length.out = 41))
for (beta in c(.25, .45, .65, 1)) for (diffusion in c(.1, 1, 3)) {
  fixture <- reference_fixture(functions, beta = beta, diffusion = diffusion)
  check(all(fixture$op$diffusion == diffusion), "Operator diffusion metadata")
  # Independent generalized-eigenvalue expression for the matrix one-step map.
  ce <- eigen(as.matrix(fixture$C), symmetric = TRUE)
  inverse_sqrt <- ce$vectors %*% diag(1 / sqrt(ce$values)) %*% t(ce$vectors)
  spectral <- eigen(inverse_sqrt %*% as.matrix(fixture$K) %*% inverse_sqrt,
                    symmetric = TRUE)
  lambda <- spectral$values
  if (beta == 1) rates <- diffusion * lambda else {
    roots <- functions$my.get.roots(4L, beta)
    z <- lambda / fixture$kappa^2
    inverse <- roots$factor * vapply(z, function(x)
      prod(1 - roots$pr_roots * x) / prod(1 - roots$pl_roots * x), numeric(1))
    rates <- diffusion * fixture$kappa^(2 * beta) / inverse
  }
  expected_map <- inverse_sqrt %*% spectral$vectors %*%
    diag(1 / (1 + fixture$dt * rates)) %*% t(spectral$vectors) %*% inverse_sqrt
  rhs <- as.numeric(fixture$C %*% fixture$initial +
    fixture$dt * fixture$spatial_load * fixture$g_average[2L])
  actual <- as.numeric(functions$my.solver.frac(fixture$op, rhs))
  check(close(actual, expected_map %*% rhs, 2e-9), "Diffusion and unscaled forcing one-step map")
  uncached <- fixture$op
  uncached$factorization <- uncached$factorizations <- NULL
  check(close(actual, as.numeric(functions$my.solver.frac(uncached, rhs))),
        "Cached and legacy-object solves agree")
  for (reference_times in grids) {
    expected <- literal_reference_error(functions, fixture, reference_times)
    for (blocks in c(1L, 2L, 7L, 100L)) {
      original <- fixture_error(functions, fixture, reference_times, "reference", blocks)
      compressed <- fixture_error(functions, fixture, reference_times, "compressed", blocks)
      check(close(original, expected), "Original reference path matches direct weighted sum")
      check(close(compressed, expected), "Compressed path preserves reference norm and weights")
    }
  }
}

# Automatic projection setup and one-mode callbacks, including dense interpolation.
fixture <- reference_fixture(functions, beta = .45, m = 1L, diffusion = 2)
reference_times <- grids[[1L]]
auto <- functions$solve_fractional_evolution_error(
  fixture$op, fixture$dt, fixture$time_seq, fixture$initial, fixture$spatial_load,
  fixture$g_average, as.matrix(fixture$Psi), fixture$reference_C, reference_times,
  error_method = "compressed", reference_basis = fixture$E,
  coefficients_at_time = fixture$coefficients, n_blocks = 3L)
check(close(auto, literal_reference_error(functions, fixture, reference_times)), "Automatic setup with dense Psi")
single <- functions$prepare_reference_error(fixture$Psi, fixture$reference_C, fixture$E[, 1L, drop = FALSE])
moments <- functions$reference_coefficient_moments(matrix(c(.2, .21, .23), 1), c(.1, .2, .3))
check(identical(dim(single$residual_gram), c(1L, 1L)), "One-mode projection dimensions")
check(close(moments$mean, sum(c(.2, .21, .23) * c(.1, .2, .3)) / .6), "One-mode weighted mean")

# Exact-representability and tiny residual tests exercise cancellation resistance.
identity_projection <- functions$prepare_reference_error(Diagonal(4), Diagonal(4),
                                                        cbind(c(1, 2, 0, -1), c(0, 1, 3, 1)))
c0 <- c(.7, -.3)
u0 <- as.numeric(identity_projection$projection %*% c0)
moments <- functions$reference_coefficient_moments(matrix(c0, 2L, 5L), rep(.2, 5L))
zero <- functions$reference_interval_error(u0, moments, identity_projection)
check(abs(zero) < 1e-28, "Exactly representable zero error")
tiny <- functions$reference_interval_error(u0 + c(1e-10, 0, 0, 0), moments, identity_projection)
check(abs(tiny - 1e-20) < 1e-25, "Tiny errors retained without subtraction of large Grams")
data <- rbind(2 + c(0, 1, 2, 4) * 1e-7, -.3 + c(1, 4, 2, 3) * 1e-7)
w <- c(.1, .2, .3, .4)
whole <- functions$reference_coefficient_moments(data, w)
merged <- functions$merge_reference_moments(
  functions$reference_coefficient_moments(data[, 1:2], w[1:2]),
  functions$reference_coefficient_moments(data[, 3:4], w[3:4]))
check(close(whole$mean, merged$mean), "Block mean merge")
check(max(abs(whole$scatter - merged$scatter)) < 1e-22, "Centered covariance merge")

# Stable exact coefficients use the full physical decay rate, with original forcing.
for (rate in c(.1, 2, 40)) {
  t <- .3
  convolution <- integrate(function(s) exp(-rate * (t - s)) * .9 * sin(1.7 * s), 0, t)$value
  check(close(functions$decayed_sine_convolution(t, rate, .9, 1.7), convolution),
        "Stable convolution equals physical integral")
}
check(identical(functions$decayed_sine_convolution(c(0, 1), 2, 1, 0), c(0, 0)), "Zero sine forcing")
check(fails(functions$my.fractional.operators.frac(fixture$K, .45, fixture$C, 4, 1, .04, diffusion = 0)),
      "Reject nonpositive diffusion")
check(fails(functions$solve_fractional_evolution_error(
  fixture$op, fixture$dt, fixture$time_seq, fixture$initial, fixture$spatial_load,
  fixture$g_average, fixture$Psi, fixture$reference_C, reference_times,
  error_method = "compressed", reference_basis = fixture$E,
  coefficients_at_time = function(t) matrix(NA_real_, 2L, length(t)))), "Reject nonfinite coefficients")

# Independent cases preserve ordering and propagate failures in both modes.
serial <- functions$run_reference_cases(c(3L, 1L, 2L), function(i) i * i, workers = 1L)
check(identical(serial, list(9L, 1L, 4L)), "Serial case order")
if (.Platform$OS.type == "unix") {
  parallel <- functions$run_reference_cases(c(3L, 1L, 2L), function(i) i * i, workers = 2L)
  check(identical(parallel, serial), "Parallel case order")
  check(fails(functions$run_reference_cases(1:2, function(i) stop("fixture failure"), 2L)),
        "Worker errors reach coordinator")
}

# Private run IDs are unique and explicit reloads enforce scientific metadata.
old_root <- Sys.getenv("NAFDEMG_PRIVATE_DIR", unset = NA_character_)
private <- tempfile("reference-test-")
Sys.setenv(NAFDEMG_PRIVATE_DIR = private)
parameters <- list(diffusion = 2, tau = .04)
run1 <- functions$reference_experiment_run("h", parameters)
run2 <- functions$reference_experiment_run("h", parameters)
check(run1 != run2, "Unique run IDs")
check(identical(run1, functions$reference_experiment_run("h", parameters, basename(run1))),
      "Matching run reload")
check(fails(functions$reference_experiment_run("h", list(diffusion = 3), basename(run1))),
      "Reject changed run parameters")
check(fails(functions$reference_experiment_run("h", parameters, "../escape")), "Reject path traversal run ID")
if (is.na(old_root)) Sys.unsetenv("NAFDEMG_PRIVATE_DIR") else Sys.setenv(NAFDEMG_PRIVATE_DIR = old_root)
unlink(private, recursive = TRUE)
cat(checks, "reference-error, factorization, diffusion, and execution checks passed.\n")
