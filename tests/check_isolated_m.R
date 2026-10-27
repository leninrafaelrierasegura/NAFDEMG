# Run in the project's renv environment: Rscript tests/check_isolated_m.R
suppressPackageStartupMessages(library(rSPDE))
suppressPackageStartupMessages(library(MetricGraph))
source("tests/reference_test_helpers.R")
functions <- load_reference_functions()
for (name in c("graphs_and_meshes", "isolated_rational_convergence")) {
  file <- tempfile(fileext = ".R")
  invisible(knitr::purl(file.path("fun", paste0(name, ".Rmd")), output = file,
                       documentation = 0L, quiet = TRUE))
  for (expression in parse(file)) if (is.call(expression) &&
      identical(expression[[1L]], as.name("<-")) && is.call(expression[[3L]]) &&
      identical(expression[[3L]][[1L]], as.name("function"))) eval(expression, functions)
  unlink(file)
}
checks <- 0L
check <- function(condition, label) {
  if (!isTRUE(condition)) stop(label)
  checks <<- checks + 1L
}
fails <- function(expression) inherits(try(force(expression), silent = TRUE), "try-error")

# A nonuniform, nonidentity-mass problem ensures no analytical discrete modes
# are being assumed in the generalized-eigen reference or norm comparison.
fem <- interval_fem(c(0, .13, .49, .7, 1))
C <- fem$C
L <- 4 * C + fem$G
spectrum <- functions$isolated_discrete_spectrum(C, L)
check(spectrum$diagnostics$relative_residual < 1e-12, "Generalized eigen residual")
check(spectrum$diagnostics$C_orthogonality_error < 1e-12, "Mass orthogonality")
check(fails(functions$isolated_discrete_spectrum(C, L, max_dimension = 4L)), "Preallocation dimension guard")
check(fails(functions$isolated_discrete_spectrum(C, L, max_dense_mib = 1e-8)), "Preallocation memory guard")
time <- seq(0, .2, length.out = 9)
initial <- c(.2, .7, -.1, .3, 1)
load <- c(.3, -.2, .1, .7, .2)
data <- list(C = C, L = L, initial = initial, spatial_load = load, kappa = 2)
for (alpha in c(.5, .9, 1.3, 2)) for (diffusion in c(.2, 1, 3)) {
  reference <- functions$isolated_exact_be_reference(C, spectrum, initial, load,
    time, alpha, diffusion, A = .8, omega = 1.7)
  # Independent nonsymmetric eigendecomposition of C^-1 L for a small matrix.
  direct <- eigen(solve(as.matrix(C), as.matrix(L)))
  check(max(abs(Im(direct$values))) < 1e-12, "Independent real spectrum")
  fractional <- direct$vectors %*% diag(direct$values^(alpha / 2)) %*% solve(direct$vectors)
  u <- initial
  dt <- diff(time)[1L]
  for (k in seq_len(length(time) - 1L)) {
    u <- as.numeric(solve(diag(length(initial)) + dt * diffusion * fractional,
      u + dt * solve(C, load) * reference$g_average[k + 1L]))
    check(max(abs(u - reference$states[, k + 1L])) < 2e-12,
          "Reference agrees with independent exact fractional BE")
  }
  for (m in c(1L, 4L)) {
    result <- functions$isolated_rational_error(data, spectrum, reference, m)
    op <- functions$my.fractional.operators.frac(L, alpha / 2, C, 4, m, dt, diffusion)
    u <- initial
    total <- 0
    for (k in seq_len(length(time) - 1L)) {
      u <- as.numeric(functions$my.solver.frac(op,
        C %*% u + dt * load * reference$g_average[k + 1L]))
      delta <- reference$states[, k + 1L] - u
      total <- total + dt * sum(delta * (C %*% delta))
    }
    check(abs(result$error - sqrt(total)) < 1e-13, "Primary norm equals literal time integral")
    check(abs(result$error - result$product_error) < 2e-11, "Sparse/product consistency")
    if (alpha == 2) check(result$error < 2e-12 && result$product_error == 0,
                          "Classical BE zero-rational-error control")
  }
}
check(fails(functions$isolated_exact_be_reference(C, spectrum, initial, load,
  time, 1, max_history_mib = 1e-8)), "History size guard")
check(fails(functions$isolated_exact_be_reference(C, spectrum, initial, load,
  time, 1, diffusion = 0)), "Reject nonpositive diffusion")
check(fails(functions$isolated_exact_be_reference(C, spectrum, initial, load,
  c(0, .03, .1), 1)), "Reject nonuniform reference step")

# Check physical forcing-load assembly independently by adaptive integration.
data <- functions$isolated_tadpole_data(.1)
graph <- functions$gets.graph.tadpole(.1)
graph$compute_fem()
PtE <- graph$mesh$PtE
expected <- numeric(nrow(data$C))
for (edge in 1:2) {
  idx <- which(PtE[, 1] == edge)
  idx <- idx[order(PtE[idx, 2])]
  x <- c(0, c(1, 2)[edge] * PtE[idx, 2], c(1, 2)[edge])
  nodes <- c(if (edge == 1) 1L else 2L, 2L + idx, 2L)
  for (j in seq_len(length(x) - 1L)) {
    a <- x[j]; b <- x[j + 1L]
    left <- integrate(function(s) 10 * sqrt(2/3) * cos(2*pi*s) * (b-s)/(b-a), a, b)$value
    right <- integrate(function(s) 10 * sqrt(2/3) * cos(2*pi*s) * (s-a)/(b-a), a, b)$value
    expected[nodes[j]] <- expected[nodes[j]] + left
    expected[nodes[j + 1L]] <- expected[nodes[j + 1L]] + right
  }
}
check(max(abs(data$spatial_load - expected)) < 1e-12, "Analytic load matches independent quadrature")

# The complete runner, refinements, output guards, and plot use real graph objects.
parameters <- functions$isolated_m_parameters()
parameters$alpha <- c(.5, .9)
parameters$scenarios$h <- c(.1, .05, .1, .05)
parameters$scenarios$target_step <- c(1e-4, 1e-4, 5e-5, 5e-5)
result <- functions$run_isolated_m_study(parameters, workers = 2L, verbose = FALSE)
check(nrow(result$errors) == 32L, "All isolated parameter pairs executed")
check(nrow(result$eigensystem_diagnostics) == 2L, "Spatial eigensystems reused")
check(nrow(result$reference_setup) == 8L, "References reused across rational orders")
check(nrow(result$zero_controls) == 4L, "Zero controls cover all scenarios")
check(result$status$numerical_checks_passed, "Experiment numerical checks")
check(length(functions$isolated_m_diagnosis(result)) >= 10L, "Diagnosis generated")
plot <- functions$isolated_m_error_plot(result)
check(inherits(plot, "ggplot"), "Plot constructed")
root <- tempfile("isolated-m-test-")
old_root <- Sys.getenv("NAFDEMG_PRIVATE_DIR", unset = NA_character_)
Sys.setenv(NAFDEMG_PRIVATE_DIR = root)
run <- functions$reference_experiment_run("m", parameters)
functions$save_isolated_m_results(result, run)
check(file.exists(file.path(run, "diagnosis.md")), "Private diagnosis saved")
check(fails(functions$save_isolated_m_results(result, run)), "Overwrite protection")
if (is.na(old_root)) Sys.unsetenv("NAFDEMG_PRIVATE_DIR") else Sys.setenv(NAFDEMG_PRIVATE_DIR = old_root)
unlink(root, recursive = TRUE)
cat(checks, "isolated-m spectral, sparse-evolution, load, refinement, and output checks passed.\n")
