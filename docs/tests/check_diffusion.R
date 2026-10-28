# Run from the repository root:
# /usr/bin/Rscript tests/check_diffusion.R
# --vanilla also works here: this script explicitly selects the project library.
project_lib <- file.path(getwd(), "renv/library/linux-ubuntu-noble/R-4.5/x86_64-pc-linux-gnu")
if (dir.exists(project_lib)) .libPaths(c(project_lib, .libPaths()))
stopifnot(getRversion() >= "4.5.0")
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(MetricGraph))

# Load function definitions only, without notebook library/plot/export setup.
load_functions <- function(path, envir) {
  for (expression in parse(path)) {
    if (is.call(expression) && identical(expression[[1]], as.name("<-")) &&
        is.call(expression[[3]]) && identical(expression[[3]][[1]], as.name("function")))
      eval(expression, envir)
  }
}
helpers <- new.env(parent = globalenv())
for (name in c("fractional_evolution", "graphs_and_meshes", "exact_solutions",
               "convergence_experiments", "convergence_settings"))
  load_functions(file.path("fun", paste0(name, ".R")), helpers)
attach(helpers, name = "diffusion-test-helpers")
near <- function(got, expected, tolerance = 1e-10) {
  stopifnot(all(is.finite(got)), max(abs(got-expected)) <= tolerance)
}
expect_error <- function(expression) {
  stopifnot(inherits(tryCatch({force(expression); NULL}, error = identity), "error"))
}

cat("Environment:", R.version.string, "; MetricGraph", as.character(packageVersion("MetricGraph")),
    "; rSPDE", as.character(packageVersion("rSPDE")), "\n")

# Independently construct generalized eigenmodes with a nonidentity,
# noncommuting mass matrix; exercise forcing and every supported order.
C <- Matrix(matrix(c(2, .2, .2, 3), 2), sparse = TRUE)
L <- Matrix(matrix(c(40, -2, -2, 160), 2), sparse = TRUE)
R <- chol(as.matrix(C))
Rinv <- solve(R)
ev <- eigen(t(Rinv) %*% as.matrix(L) %*% Rinv, symmetric = TRUE)
V <- Rinv %*% ev$vectors
dt <- 1e-4
initial <- c(1, 2)
load <- c(3, -4)
for (beta in c(.25, .45, .65, 1)) for (nu in c(.1, 1, 10)) for (m in 1:4) {
  op <- my.fractional.operators.frac(L, beta, C, 16, m, dt, diffusion = nu)
  if (beta == 1) rates <- ev$values else {
    roots <- my.get.roots(m, beta)
    z <- ev$values/16
    rates <- 16^beta / (roots$factor * vapply(z, function(x)
      prod(1-roots$pr_roots*x)/prod(1-roots$pl_roots*x), numeric(1)))
  }
  rhs <- as.numeric(C %*% initial) + dt*load
  expected <- as.numeric(V %*% ((t(V) %*% rhs)/(1+dt*nu*rates)))
  near(as.numeric(my.solver.frac(op, rhs)), expected)
  uncached <- op
  uncached$factorizations <- uncached$factorization <- NULL
  near(as.numeric(my.solver.frac(op, rhs)), as.numeric(my.solver.frac(uncached, rhs)))
}
for (nu in list(0, -1, NA_real_, Inf, c(1, 2)))
  expect_error(my.fractional.operators.frac(L, .5, C, 16, 2, dt, diffusion = nu))
expect_error(my.fractional.operators.frac(L, .5, C, 16, 5, dt))
cat("PASS: diffusion, physical forcing, cached solves, integer branch, validation\n")

# Compare with the actual main-branch constructor, not a copy of new logic.
old_text <- system2("git", c("show", "main:fun/fractional_evolution.R"), stdout = TRUE)
old_file <- tempfile(fileext = ".R")
writeLines(old_text, old_file)
old <- new.env(parent = globalenv())
load_functions(old_file, old)
for (beta in c(.25, .45, .65, 1)) for (m in 1:4) {
  new_op <- my.fractional.operators.frac(L, beta, C, 16, m, dt)
  old_op <- old$my.fractional.operators.frac(L, beta, C, 16, m, dt)
  near(as.numeric(my.solver.frac(new_op, initial)),
       as.numeric(old$my.solver.frac(old_op, initial)))
}
cat("PASS: default coefficient reproduces main-branch update\n")

for (rate in c(0, 2, 100)) for (tt in c(0, 1e-12, 1e-5, .1)) {
  expected <- tt*integrate(function(x)
    exp(-rate*tt*(1-x))*sin(pi*tt*x), 0, 1, abs.tol = 1e-30)$value
  got <- decayed_sine_convolution(tt, rate)
  stopifnot(abs(got-expected) <= 1e-12*max(abs(expected), 1e-30))
}
stopifnot(is.finite(decayed_sine_convolution(1000, 1000)))
cat("PASS: exact forcing response at zero, tiny times, and large decay\n")

spatial <- prepare_tadpole_experiment(.2, 4, c(5L, 7L), 4L)
graph <- gets.graph.tadpole(.2)
graph$compute_fem()
near(10*spatial$basis$projection[, 1], project_u0_tadpole(graph), 1e-11)
near(spatial$basis$gram+spatial$basis$residual_gram, diag(2), 1e-11)

# Independent spatial/time Gauss integration of (exact minus FEM)^2,
# using actual interval values rather than the projection decomposition.
direct_error <- function(nu, alpha, grid) {
  rule <- gauss_legendre_rule(12)
  PtE <- graph$mesh$PtE
  points <- list(); weights <- numeric(); rows <- integer(); cols <- integer(); entries <- numeric()
  npoints <- 0L
  modes <- matrix(numeric(), 0, 2)
  for (edge in 1:2) {
    ix <- which(PtE[, 1] == edge); ix <- ix[order(PtE[ix, 2])]
    x <- c(0, c(1,2)[edge]*PtE[ix,2], c(1,2)[edge])
    nodes <- c(if (edge == 1L) 1L else 2L, 2L+ix, 2L)
    for (k in seq_len(length(x)-1L)) {
      xx <- x[k]+rule$nodes*(x[k+1]-x[k])
      ee5 <- if (edge == 1L) 2*cos(3*pi*xx/2)/sqrt(3) else sin(3*pi*xx/2)/sqrt(3)
      ee7 <- cos(2*pi*xx)/sqrt(3/2)
      modes <- rbind(modes, cbind(ee5, ee7))
      ids <- npoints+seq_along(xx)
      rows <- c(rows, ids, ids)
      cols <- c(cols, rep(nodes[k], length(xx)), rep(nodes[k+1], length(xx)))
      entries <- c(entries, 1-rule$nodes, rule$nodes)
      weights <- c(weights, (x[k+1]-x[k])*rule$weights)
      npoints <- npoints+length(xx)
    }
  }
  Psi <- sparseMatrix(i = rows, j = cols, x = entries,
                      dims = c(npoints, nrow(graph$mesh$C)))
  rates <- nu*spatial$basis$eigenvalues^(alpha/2)
  u <- as.numeric(spatial$basis$projection %*% c(10,0))
  load <- as.numeric(spatial$basis$loads %*% c(0,10))
  op <- my.fractional.operators.frac(spatial$L, alpha/2, spatial$basis$C, 16, 4,
                                    grid$time_step, diffusion = nu)
  g <- average_sine_forcing(grid$time_seq)
  total <- 0
  for (k in seq_len(grid$N)) {
    u <- as.numeric(my.solver.frac(op, op$C %*% u + grid$time_step*load*g[k+1L]))
    times <- grid$time_seq[k] + rule$nodes*grid$time_step
    b <- exact_modal_coefficients(times, rates, c(10,0), c(0,10))
    residual <- modes %*% b-as.numeric(Psi %*% u)
    total <- total+grid$time_step*sum(rule$weights*colSums(residual^2*weights))
  }
  sqrt(total)
}
grid <- make_time_grid(.1, .01)
for (nu in c(.1, 1, 10)) {
  modal <- tadpole_case_error(spatial, .9, 4, nu, grid, c(10,0), c(0,10), quadrature_order = 8)
  direct <- direct_error(nu, .9, grid)
  near(modal, direct, 1e-10)
}
coarse <- tadpole_case_error(spatial, .9, 4, 1, grid, c(10,0), c(0,10))
fine <- tadpole_case_error(spatial, .9, 4, 1, grid, c(10,0), c(0,10), quadrature_order = 8)
near(coarse, fine, 1e-10)
spatial_a <- prepare_tadpole_experiment(.001, 4, c(5L,7L), 4L)
spatial_b <- prepare_tadpole_experiment(.0005, 4, c(5L,7L), 4L)
ratio <- spatial_a$basis$residual_gram[1,1]/spatial_b$basis$residual_gram[1,1]
stopifnot(abs(ratio-16) < .01, all(diag(spatial_b$basis$residual_gram) > 0))
cat("PASS: analytic loads, modal-vs-direct graph error, quadrature, fine-mesh cancellation\n")

root <- tempfile("diffusion-output-check-")
parameters <- list(diffusion = 1, h = .1)
one <- experiment_run("h", parameters, root = root)
two <- experiment_run("h", parameters, root = root)
stopifnot(one != two)
save_experiment_results(list(errors = c(1,2)), one)
near(load_experiment_results(one)$errors, c(1,2))
expect_error(save_experiment_results(list(errors = 3), one))
expect_error(experiment_run("h", list(diffusion = 2), basename(one), root))
expect_error(experiment_run("h", parameters, "..", root))
stopifnot(experiment_run("h", parameters, basename(one), root) == one)
cat("PASS: unique results, load validation, overwrite protection\n")

cached_time <- prepare_modal_time_data(grid,
  spatial$basis$eigenvalues^.45, c(10,0), c(0,10))
cached_error <- tadpole_case_error(spatial, .9, 4, 1, grid, c(10,0), c(0,10),
                                  time_data = cached_time)
near(cached_error, coarse, 1e-11)
bad_time <- cached_time
bad_time$rates <- bad_time$rates*2
expect_error(tadpole_case_error(spatial, .9, 4, 1, grid, c(10,0), c(0,10),
                                time_data = bad_time))
cache <- new_diffusion_cache(.001)
stopifnot(diffusion_cache_put(cache, "one", numeric(100)))
stopifnot(diffusion_cache_put(cache, "two", numeric(100)))
stopifnot(is.null(diffusion_cache_get(cache, "one")),
          !is.null(diffusion_cache_get(cache, "two")),
          !diffusion_cache_put(cache, "oversize", numeric(10000)))
cat("PASS: cached/streamed quadrature equivalence, cache validation and eviction\n")

# The fast evaluator must agree with independently assembled sparse FEM steps.
# These tests exercise projection norms, all graph mode families, nonzero
# physical forcing, both coefficient interpolations, and several decay scales.
for (h in c(.1, .04, .01)) {
  assembled <- prepare_tadpole_experiment(h, 4, 1:7, 4L)
  uniform <- uniform_tadpole_basis(h, 4, 1:7, 4L)
  near(as.numeric(uniform$gram), as.numeric(assembled$basis$gram), 2e-11)
  near(as.numeric(uniform$residual_gram), as.numeric(assembled$basis$residual_gram), 2e-11)
}
for (N in c(1, 2, 7, 50)) {
  rule <- discrete_uniform_rule(N, 64L)
  for (power in 0:5) near(sum(rule$weights*(rule$nodes/N)^power),
                          mean((seq_len(N)/N)^power), 1e-12)
}
rule <- discrete_uniform_rule(1e12, 64L)
near(sum(rule$weights), 1, 1e-13)
# Reduced-order Gaussian step summation vs explicit summation over all steps.
gr <- experiment_grid(.5, .5/200)
for (nu in c(1e-4, 10)) {
  a <- uniform_tadpole_error(.01, 1.3, 4, nu, gr, 4, c(5L,7L), 4L,
                             c(10,0), c(0,10), step_sum_order = 32L)
  b <- uniform_tadpole_error(.01, 1.3, 4, nu, gr, 4, c(5L,7L), 4L,
                             c(10,0), c(0,10), step_sum_order = 200L)
  near(a, b, 1e-11)
}
for (m in 1:4) for (beta in rSPDE:::m4table$beta[c(1, 20, 50)]) {
  a <- my.get.roots(m, beta, "linear")
  b <- my.get.roots(m, beta, "cubic")
  near(unlist(a), unlist(b), 1e-12)
}
for (alpha in c(.9, 1.3)) {
  linear <- max(abs(rational_modal_rates(c(38,55), alpha, 4, 1, 4, "linear")/
                      c(38,55)^(alpha/2)-1))
  cubic <- max(abs(rational_modal_rates(c(38,55), alpha, 4, 1, 4, "cubic")/
                     c(38,55)^(alpha/2)-1))
  stopifnot(cubic < linear/5)
}
for (method in c("linear", "cubic")) for (alpha in c(.5, .9, 1.3, 2))
  for (nu in c(1e-4, .7, 10)) {
    h <- .04
    assembled <- prepare_tadpole_experiment(h, 4, c(5L,7L), 4L)
    gr <- make_time_grid(.5, .01)
    for (m in c(1L,4L)) {
      direct <- tadpole_case_error(assembled, alpha, m, nu, gr, c(10,0), c(0,10),
                                   root_interpolation = method)
      fast <- uniform_tadpole_error(h, alpha, m, nu, gr, 4, c(5L,7L), 4L,
                                    c(10,0), c(0,10), root_interpolation = method)
      near(fast, direct, 2e-10)
    }
  }
# Closed recurrence checked independently at every step, including zero forcing
# frequency and very small time steps. Verify constant initial data separately.
for (dt in c(.01, 1e-11)) for (omega in c(0, pi)) {
  rates <- c(.0004, 10)
  gr <- make_time_grid(20*dt, dt)
  avg <- average_sine_forcing(gr$time_seq, omega = omega)
  # The existing forcing helper assumes nonzero frequency; sin(0*t) is zero.
  if (omega == 0) avg[] <- 0
  value <- c(10,0)
  answer <- matrix(0, 2, gr$N)
  for (k in seq_len(gr$N)) {
    value <- (value+dt*c(0,10)*avg[k+1])/(1+dt*rates)
    answer[, k] <- value
  }
  near(as.numeric(uniform_modal_coefficients(seq_len(gr$N), dt, rates,
       c(10,0), c(0,10), omega = omega)), as.numeric(answer), 1e-12)
}
for (profile in c("pilot", "full")) {
  configurations <- lapply(c("h", "tau", "m"), convergence_settings, profile = profile)
  stopifnot(identical(configurations[[1]]$shared_parameters, configurations[[2]]$shared_parameters),
            identical(configurations[[1]]$experiment_parameters$suite_id,
                      configurations[[3]]$experiment_parameters$suite_id),
            configurations[[1]]$T_final == 1, configurations[[1]]$diffusion == 1,
            max(configurations[[1]]$h_vector)/min(configurations[[1]]$h_vector) >= 100*(1-1e-12),
            max(configurations[[2]]$tau_vector)/min(configurations[[2]]$tau_vector) >= 100*(1-1e-12),
            length(configurations[[1]]$time_grid) == 3L)
}
expect_error(uniform_tadpole_basis(1, 4, c(5L,7L), 4L))
expect_error(discrete_uniform_rule(2^53))
expect_error(uniform_tadpole_error(.01, 1.3, 4, 1e8, experiment_grid(.5, 1e-5),
                                   4, c(5L,7L), 4L, c(10,0), c(0,10)))
cat("PASS: uniform modal evaluator vs sparse FEM, discrete step quadrature, closed recurrence and shared settings\n")

small_grid <- make_time_grid(.001, 1e-5)
study_args <- list(study = "m", axis_values = 1:4, alpha_vector = c(.5, .9),
  time_grids = list(small_grid), fixed_h = .04, anchor_h = .04,
  anchor_grid = small_grid, anchor_m = 4L, diffusion = .7, kappa = 4,
  modes = c(5L,7L), N_finite = 4L, initial = c(10,0), forcing = c(0,10),
  run_refinement_checks = TRUE, verbose = FALSE)
serial <- do.call(run_tadpole_study, c(study_args, list(workers = 1L)))
stopifnot(any(serial$resources$temporal_cache_hit),
          any(serial$resources$spatial_cache_hit),
          all(serial$resources$temporal_cache_mib <= 512))
thread_env <- Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"),
                         unset = NA_character_)
connections_before <- nrow(showConnections())
parallel_result <- do.call(run_tadpole_study, c(study_args,
  list(workers = 2L, worker_log = file.path(root, "parallel-test.log"))))
near(as.numeric(parallel_result$errors_projected), as.numeric(serial$errors_projected), 1e-11)
near(parallel_result$refinement_checks$refined_error, serial$refinement_checks$refined_error, 1e-11)
stopifnot(identical(parallel_result$diagnostics, serial$diagnostics),
          identical(Sys.getenv(names(thread_env), unset = NA_character_), thread_env),
          nrow(showConnections()) == connections_before)
streamed_args <- study_args
streamed_args$run_refinement_checks <- FALSE
streamed <- do.call(run_tadpole_study, c(streamed_args,
  list(workers = 1L, temporal_cache_mib = 0)))
near(as.numeric(streamed$errors_projected), as.numeric(serial$errors_projected), 1e-11)
stopifnot(!any(streamed$resources$temporal_cache_hit),
          all(streamed$resources$temporal_cache_mib == 0))
bad_args <- study_args
bad_args$axis_values <- c(1, 5)
expect_error(do.call(run_tadpole_study, c(bad_args,
  list(workers = 2L, worker_log = file.path(root, "worker-error-test.log")))))
stopifnot(nrow(showConnections()) == connections_before,
          identical(Sys.getenv(names(thread_env), unset = NA_character_), thread_env))
expect_error(diffusion_resource_plan(0L, 8L, .01))
expect_error(diffusion_resource_plan(4L, 8L, .01, memory_budget_gib = .05))
cat("PASS: serial/parallel equivalence, cache reuse, streamed fallback, worker cleanup and resource limits\n")

if ("--pilots" %in% commandArgs(trailingOnly = TRUE)) {
  for (study in c("h", "tau", "m")) {
    cat("Running notebook numerical chunks:", study, "\n")
    text <- readLines(paste0("conv_in_", study, ".Rmd"))
    starts <- which(grepl("^```\\{r", text))
    env <- new.env(parent = helpers)
    # Keep smoke-run outputs in /tmp; exercise exactly the notebook settings.
    env$experiment_run <- function(study, parameters, run_id = NULL, execution = NULL)
      helpers$experiment_run(study, parameters, run_id, root = root, execution = execution)
    for (start in starts) {
      end <- which(seq_along(text) > start & text == "```")[1]
      code <- text[(start+1L):(end-1L)]
      if (any(grepl("^# Parameters:|^# Refuse to overwrite|^# Load only this run", code)))
        eval(parse(text = code), env)
    }
    stopifnot(all(is.finite(env$errors_projected)), all(env$errors_projected > 0),
              nrow(env$results$refinement_checks) == 15L,
              all(env$results$diagnostics$positive_sampled),
              nrow(env$results$resources) > 0L,
              file.exists(file.path(env$run_dir, "resources.csv")))
    cat("PASS notebook", study, "; saved in", env$run_dir, "\n")
  }
}
cat("All diffusion checks passed. Temporary test output (removed when R exits):", root, "\n")
