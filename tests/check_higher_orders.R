# /usr/bin/Rscript --vanilla tests/check_higher_orders.R [--notebooks]
# Production tests for package m=1:4 and directly computed MATLAB m=1:8.
project_lib <- file.path(getwd(), "renv/library/linux-ubuntu-noble/R-4.5/x86_64-pc-linux-gnu")
if (dir.exists(project_lib)) .libPaths(c(project_lib, .libPaths()))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(MetricGraph))
helpers <- new.env(parent = globalenv())
for (name in c("fractional_evolution", "graphs_and_meshes", "exact_solutions",
               "convergence_experiments", "convergence_settings")) {
  for (ex in parse(file.path("fun", paste0(name, ".R"))))
    if (is.call(ex) && identical(ex[[1]], as.name("<-")) &&
        is.call(ex[[3]]) && identical(ex[[3]][[1]], as.name("function"))) eval(ex, helpers)
}
attach(helpers, name = "higher-order-test-helpers")
near <- function(a, b, tolerance = 1e-10) stopifnot(all(is.finite(c(a,b))),
                                                max(abs(a-b)) <= tolerance)
expect_error <- function(code, pattern = NULL) {
  error <- tryCatch({force(code); NULL}, error = identity)
  stopifnot(inherits(error, "error"))
  if (!is.null(pattern)) stopifnot(grepl(pattern, conditionMessage(error), fixed = TRUE))
}

near(poly.from.roots(c(-2,-3)), c(6,5,1), 0)
stopifnot(identical(rational_supported_orders("package"), 1:4),
          identical(rational_supported_orders("computed"), 1:8))
expect_error(my.get.roots(5, .25, coefficient_source = "package"), "supports m")
expect_error(my.get.roots(9, .25, coefficient_source = "computed"), "supports m")
expect_error(my.get.roots(4, .5, coefficient_source = "computed"), "unavailable at beta")
expect_error(my.get.roots(4, .25, coefficient_source = "unknown"))
expect_error(validate_fractional_roots(list(factor = 1, pr_roots = .1,
                                           pl_roots = c(-2,-1)), 1))
expect_error(validate_fractional_roots(list(factor = 1, pr_roots = -3,
                                           pl_roots = c(-2,-1)), 1), "interlace")
for (m in 1:8) for (beta in c(.25,.45,.65)) {
  a <- my.get.roots(m, beta, "linear", "computed")
  b <- my.get.roots(m, beta, "cubic", "computed")
  stopifnot(identical(a,b))
}
# Existing package coefficients are unchanged by introducing source selection.
for (m in 1:4) for (beta in c(.25,.45,.65)) {
  tab <- get(paste0("m", m, "table"), asNamespace("rSPDE"))
  r <- my.get.roots(m, beta, "linear", "package")
  near(r$factor, approx(tab$beta, tab$factor, beta)$y, 0)
  near(r$pr_roots, vapply(if (m == 1) "rc" else paste0("rc.",1:m),
                          function(col) approx(tab$beta, tab[[col]], beta)$y, numeric(1)), 0)
}
cat("PASS: source selection, package compatibility, exact-beta lookup and root validation\n")

z <- sort(unique(c(0, exp(seq(log(.001),log(1e12), length.out = 180)))))
worst_reconstruction <- 0
for (source in c("package", "computed")) for (m in rational_supported_orders(source))
  for (beta in c(.25,.45,.65)) for (strength in c(1e-16,1e-11,1e-7,1e-4,.01,1,100)) {
    roots <- my.get.roots(m, beta, "cubic", source)
    pf <- do.call(compute.partial.fraction.param,
                  c(roots, list(time_step = strength, scaling = 1)))
    stopifnot(all(pf$p < 0), all(pf$r > 0), all(is.finite(c(pf$p,pf$r))))
    inverse <- roots$factor*vapply(z, function(x)
      prod(1-roots$pr_roots*x)/prod(1-roots$pl_roots*x), numeric(1))
    expected <- inverse/(inverse+strength)
    got <- rowSums(vapply(seq_along(pf$p), function(i) pf$r[i]/(z-pf$p[i]), numeric(length(z))))
    worst_reconstruction <- max(worst_reconstruction, max(abs(got-expected)))
    near(got, expected, 1e-9)
    stopifnot(max(abs(got/expected-1)) < 1e-7)
  }
cat("PASS: positive poles/residues and reconstruction over 15 spectral decades; worst absolute discrepancy",
    format(worst_reconstruction, scientific = TRUE), "\n")

# Independent generalized eigenbasis: nonidentity/noncommuting mass matrix,
# physical forcing, several time steps and diffusion values, both sources.
C <- Matrix(matrix(c(2,.2,.2,3),2), sparse = TRUE)
L <- Matrix(matrix(c(40,-2,-2,160),2), sparse = TRUE)
Ri <- solve(chol(as.matrix(C)))
ev <- eigen(t(Ri) %*% as.matrix(L) %*% Ri, symmetric = TRUE)
V <- Ri %*% ev$vectors
initial <- c(1,2)
forcing <- c(3,-4)
worst_sparse <- 0
for (source in c("package", "computed")) for (m in rational_supported_orders(source))
  for (beta in c(.25,.45,.65)) for (dt in c(1e-11,1e-4,.01)) for (nu in c(1e-4,1,10)) {
    roots <- my.get.roots(m, beta, "cubic", source)
    x <- ev$values/16
    inverse <- roots$factor*vapply(x, function(z)
      prod(1-roots$pr_roots*z)/prod(1-roots$pl_roots*z), numeric(1))
    rates <- nu*16^beta/inverse
    rhs <- as.numeric(C %*% initial)+dt*forcing
    expected <- as.numeric(V %*% ((t(V) %*% rhs)/(1+dt*rates)))
    op <- my.fractional.operators.frac(L, beta, C, 16, m, dt,
      diffusion = nu, root_interpolation = "cubic", coefficient_source = source)
    got <- as.numeric(my.solver.frac(op, rhs))
    worst_sparse <- max(worst_sparse, max(abs(got-expected)))
    near(got, expected, 1e-9)
  }
cat("PASS: sparse solves and physical forcing for both sources; worst absolute discrepancy",
    format(worst_sparse, scientific = TRUE), "\n")

for (profile in c("pilot", "full")) {
  package <- convergence_settings("m", profile, "package")
  computed <- convergence_settings("m", profile, "computed", m_max = 8L)
  stopifnot(identical(package$m_vector,1:4), identical(computed$m_vector,1:8),
            !identical(package$experiment_parameters$suite_id, computed$experiment_parameters$suite_id))
  # Retain the legacy package problem while calibrating the computed source
  # for unit diffusion. All remaining data/conventions agree across presets.
  a <- package$shared_parameters; a[c("T_final","diffusion","coefficient_info","m_study_orders","auxiliary","error_method")] <- NULL
  b <- computed$shared_parameters; b[c("T_final","diffusion","coefficient_info","m_study_orders","auxiliary","error_method")] <- NULL
  stopifnot(identical(a,b), package$m_fixed == 4L, computed$m_fixed == 8L,
            package$T_final == .5, package$diffusion == 1e-4,
            computed$T_final == 1, computed$diffusion == 1,
            package$h_fixed == if (profile == "pilot") 1e-6 else 5e-7,
            package$time_grid$time_step == if (profile == "pilot") 1e-11 else 5e-12,
            computed$h_fixed == if (profile == "pilot") 1e-7 else 5e-8,
            computed$time_grid$N == if (profile == "pilot") 1e14 else 2e14,
            2*computed$time_grid$N <= 2^52)
  for (max_order in 4:8) {
    selected <- convergence_settings("m",profile,"computed",m_max = max_order)
    stopifnot(identical(selected$m_vector,seq_len(max_order)), selected$m_fixed == 8L,
              identical(selected$shared_parameters$m_study_orders,selected$m_vector))
    for (alpha in selected$alpha_vector) {
      evaluate <- function(h = selected$h_fixed, grid = selected$time_grid, m = selected$anchor_m)
        uniform_tadpole_error(h,alpha,m,selected$diffusion,grid,selected$kappa,
          selected$active_modes,selected$N_finite,selected$exact_initial,selected$exact_forcing,
          selected$AAA,selected$OMEGA,selected$quadrature_order,selected$step_sum_order,
          selected$root_interpolation,"computed")
      error <- evaluate()
      stopifnot(evaluate(m = 0L)/error < .01,
                abs(evaluate(h = selected$h_fixed/2)/error-1) < .01,
                abs(evaluate(grid = experiment_grid(1,selected$time_grid$time_step/2))/error-1) < .01)
    }
  }
  for (source in c("package","computed")) for (study in c("h","tau")) {
    s <- convergence_settings(study, profile, source)
    m <- convergence_settings("m",profile,source)
    stopifnot(identical(s$shared_parameters,m$shared_parameters),
              identical(s$m_vector,1:4), s$m_max == 4L,
              s$m_fixed == if (source == "package") 4L else 8L)
  }
  default <- convergence_settings("m",profile,"computed")
  stopifnot(default$h_fixed == if (profile == "pilot") 1e-5 else 5e-6,
            default$time_grid$time_step == if (profile == "pilot") 1e-10 else 5e-11,
            !identical(default$experiment_parameters$suite_id,computed$experiment_parameters$suite_id))
}
for (order in c(3,4.5,9)) expect_error(convergence_settings("m",m_max = order))
expect_error(convergence_settings("m",coefficient_source = "package",m_max = 5L))
z <- c(16+(3*pi/2)^2,16+(2*pi)^2)
for (alpha in c(.5,.9,1.3)) {
  rate <- rational_modal_rates(z, alpha, 8, 1, 4, "cubic", "computed")
  stopifnot(max(abs(rate/z^(alpha/2)-1)) < 2e-10)
}
root <- tempfile("higher-order-check-")
dir.create(root)
package <- convergence_settings("m", coefficient_source = "package")
computed <- convergence_settings("m", coefficient_source = "computed")
run <- experiment_run("m", package$experiment_parameters, root = root)
expect_error(experiment_run("m", computed$experiment_parameters, basename(run), root),
             "parameters do not match")
cat("PASS: independent m range/order, unit-diffusion calibration (<1% floor for m_max=4:8), legacy package and run separation\n")

# Ordinary recurrence formulas evaluated independently with 100-digit Decimal
# arithmetic. Near 1e-14, subtracting binary64 solution coefficients is not an
# adequate reference for this cancellation-sensitive error calculation.
reference <- read.csv(text = system2("python3", "tests/reference_modal_difference.py", stdout = TRUE))
for (i in seq_len(nrow(reference))) {
  v <- reference[i, ]
  got <- uniform_modal_difference(v$k, v$fraction, v$time_step, v$exact_rate,
    v$numerical_rate, v$initial, v$forcing, v$A, v$omega)
  near(as.numeric(got), v$difference, max(1e-29, 5e-13*abs(v$difference)))
}
for (dt in c(.1,.01,1e-6)) for (fraction in c(0,.37,1)) for (omega in c(0,-pi,pi)) {
  k <- c(1,2,5)
  exact <- exact_modal_coefficients((k-1+fraction)*dt,c(.17,.63),c(3,-2),c(-4,6),.7,omega)
  numerical <- uniform_modal_coefficients(k,dt,c(.18,.59),c(3,-2),c(-4,6),.7,omega)
  stable <- uniform_modal_difference(k,fraction,dt,c(.17,.63),c(.18,.59),c(3,-2),c(-4,6),.7,omega)
  near(stable,exact-numerical,4e-15)
}
for (rates in list(c(2000,10),c(10,2000))) {
  stable <- uniform_modal_difference(1e6,.37,1e-6,rates[1],rates[2],7,3)
  direct <- exact_modal_coefficients((1e6-1+.37)*1e-6,rates[1],7,3)-
    uniform_modal_coefficients(1e6,1e-6,rates[2],7,3)
  near(stable,direct,2e-15)
}
cat("PASS: stable modal differences against independent 100-digit recurrence and coarse-grid formulas\n")

grid <- make_time_grid(.01,1e-4)
args <- list(study = "m", axis_values = 1:8, alpha_vector = c(.5,.9,1.3),
  time_grids = list(grid), fixed_h = .04, anchor_h = .04, anchor_grid = grid,
  anchor_m = 8L, diffusion = .7, kappa = 4, modes = c(5L,7L), N_finite = 4L,
  initial = c(10,0), forcing = c(0,10), run_refinement_checks = FALSE,
  root_interpolation = "cubic", coefficient_source = "computed", verbose = FALSE)
serial <- do.call(run_tadpole_study,c(args,list(workers = 1L)))
parallel_result <- do.call(run_tadpole_study,c(args,list(workers = 2L,
  worker_log = file.path(root,"computed-workers.log"))))
near(as.numeric(serial$errors_projected),as.numeric(parallel_result$errors_projected),1e-11)
stopifnot(all(parallel_result$diagnostics$coefficient_source == "computed"),
          identical(serial$coefficient_info,rational_coefficient_info("computed")))
modal <- do.call(run_tadpole_study,c(args,list(workers = 1L,solver_backend = "uniform-modal")))
near(as.numeric(serial$errors_projected),as.numeric(modal$errors_projected),1e-9)
unit_args <- args
unit_args$axis_values <- c(1L,4L,8L)
unit_args$time_grids <- list(make_time_grid(1,.002))
unit_args$anchor_grid <- unit_args$time_grids[[1]]
unit_args$diffusion <- 1
unit_direct <- do.call(run_tadpole_study,c(unit_args,list(workers = 1L)))
unit_modal <- do.call(run_tadpole_study,c(unit_args,list(workers = 1L,solver_backend = "uniform-modal")))
near(as.numeric(unit_direct$errors_projected),as.numeric(unit_modal$errors_projected),1e-10)
stopifnot(max(abs(unit_direct$errors_projected/unit_modal$errors_projected-1)) < 1e-8)
cat("PASS: T=1, diffusion=1 closed recurrence against direct sparse stepping for m=1,4,8 and all alpha\n")
args$alpha_vector <- 1.1
expect_error(do.call(run_tadpole_study,c(args,list(workers = 2L))), "unavailable at beta")
cat("PASS: higher orders through direct/modal backends and parallel workers; unsupported beta rejected early\n")

if ("--notebooks" %in% commandArgs(trailingOnly = TRUE)) {
  for (study in c("h","tau","m")) {
    lines <- readLines(paste0("conv_in_",study,".Rmd"))
    env <- new.env(parent = helpers)
    env$convergence_settings <- function(study)
      helpers$convergence_settings(study,coefficient_source = "computed")
    env$experiment_run <- function(study,parameters,run_id=NULL,execution=NULL)
      helpers$experiment_run(study,parameters,run_id,root=root,execution=execution)
    for (start in which(grepl("^```\\{r",lines))) {
      end <- which(seq_along(lines)>start & lines=="```")[1]
      code <- lines[(start+1):(end-1)]
      if (any(grepl("^# Parameters:|^# Refuse to overwrite|^# Load only this run",code)))
        eval(parse(text=code),env)
    }
    stopifnot(all(is.finite(env$errors_projected)), all(env$errors_projected>0),
              all(env$results$diagnostics$coefficient_source=="computed"),
              identical(env$results$parameters$shared$coefficient_info,
                        env$results$coefficient_info))
    if (study=="m") {
      stopifnot(nrow(env$errors_projected)==env$m_max,env$m_fixed==8L)
      checks <- env$results$refinement_checks
      floor <- checks$check=="spatial-temporal floor / total error"
      stopifnot(all(checks$required[floor]),all(checks$passed[checks$required]),
                all(checks$relative_change[floor] < .01),
                all(apply(env$errors_projected,2,function(e) all(diff(e)<0))))
    }
    cat("PASS computed-source notebook:",study,"\n")
  }
}
cat("All higher-order production checks passed. Temporary results:",root,"\n")
