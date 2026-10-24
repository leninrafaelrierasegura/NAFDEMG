# Investigate m=1,...,8 without changing solvers, coefficient tables or notebooks.
# /usr/bin/Rscript --vanilla tests/investigate_rational_orders.R [--full]
# --full includes every fractional-beta knot in both existing table sets.
project_lib <- file.path(getwd(), "renv/library/linux-ubuntu-noble/R-4.5/x86_64-pc-linux-gnu")
if (dir.exists(project_lib)) .libPaths(c(project_lib, .libPaths()))
stopifnot(getRversion() >= "4.5.0")
suppressPackageStartupMessages(library(Matrix))

input_files <- c("matlab/chebfun.m", "matlab/chebfun_tables.R",
  "matlab/chebfun_tables.RDS", "data_files/chebfun_tables.RDS",
  sprintf("matlab/m%d.bin", 1:8), sprintf("matlab/m%dinfo.txt", 1:8),
  paste0("conv_in_", c("h", "tau", "m"), ".Rmd"),
  "fun/convergence_settings.Rmd", "fun/fractional_evolution.Rmd")
input_hashes <- tools::md5sum(input_files)

load_functions <- function(env) {
  for (name in c("fractional_evolution", "graphs_and_meshes", "exact_solutions",
                 "convergence_experiments", "convergence_settings")) {
    for (ex in parse(file.path("fun", paste0(name, ".R")))) {
      if (is.call(ex) && identical(ex[[1]], as.name("<-")) &&
          is.call(ex[[3]]) && identical(ex[[3]][[1]], as.name("function"))) eval(ex, env)
    }
  }
}
production <- new.env(parent = globalenv())
load_functions(production)
# Freeze the pre-fix algorithms so the investigation remains reproducible after
# the production implementation gains the numerical fixes.
production$poly.from.roots <- function(roots) {
  out <- 1
  for (r in roots) out <- convolve(out, c(1,-r), type = "open")
  out
}
legacy_partial_fractions <- function(factor, pr_roots, pl_roots, time_step, scaling) {
  numerator <- c(0, poly.from.roots(pr_roots))
  denominator <- numerator+time_step*scaling/factor*poly.from.roots(pl_roots)
  poles <- polyroot(rev(denominator))
  stopifnot(all(abs(Im(poles)) <= 1e-10*pmax(1, abs(Re(poles)))), all(Re(poles) < 0))
  poles <- Re(poles)
  evaluate <- function(coef, x) Reduce(function(value, a) value*x+a, coef, init = 0)
  residues <- vapply(seq_along(poles), function(j)
    evaluate(numerator, poles[j])/(denominator[1]*prod(poles[j]-poles[-j])), numeric(1))
  list(r = residues, p = poles, k = numeric())
}
production$compute.partial.fraction.param <- legacy_partial_fractions
environment(production$compute.partial.fraction.param) <- production
guard_message <- tryCatch({production$my.get.roots(5L, .25); "No rejection"},
                          error = conditionMessage)

read_binary_table <- function(m) {
  info <- sprintf("matlab/m%dinfo.txt", m)
  path <- sprintf("matlab/m%d.bin", m)
  dims <- scan(info, what = integer(), quiet = TRUE)
  stopifnot(length(dims) == 2L, dims[1] >= 2L, dims[2] == 2L*m+3L,
            file.info(path)$size == 8*prod(dims))
  con <- file(path, "rb")
  on.exit(close(con))
  values <- readBin(con, "double", n = prod(dims), size = 8, endian = "little")
  stopifnot(length(values) == prod(dims), all(is.finite(values)))
  out <- matrix(values, nrow = dims[1], ncol = dims[2])
  colnames(out) <- c("beta", "factor", if (m == 1) "rc" else paste0("rc.", 1:m),
                     paste0("rb.", 1:(m+1)))
  as.data.frame(out)
}
binary_tables <- lapply(1:8, read_binary_table)
matlab_tables <- readRDS("matlab/chebfun_tables.RDS")
stopifnot(length(matlab_tables) == 8L,
          all(vapply(1:8, function(m) identical(matlab_tables[[m]], binary_tables[[m]]), logical(1))))
sources <- list(matlab200 = matlab_tables,
                legacy1000 = readRDS("data_files/chebfun_tables.RDS"))

interpolate_roots <- function(tables, m, beta, method = "cubic") {
  mt <- tables[[m]]
  at <- function(col) if (method == "linear") approx(mt$beta, mt[[col]], beta)$y else
    splinefun(mt$beta, mt[[col]], method = "fmm")(beta)
  list(factor = at("factor"),
       pr_roots = vapply(if (m == 1) "rc" else paste0("rc.", 1:m), at, numeric(1)),
       pl_roots = vapply(paste0("rb.", 1:(m+1)), at, numeric(1)))
}

# Preserve the existing convention exactly: coefficients in descending powers
# of prod(1-root*z), not the monic polynomial prod(z-root). No FFT is used.
direct_product <- function(roots) {
  out <- 1
  for (r in roots) out <- c(-r*out, 0)+c(0, out)
  out
}
stopifnot(identical(direct_product(c(-2,-3)), c(6,5,1)))
candidate <- new.env(parent = globalenv())
load_functions(candidate)
candidate$poly.from.roots <- direct_product
candidate$compute.partial.fraction.param <- legacy_partial_fractions
environment(candidate$compute.partial.fraction.param) <- candidate
factored <- new.env(parent = globalenv())
load_functions(factored)
factored$poly.from.roots <- direct_product
factored$compute.partial.fraction.param <- function(factor, pr_roots, pl_roots,
                                                    time_step, scaling) {
  out <- candidate$compute.partial.fraction.param(factor, pr_roots, pl_roots,
                                                 time_step, scaling)
  derivative <- function(roots, x) sum(vapply(seq_along(roots), function(i)
    -roots[i]*prod(1-roots[-i]*x), numeric(1)))
  strength <- time_step*scaling/factor
  # At a pole, P=-strength*Q. This avoids evaluating a nearly zero P.
  out$r <- vapply(out$p, function(x)
    -strength*prod(1-pl_roots*x)/
      (derivative(pr_roots, x)+strength*derivative(pl_roots, x)), numeric(1))
  out
}

base <- file.path("data_files", "diffusion_experiments", "rational_order_diagnostics")
dir.create(base, recursive = TRUE, showWarnings = FALSE)
run_dir <- tempfile(paste0(format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"), "-"), base)
stopifnot(dir.create(run_dir))
run_dir <- normalizePath(run_dir)
write.csv(data.frame(file = names(input_hashes), md5 = unname(input_hashes)),
          file.path(run_dir, "input_hashes.csv"), row.names = FALSE)

audit <- list()
for (source in names(sources)) for (m in 1:8) {
  mt <- sources[[source]][[m]]
  fractional <- mt[mt$beta >= .25 & mt$beta < 1, , drop = FALSE]
  audit[[length(audit)+1L]] <- data.frame(source = source, m = m, rows = nrow(mt),
    fractional_rows = nrow(fractional), finite = all(is.finite(as.matrix(mt))),
    beta_increasing = all(diff(mt$beta) > 0),
    fractional_roots_negative = all(as.matrix(fractional[, -(1:2)]) < 0),
    fractional_factor_positive = all(fractional$factor > 0))
}
audit <- do.call(rbind, audit)
stopifnot(all(audit$finite & audit$beta_increasing & audit$fractional_roots_negative &
                audit$fractional_factor_positive))
write.csv(audit, file.path(run_dir, "table_audit.csv"), row.names = FALSE)

# Compare the SAME rational update in two independent representations.
# This is a reconstruction test, not an accuracy claim for z^(-beta) on the
# entire sampled spectrum (which can extend beyond the generator's interval).
z <- sort(unique(c(exp(seq(log(1), log(1e8), length.out = 150)),
                    (16+c((3*pi/2)^2, (2*pi)^2))/16)))
strengths <- c(1e-16, 1e-11, 1e-7, 1e-4, .01, 1, 100)
full <- "--full" %in% commandArgs(trailingOnly = TRUE)
rows <- failures <- list()
for (source in names(sources)) for (m in 1:8) {
  mt <- sources[[source]][[m]]
  betas <- c(.25, .45, .65, .95, .99)
  if (full) betas <- c(betas, mt$beta[mt$beta >= .25 & mt$beta < 1])
  betas <- sort(unique(betas))
  cat("Testing", source, "m =", m, "at", length(betas), "beta values\n")
  for (method in c("fft", "direct", "direct_factored_residues")) {
    env <- switch(method, fft = production, direct = candidate,
                  direct_factored_residues = factored)
    maximum_absolute <- maximum_relative <- 0
    failed <- count <- 0L
    worst_beta <- worst_strength <- NA_real_
    for (beta in betas) for (strength in strengths) {
      count <- count+1L
      roots <- interpolate_roots(sources[[source]], m, beta)
      pf <- tryCatch(do.call(env$compute.partial.fraction.param,
        c(roots, list(time_step = strength, scaling = 1))), error = identity)
      if (inherits(pf, "error")) {
        failed <- failed+1L
        failures[[length(failures)+1L]] <- data.frame(source = source, m = m,
          method = method, beta = beta, strength = strength, error = conditionMessage(pf))
        next
      }
      inverse <- roots$factor*vapply(z, function(x)
        prod(1-roots$pr_roots*x)/prod(1-roots$pl_roots*x), numeric(1))
      target <- inverse/(inverse+strength)
      got <- rowSums(vapply(seq_along(pf$p), function(i) pf$r[i]/(z-pf$p[i]), numeric(length(z))))
      absolute <- max(abs(got-target))
      relative <- max(abs(got/target-1))
      if (!all(is.finite(c(got, absolute, relative))) || any(pf$r <= 0)) {
        failed <- failed+1L
        failures[[length(failures)+1L]] <- data.frame(source = source, m = m,
          method = method, beta = beta, strength = strength, error = "Nonfinite reconstruction or nonpositive residue")
        next
      }
      if (absolute > maximum_absolute) {
        maximum_absolute <- absolute
        worst_beta <- beta
        worst_strength <- strength
      }
      maximum_relative <- max(maximum_relative, relative)
    }
    rows[[length(rows)+1L]] <- data.frame(source = source, m = m, method = method,
      cases = count, exceptions_or_invalid = failed, max_absolute_error = maximum_absolute,
      max_relative_error = maximum_relative, worst_beta = worst_beta,
      worst_strength = worst_strength)
  }
}
stability <- do.call(rbind, rows)
write.csv(stability, file.path(run_dir, "partial_fraction_stability.csv"), row.names = FALSE)
if (length(failures)) write.csv(do.call(rbind, failures), file.path(run_dir, "failures.csv"), row.names = FALSE)

# Separate approximation quality at the active eigenvalues from interpolation
# error. Nearby beta-knot results are diagnostics for a DIFFERENT beta, not
# replacements for the requested experiment parameters.
active <- list()
active_z <- (16+c((3*pi/2)^2, (2*pi)^2))/16
rate_error <- function(roots, beta) {
  R <- roots$factor*vapply(active_z, function(x)
    prod(1-roots$pr_roots*x)/prod(1-roots$pl_roots*x), numeric(1))
  max(abs(1/(R*active_z^beta)-1))
}
for (source in names(sources)) for (beta in c(.25, .45, .65)) for (m in 1:8) {
  mt <- sources[[source]][[m]]
  knot <- mt$beta[which.min(abs(mt$beta-beta))]
  active[[length(active)+1L]] <- data.frame(source = source, m = m, beta = beta,
    linear_rate_error = rate_error(interpolate_roots(sources[[source]], m, beta, "linear"), beta),
    cubic_rate_error = rate_error(interpolate_roots(sources[[source]], m, beta, "cubic"), beta),
    nearby_knot_beta = knot,
    nearby_knot_rate_error = rate_error(interpolate_roots(sources[[source]], m, knot), knot))
}
write.csv(do.call(rbind, active), file.path(run_dir, "active_rate_accuracy.csv"), row.names = FALSE)

# Exercise the actual sparse operator and solve with diagnostic-only overrides.
# No production function or table file is edited or replaced.
for (env in list(candidate, factored)) {
  env$.tables <- matlab_tables
  env$my.get.roots <- function(m, beta, interpolation = c("linear", "cubic"))
    interpolate_roots(.tables, m, beta, match.arg(interpolation))
  environment(env$my.get.roots) <- env
}
C <- Matrix(matrix(c(2,.2,.2,3), 2), sparse = TRUE)
R <- chol(as.matrix(C))
V <- solve(R)
L <- Matrix(16*t(R) %*% diag(c(1,1e4)) %*% R, sparse = TRUE)
initial <- c(2,-.7)
u <- as.numeric(V %*% initial)
sparse_checks <- list()
for (method in c("direct", "direct_factored_residues"))
  for (m in 5:8) for (beta in c(.25,.45,.65)) for (strength in c(1e-15,1e-4,.01)) {
  env <- if (method == "direct") candidate else factored
  roots <- env$my.get.roots(m, beta, "cubic")
  inverse <- roots$factor*vapply(c(1,1e4), function(x)
    prod(1-roots$pr_roots*x)/prod(1-roots$pl_roots*x), numeric(1))
  expected <- as.numeric(V %*% (initial*inverse/(inverse+strength)))
  op <- env$my.fractional.operators.frac(L, beta, C, 16, m,
    time_step = strength/16^beta, root_interpolation = "cubic")
  got <- as.numeric(env$my.solver.frac(op, C %*% u))
  error <- max(abs(got-expected))
  sparse_checks[[length(sparse_checks)+1L]] <- data.frame(method = method, m = m, beta = beta,
    strength = strength, max_absolute_error = error)
}
sparse_checks <- do.call(rbind, sparse_checks)
write.csv(sparse_checks, file.path(run_dir, "sparse_solver_checks.csv"), row.names = FALSE)
stopifnot(max(sparse_checks$max_absolute_error) < 1e-9)
stopifnot(identical(input_hashes, tools::md5sum(input_files)))
capture.output(sessionInfo(), file = file.path(run_dir, "sessionInfo.txt"))

report <- c("# Investigation of rational orders above four", "",
  paste("Run:", format(Sys.time(), tz = "UTC")), "",
  "## Scope", "",
  "The production solver, notebooks, physical parameters and original tables were not changed.",
  "The candidate functions exist only in an isolated R environment inside this diagnostic.",
  "Input hashes were checked again after the investigation.", "",
  "## Findings", "",
  paste("The current production m=5 rejection is:", guard_message),
  "The MATLAB folder already exports m=1,...,8. Its RDS matches all binary files exactly.",
  "The checked-in generator currently loops over 1:4; it does not reproduce all eight exports as written.",
  "The data_files RDS contains a distinct, denser 1000-row table set; it is not the MATLAB 200-row RDS.",
  "Both exported sets have finite entries, increasing beta knots, negative roots and positive factors for 0.25<=beta<1.", "",
  "FFT polynomial convolution loses relative accuracy in small coefficients at higher orders.",
  "An equivalent direct product, preserving prod(1-root*z), substantially improves partial-fraction reconstruction.",
  "The direct_factored_residues prototype also avoids evaluating a nearly zero numerator at the computed pole; inspect sign checks as well as reconstruction accuracy.",
  "See partial_fraction_stability.csv: reconstruction compares the same rational function in factored and partial-fraction forms, not the fractional-power target.",
  paste("Largest sparse-solver discrepancy with the direct-product candidate:",
        format(max(sparse_checks$max_absolute_error), scientific = TRUE)), "",
  "High-order approximation quality is a separate issue. Interpolating roots and factor can produce a floor or nonmonotonic m convergence away from beta knots.",
  "See active_rate_accuracy.csv. Nearby-knot results use another beta and must not replace the requested alpha.", "",
  "## Remaining work before enabling higher orders", "",
  "1. Use the existing m=5,...,8 tables through an explicit, validated loading path.",
  "2. Replace FFT coefficient construction and add high-order reconstruction/solve checks.",
  "3. Resolve or quantify interpolation errors; increasing m alone does not guarantee better results.",
  "4. Audit/regenerate MATLAB coefficients with real/imaginary, pole, degree and approximation-quality checks before export.",
  "5. Only then increase a study's m range; keep the other physical and numerical settings fixed during that experiment.", "",
  "This R diagnostic does not rerun MATLAB. Use tests/audit_rational_generation.m separately to audit generation at exact beta values.",
  "The existing real-only binary files cannot establish whether imaginary parts existed before export.",
  "No claim of support beyond m=8 or straight root-exponential convergence is established by this diagnostic.")
writeLines(report, file.path(run_dir, "findings.md"))
cat("\nPartial-fraction reconstruction summary:\n")
print(stability, row.names = FALSE, digits = 5)
cat("\nAll sparse candidate checks passed; original inputs unchanged.\nReport:",
    file.path(run_dir, "findings.md"), "\n")
