# Project environment: Rscript tests/check_computed_isolated_m.R
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
file <- functions$isolated_default_coefficient_file()
tab <- functions$isolated_read_computed_coefficients(file)
check(nrow(tab) == 264L && identical(sort(unique(tab$m)), 1:8), "Production artifact orders and rows")
parameters <- functions$isolated_m_parameters("computed")
check(identical(parameters$m, 1:8), "Computed source defaults to orders 1--8")
check(identical(parameters$coefficient_source, "computed"), "Computed source recorded")
check(identical(parameters$coefficient_file, normalizePath(file)), "Branch-local file recorded")
check(identical(parameters$coefficient_sha256, digest::digest(file = file, algo = "sha256")), "Artifact fingerprint")
blob <- system2("git", c("hash-object", "tests/audit_rational_generation.m"), stdout = TRUE)
check(identical(blob, parameters$coefficient_generator$git_blob), "Copied generator matches verified origin")
package <- functions$isolated_m_parameters("package")
check(identical(package$m, 1:4), "Package default order range retained")
check(identical(package$coefficient_source, "original rSPDE package tables with linear interpolation"),
      "Historical package metadata retained")
check(fails(functions$isolated_fractional_roots(5, .25, "package")), "Package source refuses higher orders")
check(fails(functions$isolated_fractional_roots(9, .25, "computed")), "Computed source refuses unsupported orders")
check(fails(functions$isolated_fractional_roots(4, .5, "computed")), "No beta interpolation or package fallback")
check(fails(functions$isolated_read_computed_coefficients(tempfile())), "Missing artifact")

fixture <- reference_fixture(functions)
data <- list(C = fixture$C, L = fixture$K, initial = fixture$initial,
             spatial_load = fixture$spatial_load, kappa = fixture$kappa)
spectrum <- functions$isolated_discrete_spectrum(data$C, data$L)
for (alpha in c(.5, .9, 1.3)) {
  reference <- functions$isolated_exact_be_reference(data$C, spectrum, data$initial,
    data$spatial_load, fixture$time_seq, alpha, diffusion = 2, A = .9, omega = 1.7)
  for (m in 1:8) {
    roots <- functions$isolated_fractional_roots(m, alpha / 2, "computed")
    selected <- tab[tab$m == m & abs(tab$beta - alpha / 2) < 1e-14, ]
    check(identical(roots$pr_roots,
      selected$value[match(seq_len(m), selected$index + ifelse(selected$kind == "rc", 0, 1000))]),
      "CSV numerator index selection")
    rates <- functions$isolated_rational_rates(spectrum$values, alpha, m, data$kappa, 2,
                                               coefficient_source = "computed")
    z <- spectrum$values / data$kappa^2
    inverse <- roots$factor * vapply(z, function(x)
      prod(1 - roots$pr_roots * x) / prod(1 - roots$pl_roots * x), numeric(1))
    check(max(abs(rates / (2 * data$kappa^alpha / inverse) - 1)) < 1e-12,
          "Computed rates match independently evaluated root product")
    op <- functions$my.fractional.operators.frac(data$L, alpha / 2, data$C,
      data$kappa^2, m, reference$dt, diffusion = 2, fractional_roots = roots)
    rhs <- as.numeric(data$C %*% data$initial +
      reference$dt * data$spatial_load * reference$g_average[2L])
    expected <- spectrum$vectors %*% ((1 / (1 + reference$dt * rates)) *
      as.numeric(crossprod(spectrum$vectors, rhs)))
    check(max(abs(as.numeric(functions$my.solver.frac(op, rhs)) - expected)) < 2e-9,
          "Higher-order sparse update matches numerical spectral map")
    result <- functions$isolated_rational_error(data, spectrum, reference, m,
      coefficient_source = "computed", coefficient_table = tab)
    check(abs(result$error - result$product_error) < 2e-10, "Full computed sparse/product trajectory consistency")
    check(result$coefficient_source == "computed" && result$coefficient_beta_exact &&
      is.na(result$beta_exact_in_package_table), "Exact-beta source diagnostics")
  }
}

# Corrupt/incomplete inputs must fail explicitly, never silently select a source.
temporary <- tempfile(fileext = ".csv")
write.csv(rbind(tab, tab[1L, ]), temporary, row.names = FALSE)
check(fails(functions$isolated_read_computed_coefficients(temporary)), "Duplicate keys rejected")
write.csv(tab[!(tab$m == 8 & tab$beta == .25 & tab$kind == "rc" & tab$index == 8), ],
          temporary, row.names = FALSE)
check(fails(functions$isolated_fractional_roots(8, .25, "computed", temporary)), "Incomplete block rejected")
bad <- tab
bad$value[bad$m == 8 & bad$beta == .25 & bad$kind == "rc" & bad$index == 1] <- 1
write.csv(bad, temporary, row.names = FALSE)
check(fails(functions$isolated_fractional_roots(8, .25, "computed", temporary)), "Positive root rejected")
roots <- functions$isolated_fractional_roots(4, .25, "computed")
roots$pr_roots[1L] <- min(roots$pl_roots) - 1
check(fails(functions$isolated_validate_roots(roots, 4)), "Noninterlacing roots rejected")
invisible(file.copy(file, temporary, overwrite = TRUE))
changed <- functions$isolated_m_parameters("computed", temporary)
cat("\n", file = temporary, append = TRUE)
check(fails(functions$run_isolated_m_study(changed, workers = 1L, verbose = FALSE)),
      "Changed artifact rejected before dense setup")
unlink(temporary)
bad <- parameters
bad$alpha <- 1
check(fails(functions$run_isolated_m_study(bad, workers = 1L, verbose = FALSE)),
      "Unavailable beta rejected before dense setup")

# Exercise all eight orders through runner/diagnosis on a small real mesh.
small <- parameters
small$alpha <- .5
small$scenarios <- small$scenarios[1L, ]
small$scenarios$h <- .1
small$scenarios$target_step <- 1e-4
result <- functions$run_isolated_m_study(small, workers = 2L, verbose = FALSE)
check(nrow(result$errors) == 8L && identical(sort(result$errors$m), 1:8), "Runner executes computed orders 1--8")
check(result$status$numerical_checks_passed, "Computed runner numerical checks")
check(all(c("resolved_order_count", "resolved_only_slope") %in% names(result$slopes)),
      "High-order resolution-aware fit diagnostics")
check(any(grepl("Coefficient source: computed", functions$isolated_m_diagnosis(result), fixed = TRUE)),
      "Computed diagnosis source label")
cat(checks, "computed-source, higher-order sparse, provenance, and metadata checks passed.\n")
