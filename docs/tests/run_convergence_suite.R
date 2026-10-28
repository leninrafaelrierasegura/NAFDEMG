# Run the actual numerical chunks of all three notebooks, with persistent local
# output and a common settings fingerprint. No package installs or paper edits.
# /usr/bin/Rscript --vanilla tests/run_convergence_suite.R [--full] [--package|--computed]
project_lib <- file.path(getwd(), "renv/library/linux-ubuntu-noble/R-4.5/x86_64-pc-linux-gnu")
if (dir.exists(project_lib)) .libPaths(c(project_lib, .libPaths()))
stopifnot(getRversion() >= "4.5.0")
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(MetricGraph))
helpers <- new.env(parent = globalenv())
for (name in c("fractional_evolution", "graphs_and_meshes", "exact_solutions",
               "convergence_experiments", "convergence_settings")) {
  source_file <- file.path("fun", paste0(name, ".Rmd"))
  stopifnot(identical(readLines(file.path("fun", paste0(name, ".hash"))),
                      digest::digest(file = source_file, algo = "md5")))
  for (ex in parse(file.path("fun", paste0(name, ".R"))))
    if (is.call(ex) && identical(ex[[1]], as.name("<-")) &&
        is.call(ex[[3]]) && identical(ex[[3]][[1]], as.name("function"))) eval(ex, helpers)
}
profile <- if ("--full" %in% commandArgs(trailingOnly = TRUE)) "full" else "pilot"
source_flags <- intersect(commandArgs(trailingOnly = TRUE),c("--package","--computed"))
stopifnot(length(source_flags) <= 1L)
coefficient_source <- if (length(source_flags)) sub("^--","",source_flags) else
  formals(helpers$convergence_settings)$coefficient_source
tables <- list()
suite_id <- NULL
paths <- character()
for (study in c("h", "tau", "m")) {
  lines <- readLines(paste0("conv_in_", study, ".Rmd"))
  env <- new.env(parent = helpers)
  env$convergence_settings <- function(study)
    helpers$convergence_settings(study, profile = profile, coefficient_source = coefficient_source)
  for (start in which(grepl("^```\\{r", lines))) {
    end <- which(seq_along(lines) > start & lines == "```")[1]
    code <- lines[(start+1):(end-1)]
    if (any(grepl("^# Parameters:|^# Refuse to overwrite|^# Load only this run", code))) {
      eval(parse(text = code), env)
    }
  }
  id <- env$experiment_parameters$suite_id
  if (is.null(suite_id)) suite_id <- id else stopifnot(identical(suite_id, id))
  checks <- env$results$refinement_checks
  stopifnot(all(is.finite(env$errors_projected)), all(env$errors_projected > 0),
            all(checks$passed[checks$required]))
  transformed_x <- if (study == "m") sqrt(env$axis_values) else log(env$axis_values)
  observed <- apply(env$errors_projected, 2, function(error)
    unname(coef(lm(log(error) ~ transformed_x))[2]))
  expected <- switch(study, h = rep(2, length(observed)),
                     tau = rep(1, length(observed)), m = -2*pi*sqrt(1-env$alpha_vector/2))
  monotone <- apply(env$errors_projected, 2, function(error)
    if (study == "m") all(diff(error) < 0) else all(diff(error) > 0))
  rates <- data.frame(study = study, profile = profile, alpha = env$alpha_vector,
                       expected = expected, observed = observed, monotone = monotone)
  # The m theorem is an upper-bound benchmark, not a finite-order slope equality.
  stopifnot(all(monotone))
  if (study == "h") stopifnot(all(abs(observed-2) < .05))
  if (study == "tau") stopifnot(all(abs(observed-1) < .05))
  write.csv(rates, file.path(env$run_dir, "observed_rates.csv"), row.names = FALSE)
  error_table <- do.call(rbind, lapply(seq_along(env$alpha_vector), function(j)
    data.frame(axis = env$axis_values, alpha = env$alpha_vector[j],
               error = env$errors_projected[, j])))
  write.csv(error_table, file.path(env$run_dir, "errors.csv"), row.names = FALSE)
  png(file.path(env$run_dir, "calibration_curve.png"), width = 1100, height = 800, res = 140)
  if (study == "m") {
    matplot(sqrt(env$axis_values), log(env$errors_projected), type = "b", pch = 1:3,
            lty = 1, xaxt = "n", xlab = "m (sqrt scale)", ylab = "ln(space-time L2 error)",
            main = paste("Convergence in m;", profile, "; T =", env$T_final))
    axis(1, at = sqrt(env$axis_values), labels = env$axis_values)
  } else {
    matplot(env$axis_values, env$errors_projected, type = "b", log = "xy",
            pch = 1:3, lty = 1, xlab = study, ylab = "space-time L2 error",
            main = paste("Convergence in", study, ";", profile, "; T =", env$T_final))
  }
  legend(if (study == "m") "bottomleft" else "topleft",
         legend = paste("alpha =", env$alpha_vector), col = 1:3,
         pch = 1:3, lty = 1, bty = "n")
  dev.off()
  tables[[study]] <- rates
  paths[study] <- env$run_dir
}
cat("\nShared suite:", suite_id, "\n")
print(do.call(rbind, tables), row.names = FALSE)
cat("\nAll required accuracy checks passed. Results stay local:\n")
print(paths)
