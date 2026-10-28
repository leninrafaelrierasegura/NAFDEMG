# From the repository root, with the project's R:
# /usr/bin/Rscript --vanilla tests/benchmark_diffusion.R --steps=5000 --workers=1,4,8
# Uses full-profile h meshes and shortened time grids. Saves private
# results and resource measurements; does not execute full production runs.
project_lib <- file.path(getwd(), "renv/library/linux-ubuntu-noble/R-4.5/x86_64-pc-linux-gnu")
if (dir.exists(project_lib)) .libPaths(c(project_lib, .libPaths()))
suppressPackageStartupMessages(library(MetricGraph))
suppressPackageStartupMessages(library(Matrix))
helpers <- new.env(parent = globalenv())
for (name in c("fractional_evolution", "graphs_and_meshes", "exact_solutions",
               "convergence_experiments")) {
  for (expression in parse(file.path("fun", paste0(name, ".R")))) {
    if (is.call(expression) && identical(expression[[1]], as.name("<-")) &&
        is.call(expression[[3]]) && identical(expression[[3]][[1]], as.name("function")))
      eval(expression, helpers)
  }
}
args <- commandArgs(trailingOnly = TRUE)
argument <- function(name, default) {
  value <- args[startsWith(args, paste0("--", name, "="))]
  if (!length(value)) default else sub(paste0("^--", name, "="), "", value[1])
}
n_steps <- as.integer(argument("steps", "5000"))
worker_counts <- as.integer(strsplit(argument("workers", "1,4,8"), ",", fixed = TRUE)[[1]])
stopifnot(is.finite(n_steps), n_steps >= 2L, all(is.finite(worker_counts)),
          all(worker_counts >= 1L), worker_counts[1] == 1L)
T_final <- 1e-3
grid <- list(N = n_steps, time_step = T_final/n_steps,
             time_seq = seq(0, T_final, length.out = n_steps+1L))
parameters <- list(diffusion = 1, study = "h benchmark", steps = n_steps,
  h = 10^seq(log10(.0009), log10(.01), length.out = 4),
  alpha = c(.5,.9,1.3), m = 4L, workers = worker_counts,
  temporal_cache_mib = 512, production_time_grid = FALSE)
run_dir <- helpers$experiment_run("h", parameters,
  root = file.path(helpers$private_project_dir(), "diffusion_experiments", "benchmarks"))
cat("Benchmark directory:", run_dir, "\n")
baseline <- NULL
summary <- list()
for (workers in worker_counts) {
  gc(verbose = FALSE)
  cat("Benchmark requested workers:", workers, "\n")
  result <- helpers$run_tadpole_study(
    study = "h", axis_values = parameters$h, alpha_vector = parameters$alpha,
    time_grids = list(grid), fixed_m = 4L, anchor_h = min(parameters$h),
    anchor_grid = grid, anchor_m = 4L, diffusion = 1, kappa = 4,
    modes = c(5L,7L), N_finite = 4L, initial = c(10,0), forcing = c(0,10),
    run_refinement_checks = FALSE, workers = workers,
    worker_log = file.path(run_dir, paste0("workers-", workers, ".log")))
  if (is.null(baseline)) baseline <- result
  discrepancy <- max(abs(result$errors_projected-baseline$errors_projected))
  stopifnot(discrepancy < 1e-11)
  summary[[length(summary)+1L]] <- data.frame(
    workers_requested = workers, workers_used = result$execution$workers_used,
    elapsed_seconds = result$execution$elapsed_seconds,
    speedup = baseline$execution$elapsed_seconds/result$execution$elapsed_seconds,
    sum_worker_peak_rss_mib = result$execution$sum_worker_peak_rss_mib,
    temporal_cache_hits = sum(result$resources$temporal_cache_hit),
    max_error_difference = discrepancy)
  saveRDS(result, file.path(run_dir, paste0("workers-", workers, ".rds")))
  write.csv(result$resources, file.path(run_dir, paste0("resources-", workers, ".csv")),
            row.names = FALSE)
  write.csv(do.call(rbind, summary), file.path(run_dir, "summary.csv"), row.names = FALSE)
  print(tail(summary, 1)[[1]])
}

# Measure the actual memory of a full h time-data entry, without running a
# million-step solve. This cache can be reused across h or m on a worker.
large_grid <- list(N = 1000000L, time_step = 1e-9,
  time_seq = seq(0, 1e-3, length.out = 1000001L))
rates <- (16+c((3*pi/2)^2, (4*pi/2)^2))^.45
elapsed <- system.time(time_data <- helpers$prepare_modal_time_data(
  large_grid, rates, c(10,0), c(0,10)))[["elapsed"]]
cache_example <- list(steps = large_grid$N,
  object_mib = as.numeric(object.size(time_data))/1024^2, elapsed_seconds = elapsed)
saveRDS(cache_example, file.path(run_dir, "million-step-cache-memory.rds"))
print(cache_example)
cat("Benchmark finished; all worker counts matched the sequential errors.\n")
cat("Saved summary:", file.path(run_dir, "summary.csv"), "\n")
