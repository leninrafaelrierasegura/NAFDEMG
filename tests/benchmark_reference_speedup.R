# Representative local benchmark; writes no experiment files or notifications.
# Run: Rscript --vanilla tests/benchmark_reference_speedup.R
source("tests/reference_test_helpers.R")
functions <- load_reference_functions()
fixture <- reference_fixture(functions, nodes = seq(0, 1, length.out = 51),
  reference_nodes = seq(0, 1, length.out = 2001),
  time_seq = seq(0, .2, length.out = 51), diffusion = 2)
times <- seq(0, .2, length.out = 10001)
methods <- c("reference", "compressed")
answers <- elapsed <- setNames(numeric(2), methods)
for (method in methods) {
  measurements <- replicate(3L, {
    gc()
    timing <- system.time(value <- fixture_error(functions, fixture, times, method, 40L))
    answers[method] <<- value
    timing[["elapsed"]]
  })
  elapsed[method] <- median(measurements)
}
stopifnot(abs(diff(answers)) <= 2e-11 * max(1, answers))
print(data.frame(method = methods, error = unname(answers), median_seconds = unname(elapsed)), row.names = FALSE)
cat("Measured end-to-end ratio:", elapsed["reference"] / elapsed["compressed"], "\n")
cat("Includes ordinary time stepping in both cases; projection setup is shared and excluded.\n")
