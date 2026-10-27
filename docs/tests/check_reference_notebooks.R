# Exercise each notebook's real setup/computation/load chunks at small sizes.
# Use the project environment: Rscript tests/check_reference_notebooks.R
suppressPackageStartupMessages(library(rSPDE))
suppressPackageStartupMessages(library(MetricGraph))
source("tests/reference_test_helpers.R")
functions <- load_reference_functions()
extracted <- tempfile(fileext = ".R")
knitr::purl("fun/graphs_and_meshes.Rmd", output = extracted, documentation = 0L, quiet = TRUE)
for (expression in parse(extracted)) {
  if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
      is.call(expression[[3L]]) && identical(expression[[3L]][[1L]], as.name("function")))
    eval(expression, envir = functions)
}
unlink(extracted)

chunks <- function(path) {
  lines <- readLines(path)
  starts <- which(grepl("^```\\{r", lines))
  ends <- which(lines == "```")
  lapply(starts, function(start) paste(lines[seq.int(start + 1L, min(ends[ends > start]) - 1L)], collapse = "\n"))
}

root <- tempfile("reference-notebook-smoke-")
old_root <- Sys.getenv("NAFDEMG_PRIVATE_DIR", unset = NA_character_)
Sys.setenv(NAFDEMG_PRIVATE_DIR = root)
answers <- list()
for (method in c("reference", "compressed")) for (study in c("h", "m", "tau")) {
  environment <- new.env(parent = functions)
  overrides <- list(T_final = .02, diffusion = 2, workers = if (method == "compressed") 2L else 1L,
    error_method = method, h_vector = c(.2, .1), alpha_vector = c(.5, .9),
    m_vector = 1:2, m_fixed = 2L, h_fixed = .1, time_step_target = .004,
    tau_targets = c(.004, .002), overkill_h = .04, overkill_time_step = .0007,
    error_verbose = FALSE, error_n_blocks = 3L)
  started <- FALSE
  for (chunk in chunks(paste0("conv_in_", study, ".Rmd"))) {
    if (grepl("^# Parameters", trimws(chunk))) started <- TRUE
    if (!started) next
    if (grepl("observed_rates <-", chunk, fixed = TRUE)) break
    for (expression in parse(text = chunk)) {
      if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
          is.symbol(expression[[2L]]) && as.character(expression[[2L]]) %in% names(overrides)) {
        assign(as.character(expression[[2L]]), overrides[[as.character(expression[[2L]])]], envir = environment)
      } else if (!(is.call(expression) && identical(expression[[1L]], as.name("print"))))
        eval(expression, envir = environment)
    }
  }
  stopifnot(identical(dim(environment$errors_projected), c(2L, 2L)),
            all(is.finite(environment$errors_projected)),
            environment$study_parameters$diffusion == 2,
            file.exists(environment$result_file))
  answers[[paste(study, method)]] <- environment$errors_projected
}
for (study in c("h", "m", "tau")) {
  stopifnot(max(abs(answers[[paste(study, "reference")]] -
                    answers[[paste(study, "compressed")]])) < 2e-11)
}
if (is.na(old_root)) Sys.unsetenv("NAFDEMG_PRIVATE_DIR") else Sys.setenv(NAFDEMG_PRIVATE_DIR = old_root)
unlink(root, recursive = TRUE)
cat("All three notebooks passed reference/compressed and serial/parallel smoke comparisons.\n")
