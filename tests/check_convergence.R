# From the project root:
# /usr/bin/Rscript --vanilla tests/check_convergence.R [--save] [--plots]
# --save executes the notebooks' original fixed-name saves; --plots also exports figures.
# Frozen reference: diffusion-parameter commit 019e6c3, computed pilot,
# T=1, diffusion=1, alpha=(0.5,0.9,1.3), h/tau=1e-4:1e-2, m=1:4.
project_lib <- file.path(getwd(), "renv/library/linux-ubuntu-noble/R-4.5/x86_64-pc-linux-gnu")
if (dir.exists(project_lib)) .libPaths(c(project_lib, .libPaths()))
suppressPackageStartupMessages(library(Matrix))
helpers <- new.env(parent = globalenv())
for (name in c("fractional_evolution", "exact_solutions", "tadpole_convergence",
               "plotting", "metafunctions")) {
  for (ex in parse(file.path("fun", paste0(name, ".R"))))
    if (is.call(ex) && identical(ex[[1]], as.name("<-")) &&
        is.call(ex[[3]]) && identical(ex[[3]][[1]], as.name("function"))) eval(ex, helpers)
}
plots <- "--plots" %in% commandArgs(trailingOnly = TRUE)
save_outputs <- plots || "--save" %in% commandArgs(trailingOnly = TRUE)
if (plots) {
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(reshape2))
}
reference <- read.csv("tests/fixtures/convergence_reference.csv")
shared <- NULL
for (study in c("h", "tau", "m")) {
  env <- new.env(parent = helpers)
  lines <- readLines(paste0("conv_in_", study, ".Rmd"))
  for (start in which(grepl("^```\\{r", lines))) {
    end <- which(seq_along(lines) > start & lines == "```")[1]
    code <- lines[(start+1):(end-1)]
    if (any(grepl("^# Parameters:|^# Compute", code))) {
      for (ex in parse(text = code))
        if (save_outputs || !is.call(ex) || !identical(ex[[1]], as.name("save"))) eval(ex, env)
    } else if (save_outputs && any(grepl("^# Verify", code))) {
      eval(parse(text = code), env)
    } else if (plots && any(grepl("^observed_rates <-", code))) {
      eval(parse(text = code), env)
    }
  }
  if (is.null(shared)) shared <- env$convergence_parameters else
    stopifnot(identical(shared, env$convergence_parameters))
  axis <- switch(study, h = env$h_vector, tau = env$tau_vector, m = env$m_vector)
  stopifnot(length(env$alpha_vector) == 3L, length(axis) >= 4L,
            is.list(env$convergence_parameters),
            identical(dim(env$errors_projected),c(length(axis),3L)))
  for (j in seq_along(env$alpha_vector)) {
    truth <- reference[reference$study == study & reference$alpha == env$alpha_vector[j], ]
    stopifnot(length(axis) == nrow(truth),
              max(abs(axis/truth$axis-1)) < 1e-12,
              max(abs(env$errors_projected[,j]/truth$error-1)) < 1e-10)
  }
  predictor <- if (study == "m") sqrt(axis) else log(axis)
  rate <- apply(env$errors_projected, 2, function(error)
    unname(coef(lm(log(error) ~ predictor))[2]))
  if (study == "h") stopifnot(all(abs(rate-2) < .02))
  if (study == "tau") stopifnot(all(abs(rate-1) < .02))
  stopifnot(all(apply(env$errors_projected, 2, function(e)
    if (study == "m") all(diff(e)<0) else all(diff(e)>0))))
  h <- if (study == "h") min(axis) else env$h_fixed
  tau <- if (study == "tau") min(axis) else env$time_step_fixed
  m <- if (study == "m") max(axis) else env$m_fixed
  for (alpha in env$alpha_vector) {
    evaluate <- function(hh = h, tt = tau, mm = m, q = 4L, qs = 64L)
      helpers$tadpole_convergence_error(hh,tt,mm,alpha,shared,q,qs)
    error <- evaluate()
    changes <- c(tau = abs(evaluate(tt=tau/2)/error-1),
      h = abs(evaluate(hh=h/2)/error-1), quadrature = abs(evaluate(q=8L)/error-1),
      step_sum = abs(evaluate(qs=128L)/error-1),
      rational = if (study == "m") evaluate(mm=0L)/error else abs(evaluate(mm=0L)/error-1))
    required <- if (study == "m") names(changes) else setdiff(names(changes),study)
    stopifnot(all(changes[required] < .05))
    if (study == "m") stopifnot(changes["rational"] < .01)
  }
  cat("PASS:",study,"matches preserved results; accuracy checks pass; slopes",format(rate,digits=7),"\n")
}

# Independent sparse stepping with a nonidentity mass matrix and physical load.
# This verifies the general operator, not just the uniform-tadpole shortcut.
C <- Matrix(matrix(c(2,.2,.2,3),2), sparse=TRUE)
L <- Matrix(matrix(c(40,-2,-2,160),2), sparse=TRUE)
Ri <- solve(chol(as.matrix(C)))
ev <- eigen(t(Ri)%*%as.matrix(L)%*%Ri,symmetric=TRUE)
V <- Ri%*%ev$vectors
dt <- .002
for (m in c(1L,4L,8L)) for (alpha in c(.5,.9,1.3)) {
  initial <- c(1,2)
  forcing <- c(3,-4)
  op <- helpers$my.fractional.operators.frac(L,alpha/2,C,16,m,dt,
                                           diffusion=1,coefficient_source="computed")
  numerical <- initial
  for (k in 1:500) {
    average <- sin(pi*(k-.5)*dt)*sin(pi*dt/2)/(pi*dt/2)
    numerical <- as.numeric(helpers$my.solver.frac(op,C%*%numerical+dt*forcing*average))
  }
  rates <- helpers$rational_modal_rates(ev$values,alpha,m,1,4)
  expected <- as.numeric(V%*%helpers$uniform_modal_coefficients(500,dt,rates,
    as.numeric(t(V)%*%C%*%initial),as.numeric(t(V)%*%forcing)))
  stopifnot(max(abs(numerical-expected)) < 1e-9)
}
cat("PASS: sparse physical-forcing updates for m=1,4,8 at T=1, diffusion=1.\n")
if (save_outputs) cat("Saved the notebooks' original fixed-name outputs under data_files/.\n")
cat("All convergence checks passed.\n")
