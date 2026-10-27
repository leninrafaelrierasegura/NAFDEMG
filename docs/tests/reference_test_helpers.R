suppressPackageStartupMessages(library(Matrix))

load_reference_functions <- function() {
  environment <- new.env(parent = globalenv())
  for (name in c("fractional_evolution", "exact_solutions", "reference_experiments")) {
    extracted <- tempfile(fileext = ".R")
    knitr::purl(file.path("fun", paste0(name, ".Rmd")), output = extracted,
                documentation = 0L, quiet = TRUE)
    for (expression in parse(extracted)) {
      if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
          is.call(expression[[3L]]) && identical(expression[[3L]][[1L]], as.name("function")))
        eval(expression, envir = environment)
    }
    unlink(extracted)
  }
  environment
}

interval_fem <- function(nodes) {
  count <- length(nodes)
  C <- G <- matrix(0, count, count)
  for (i in seq_len(count - 1L)) {
    indices <- i:(i + 1L)
    h <- nodes[i + 1L] - nodes[i]
    C[indices, indices] <- C[indices, indices] + h / 6 * matrix(c(2, 1, 1, 2), 2)
    G[indices, indices] <- G[indices, indices] + 1 / h * matrix(c(1, -1, -1, 1), 2)
  }
  list(C = as(C, "sparseMatrix"), G = as(G, "sparseMatrix"))
}

reference_fixture <- function(functions, nodes = c(0, .13, .49, .7, 1),
                              reference_nodes = sort(c(seq(0, 1, length.out = 23), .19, .81)),
                              beta = .45, m = 4L, diffusion = 1,
                              time_seq = seq(0, .2, length.out = 6)) {
  fem <- interval_fem(nodes)
  ref <- interval_fem(reference_nodes)
  identity <- diag(length(nodes))
  Psi <- sapply(seq_along(nodes), function(j)
    approx(nodes, identity[, j], xout = reference_nodes)$y)
  Psi <- Matrix(Psi, sparse = TRUE)
  E <- cbind(cos(pi * reference_nodes), cos(2 * pi * reference_nodes))
  initial <- c(.8, -.3)
  forcing <- c(1.2, .7)
  kappa <- 2
  K <- kappa^2 * fem$C + fem$G
  spatial <- as.matrix(crossprod(Psi, ref$C %*% E))
  dt <- time_seq[2L] - time_seq[1L]
  rates <- diffusion * (kappa^2 + c(pi, 2 * pi)^2)^beta
  coefficients <- function(t) functions$exact_modal_coefficients(
    t, rates, initial, forcing, A = .9, omega = 1.7)
  list(C = fem$C, K = K, beta = beta, diffusion = diffusion, kappa = kappa,
       op = functions$my.fractional.operators.frac(K, beta, fem$C, kappa^2, m, dt,
                                                 diffusion = diffusion),
       dt = dt, time_seq = time_seq,
       initial = as.numeric(solve(fem$C, spatial %*% initial)),
       spatial_load = as.numeric(spatial %*% forcing),
       g_average = functions$average_sine_forcing(time_seq, A = .9, omega = 1.7),
       Psi = Psi, reference_C = ref$C, E = E,
       coefficients = coefficients, exact = function(t) E %*% coefficients(t),
       projection = functions$prepare_reference_error(Psi, ref$C, E))
}

fixture_error <- function(functions, fixture, reference_times, method, blocks) {
  functions$solve_fractional_evolution_error(
    fixture$op, fixture$dt, fixture$time_seq, fixture$initial,
    fixture$spatial_load, fixture$g_average, fixture$Psi, fixture$reference_C,
    reference_times, fixture$exact, n_blocks = blocks, error_method = method,
    coefficients_at_time = fixture$coefficients, error_projection = fixture$projection)
}

literal_reference_error <- function(functions, fixture, reference_times) {
  u <- fixture$initial
  states <- matrix(0, length(u), length(fixture$time_seq))
  states[, 1L] <- u
  for (k in seq_len(ncol(states) - 1L)) {
    rhs <- fixture$C %*% u + fixture$dt * fixture$spatial_load * fixture$g_average[k + 1L]
    u <- as.numeric(functions$my.solver.frac(fixture$op, rhs))
    states[, k + 1L] <- u
  }
  total <- 0
  for (j in 2:length(reference_times)) {
    k <- findInterval(reference_times[j], fixture$time_seq, left.open = TRUE)
    delta <- as.numeric(fixture$exact(reference_times[j]) - fixture$Psi %*% states[, k + 1L])
    total <- total + (reference_times[j] - reference_times[j - 1L]) *
      sum(delta * (fixture$reference_C %*% delta))
  }
  sqrt(total)
}
