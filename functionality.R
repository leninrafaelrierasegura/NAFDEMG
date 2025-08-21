## -----------------------------------------------------------------------------
# remotes::install_github("davidbolin/rspde", ref = "devel")
# remotes::install_github("davidbolin/metricgraph", ref = "devel")
library(rSPDE)
library(MetricGraph)
library(grateful)

library(ggplot2)
library(reshape2)
library(plotly)


## -----------------------------------------------------------------------------
# Function to compute the roots and factor for the rational approximation
my.get.roots <- function(m, # rational order, m = 1, 2, 3, or 4
                         beta # smoothness parameter, beta = alpha/2 with alpha between 0.5 and 2
                         ) {
  m1table <- rSPDE:::m1table
  m2table <- rSPDE:::m2table
  m3table <- rSPDE:::m3table
  m4table <- rSPDE:::m4table
  mt <- get(paste0("m", m, "table"))
  rb <- rep(0, m + 1)
  rc <- rep(0, m)
  if(m == 1) {
    rc = approx(mt$beta, mt[[paste0("rc")]], beta)$y
  } else {
    rc = sapply(1:m, function(i) {
      approx(mt$beta, mt[[paste0("rc.", i)]], beta)$y
    })
  }
  rb = sapply(1:(m+1), function(i) {
    approx(mt$beta, mt[[paste0("rb.", i)]], xout = beta)$y
  })
  factor = approx(mt$beta, mt$factor, xout = beta)$y
  return(list(pl_roots = rb, # roots \{r_{2j}\}_{j=1}^{m+1}
              pr_roots = rc, # roots \{r_{1i}\}_{i=1}^m
              factor = factor # this is c_m/b_{m+1}
              ))
}


## -----------------------------------------------------------------------------
# Function to compute polynomial coefficients from roots
poly.from.roots <- function(roots) {
  coef <- 1
  for (r in roots) {coef <- convolve(coef, c(1, -r), type = "open")}
  return(coef) # returned in increasing order like a+bx+cx^2+...
}


## -----------------------------------------------------------------------------
# Function to compute the parameters for the partial fraction decomposition
compute.partial.fraction.param <- function(factor, # c_m/b_{m+1}
                                           pr_roots, # roots \{r_{1i}\}_{i=1}^m
                                           pl_roots, # roots \{r_{2j}\}_{j=1}^{m+1}
                                           time_step, # \tau
                                           scaling # \kappa^{2\beta}
                                           ) {
  pr_coef <- c(0, poly.from.roots(pr_roots)) 
  pl_coef <- poly.from.roots(pl_roots) 
  factor_pr_coef <- pr_coef
  pr_plus_pl_coef <- factor_pr_coef + ((scaling * time_step)/factor) * pl_coef
  res <- gsignal::residue(factor_pr_coef, pr_plus_pl_coef)
  return(list(r = res$r, # residues \{a_k\}_{k=1}^{m+1}
              p = res$p, # poles \{p_k\}_{k=1}^{m+1}
              k = res$k # remainder r
              )) 
}


## -----------------------------------------------------------------------------
# Function to compute the fractional operator
my.fractional.operators.frac <- function(L, # Laplacian matrix
                                         beta, # smoothness parameter beta
                                         C, # mass matrix (not lumped)
                                         scale.factor, # scaling parameter = kappa^2
                                         m = 1, # rational order, m = 1, 2, 3, or 4
                                         time_step # time step = tau
                                         ) {
  I <- Matrix::Diagonal(dim(C)[1])
  L <- L / scale.factor 
  if(beta == 1){
    L <- L * scale.factor^beta
    return(list(C = C, # mass matrix
                L = L, # Laplacian matrix scaled
                m = m, # rational order
                beta = beta, # smoothness parameter
                LHS = C + time_step * L # left-hand side of the linear system
                ))
  } else {
    scaling <- scale.factor^beta
    roots <- my.get.roots(m, beta)
    poles_rs_k <- compute.partial.fraction.param(roots$factor, roots$pr_roots, roots$pl_roots, time_step, scaling)

    partial_fraction_terms <- list()
    for (i in 1:(m+1)) {
      # Here is where the terms in the sum in eq 12 are computed
      partial_fraction_terms[[i]] <- (L - poles_rs_k$p[i] * C)/poles_rs_k$r[i]
      }
    return(list(C = C, # mass matrix
                L = L, # Laplacian matrix scaled
                m = m, # rational order
                beta = beta, # smoothness parameter
                partial_fraction_terms = partial_fraction_terms # partial fraction terms
                ))
  }
}


## -----------------------------------------------------------------------------
# Function to solve the iteration
my.solver.frac <- function(obj, # object returned by my.fractional.operators.frac()
                           v # vector to be solved for
                           ){
  beta <- obj$beta
  m <- obj$m
  if (beta == 1){
    return(solve(obj$LHS, v) # solve the linear system directly for beta = 1
           )
  } else {
    partial_fraction_terms <- obj$partial_fraction_terms
    output <- v*0
    for (i in 1:(m+1)) {output <- output + solve(partial_fraction_terms[[i]], v)}
    return(output # solve the linear system using the partial fraction decomposition
           )
  }
}


## -----------------------------------------------------------------------------
solve_fractional_evolution <- function(my_op_frac, time_step, time_seq, val_at_0, RHST) {
  CC <- my_op_frac$C
  SOL <- matrix(NA, nrow = nrow(CC), ncol = length(time_seq))
  SOL[, 1] <- val_at_0
  for (k in 1:(length(time_seq) - 1)) {
    rhs <- CC %*% SOL[, k] + time_step * RHST[, k + 1]
    SOL[, k + 1] <- as.matrix(my.solver.frac(my_op_frac, rhs))
  }
  return(SOL)
}


## -----------------------------------------------------------------------------
# Function to build a tadpole graph and create a mesh
gets.graph.tadpole <- function(h){
  edge1 <- rbind(c(0,0),c(1,0))
  theta <- seq(from=-pi,to=pi,length.out = 10000)
  edge2 <- cbind(1+1/pi+cos(theta)/pi,sin(theta)/pi)
  edges <- list(edge1, edge2)
  graph <- metric_graph$new(edges = edges, verbose = 0)
  graph$set_manual_edge_lengths(edge_lengths = c(1,2))
  graph$build_mesh(h = h)
  return(graph)
}


## -----------------------------------------------------------------------------
# Function to compute the eigenfunctions of the tadpole graph
tadpole.eig <- function(k,graph){
x1 <- c(0,graph$get_edge_lengths()[1]*graph$mesh$PtE[graph$mesh$PtE[,1]==1,2]) 
x2 <- c(0,graph$get_edge_lengths()[2]*graph$mesh$PtE[graph$mesh$PtE[,1]==2,2]) 

if(k==0){ 
  f.e1 <- rep(1,length(x1)) 
  f.e2 <- rep(1,length(x2)) 
  f1 = c(f.e1[1],f.e2[1],f.e1[-1], f.e2[-1]) 
  f = list(phi=f1/sqrt(3)) 
  
} else {
  f.e1 <- -2*sin(pi*k*1/2)*cos(pi*k*x1/2) 
  f.e2 <- sin(pi*k*x2/2)                  
  
  f1 = c(f.e1[1],f.e2[1],f.e1[-1], f.e2[-1]) 
  
  if((k %% 2)==1){ 
    f = list(phi=f1/sqrt(3)) 
  } else { 
    f.e1 <- (-1)^{k/2}*cos(pi*k*x1/2)
    f.e2 <- cos(pi*k*x2/2)
    f2 = c(f.e1[1],f.e2[1],f.e1[-1],f.e2[-1]) 
    f <- list(phi=f1,psi=f2/sqrt(3/2))
  }
}
return(f)
}


## -----------------------------------------------------------------------------
# Function to compute the eigenpairs of the tadpole graph
gets.eigen.params <- function(N_finite = 4, kappa = 1, alpha = 0.5, graph){
  EIGENVAL <- NULL
  EIGENVAL_ALPHA <- NULL
  EIGENVAL_MINUS_ALPHA <- NULL
  EIGENFUN <- NULL
  INDEX <- NULL
  for (j in 0:N_finite) {
    lambda_j <- kappa^2 + (j*pi/2)^2
    lambda_j_alpha_half <- lambda_j^(alpha/2)
    lambda_j_minus_alpha_half <- lambda_j^(-alpha/2)
    e_j <- tadpole.eig(j,graph)$phi
    EIGENVAL <- c(EIGENVAL, lambda_j)
    EIGENVAL_ALPHA <- c(EIGENVAL_ALPHA, lambda_j_alpha_half)  
    EIGENVAL_MINUS_ALPHA <- c(EIGENVAL_MINUS_ALPHA, lambda_j_minus_alpha_half)
    EIGENFUN <- cbind(EIGENFUN, e_j)
    INDEX <- c(INDEX, j)
    if (j>0 && (j %% 2 == 0)) {
      lambda_j <- kappa^2 + (j*pi/2)^2
      lambda_j_alpha_half <- lambda_j^(alpha/2)
      lambda_j_minus_alpha_half <- lambda_j^(-alpha/2)
      e_j <- tadpole.eig(j,graph)$psi
      EIGENVAL <- c(EIGENVAL, lambda_j)
      EIGENVAL_ALPHA <- c(EIGENVAL_ALPHA, lambda_j_alpha_half)    
      EIGENVAL_MINUS_ALPHA <- c(EIGENVAL_MINUS_ALPHA, lambda_j_minus_alpha_half)
      EIGENFUN <- cbind(EIGENFUN, e_j)
      INDEX <- c(INDEX, j+0.1)
      }
    }
  return(list(EIGENVAL = EIGENVAL,
              EIGENVAL_ALPHA = EIGENVAL_ALPHA, 
              EIGENVAL_MINUS_ALPHA = EIGENVAL_MINUS_ALPHA,
              EIGENFUN = EIGENFUN,
              INDEX = INDEX))
}


## -----------------------------------------------------------------------------
# Function to construct a piecewise constant projection of approximated values
construct_piecewise_projection <- function(projected_U_approx, time_seq, overkill_time_seq) {
  projected_U_piecewise <- matrix(NA, nrow = nrow(projected_U_approx), ncol = length(overkill_time_seq))
  
  # Assign value at t = 0
  projected_U_piecewise[, which(overkill_time_seq == 0)] <- projected_U_approx[, 1]
  
  # Assign values for intervals (t_{k-1}, t_k]
  for (k in 2:length(time_seq)) {
    idxs <- which(overkill_time_seq > time_seq[k - 1] & overkill_time_seq <= time_seq[k])
    projected_U_piecewise[, idxs] <- projected_U_approx[, k]
  }
  
  return(projected_U_piecewise)
}


## -----------------------------------------------------------------------------
loglog_line_equation <- function(x1, y1, slope) {
  b <- log10(y1 / (x1 ^ slope))
  
  function(x) {
    (x ^ slope) * (10 ^ b)
  }
}
exp_line_equation <- function(x1, y1, slope) {
  lnC <- log(y1) - slope * x1
  
  function(x) {
    exp(lnC + slope * x)
  }
}
compute_guiding_lines <- function(x_axis_vector, errors, theoretical_rates, line_equation_fun) {
  guiding_lines <- matrix(NA, nrow = length(x_axis_vector), ncol = length(theoretical_rates))
  
  for (j in seq_along(theoretical_rates)) {
    guiding_lines_aux <- matrix(NA, nrow = length(x_axis_vector), ncol = length(x_axis_vector))
    
    for (k in seq_along(x_axis_vector)) {
      point_x1 <- x_axis_vector[k]
      point_y1 <- errors[k, j]
      slope <- theoretical_rates[j]
      
      line <- line_equation_fun(x1 = point_x1, y1 = point_y1, slope = slope)
      guiding_lines_aux[, k] <- line(x_axis_vector)
    }
    
    guiding_lines[, j] <- rowMeans(guiding_lines_aux)
  }
  
  return(guiding_lines)
}


## -----------------------------------------------------------------------------
# Functions to compute the exact solution to the fractional diffusion equation
g_linear <- function(r, A, lambda_j_alpha_half) {
  return(A * exp(-lambda_j_alpha_half * r))
  }
G_linear <- function(t, A) {
  return(A * t)
  }
g_exp <- function(r, A, mu) {
  return(A * exp(mu * r))
  }
G_exp <- function(t, A, lambda_j_alpha_half, mu) {
  exponent <- lambda_j_alpha_half + mu
  return(A * (exp(exponent * t) - 1) / exponent)
  }
g_poly <- function(r, A, n) {
  return(A * r^n)
}
G_poly <- function(t, A, lambda_j_alpha_half, n) {
  t <- as.vector(t)
  k_vals <- 0:n
  sum_term <- sapply(t, function(tt) {
    sum(((-lambda_j_alpha_half * tt)^k_vals) / factorial(k_vals))
  })
  coeff <- ((-1)^(n + 1)) * factorial(n) / (lambda_j_alpha_half^(n + 1))
  return(A * coeff * (1 - exp(lambda_j_alpha_half * t) * sum_term))
}
g_sin <- function(r, A, omega) {
  return(A * sin(omega * r))
}
G_sin <- function(t, A, lambda_j_alpha_half, omega) {
  denom <- lambda_j_alpha_half^2 + omega^2
  numerator <- exp(lambda_j_alpha_half * t) * (lambda_j_alpha_half * sin(omega * t) - omega * cos(omega * t)) + omega
  return(A * numerator / denom)
}
g_cos <- function(r, A, theta) {
  return(A * cos(theta * r)) 
}
G_cos <- function(t, A, lambda_j_alpha_half, theta) {
  denom <- lambda_j_alpha_half^2 + theta^2
  numerator <- exp(lambda_j_alpha_half * t) * (lambda_j_alpha_half * cos(theta * t) + theta * sin(theta * t)) - lambda_j_alpha_half
  return(A * numerator / denom)
}


## -----------------------------------------------------------------------------
reversecolumns <- function(mat) {
  return(mat[, rev(seq_len(ncol(mat)))])
}


## -----------------------------------------------------------------------------
# helper: measure change relative to the size of the previous iterate 
change_comparer <- function(X_new, X_old, time_step, time_seq, weights, relative = TRUE) {
  num <- sqrt(as.double(t(weights) %*% ((X_new - X_old)^2) %*% rep(time_step, length(time_seq))))
  if (!relative) {
    return(num)
    }
  den <- sqrt(as.double(t(weights) %*% (X_new^2) %*% rep(time_step, length(time_seq))))
  if (den < .Machine$double.eps) {
    return(ifelse(num < .Machine$double.eps, 0, num))
  } else {
    return(num / den)
  }
}


## -----------------------------------------------------------------------------
# Coupled solver with multi-criteria convergence
solve_coupled_system_multi_tol <- function(
  my_op_frac,           # operator
  time_step,            # tau
  time_seq,             # vector of times 
  u_0,                  # initial state U^0 
  F_proj,               # matrix of F 
  Z_ini,
  V_d,                  # matrix of 
  u_d,
  Psi,                  # Psi matrix
  R,                    # R matrix
  A, B,                 # lower/upper bounds (vector or matrix broadcastable to time grid)
  mu,                   # positive scalar
  weights,
  tol = 1e-8,           # scalar or named list: list(Z=..., U=..., P=...)
  maxit = 200,
  verbose = FALSE,
  true_sol
) {

  if (is.numeric(tol) && length(tol) == 1) {
    tol_list <- list(Z = tol, U = tol, P = tol)
  } else if (is.list(tol)) {
    tol_list <- modifyList(list(Z = 1e-8, U = 1e-8, P = 1e-8), tol)
  } else stop("tol must be scalar or list(Z=...,U=...,P=...)")

  it <- 0
  converged <- FALSE
  
  rel_history <- data.frame(iter = integer(0), variable = character(0), value = numeric(0))
  abs_history <- data.frame(iter = integer(0), variable = character(0), value = numeric(0))
  min_history <- data.frame(iter = integer(0), variable = character(0), value = numeric(0))

  Z_list <- list()
  U_list <- list()
  P_list <- list()
  
  z_prev <- Z_ini
  Z_mat <- R %*% Psi %*% z_prev
  U_prev <- F_proj*0
  P_prev <- F_proj*0

  repeat {
    it <- it + 1

    U_mat <- solve_fractional_evolution(my_op_frac, time_step, time_seq, val_at_0 = u_0, RHST = F_proj + Z_mat)
    V_mat <- reversecolumns(R %*% Psi %*% U_mat)
    Q_mat <- solve_fractional_evolution(my_op_frac, time_step, time_seq, val_at_0 = u_0 * 0, RHST = V_mat - V_d)
    P_mat <- reversecolumns(Q_mat)
    z_new <- pmax(A, pmin(B, - P_mat / mu))
    Z_mat <- R %*% Psi %*% z_new
    
    # relative changes
    rel_changes_Z <- change_comparer(z_new, z_prev, time_step, time_seq, weights, relative = TRUE)  
    rel_changes_U <- change_comparer(U_mat, U_prev, time_step, time_seq, weights, relative = TRUE)
    rel_changes_P <- change_comparer(P_mat, P_prev, time_step, time_seq, weights, relative = TRUE)
    abs_changes_Z <- change_comparer(z_new, true_sol$z_bar, time_step, time_seq, weights, relative = FALSE)
    abs_changes_U <- change_comparer(U_mat, true_sol$u_bar, time_step, time_seq, weights, relative = FALSE)
    abs_changes_P <- change_comparer(P_mat, true_sol$p_bar, time_step, time_seq, weights, relative = FALSE)
    min_change <- 0.5 * (as.double(t(weights) %*% ((U_mat - u_d)^2 + mu * z_new^2) %*% rep(time_step, length(time_seq))))

    rel_history <- rbind(rel_history,
      data.frame(iter = it, variable = "Z", value = rel_changes_Z),
      data.frame(iter = it, variable = "U", value = rel_changes_U),
      data.frame(iter = it, variable = "P", value = rel_changes_P))
    abs_history <- rbind(abs_history,
      data.frame(iter = it, variable = "Z", value = abs_changes_Z),
      data.frame(iter = it, variable = "U", value = abs_changes_U),
      data.frame(iter = it, variable = "P", value = abs_changes_P))
    min_history <- rbind(min_history,
      data.frame(iter = it, variable = "min", value = min_change))
    
    if (verbose) {message(sprintf("iter %3d: rel(Z) = %.3e, rel(U) = %.3e, rel(P) = %.3e", it, rel_changes_Z, rel_changes_U, rel_changes_P))}

    # update stored previous iterates
    z_prev <- z_new
    U_prev <- U_mat
    P_prev <- P_mat
    
    Z_list[[paste0("iteration ", it)]] <- z_new
    U_list[[paste0("iteration ",it)]] <- U_mat
    P_list[[paste0("iteration ",it)]] <- P_mat

    # convergence check: require all rel_changes <= respective tol
    cond_Z <- rel_changes_Z <= tol_list$Z
    cond_U <- rel_changes_U <= tol_list$U
    cond_P <- rel_changes_P <= tol_list$P

    if ((cond_Z && cond_U && cond_P) || it >= maxit) {
      converged <- (cond_Z && cond_U && cond_P)
      break
    }
  }

  if (verbose && !converged) {
    message(sprintf(
      "Stopped at maxit=%d; rel_changes: Z = %.3e (tol %.3e), U = %.3e (tol %.3e), P = %.3e (tol %.3e)",
      it, rel_changes_Z, tol_list$Z, rel_changes_U, tol_list$U, rel_changes_P, tol_list$P
    ))
  }

  return(list(U = U_mat,  # solution U
              Z = z_new,  # solution z
              P = P_mat, # solution P
              iterations = it,
              converged = converged,
              tol_list = tol_list,
              rel_history = rel_history,
              abs_history = abs_history,
              min_history = min_history,
              Z_list = Z_list,
              U_list = U_list,
              P_list = P_list))
}


## -----------------------------------------------------------------------------
plot_convergence_history <- function(history_df, tol_list = NULL, type = "relative") {
  if (type == "relative"){
    text_title <- "|X_{iter} - X_{iter-1}| / |X_{iter}|"
  } else if (type == "absolute") {
    text_title <- "|X_{exact} - X_{iter}|"
  } else if (type == "minimum") {
    text_title <- "J(U_{iter},z_{iter})"
  }

  p <- ggplot(history_df, aes(x = iter, y = value, color = variable)) +
    geom_line() +
    geom_point(size = 1.5) +
    scale_y_log10() +
    labs(
      title = text_title,
      x = "Iteration",
      y = "Error",
      color = "Quantity"
    ) +
    theme_minimal()
  
  # Add tolerance lines if provided
  if (!is.null(tol_list)) {
    tol_df <- data.frame(
      variable = names(tol_list),
      tol = unlist(tol_list)
    )
    p <- p + geom_hline(
      data = tol_df,
      aes(yintercept = tol, color = variable),
      linetype = "dashed"
    )
  }
  
  return(plotly::ggplotly(p))
}


## -----------------------------------------------------------------------------
# Function to order the vertices for plotting
plotting.order <- function(v, graph){
  edge_number <- graph$mesh$VtE[, 1]
  pos <- sum(edge_number == 1)+1
  return(c(v[1], v[3:pos], v[2], v[(pos+1):length(v)], v[2]))
}


## -----------------------------------------------------------------------------
# Function to set the scene for 3D plots
global.scene.setter <- function(x_range, y_range, z_range, z_aspectratio = 4) {
  
  return(list(xaxis = list(title = "x", range = x_range),
              yaxis = list(title = "y", range = y_range),
              zaxis = list(title = "z", range = z_range),
              aspectratio = list(x = 2*(1+2/pi), 
                                 y = 2*(2/pi), 
                                 z = z_aspectratio*(2/pi)),
              camera = list(eye = list(x = (1+2/pi)/2, 
                                       y = 4, 
                                       z = 2),
                            center = list(x = (1+2/pi)/2, 
                                          y = 0, 
                                          z = 0))))
}


## -----------------------------------------------------------------------------
# Function to plot in 3D
graph.plotter.3d.old <- function(graph, time_seq, frame_val_to_display, ...) {
  U_list <- list(...)
  U_names <- sapply(substitute(list(...))[-1], deparse)

  # Spatial coordinates
  x <- plotting.order(graph$mesh$V[, 1], graph)
  y <- plotting.order(graph$mesh$V[, 2], graph)
  weights <- graph$mesh$weights

  # Apply plotting.order to each U
  U_list <- lapply(U_list, function(U) apply(U, 2, plotting.order, graph = graph))
  n_vars <- length(U_list)

  # Create plot_data frame with time and position replicated
  n_time <- ncol(U_list[[1]])
  base_data <- data.frame(
    x = rep(x, times = n_time),
    y = rep(y, times = n_time),
    the_graph = 0,
    frame = rep(time_seq, each = length(x))
  )

  # Add U columns to plot_data
  for (i in seq_along(U_list)) {
    base_data[[paste0("u", i)]] <- as.vector(U_list[[i]])
  }

  plot_data <- base_data

  # Generate vertical lines
  vertical_lines_list <- lapply(seq_along(U_list), function(i) {
    do.call(rbind, lapply(time_seq, function(t) {
      idx <- which(plot_data$frame == t)
      z_vals <- plot_data[[paste0("u", i)]][idx]
      data.frame(
        x = rep(plot_data$x[idx], each = 3),
        y = rep(plot_data$y[idx], each = 3),
        z = as.vector(t(cbind(0, z_vals, NA))),
        frame = rep(t, each = length(idx) * 3)
      )
    }))
  })

  # Set axis ranges
  z_range <- range(unlist(U_list))
  x_range <- range(x)
  y_range <- range(y)

  # Create plot
  p <- plot_ly(plot_data, frame = ~frame) %>%
    add_trace(x = ~x, y = ~y, z = ~the_graph, type = "scatter3d", mode = "lines",
              name = "", showlegend = FALSE,
              line = list(color = "black", width = 3))

  # Add traces for each variable
  colors <- rev(viridisLite::viridis(n_vars)) #RColorBrewer::brewer.pal(min(n_vars, 8), "Set1")
  for (i in seq_along(U_list)) {
    p <- add_trace(p,
      x = ~x, y = ~y, z = as.formula(paste0("~u", i)),
      type = "scatter3d", mode = "lines", name = U_names[i],
      line = list(color = colors[i], width = 3))
  }

  # Add vertical lines
  for (i in seq_along(vertical_lines_list)) {
    p <- add_trace(p,
      data = vertical_lines_list[[i]],
      x = ~x, y = ~y, z = ~z, frame = ~frame,
      type = "scatter3d", mode = "lines",
      line = list(color = "gray", width = 0.5),
      name = "Vertical lines",
      showlegend = FALSE)
  }
  frame_name <- deparse(substitute(frame_val_to_display))
  # Layout and animation controls
  p <- p %>%
    layout(
      scene = global.scene.setter(x_range, y_range, z_range),
      updatemenus = list(list(type = "buttons", showactive = FALSE,
                              buttons = list(
                                list(label = "Play", method = "animate",
                                     args = list(NULL, list(frame = list(duration = 2000 / length(time_seq), redraw = TRUE), fromcurrent = TRUE))),
                                list(label = "Pause", method = "animate",
                                     args = list(NULL, list(mode = "immediate", frame = list(duration = 0), redraw = FALSE)))
                              )
      )),
      title = paste0(frame_name,": ", formatC(frame_val_to_display[1], format = "f", digits = 4))
    ) %>%
    plotly_build()

  for (i in seq_along(p$x$frames)) {
    p$x$frames[[i]]$layout <- list(title = paste0(frame_name,": ", formatC(frame_val_to_display[i], format = "f", digits = 4)))
  }

  return(p)
}


## -----------------------------------------------------------------------------
graph.plotter.3d <- function(graph, time_seq, frame_val_to_display, U_list) {
  U_names <- names(U_list) 
  # Spatial coordinates
  x <- plotting.order(graph$mesh$V[, 1], graph)
  y <- plotting.order(graph$mesh$V[, 2], graph)
  weights <- graph$mesh$weights

  # Apply plotting.order to each U
  U_list <- lapply(U_list, function(U) apply(U, 2, plotting.order, graph = graph))
  n_vars <- length(U_list)
  
  # Create plot_data frame with time and position replicated
  n_time <- ncol(U_list[[1]])
  base_data <- data.frame(
    x = rep(x, times = n_time),
    y = rep(y, times = n_time),
    the_graph = 0,
    frame = rep(time_seq, each = length(x))
  )

  # Add U columns to plot_data
  for (i in seq_along(U_list)) {
    base_data[[paste0("u", i)]] <- as.vector(U_list[[i]])
  }

  plot_data <- base_data

  # Generate vertical lines
  vertical_lines_list <- lapply(seq_along(U_list), function(i) {
    do.call(rbind, lapply(time_seq, function(t) {
      idx <- which(plot_data$frame == t)
      z_vals <- plot_data[[paste0("u", i)]][idx]
      data.frame(
        x = rep(plot_data$x[idx], each = 3),
        y = rep(plot_data$y[idx], each = 3),
        z = as.vector(t(cbind(0, z_vals, NA))),
        frame = rep(t, each = length(idx) * 3)
      )
    }))
  })

  # Set axis ranges
  z_range <- range(unlist(U_list))
  x_range <- range(x)
  y_range <- range(y)

  # Create plot
  p <- plot_ly(plot_data, frame = ~frame) %>%
    add_trace(x = ~x, y = ~y, z = ~the_graph, type = "scatter3d", mode = "lines",
              name = "", showlegend = FALSE,
              line = list(color = "black", width = 3))

  if (n_vars == 2) {
    colors <- RColorBrewer::brewer.pal(min(n_vars, 8), "Set1") 
    } else {
    colors <- rev(viridisLite::viridis(n_vars)) 
  }
  # RColorBrewer::brewer.pal(min(n_vars, 8), "Set1")
  for (i in seq_along(U_list)) {
    p <- add_trace(p,
      x = ~x, y = ~y, z = as.formula(paste0("~u", i)),
      type = "scatter3d", mode = "lines", name = U_names[i],
      line = list(color = colors[i], width = 3))
  }

  # Add vertical lines
  for (i in seq_along(vertical_lines_list)) {
    p <- add_trace(p,
      data = vertical_lines_list[[i]],
      x = ~x, y = ~y, z = ~z, frame = ~frame,
      type = "scatter3d", mode = "lines",
      line = list(color = "gray", width = 0.5),
      name = "Vertical lines",
      showlegend = FALSE)
  }
  frame_name <- deparse(substitute(frame_val_to_display))
  # Layout and animation controls
  p <- p %>%
    layout(
      scene = global.scene.setter(x_range, y_range, z_range),
      updatemenus = list(list(type = "buttons", showactive = FALSE,
                              buttons = list(
                                list(label = "Play", method = "animate",
                                     args = list(NULL, list(frame = list(duration = 2000 / length(time_seq), redraw = TRUE), fromcurrent = TRUE))),
                                list(label = "Pause", method = "animate",
                                     args = list(NULL, list(mode = "immediate", frame = list(duration = 0), redraw = FALSE)))
                              )
      )),
      title = paste0(frame_name,": ", formatC(frame_val_to_display[1], format = "f", digits = 4))
    ) %>%
    plotly_build()

  for (i in seq_along(p$x$frames)) {
    p$x$frames[[i]]$layout <- list(title = paste0(frame_name,": ", formatC(frame_val_to_display[i], format = "f", digits = 4)))
  }

  return(p)
}


## -----------------------------------------------------------------------------
# Function to plot the error at each time step
error.at.each.time.plotter <- function(graph, U_true, U_approx, time_seq, time_step) {
  weights <- graph$mesh$weights
  error_at_each_time <- t(weights) %*% (U_true - U_approx)^2
  error <- sqrt(as.double(t(weights) %*% (U_true - U_approx)^2 %*% rep(time_step, ncol(U_true))))
  p <- plot_ly() %>% 
  add_trace(
  x = ~time_seq, y = ~error_at_each_time, type = 'scatter', mode = 'lines+markers',
  line = list(color = 'blue', width = 2),
  marker = list(color = 'blue', size = 4),
  name = "",
  showlegend = TRUE
) %>% 
  layout(
  title = paste0("Error at Each Time Step (Total error = ", formatC(error, format = "f", digits = 9), ")"),
  xaxis = list(title = "t"),
  yaxis = list(title = "Error"),
  legend = list(x = 0.1, y = 0.9)
)
  return(p)
}


## -----------------------------------------------------------------------------
# Function to plot the 3D comparison of U_true and U_approx
graph.plotter.3d.comparer <- function(graph, U_true, U_approx, time_seq) {
  x <- graph$mesh$V[, 1]; y <- graph$mesh$V[, 2]
  x <- plotting.order(x, graph); y <- plotting.order(y, graph)

  U_true <- apply(U_true, 2, plotting.order, graph = graph)
  U_approx <- apply(U_approx, 2, plotting.order, graph = graph)
  n_times <- length(time_seq)
  
  x_range <- range(x); y_range <- range(y); z_range <- range(c(U_true, U_approx))
  
  # Normalize time_seq
  time_normalized <- (time_seq - min(time_seq)) / (max(time_seq) - min(time_seq))
  blues <- colorRampPalette(c("lightblue", "blue"))(n_times)
  reds <- colorRampPalette(c("mistyrose", "red"))(n_times)
  
  # Accurate colorscales
  colorscale_greens <- Map(function(t, col) list(t, col), time_normalized, blues)
  colorscale_reds <- Map(function(t, col) list(t, col), time_normalized, reds)
  
  p <- plot_ly()
  
  # Static black graph structure
  p <- p %>%
    add_trace(x = x, y = y, z = rep(0, length(x)),
              type = "scatter3d", mode = "lines",
              line = list(color = "black", width = 4),
              name = "Graph", showlegend = FALSE)
  
  # U_true traces (green)
  for (i in seq_len(n_times)) {
    z <- U_true[, i]
    p <- add_trace(
      p,
      type = "scatter3d",
      mode = "lines",
      x = x, y = y, z = z,
      line = list(color = blues[i], width = 4),
      showlegend = FALSE,
      scene = "scene"
    )
  }
  
  # U_approx traces (dashed red)
  for (i in seq_len(n_times)) {
    z <- U_approx[, i]
    p <- add_trace(
      p,
      type = "scatter3d",
      mode = "lines",
      x = x, y = y, z = z,
      line = list(color = reds[i], width = 4, dash = "dot"),
      showlegend = FALSE,
      scene = "scene"
    )
  }
  
  # Dummy green colorbar (True) – with ticks
  p <- add_trace(
    p,
    type = "heatmap",
    z = matrix(time_seq, nrow = 1),
    showscale = TRUE,
    colorscale = colorscale_greens,
    colorbar = list(
      title = list(font = list(size = 12, color = "black"), text = "Time", side = "top"),
      len = 0.9,
      thickness = 15,
      x = 1.02,
      xanchor = "left",
      y = 0.5,
      yanchor = "middle",
      tickvals = NULL,   # hide tick values
      ticktext = NULL,
      ticks = ""         # also hides tick marks
    ),
    x = matrix(time_seq, nrow = 1),
    y = matrix(1, nrow = 1),
    hoverinfo = "skip",
    opacity = 0
  )

# Dummy red colorbar (Approx) – no ticks
  p <- add_trace(
    p,
    type = "heatmap",
    z = matrix(time_seq, nrow = 1),
    showscale = TRUE,
    colorscale = colorscale_reds,
    colorbar = list(
      title = list(font = list(size = 12, color = "black"), text = ".", side = "top"),
      len = 0.9,
      thickness = 15,
      x = 1.05,
      xanchor = "left",
      y = 0.5,
      yanchor = "middle"
    ),
    x = matrix(time_seq, nrow = 1),
    y = matrix(1, nrow = 1),
    hoverinfo = "skip",
    opacity = 0
  )
  p <- p %>%
    add_trace(x = x, y = y, z = rep(0, length(x)),
              type = "scatter3d", mode = "lines",
              line = list(color = "black", width = 4),
              name = "Graph", showlegend = FALSE)
  p <- layout(p,
            scene = global.scene.setter(x_range, y_range, z_range),
            xaxis = list(visible = FALSE),
            yaxis = list(visible = FALSE),
            annotations = list(
  list(
    text = "Exact",
    x = 1.045,
    y = 0.5,
    xref = "paper",
    yref = "paper",
    showarrow = FALSE,
    font = list(size = 12, color = "black"),
    textangle = -90
  ),
  list(
    text = "Approx",
    x = 1.075,
    y = 0.5,
    xref = "paper",
    yref = "paper",
    showarrow = FALSE,
    font = list(size = 12, color = "black"),
    textangle = -90
  )
)

)

  
  return(p)
}


## -----------------------------------------------------------------------------
# Function to plot a single 3D line for 
graph.plotter.3d.single <- function(graph, U_true, time_seq) {
  x <- graph$mesh$V[, 1]; y <- graph$mesh$V[, 2]
  x <- plotting.order(x, graph); y <- plotting.order(y, graph)

  U_true <- apply(U_true, 2, plotting.order, graph = graph)
  n_times <- length(time_seq)
  
  x_range <- range(x); y_range <- range(y); z_range <- range(U_true)
  z_range[1] <- z_range[1] - 10^-6
  viridis_colors <- viridisLite::viridis(100)
  
  # Normalize time_seq
  time_normalized <- (time_seq - min(time_seq)) / (max(time_seq) - min(time_seq))
  #greens <- colorRampPalette(c("palegreen", "darkgreen"))(n_times)
  greens <- colorRampPalette(c(viridis_colors[1], viridis_colors[50],  viridis_colors[100]))(n_times)
  # Accurate colorscales
  colorscale_greens <- Map(function(t, col) list(t, col), time_normalized, greens)
  
  p <- plot_ly()
  
  # Add the 3D lines with fading green color
  for (i in seq_len(n_times)) {
    z <- U_true[, i]
    
    p <- add_trace(
      p,
      type = "scatter3d",
      mode = "lines",
      x = x,
      y = y,
      z = z,
      line = list(color = greens[i], width = 2),
      showlegend = FALSE,
      scene = "scene"
    )
  }
  p <- p %>%
    add_trace(x = x, y = y, z = rep(0, length(x)),
              type = "scatter3d", mode = "lines",
              line = list(color = "black", width = 5),
              name = "Graph", showlegend = FALSE)
  # Add dummy heatmap to show colorbar (not part of scene)
  p <- add_trace(
    p,
    type = "heatmap",
    z = matrix(time_seq, nrow = 1),
    showscale = TRUE,
    colorscale = colorscale_greens,
    colorbar = list(
    title = list(font = list(size = 12, color = "black"), text = "Time", side = "top"),
    len = 0.9,         # height (0 to 1)
    thickness = 15,     # width in pixels
    x = 1.02,           # shift it slightly right of the plot
    xanchor = "left",
    y = 0.5,
    yanchor = "middle"),
    x = matrix(time_seq, nrow = 1),
    y = matrix(1, nrow = 1),
    hoverinfo = "skip",
    opacity = 0
  )
  
  p <- layout(p,
              scene = global.scene.setter(x_range, y_range, z_range),
              xaxis = list(visible = FALSE),
              yaxis = list(visible = FALSE)
  )
  
  return(p)
}


## -----------------------------------------------------------------------------
# Function to plot the error convergence
error.convergence.plotter <- function(x_axis_vector, 
                                      alpha_vector, 
                                      errors, 
                                      theoretical_rates, 
                                      observed_rates,
                                      line_equation_fun,
                                      fig_title,
                                      x_axis_label,
                                      apply_sqrt = FALSE) {
  
  x_vec <- if (apply_sqrt) sqrt(x_axis_vector) else x_axis_vector
  
  guiding_lines <- compute_guiding_lines(x_axis_vector = x_vec, 
                                         errors = errors, 
                                         theoretical_rates = theoretical_rates, 
                                         line_equation_fun = line_equation_fun)
  
  default_colors <- scales::hue_pal()(length(alpha_vector))
  
  plot_lines <- lapply(1:ncol(guiding_lines), function(i) {
    geom_line(
      data = data.frame(x = x_vec, y = guiding_lines[, i]),
      aes(x = x, y = y),
      color = default_colors[i],
      linetype = "dashed",
      show.legend = FALSE
    )
  })
  
  df <- as.data.frame(cbind(x_vec, errors))
  colnames(df) <- c("x_axis_vector", alpha_vector)
  df_melted <- melt(df, id.vars = "x_axis_vector", variable.name = "column", value.name = "value")
  
  custom_labels <- paste0(formatC(alpha_vector, format = "f", digits = 2), 
                          " | ", 
                          formatC(theoretical_rates, format = "f", digits = 4), 
                          " | ", 
                          formatC(observed_rates, format = "f", digits = 4))
  
  df_melted$column <- factor(df_melted$column, levels = alpha_vector, labels = custom_labels)

  p <- ggplot() +
    geom_line(data = df_melted, aes(x = x_axis_vector, y = value, color = column)) +
    geom_point(data = df_melted, aes(x = x_axis_vector, y = value, color = column)) +
    plot_lines +
    labs(
      title = fig_title,
      x = x_axis_label,
      y = expression(Error),
      color = "          α  | theo  | obs"
    ) +
    (if (apply_sqrt) {
      scale_x_continuous(breaks = x_vec, labels = round(x_axis_vector, 4))
    } else {
      scale_x_log10(breaks = x_axis_vector, labels = round(x_axis_vector, 4))
    }) +
    (if (apply_sqrt) {
      scale_y_continuous(trans = "log", labels = scales::scientific_format())
    } else {
      scale_y_log10(labels = scales::scientific_format())
    }) +
    theme_minimal() +
    theme(text = element_text(family = "Palatino"),
          legend.position = "bottom",
          legend.direction = "vertical",
          plot.margin = margin(0, 0, 0, 0),
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))
  
  return(p)
}



## -----------------------------------------------------------------------------
graph.plotter.3d.static <- function(graph, z_list) {
  x <- plotting.order(graph$mesh$V[, 1], graph)
  y <- plotting.order(graph$mesh$V[, 2], graph)
  U_names <- names(z_list)
  n_vars <- length(z_list)
  z_list <- lapply(z_list, function(z) plotting.order(z, graph))

  # Axis ranges
  z_range <- range(unlist(z_list))
  x_range <- range(x)
  y_range <- range(y)

  if (n_vars == 2) {
    colors <- RColorBrewer::brewer.pal(min(n_vars, 8), "Set1") 
    } else {
    colors <- rev(viridisLite::viridis(n_vars)) 
  }
  p <- plot_ly()

  for (i in seq_along(z_list)) {
    z <- z_list[[i]]

    # Main 3D curve
    p <- add_trace(
      p,
      x = x, y = y, z = z,
      type = "scatter3d", mode = "lines",
      line = list(color = colors[i], width = 3),
      name = U_names[i], showlegend = TRUE
    )

    # Efficient vertical lines: one trace with breaks (NA)
    x_vert <- rep(x, each = 3)
    y_vert <- rep(y, each = 3)
    z_vert <- unlist(lapply(z, function(zj) c(0, zj, NA)))

    p <- add_trace(
      p,
      x = x_vert, y = y_vert, z = z_vert,
      type = "scatter3d", mode = "lines",
      line = list(color = "gray", width = 0.5),
      showlegend = FALSE
    )
  }
  p <- p %>% add_trace(x = x, y = y, z = x*0, type = "scatter3d", mode = "lines",
              line = list(color = "black", width = 3),
              name = "thegraph", showlegend = FALSE) %>%
    layout(scene = global.scene.setter(x_range, y_range, z_range))
  return(p)
}



## -----------------------------------------------------------------------------
graph.plotter.3d.two.meshes.time <- function(graph_finer, graph_coarser, 
                                             time_seq, frame_val_to_display,
                                             fs_finer = list(), fs_coarser = list()) {
  # Spatial coordinates (ordered for plotting)
  x_finer <- plotting.order(graph_finer$mesh$V[, 1], graph_finer)
  y_finer <- plotting.order(graph_finer$mesh$V[, 2], graph_finer)
  x_coarser <- plotting.order(graph_coarser$mesh$V[, 1], graph_coarser)
  y_coarser <- plotting.order(graph_coarser$mesh$V[, 2], graph_coarser)
  
  n_time <- if (length(fs_finer) > 0) ncol(fs_finer[[1]]) else ncol(fs_coarser[[1]])

  # Helper: make dataframe from one function
  make_df <- function(f_mat, graph, x, y, mesh_name) {
    z <- apply(f_mat, 2, plotting.order, graph = graph)
    data.frame(
      x = rep(x, times = n_time),
      y = rep(y, times = n_time),
      z = as.vector(z),
      frame = rep(time_seq, each = length(x)),
      mesh = mesh_name
    )
  }
  
  # Build data for finer functions
  data_finer_list <- lapply(names(fs_finer), function(nm) {
    make_df(fs_finer[[nm]], graph_finer, x_finer, y_finer, nm)
  })
  
  # Build data for coarser functions
  data_coarser_list <- lapply(names(fs_coarser), function(nm) {
    make_df(fs_coarser[[nm]], graph_coarser, x_coarser, y_coarser, nm)
  })
  
  # Combine
  all_data <- c(data_finer_list, data_coarser_list)
  
  # Baseline graph (on finer mesh for consistency)
  data_graph <- data.frame(
    x = rep(x_finer, times = n_time),
    y = rep(y_finer, times = n_time),
    z = 0,
    frame = rep(time_seq, each = length(x_finer)),
    mesh = "Graph"
  )
  
# --------- Vertical lines helper ----------
vertical_lines <- function(x, y, z, frame_vals, mesh_name) {
  do.call(rbind, lapply(seq_along(frame_vals), function(i) {
    idx <- ((i - 1) * length(x) + 1):(i * length(x))
    data.frame(
      x = rep(x, each = 3),
      y = rep(y, each = 3),
      z = as.vector(t(cbind(0, z[idx], NA))),
      frame = rep(frame_vals[i], each = length(x) * 3),
      mesh = mesh_name
    )
  }))
}

# --------- Compute vertical lines per mesh using max absolute value ---------
make_vertical_from_list <- function(data_list, x, y, mesh_name) {
  if (length(data_list) == 0) return(NULL)
  
  # Reshape each function's z back to matrix: (nodes × time)
  z_mats <- lapply(data_list, function(df) {
    matrix(df$z, nrow = length(x), ncol = length(time_seq))
  })
  
  # Stack into 3D array: (nodes × time × functions)
  arr <- array(unlist(z_mats), dim = c(length(x), length(time_seq), length(z_mats)))
  
  # For each node × time, select the entry with largest absolute value (keep sign)
  idx <- apply(arr, c(1, 2), function(v) which.max(abs(v)))
  z_signed_max <- mapply(function(i, j) arr[i, j, idx[i, j]],
                         rep(1:length(x), times = length(time_seq)),
                         rep(1:length(time_seq), each = length(x)))
  
  # Flatten back into long vector
  z_signed_max <- as.vector(z_signed_max)
  
  vertical_lines(x, y, z_signed_max, time_seq, mesh_name)
}


vertical_finer   <- make_vertical_from_list(data_finer_list,   x_finer,   y_finer,   "finer")
vertical_coarser <- make_vertical_from_list(data_coarser_list, x_coarser, y_coarser, "coarser")

  
  # Compute ranges
  all_z <- unlist(lapply(all_data, function(df) df$z))
  x_range <- range(c(x_finer, x_coarser))
  y_range <- range(c(y_finer, y_coarser))
  z_range <- range(all_z)
  
  # --------- Plotly object ----------
  p <- plot_ly(frame = ~frame)
  
  # Add traces for finer + coarser (looping automatically with names)
  for (df in all_data) {
    p <- p %>%
      add_trace(data = df,
                x = ~x, y = ~y, z = ~z,
                type = "scatter3d", mode = "lines",
                line = list(width = 3),
                name = unique(df$mesh))
  }
  
  # Add baseline
  p <- p %>%
    add_trace(data = data_graph,
              x = ~x, y = ~y, z = ~z,
              type = "scatter3d", mode = "lines",
              line = list(color = "black", width = 2),
              name = "Graph", showlegend = FALSE)
  
# Add verticals (one per mesh, envelope of all functions)
if (!is.null(vertical_finer)) {
  p <- p %>%
    add_trace(data = vertical_finer,
              x = ~x, y = ~y, z = ~z,
              type = "scatter3d", mode = "lines",
              line = list(color = "gray", width = 0.5),
              name = "Vertical finer", showlegend = FALSE)
}
if (!is.null(vertical_coarser)) {
  p <- p %>%
    add_trace(data = vertical_coarser,
              x = ~x, y = ~y, z = ~z,
              type = "scatter3d", mode = "lines",
              line = list(color = "gray", width = 0.5),
              name = "Vertical coarser", showlegend = FALSE)
}

  
  frame_name <- deparse(substitute(frame_val_to_display))
  
  p <- p %>%
    layout(
      scene = global.scene.setter(x_range, y_range, z_range),
      updatemenus = list(list(
        type = "buttons", showactive = FALSE,
        buttons = list(
          list(label = "Play", method = "animate",
               args = list(NULL, list(frame = list(duration = 2000 / length(time_seq), redraw = TRUE),
                                      fromcurrent = TRUE))),
          list(label = "Pause", method = "animate",
               args = list(NULL, list(mode = "immediate", 
                                      frame = list(duration = 0), redraw = FALSE)))
        )
      )),
      title = paste0(frame_name, ": ", formatC(frame_val_to_display[1], format = "f", digits = 4))
    ) %>%
    plotly_build()
  
  # Update frame titles
  for (i in seq_along(p$x$frames)) {
    p$x$frames[[i]]$layout <- list(
      title = paste0(frame_name, ": ", formatC(frame_val_to_display[i], format = "f", digits = 4))
    )
  }
  
  return(p)
}


