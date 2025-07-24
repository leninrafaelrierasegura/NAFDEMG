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
  
  C <- Matrix::Diagonal(dim(C)[1], rowSums(C)) 
  Ci <- Matrix::Diagonal(dim(C)[1], 1 / rowSums(C)) 
  I <- Matrix::Diagonal(dim(C)[1])
  L <- L / scale.factor 
  LCi <- L %*% Ci
  if(beta == 1){
    L <- L * scale.factor^beta
    return(list(Ci = Ci, # inverse of lumped mass matrix
                C = C, # lumped mass matrix
                LCi = LCi, # Laplacian matrix times inverse of lumped mass matrix
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
      # Here is where the terms in the sum in eq 11 are computed
      partial_fraction_terms[[i]] <- (LCi - poles_rs_k$p[i] * I)/poles_rs_k$r[i]
      }
    partial_fraction_terms[[m+2]] <- ifelse(is.null(poles_rs_k$k), 0, poles_rs_k$k) * I
    return(list(Ci = Ci, # inverse of lumped mass matrix
                C = C, # lumped mass matrix
                LCi = LCi, # Laplacian matrix times inverse of lumped mass matrix
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
  C <- obj$C
  Ci <- obj$Ci
  if (beta == 1){
    return(solve(obj$LHS, v) # solve the linear system directly for beta = 1
           )
  } else {
    partial_fraction_terms <- obj$partial_fraction_terms
    output <- partial_fraction_terms[[m+2]] %*% v
    for (i in 1:(m+1)) {output <- output + solve(partial_fraction_terms[[i]], v)}
    return(Ci %*% output # solve the linear system using the partial fraction decomposition
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
graph.plotter.3d <- function(graph, time_seq, frame_val_to_display, ...) {
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
  colors <- RColorBrewer::brewer.pal(min(n_vars, 8), "Set1")
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

