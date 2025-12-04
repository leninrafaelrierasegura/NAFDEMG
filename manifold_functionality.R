## ---------------------------------------------------------------------------------------
library(plotly)
library(akima)
library(orthopolynom)
library(Matrix)
library(pracma) 


## ----setup, include = FALSE-------------------------------------------------------------
# to install in terminal
# conda activate fenicsenv
# conda install -c conda-forge matplotlib plotly ipywidgets
library(reticulate)
# VERY IMPORTANT: tell reticulate to use your conda env
use_condaenv("fenicsenv", required = TRUE)
py_config()


## ---------------------------------------------------------------------------------------
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


## ---------------------------------------------------------------------------------------
# Function to compute polynomial coefficients from roots
poly.from.roots <- function(roots) {
  coef <- 1
  for (r in roots) {coef <- convolve(coef, c(1, -r), type = "open")}
  return(coef) # returned in increasing order like a+bx+cx^2+...
}


## ---------------------------------------------------------------------------------------
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


## ---------------------------------------------------------------------------------------
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
      partial_fraction_terms[[i]] <- (L - poles_rs_k$p[i] * C)#/poles_rs_k$r[i]
      }
    return(list(C = C, # mass matrix
                L = L, # Laplacian matrix scaled
                m = m, # rational order
                beta = beta, # smoothness parameter
                partial_fraction_terms = partial_fraction_terms, # partial fraction terms
                residues = poles_rs_k$r # residues \{a_k\}_{k=1}^{m+1}
                ))
  }
}


## ---------------------------------------------------------------------------------------
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
    residues <- obj$residues
    output <- v*0
    for (i in 1:(m+1)) {output <- output + residues[i] * solve(partial_fraction_terms[[i]], v)}
    return(output # solve the linear system using the partial fraction decomposition
           )
  }
}


## ---------------------------------------------------------------------------------------
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


## ---------------------------------------------------------------------------------------
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


## ---------------------------------------------------------------------------------------
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
g_cos <- function(r, A, omega) {
  return(A * cos(omega * r)) 
}
G_cos <- function(t, A, lambda_j_alpha_half, omega) {
  denom <- lambda_j_alpha_half^2 + omega^2
  numerator <- exp(lambda_j_alpha_half * t) * (lambda_j_alpha_half * cos(omega * t) + omega * sin(omega * t)) - lambda_j_alpha_half
  return(A * numerator / denom)
}


## ---------------------------------------------------------------------------------------
# Call the Python function
fem <- py_run_string("
from dolfin import *
import numpy as np
import scipy.sparse

def fem_shifted_laplacian_sparse(a, b, nx, ny):
    nx = int(nx)
    ny = int(ny)
    mesh = RectangleMesh(Point(0,0), Point(a,b), nx, ny)
    V = FunctionSpace(mesh, 'Lagrange', 1)
    u = TrialFunction(V)
    v = TestFunction(V)
    G = inner(grad(u), grad(v))*dx
    C = u * v * dx
    G = assemble(G)
    C = assemble(C)
    G = scipy.sparse.csr_matrix(G.array())
    C = scipy.sparse.csr_matrix(C.array())
    mesh_loc = mesh.coordinates()
    x = np.linspace(0, a, nx+1)
    y = np.linspace(0, b, ny+1)
    return {
        'G': G,
        'C': C,
        'mesh_loc': mesh_loc,
        'x': x,
        'y': y,
        'mesh': mesh,
        'V': V
    }

def transfer_solution_to_fine(V_coarse, U_coarse_matrix, a, b, nx_fine, ny_fine):
    nx_fine = int(nx_fine)
    ny_fine = int(ny_fine)
    
    fine_mesh = RectangleMesh(Point(0,0), Point(a,b), nx_fine, ny_fine)
    V_fine = FunctionSpace(fine_mesh, 'Lagrange', 1)
    mesh_loc = fine_mesh.coordinates()
    x_fine = np.linspace(0, a, nx_fine+1)
    y_fine = np.linspace(0, b, ny_fine+1)
    
    U_fine_list = []
    
    # Loop over time steps / columns of U_coarse_matrix
    for j in range(U_coarse_matrix.shape[1]):
        u_coarse = Function(V_coarse)
        u_coarse.vector()[:] = U_coarse_matrix[:, j]
        
        # Interpolate to fine mesh, if you prefer projection, use project(u_coarse, V_fine)
        u_fine = interpolate(u_coarse, V_fine)
        
        # Store as matrix compatible with x_fine, y_fine
        U_fine_list.append(np.reshape(u_fine.vector().get_local(), (nx_fine+1, ny_fine+1), order='F'))
    
    return {
    'U_fine_list': U_fine_list,
    'mesh_loc': mesh_loc,
    'x_fine': x_fine, 
    'y_fine': y_fine}
")


## ---------------------------------------------------------------------------------------
# Function to find globe for a given h via interpolation
globe_from_h <- function(h) {
  readed_data <- readRDS(
    file = here::here("data_files/h_min_max_vs_globe.RDS")
  )
  
  globe_vector <- readed_data$globe_vector
  h_min_vector <- readed_data$h_min_vector
  h_max_vector <- readed_data$h_max_vector
  
  ord <- order(h_max_vector)
  h_sorted <- h_max_vector[ord]
  globe_sorted <- globe_vector[ord]
  
  return(round(approx(x = h_sorted, y = globe_sorted, xout = h)$y))
}


## ---------------------------------------------------------------------------------------
from_matrix_to_list <- function(M, nrow, ncol) {
  return(lapply(1:ncol(M), function(j) matrix(M[, j], nrow = nrow, ncol = ncol, byrow = FALSE)))
}
from_list_to_matrix <- function(L) {
  return(do.call(cbind, lapply(L, function(M) as.vector(M))))
}


## ---------------------------------------------------------------------------------------
quad_weights <- function(z) {
  n <- length(z)
  w <- numeric(n)
  w[1] <- (z[2] - z[1]) / 2
  w[n] <- (z[n] - z[n-1]) / 2
  if (n > 2) {
    for (i in 2:(n-1)) {
      w[i] <- (z[i+1] - z[i-1]) / 2
    }
  }
  return(w)
}


## ---------------------------------------------------------------------------------------
compute_true_eigen_rectangle <- function(a = 1,
                                         b = 1,
                                         loc,
                                         kappa = 1, 
                                         m_max = 3, 
                                         n_max = 3) {
  loc_x <- loc[,1]
  loc_y <- loc[,2]
  N <- (m_max+1)*(n_max+1)
  # Prepare lists
  eigvals <- numeric(N)
  m_vector <- numeric(N)
  n_vector <- numeric(N)
  eigfuncs <- matrix(0, nrow = nrow(loc), ncol = N)
  
  i <- 0 
  # Loop over mode indices
  for (m in 0:m_max) {
    for (n in 0:n_max) {
      lambda_mn <- kappa^2 + pi^2 * ((m^2 / a^2) + (n^2 / b^2))
      eigvals[i + 1] <- lambda_mn
      # Evaluate eigenfunction on mesh grid
      phi_mn <-  cos(m*pi*loc_x/a) * cos(n*pi*loc_y/b)
      eigfuncs[, i + 1] <- phi_mn
      m_vector[i + 1] <- m
      n_vector[i + 1] <- n
      i <- i + 1
    }
  }
  # Sort eigenvalues and corresponding eigenfunctions
  idx <- order(eigvals)
  eigvals <- eigvals[idx]
  eigfuncs <- eigfuncs[, idx]
  m_vector <- m_vector[idx]
  n_vector <- n_vector[idx]
  
  return(list(eigvals = eigvals,
              eigfuncs = eigfuncs,
              m_vector = m_vector,
              n_vector = n_vector))
}

compute_true_eigen_sphere <- function(mesh, 
                                      kappa,
                                      L_max,
                                      rot.inv){
  eigfucs <- fmesher::fm_raw_basis(mesh = mesh, 
                                   type = "sph.harm",
                                   n = L_max, 
                                   rot.inv = rot.inv)
  eigvals <- calculate_laplace_beltrami_eigenvalues(kappa = kappa, 
                                                    L_max = L_max, 
                                                    rot.inv = rot.inv)
  return(list(eigvals = eigvals,
              eigfuncs = eigfucs))
}


## ---------------------------------------------------------------------------------------
# Function to calculate real spherical harmonics R_lm(theta, phi)
calculate_real_spherical_harmonics_vectorized <- function(l, m, theta, phi) {
  # theta and phi are vectors of length N (number of points)
  N <- length(theta)
  Y_lm_real_values <- numeric(N)
  abs_m <- abs(m)

  # Normalization constant (can be calculated once outside the loop over N points)
  norm_const <- sqrt((2*l + 1) / (4 * pi) * factorial(l - abs_m) / factorial(l + abs_m))

  # Iterate over all N points to apply the pracma::legendre function correctly
  for (i in 1:N) {
    # pracma::legendre returns a 1x(l+1) matrix for a single input cos(theta[i])
    Plm_val <- pracma::legendre(l, cos(theta[i]))[abs_m + 1]
    
    # Calculate the Y_lm value for this specific point
    if (m > 0) {
      Y_lm_real_values[i] <- sqrt(2) * norm_const * Plm_val * cos(m * phi[i])
    } else if (m < 0) {
      Y_lm_real_values[i] <- sqrt(2) * norm_const * Plm_val * sin(abs_m * phi[i])
    } else { # m = 0
      Y_lm_real_values[i] <- norm_const * Plm_val
    }
  }
  
  return(Y_lm_real_values)
}
calculate_laplace_beltrami_eigens_manual <- function(loc, kappa = 0, L_max = 5) {
  
  # Convert Cartesian to Spherical (Physics convention: theta inclination, phi azimuth)
  theta <- acos(loc[, 3])
  phi <- atan2(loc[, 2], loc[, 1])
  phi <- ifelse(phi < 0, phi + 2 * pi, phi) 

  eigenvalues <- numeric(0)
  eigenfunctions_matrix <- matrix(nrow = nrow(loc), ncol = 0)
  col_names <- character(0)
  
  for (l in 0:L_max) {
    # The eigenvalue for the Laplace-Beltrami is actually just l*(l+1). 
    # Your definition "kappa^2 + l*(l+1)" might be for a modified Helmholtz equation/operator.
    # Sticking to your provided calculation:
    lambda_l <- kappa^2 + l * (l + 1)
    
    for (m in -l:l) {
      # Use the new vectorized function:
      Y_lm_values <- calculate_real_spherical_harmonics_vectorized(l, m, theta, phi)
      
      eigenvalues <- c(eigenvalues, lambda_l)
      eigenfunctions_matrix <- cbind(eigenfunctions_matrix, Y_lm_values)
      col_names <- c(col_names, paste0("Y_", l, "_", m, "_k_", kappa))
    }
  }
  
  colnames(eigenfunctions_matrix) <- col_names
  
  return(list(
    eigenvalues = eigenvalues,
    eigenfunctions = eigenfunctions_matrix
  ))
}


calculate_laplace_beltrami_eigenvalues <- function(kappa = 0, L_max = 5, rot.inv = FALSE) {
  
  eigenvalues <- numeric(0)
  
  for (l in 0:L_max) {
    lambda_l <- kappa^2 + l * (l + 1)
    
    if (rot.inv) {
      # Only one eigenvalue per l (rotational invariance)
      eigenvalues <- c(eigenvalues, lambda_l)
    } else {
      # Full multiplicity: (m = -l,...,l)
      multiplicity <- 2 * l + 1
      eigenvalues <- c(eigenvalues, rep(lambda_l, multiplicity))
    }
  }
  
  return(eigenvalues)
}




## ---------------------------------------------------------------------------------------
global.scene.setter <- function(x_range, y_range, z_range) {
  
  return(list(xaxis = list(title = "x", range = x_range),
              yaxis = list(title = "y", range = y_range),
              zaxis = list(title = "z", range = z_range),
              aspectratio = list(x = 1, 
                                 y = 1, 
                                 z = 1),
              camera = list(eye = list(x = 1.5, 
                                       y = 1.5, 
                                       z = 1.5),
                            center = list(x = 0, 
                                          y = 0, 
                                          z = 0))))
}


## ---------------------------------------------------------------------------------------
simple3d_plotter <- function(mesh_loc, U_0) {
interp_res <- interp(
  x = mesh_loc[,2],
  y = mesh_loc[,1],
  z = U_0,
  duplicate = "mean"
)

plot_ly(
  x = interp_res$y,
  y = interp_res$x,
  z = interp_res$z,
  type = "surface"
)}


## ---------------------------------------------------------------------------------------
plot_3d_slider <- function(loc, nx, ny, eigvals, eigfuncs) {
  eigfuncs <- from_matrix_to_list(eigfuncs, nx+1, ny+1)
  colorscale <- "Viridis"
  # Build frames (each frame *must* include the colorscale)
  frames <- lapply(seq_along(eigfuncs), function(i) {
    list(
      name = as.character(i),
      data = list(list(
        z = eigfuncs[[i]],
        type = "surface",
        colorscale = colorscale,
        showscale = TRUE
      ))
    )
  })
  
  X <- matrix(loc[,1], nx+1,ny+1)
  Y <- matrix(loc[,2], nx+1,ny+1)
  # Initial plot
  p <- plot_ly(
    x = X,
    y = Y,
    z = eigfuncs[[1]],
    type = "surface",
    colorscale = colorscale,
    showscale = TRUE,
    frame = "1"
  )
  
  p$x$frames <- frames
  
  z_range <- range(unlist(eigfuncs))
  x_range <- range(x)
  y_range <- range(y)
  
  frame_name <- deparse(substitute(eigvals))
  
  # Layout + slider
  p <- p %>% layout(
    title = paste0(frame_name, ": ", eigvals[1]),
    scene = global.scene.setter(x_range, y_range, z_range),
    sliders = list(
      list(
        active = 0,
        currentvalue = list(prefix = "Frame: "),
        pad = list(t = 50),
        steps = lapply(seq_along(eigfuncs), function(i) {
          list(
            label = as.character(i),
            method = "animate",
            args = list(list(as.character(i)),
                        list(mode = "immediate",
                             frame = list(duration = 300, redraw = TRUE),
                             transition = list(duration = 0)))
          )
        })
      )
    ),
    updatemenus = list(
      list(
        type = "buttons",
        showactive = FALSE,
        y = 1,
        x = 1.15,
        xanchor = "right",
        yanchor = "top",
        buttons = list(
          list(label = "Play",
               method = "animate",
               args = list(NULL, list(frame = list(duration = 300, redraw = TRUE),
                                      fromcurrent = TRUE, mode = "immediate"))),
          list(label = "Pause",
               method = "animate",
               args = list(NULL, list(frame = list(duration = 0, redraw = FALSE),
                                      mode = "immediate")))
        )
      )
    )
  ) %>% plotly_build()
  for (i in seq_along(p$x$frames)) {
    p$x$frames[[i]]$layout <- list(title = paste0(frame_name, ": ", eigvals[i]))
  }
  return(p)
}


## ---------------------------------------------------------------------------------------
plot_3d_slider_scatter <- function(loc, eigvals, eigfuncs) {

  x <- loc[,1]
  y <- loc[,2]

  colorscale <- "Viridis"

  # Frames ----------------------------------------------------------------
  frames <- lapply(seq_len(ncol(eigfuncs)), function(i) {
    list(
      name = as.character(i),
      data = list(list(
        x = x,
        y = y,
        z = eigfuncs[,i],
        type = "scatter3d",
        mode = "markers",
        marker = list(
          size = 5,
          color = eigfuncs[,i],
          colorscale = colorscale,
          showscale = TRUE
        )
      ))
    )
  })

  # Initial Plot ----------------------------------------------------------
  p <- plot_ly(
    x = x,
    y = y,
    z = eigfuncs[,1],
    type = "scatter3d",
    mode = "markers",
    marker = list(
      size = 5,
      color = eigfuncs[,1],
      colorscale = colorscale,
      showscale = TRUE
    ),
    frame = "1"
  )

  p$x$frames <- frames

  z_range <- range(eigfuncs)
  x_range <- range(x)
  y_range <- range(y)

  frame_name <- deparse(substitute(eigvals))

  # Layout + slider -------------------------------------------------------
  p <- p %>% layout(
    title = paste0(frame_name, ": ", eigvals[1]),
    scene = global.scene.setter(x_range, y_range, z_range),
    sliders = list(
      list(
        active = 0,
        currentvalue = list(prefix = "Frame: "),
        pad = list(t = 50),
        steps = lapply(seq_len(ncol(eigfuncs)), function(i) {
          list(
            label = as.character(i),
            method = "animate",
            args = list(
              list(as.character(i)),
              list(frame = list(duration = 300, redraw = TRUE),
                   mode = "immediate")
            )
          )
        })
      )
    ),
    updatemenus = list(
      list(
        type = "buttons",
        showactive = FALSE,
        y = 1,
        x = 1.15,
        xanchor = "right",
        yanchor = "top",
        buttons = list(
          list(label = "Play",
               method = "animate",
               args = list(
                 NULL,
                 list(frame = list(duration = 300, redraw = TRUE),
                      fromcurrent = TRUE)
               )),
          list(label = "Pause",
               method = "animate",
               args = list(NULL, list(frame = list(duration = 0))))
        )
      )
    )
  ) %>% plotly_build()

  # Update title in each frame --------------------------------------------
  for (i in seq_len(ncol(eigfuncs))) {
    p$x$frames[[i]]$layout <- list(
      title = paste0(frame_name, ": ", eigvals[i])
    )
  }

  return(p)
}

# x <- mesh_fine$loc[,1]
# y <- mesh_fine$loc[,2]
# z <- eigfuncs[,9]
# 
# plotly::plot_ly(
#   x = x,
#   y = y,
#   z = z,
#   type = "scatter3d",
#   mode = "markers",
#   marker = list(size = 5, color = z, colorscale = "Viridis")
# )



## ---------------------------------------------------------------------------------------
plot_3d_slider_sphere <- function(mesh, eigvals, eigfuncs, fixed_colorscale = TRUE) {
  
  colorscale = "Viridis"
  opacity = 1
  
  x <- mesh$loc[, 1]
  y <- mesh$loc[, 2]
  z <- mesh$loc[, 3]
  tri <- mesh$graph$tv
  
  # Global intensity limits if fixed
  if (fixed_colorscale) {
    cmin <- min(eigfuncs)
    cmax <- max(eigfuncs)
  } else {
    cmin <- NULL
    cmax <- NULL
  }
  
  # Create frames
  frames <- lapply(seq_len(ncol(eigfuncs)), function(i) {
    fvals <- eigfuncs[,i]
    list(
      name = as.character(i),
      data = list(list(
        x = x,
        y = y,
        z = z,
        i = tri[,1] - 1,
        j = tri[,2] - 1,
        k = tri[,3] - 1,
        type = "mesh3d",
        intensity = fvals,
        colorscale = colorscale,
        opacity = opacity,
        flatshading = TRUE,
        cmin = cmin,
        cmax = cmax,
        text = paste0(
          "x: ", sprintf("%.3f", x), "<br>",
          "y: ", sprintf("%.3f", y), "<br>",
          "z: ", sprintf("%.3f", z), "<br>",
          "f: ", sprintf("%.5f", fvals)
        ),
        hoverinfo = "text"
      ))
    )
  })
  
  # Initial plot
  fvals0 <- eigfuncs[,1]
  
  p <- plot_ly(
    x = x, y = y, z = z,
    i = tri[,1] - 1,
    j = tri[,2] - 1,
    k = tri[,3] - 1,
    type = "mesh3d",
    intensity = fvals0,
    colorscale = colorscale,
    opacity = opacity,
    flatshading = TRUE,
    cmin = cmin,
    cmax = cmax,
    text = paste0(
      "x: ", sprintf("%.3f", x), "<br>",
      "y: ", sprintf("%.3f", y), "<br>",
      "z: ", sprintf("%.3f", z), "<br>",
      "f: ", sprintf("%.5f", fvals0)
    ),
    hoverinfo = "text",
    frame = "1"
  )
  
  frame_name <- deparse(substitute(eigvals))
  
  p$x$frames <- frames
  
  # Layout + slider + buttons
  p <- p %>% layout(
    title = paste0(frame_name, ": ", eigvals[1]),
    sliders = list(
      list(
        active = 0,
        currentvalue = list(prefix = "Frame: "),
        pad = list(t = 50),
        steps = lapply(seq_len(ncol(eigfuncs)), function(i) {
          list(
            label = as.character(i),
            method = "animate",
            args = list(list(as.character(i)),
                        list(mode = "immediate",
                             frame = list(duration = 300, redraw = TRUE),
                             transition = list(duration = 0)))
          )
        })
      )
    ),
    updatemenus = list(
      list(
        type = "buttons",
        showactive = FALSE,
        y = 1,
        x = 1.15,
        xanchor = "right",
        yanchor = "top",
        buttons = list(
          list(label = "Play",
               method = "animate",
               args = list(NULL, list(frame = list(duration = 300, redraw = TRUE),
                                      fromcurrent = TRUE, mode = "immediate"))),
          list(label = "Pause",
               method = "animate",
               args = list(NULL, list(frame = list(duration = 0, redraw = FALSE),
                                      mode = "immediate")))
        )
      )
    )
  ) %>% plotly_build()
  
  # Add titles for frames
  for (i in seq_len(ncol(eigfuncs))) {
    p$x$frames[[i]]$layout <- list(
      title = paste0(frame_name, ": ", eigvals[i])
    )
  }
  
  return(p)
}



## ---------------------------------------------------------------------------------------
plot_3d_slider_sphere_scatter <- function(mesh, eigvals, eigfuncs,
                                          fixed_colorscale = TRUE) {
  
  colorscale = "Viridis"
  
  x <- mesh$loc[,1]
  y <- mesh$loc[,2]
  z <- mesh$loc[,3]
  
  # Compute global color limits if fixed
  if (fixed_colorscale) {
    cmin <- min(eigfuncs)
    cmax <- max(eigfuncs)
  } else {
    cmin <- NULL
    cmax <- NULL
  }
  
  # Create frames
  frames <- lapply(seq_len(ncol(eigfuncs)), function(i) {
    
    fvals <- eigfuncs[,i]
    
    list(
      name = as.character(i),
      data = list(list(
        x = x,
        y = y,
        z = z,
        type = "scatter3d",
        mode = "markers",
        marker = list(
          size = 5,
          color = fvals,
          colorscale = colorscale,
          showscale = TRUE,
          cmin = cmin,
          cmax = cmax
        ),
        text = paste0(
          "x: ", sprintf("%.3f", x), "<br>",
          "y: ", sprintf("%.3f", y), "<br>",
          "z: ", sprintf("%.3f", z), "<br>",
          "f: ", sprintf("%.5f", fvals)
        ),
        hoverinfo = "text"
      ))
    )
  })
  
  # Initial plot (frame 1)
  fvals0 <- eigfuncs[,1]
  
  p <- plot_ly(
    x = x, y = y, z = z,
    type = "scatter3d",
    mode = "markers",
    marker = list(
      size = 5,
      color = fvals0,
      colorscale = colorscale,
      showscale = TRUE,
      cmin = cmin,
      cmax = cmax
    ),
    text = paste0(
      "x: ", sprintf("%.3f", x), "<br>",
      "y: ", sprintf("%.3f", y), "<br>",
      "z: ", sprintf("%.3f", z), "<br>",
      "f: ", sprintf("%.5f", fvals0)
    ),
    hoverinfo = "text",
    frame = "1"
  )
  
  frame_name <- deparse(substitute(eigvals))
  
  p$x$frames <- frames
  
  # Slider + play/pause
  p <- p %>% layout(
    title = paste0(frame_name, ": ", eigvals[1]),
    sliders = list(
      list(
        active = 0,
        currentvalue = list(prefix = "Mode: "),
        pad = list(t = 50),
        steps = lapply(seq_len(ncol(eigfuncs)), function(i) {
          list(
            label = as.character(i),
            method = "animate",
            args = list(list(as.character(i)),
                        list(mode = "immediate",
                             frame = list(duration = 300, redraw = TRUE),
                             transition = list(duration = 0)))
          )
        })
      )
    ),
    updatemenus = list(
      list(
        type = "buttons",
        showactive = FALSE,
        y = 1,
        x = 1.15,
        xanchor = "right",
        yanchor = "top",
        buttons = list(
          list(
            label = "Play",
            method = "animate",
            args = list(NULL,
                        list(frame = list(duration = 300, redraw = TRUE),
                             fromcurrent = TRUE,
                             mode = "immediate"))
          ),
          list(
            label = "Pause",
            method = "animate",
            args = list(NULL,
                        list(frame = list(duration = 0, redraw = FALSE),
                             mode = "immediate"))
          )
        )
      )
    )
  ) %>% plotly_build()
  
  # Update title per frame
  for (i in seq_len(ncol(eigfuncs))) {
    p$x$frames[[i]]$layout <- list(
      title = paste0(frame_name, ": ", eigvals[i])
    )
  }
  
  return(p)
}

