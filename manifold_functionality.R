## -----------------------------------------------------------------------------
library(plotly)


## ----setup, include = FALSE---------------------------------------------------
# to install in terminal
# conda activate fenicsenv
# conda install -c conda-forge matplotlib plotly ipywidgets
library(reticulate)
# VERY IMPORTANT: tell reticulate to use your conda env
use_condaenv("fenicsenv", required = TRUE)
py_config()


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
    residues <- obj$residues
    output <- v*0
    for (i in 1:(m+1)) {output <- output + residues[i] * solve(partial_fraction_terms[[i]], v)}
    return(as.vector(output) # solve the linear system using the partial fraction decomposition
           )
  }
}


## -----------------------------------------------------------------------------
solve_fractional_evolution <- function(my_op_frac, time_step, time_seq, val_at_0) {
  CC <- my_op_frac$C
  SOL <- matrix(NA, nrow = nrow(CC), ncol = length(time_seq))
  SOL[, 1] <- val_at_0
  for (k in 1:(length(time_seq) - 1)) {
    rhs <- as.vector(CC %*% SOL[, k])
    SOL[, k + 1] <- my.solver.frac(my_op_frac, rhs)
  }
  return(SOL)
}


## -----------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------
from_matrix_to_list <- function(M, nrow, ncol) {
  return(lapply(1:ncol(M), function(j) matrix(M[, j], nrow = nrow, ncol = ncol, byrow = FALSE)))
}
from_list_to_matrix <- function(L) {
  return(do.call(cbind, lapply(L, function(M) as.vector(M))))
}


## -----------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------
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
      i <- i + 1
    }
  }
  # Sort eigenvalues and corresponding eigenfunctions
  idx <- order(eigvals)
  eigvals <- eigvals[idx]
  eigfuncs <- eigfuncs[, idx]
  
  return(list(eigvals = eigvals,
              eigfuncs = eigfuncs))
}


## -----------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------
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

