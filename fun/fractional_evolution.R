## ---------------------------------------------------------------------------------------------------------------
# remotes::install_github("davidbolin/rspde", ref = "devel")
# remotes::install_github("davidbolin/metricgraph", ref = "devel")
library(rSPDE)
library(MetricGraph)
library(grateful)

library(ggplot2)
library(reshape2)
library(plotly)


## ---------------------------------------------------------------------------------------------------------------
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


## ---------------------------------------------------------------------------------------------------------------
# Function to compute polynomial coefficients from roots
poly.from.roots <- function(roots) {
  coef <- 1
  for (r in roots) {coef <- convolve(coef, c(1, -r), type = "open")}
  return(coef) # returned in increasing order like a+bx+cx^2+...
}


## ---------------------------------------------------------------------------------------------------------------
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
  #res <- gsignal::residue(factor_pr_coef, pr_plus_pl_coef)
  # Coefficients passed here are in descending powers.
  # polyroot() expects ascending powers.
  p <- polyroot(rev(pr_plus_pl_coef))
  
  # Verify that these poles are numerically real and negative.
  stopifnot(
    all(abs(Im(p)) <= 1e-10 * pmax(1, abs(Re(p)))),
    all(Re(p) < 0)
  )
  p <- Re(p)
  
  eval_poly <- function(coef, x) {
    Reduce(function(value, a) value * x + a, coef, init = 0)
  }
  
  r <- vapply(seq_along(p), function(j) {
    eval_poly(factor_pr_coef, p[j]) /
      (pr_plus_pl_coef[1] * prod(p[j] - p[-j]))
  }, numeric(1))
  
  res <- list(r = r, p = p, k = numeric(0))
  return(list(r = res$r, # residues \{a_k\}_{k=1}^{m+1}
              p = res$p, # poles \{p_k\}_{k=1}^{m+1}
              k = res$k # remainder r
              ))
}


## ---------------------------------------------------------------------------------------------------------------
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
    poles_rs_k <- compute.partial.fraction.param(factor = roots$factor,
                                                 pr_roots = roots$pr_roots,
                                                 pl_roots = roots$pl_roots,
                                                 time_step = time_step, 
                                                 scaling = scaling)

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


## ---------------------------------------------------------------------------------------------------------------
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

# my.solver.frac <- function(obj, v) {
#   beta <- obj$beta
#   m <- obj$m
#   
#   if (beta == 1) {
#     return(solve(obj$LHS, v))
#   } 
#   
#   partial_fraction_terms <- obj$partial_fraction_terms
#   residues <- obj$residues
#   
#   # --- PARALLEL PART ---
#   results <- parallel::mclapply(
#     1:(m+1),
#     function(i) residues[i] * solve(partial_fraction_terms[[i]], v),
#     mc.cores = parallel::detectCores()  
#   )
#   
#   return(Reduce(`+`, results))
# }


## ---------------------------------------------------------------------------------------------------------------
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


## ---------------------------------------------------------------------------------------------------------------
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

