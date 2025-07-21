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
poly.from.roots <- function(roots) {
  coef <- 1
  for (r in roots) {coef <- convolve(coef, c(1, -r), type = "open")}
  return(coef) # returned in increasing order like a+bx+cx^2+...
}


## -----------------------------------------------------------------------------
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
      # Here is where the terms in the sum eq 11 are computed
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
# Function to compute the eigenfunctions 
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
# Function to order the vertices for plotting
plotting.order <- function(v, graph){
  edge_number <- graph$mesh$VtE[, 1]
  pos <- sum(edge_number == 1)+1
  return(c(v[1], v[3:pos], v[2], v[(pos+1):length(v)], v[2]))
}


## -----------------------------------------------------------------------------
grateful::cite_packages(output = "paragraph", out.dir = ".")

