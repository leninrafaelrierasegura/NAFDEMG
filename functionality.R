# This file contains functions that are used in the vignettes
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# For each rational order m (1,2,3,4) and smoothness parameter beta (= alpha/2 with alpha between 0.5 and 2), Function my.get.roots() (adapted from the rSPDE package) returns c_m/b_{m+1} and the roots of p_\ell and p_r. 
# See David Bolin & Kristin Kirchner (2020) The Rational SPDE Approach for Gaussian Random Fields With General Smoothness, Journal of Computational and Graphical Statistics, 29:2, 274-285, DOI: 10.1080/10618600.2019.1665537
# for the definition of c_m/b_{m+1} and p_\ell and p_r.
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
  return(list(rb = rb, # roots of p_\ell
              rc = rc, # roots of p_r
              factor = factor # this is c_m/b_{m+1}
              ))
}
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Function poly.from.roots() computes the coefficients of a polynomial from its roots.
poly.from.roots <- function(roots) {
  coef <- 1
  for (r in roots) {coef <- convolve(coef, c(1, -r), type = "open")}
  return(coef) # returned in increasing order like a+bx+cx^2+...
}
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Function compute.partial.fraction.param() computes the parameters for the partial fraction decomposition of the rational function
# \dfrac{p_r(x)}{p_r(x)+\tau \kappa^{2\beta}p_\ell(x) } = \sum_{k=1}^{m+1}a_k(x-p_k)^{-1} + r
compute.partial.fraction.param <- function(factor, pr_roots, pl_roots, time_step, scaling) {
  pr_coef <- c(0, poly.from.roots(pr_roots)) 
  pl_coef <- poly.from.roots(pl_roots) 
  factor_pr_coef <- pr_coef
  pr_plus_pl_coef <- factor_pr_coef + ((scaling * time_step)/factor) * pl_coef
  res <- gsignal::residue(factor_pr_coef, pr_plus_pl_coef)
  return(list(r = res$r, p = res$p, k = res$k)) 
}