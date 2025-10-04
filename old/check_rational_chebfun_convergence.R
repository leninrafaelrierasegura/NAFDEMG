source(here::here("control_functionality.R"))
exp_x <- function(x, beta){
  return(x^(beta-1))
}
rat_x <- function(x, pr_roots, pl_roots, factor) {
  num <- apply(outer(x, pr_roots, "-"), 1, prod)
  den <- apply(outer(x, pl_roots, "-"), 1, prod)
  return(factor * num / den)
}
m_vector <- 1:8
beta_vector_index <- c(1, 10, 100, 200, 400)
beta_vector <- rep(NA, length(beta_vector_index))
errors <- matrix(NA, nrow = length(m_vector), ncol = length(beta_vector_index))
for (j in 1:length(beta_vector_index)) {
  beta <- readRDS("data_files/chebfun_tables.RDS")[[1]]$beta[beta_vector_index[j]]
  beta_vector[j] <- beta
  for (i in 1:length(m_vector)) {
    m <- m_vector[i]
    
    res <- my.get.roots(m = m, beta = beta)
    
    factor <- res$factor
    pr_roots <- res$pr_roots
    pl_roots <- res$pl_roots
    
    h <- exp(-pi * sqrt(m / (1-beta)))
    upper_x <- 1
    lower_x <- 10^(-(5+m)/2) + 0
    x <- seq(lower_x, upper_x, length.out = ((upper_x - lower_x) / h + 1))
    
    f <- exp_x(x, beta) 
    g <- rat_x(x, pr_roots, pl_roots, factor) 
    
    errors[i, j] <- max(abs(f-g))
  }
}


observed_rates <- rep(NA, length(beta_vector_index))
for (u in 1:length(beta_vector_index)) {
  observed_rates[u] <- coef(lm(log(errors[, u]) ~ sqrt(m_vector)))[2]
}

theoretical_rates <- - 4 * pi * sqrt(1-beta_vector)

p_m <- error.convergence.plotter(x_axis_vector = m_vector, 
                                 beta_vector, 
                                 errors, 
                                 theoretical_rates, 
                                 observed_rates,
                                 line_equation_fun = exp_line_equation,
                                 fig_title = expression("Convergence in " * italic(m)),
                                 x_axis_label = expression(italic(m)),
                                 apply_sqrt = TRUE)
p_m






