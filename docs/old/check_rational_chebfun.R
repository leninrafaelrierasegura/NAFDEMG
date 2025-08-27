
m <- 4
beta <- readRDS("data_files/chebfun_tables.RDS")[[1]]$beta[1]

res <- my.get.roots(m = m, beta = beta)

factor <- res$factor
pr_roots <- res$pr_roots
pl_roots <- res$pl_roots

kappa <- 15
scaling <- (kappa^2)^beta
time_step <- 0.001

pr_coef <- poly.from.roots(pr_roots)
pl_coef <- poly.from.roots(pl_roots)

pr_plus_pl_coef <- c(0, pr_coef) + ((scaling * time_step)/factor) * pl_coef
res <- gsignal::residue(pr_coef, pr_plus_pl_coef)
res

den_roots <- Re(polyroot(rev(pr_plus_pl_coef)))
den_roots

num_vals <- polyval(pr_coef, den_roots)
den_deriv <- Re(polyval(polyder(pr_plus_pl_coef), den_roots))
residues <- num_vals / den_deriv
Re(residues)



g_x <- function(x) {
  term_matrix <- outer(x, den_roots, "-") # matrix of (x - p[i])
  term_values <- sweep(1 / term_matrix, 2, Re(residues), "*") # multiply columns by r[i]
  return(rowSums(term_values))
}

f_x <- function(x) {
  up <- sapply(x, function(xx) sum(pr_coef * xx^(rev(seq_along(pr_coef))-1)))
  down <- sapply(x, function(xx) sum(pr_plus_pl_coef * xx^(rev(seq_along(pr_plus_pl_coef))-1)))
  up / down
}

exp_x <- function(x){
  return(x^(beta-1))
}

# rat_x <- function(x){
#   factor <- res$factor
#   up_roots <- res$pr_roots
#   down_roots <- res$pl_roots
#   rat <- factor * prod((x - up_roots)) / prod((x - down_roots))
#   return(rat)
# }

rat_x <- function(x) {
  num <- apply(outer(x, pr_roots, "-"), 1, prod)
  den <- apply(outer(x, pl_roots, "-"), 1, prod)
  factor * num / den
}


x <- seq(0.01, 1, by = 0.01)


f <- f_x(x)
g <- g_x(x)


sum((f - g)^2)  # Check if the two polynomials are equal
df <- data.frame(x = x, f = f, g = g)

p <- ggplot(df, aes(x = x)) +
  geom_line(aes(y = f, color = "p/q"), size = 1) +
  geom_point(aes(y = f, color = "p/q"), size = 1.5) +
  geom_line(aes(y = g, color = "partialfraction p/q"), size = 1, linetype = "dashed") +
  geom_point(aes(y = g, color = "partialfraction p/q"), size = 1.5) +
  labs(
    title = "Polynomials Comparison",
    x = "x",
    y = "f(x) and g(x)",
    color = "Function"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
plotly::ggplotly(p)