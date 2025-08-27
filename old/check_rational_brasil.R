library(pracma)

real_residue <- function(b, a, tol = 1e-8) {
  # Compute residues
  res <- gsignal::residue(b, a)
  r <- res$r
  p <- res$p
  k <- res$k  # direct term (polynomial) if any
  
  real_terms <- list()
  used <- rep(FALSE, length(r))
  
  for (i in seq_along(r)) {
    if (used[i]) next
    
    if (abs(Im(p[i])) < tol) {
      # Real pole
      real_terms[[length(real_terms) + 1]] <- list(
        numerator = Re(r[i]),
        denominator = c(1, -Re(p[i]))  # s - p_i
      )
      used[i] <- TRUE
    } else {
      # Complex pole: find its conjugate
      conj_idx <- which(abs(p - Conj(p[i])) < tol & !used)
      if (length(conj_idx) == 0) stop("Conjugate not found")
      j <- conj_idx[1]
      
      # Form real quadratic numerator / denominator
      x <- Re(p[i])
      y <- Im(p[i])
      r1 <- r[i]
      r2 <- r[j]
      
      # Real numerator: A s + B
      A <- 2 * Re(r1)
      B <- -2 * Re(r1) * x + 2 * Im(r1) * y
      
      real_terms[[length(real_terms) + 1]] <- list(
        numerator = c(A, B),
        denominator = c(1, -2*x, x^2 + y^2)  # (s - x)^2 + y^2
      )
      used[c(i,j)] <- TRUE
    }
  }
  
  return(list(real_terms = real_terms, direct = k))
}

alpha <- 0.6
# to check if the rational approximation is correct
res <- rSPDE:::interp_rational_coefficients(order = 20, 
                                            type_rational_approx = "brasil", 
                                            type_interp = "spline", 
                                            alpha = alpha)
res
pol <- gsignal::rresidue(res$r, res$p, res$k, tol = 0.00001)
pol

g_x <- function(x) {
  term_matrix <- outer(x, res$p, "-") # matrix of (x - p[i])
  term_values <- sweep(1 / term_matrix, 2, res$r, "*") # multiply columns by r[i]
  return(rowSums(term_values) + res$k)
}

f_x <- function(x) {
  up <- sapply(x, function(xx) sum(pol$b * xx^(rev(seq_along(pol$b))-1)))
  down <- sapply(x, function(xx) sum(pol$a * xx^(rev(seq_along(pol$a))-1)))
  up / down
}

x <- seq(0, 1, by = 0.01)
sum((f_x(x) - g_x(x))^2)  # Check if the two polynomials are equal
df <- data.frame(x = x, f = f_x(x), g = g_x(x))

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

kappa <- 4
time_step <- 0.001
scaling <- (kappa^2)^alpha
out <- gsignal::residue(pol$b, pol$b + scaling * time_step * pol$a)
out

out2 <- real_residue(pol$b, pol$b + scaling * time_step * pol$a)

roots_of_a <- polyroot(rev(pol$b + scaling * time_step * pol$a))
roots_of_a

G_x <- function(x) {
  term_matrix <- outer(x, out$p, "-") # matrix of (x - p[i])
  term_values <- sweep(1 / term_matrix, 2, out$r, "*") # multiply columns by r[i]
  return(rowSums(term_values) + out$k)
}

F_x <- function(x) {
  up <- sapply(x, function(xx) sum(pol$b * xx^(rev(seq_along(pol$b))-1)))
  down <- sapply(x, function(xx) sum(pol$a * xx^(rev(seq_along(pol$a))-1)))
  return(up/(up + time_step * scaling * down))
}


sum((F_x(x) - G_x(x))^2)  # Check if the two polynomials are equal
df <- data.frame(x = x, f = F_x(x), g = G_x(x))

p <- ggplot(df, aes(x = x)) +
  geom_line(aes(y = f, color = "p/(p+tau kappa ^ alpha q)"), size = 1) +
  geom_point(aes(y = f, color = "p/(p+tau kappa ^ alpha q)"), size = 1.5) +
  geom_line(aes(y = g, color = "partial fraction p/(p+tau kappa ^ alpha q)"), size = 1, linetype = "dashed") +
  geom_point(aes(y = g, color = "partial fraction p/(p+tau kappa ^ alpha q)"), size = 1.5) +
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
