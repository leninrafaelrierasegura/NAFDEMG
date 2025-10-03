source(here::here("control_functionality.R"))
library(pracma)
m_vector <- 1:4
errors <- rep(0, length(m_vector))
for (i in 1:length(m_vector)) {
m <- m_vector[i]
beta <- readRDS("data_files/chebfun_tables.RDS")[[1]]$beta[100]

res <- my.get.roots(m = m, beta = beta)

factor <- res$factor
pr_roots <- res$pr_roots
pl_roots <- res$pl_roots

exp_x <- function(x){
  return(x^(beta-1))
}

rat_x <- function(x) {
  num <- apply(outer(x, pr_roots, "-"), 1, prod)
  den <- apply(outer(x, pl_roots, "-"), 1, prod)
  factor * num / den
}

h <- pi * sqrt(m / (1-beta))
upper_x <- 1
lower_x <- 10^(-(5+m)/2)+0
x <- seq(lower_x, upper_x, length.out = ((upper_x - 0) / h + 1))


f <- exp_x(x) #f_x(x)
g <- rat_x(x) #g_x(x)

errors[i] <- max(abs(f-g))
}


coef(lm(log(errors) ~ sqrt(m_vector)))[2]

- 2 * pi * sqrt(1-beta)











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