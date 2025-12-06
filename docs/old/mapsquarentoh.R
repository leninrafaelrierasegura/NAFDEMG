library(fmesher)
library(Matrix)
library(ggplot2)
library(tidyr)
library(dplyr)

a <- 100
b <- 100
n_vector <-10:100
h_min_vector <- numeric(length(n_vector))
h_max_vector <- numeric(length(n_vector))
another_h_vector <- numeric(length(n_vector))
for (i in 1:length(n_vector)) {
  n <- n_vector[i]
  mesh <- fm_rcdt_2d(
    lattice = fm_lattice_2d(
      x = seq(0, a, length.out = n + 1),
      y = seq(0, b, length.out = n + 1)
    ), extend = FALSE)
  C_ii <- diag(fm_fem(mesh, order = 2)$c0)
  h_min_vector[i] <- mean(sqrt(C_ii))
  another_h_vector[i] <- sqrt(mean(C_ii))
  h_max_vector[i] <- sqrt(a*b)/(n+1)
  print(i)
}


saveRDS(
  list(
    n_vector = n_vector,
    h_min_vector = h_min_vector,
    h_max_vector = h_max_vector,
    another_h_vector = another_h_vector
  ),
  file = here::here("data_files/h_min_max_vs_n.RDS")
)

readed_data <- readRDS(
  file = here::here("data_files/h_min_max_vs_n.RDS")
)

n_vector <- readed_data$n_vector
h_min_vector <- readed_data$h_min_vector
h_max_vector <- readed_data$h_max_vector
another_h_vector <- readed_data$another_h_vector

# Create a data frame
df <- data.frame(
  n = n_vector,
  h_min = h_min_vector,
  h_max = h_max_vector,
  another_h = another_h_vector
)

# Convert to long format for ggplot
df_long <- df %>%
  pivot_longer(cols = c(h_min, h_max, another_h),
               names_to = "type",
               values_to = "h_value")

# Plot
p <- ggplot(df_long, aes(x = n, y = h_value, color = type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(x = "n", y = "h value", color = "Type") +
  theme_minimal()

plotly::ggplotly(p)

closest_n <- function(h) {
  idx <- which.min(abs(h_max_vector - h))  # index of closest value
  n_vector[idx]  # return corresponding n
}

# Example usage
closest_n(6e-3)

ord <- order(h_max_vector)
h_sorted <- h_max_vector[ord]
n_sorted <- n_vector[ord]

# Function to find n for a given h via interpolation
n_from_h <- function(h) {
  readed_data <- readRDS(
    file = here::here("data_files/h_min_max_vs_n.RDS")
  )
  
  n_vector <- readed_data$n_vector
  h_min_vector <- readed_data$h_min_vector
  h_max_vector <- readed_data$h_max_vector
  
  ord <- order(h_max_vector)
  h_sorted <- h_max_vector[ord]
  n_sorted <- n_vector[ord]
  
  return(round(approx(x = h_sorted, y = n_sorted, xout = h)$y))
}

# Example usage
n_from_h(6e-3)

ord <- order(h_max_vector)
h_sorted <- h_max_vector[ord]
n_sorted <- n_vector[ord]

# Create a smooth spline function: n as function of h_max
spline_func <- splinefun(x = h_sorted, y = n_sorted, method = "monoH.FC")

# Function to get n from h using spline
n_from_h_spline <- function(h) {
  spline_func(h)
}
# Example usage
n_from_h_spline(6e-3)