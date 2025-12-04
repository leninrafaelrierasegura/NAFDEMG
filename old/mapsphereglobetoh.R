library(fmesher)
library(Matrix)
library(ggplot2)
library(tidyr)
library(dplyr)

# globe_vector <-1:100
# h_min_vector <- numeric(length(globe_vector))
# h_max_vector <- numeric(length(globe_vector))
# for (i in 1:length(globe_vector)) {
#   globe <- globe_vector[i]
#   rangeof_h <- range(diag(fm_fem(fm_rcdt_2d(globe = globe), order = 2)$c0))
#   h_min_vector[i] <- rangeof_h[1]
#   h_max_vector[i] <- rangeof_h[2]
#   print(i)
# }
# 
# 
# saveRDS(
#   list(
#     globe_vector = globe_vector,
#     h_min_vector = h_min_vector,
#     h_max_vector = h_max_vector
#   ),
#   file = here::here("data_files/h_min_max_vs_globe.RDS")
# )

readed_data <- readRDS(
  file = here::here("data_files/h_min_max_vs_globe.RDS")
)

globe_vector <- readed_data$globe_vector
h_min_vector <- readed_data$h_min_vector
h_max_vector <- readed_data$h_max_vector

# Create a data frame
df <- data.frame(
  globe = globe_vector,
  h_min = h_min_vector,
  h_max = h_max_vector
)

# Convert to long format for ggplot
df_long <- df %>%
  pivot_longer(cols = c(h_min, h_max),
               names_to = "type",
               values_to = "h_value")

# Plot
p <- ggplot(df_long, aes(x = globe, y = h_value, color = type)) +
  geom_line(linewidth = 1) +
  labs(x = "Globe", y = "h value", color = "Type") +
  theme_minimal()

plotly::ggplotly(p)

closest_globe <- function(h) {
  idx <- which.min(abs(h_max_vector - h))  # index of closest value
  globe_vector[idx]  # return corresponding globe
}

# Example usage
closest_globe(6e-3)

ord <- order(h_max_vector)
h_sorted <- h_max_vector[ord]
globe_sorted <- globe_vector[ord]

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

# Example usage
globe_from_h(6e-3)

ord <- order(h_max_vector)
h_sorted <- h_max_vector[ord]
globe_sorted <- globe_vector[ord]

# Create a smooth spline function: globe as function of h_max
spline_func <- splinefun(x = h_sorted, y = globe_sorted, method = "monoH.FC")

# Function to get globe from h using spline
globe_from_h_spline <- function(h) {
  spline_func(h)
}
# Example usage
globe_from_h_spline(6e-3)