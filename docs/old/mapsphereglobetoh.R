library(fmesher)
library(Matrix)
library(ggplot2)
library(tidyr)
library(dplyr)

# globe_vector <- c(1:100, seq(110, 200, by = 10), seq(250, 500, by =50), seq(600, 1000, by =100))
# h_min_vector <- numeric(length(globe_vector))
# h_max_vector <- numeric(length(globe_vector))
# for (i in 1:length(globe_vector)) {
#   globe <- globe_vector[i]
#   C_ii <- diag(fm_fem(fm_rcdt_2d(globe = globe), order = 2)$c0)
#   h_min_vector[i] <- sqrt(mean(C_ii))
#   rm(C_ii)
#   h_max_vector[i] <- 1/globe
#   print(i)
#   if(i>114){
#   saveRDS(
#     list(
#       globe_vector = globe_vector,
#       h_min_vector = h_min_vector,
#       h_max_vector = h_max_vector
#     ),
#     file = here::here("data_files/h_min_max_vs_globe.RDS")
#   )
#   }
# }
# 



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
  geom_point(size = 2) +
  labs(x = "Globe", y = "h value", color = "Type") +
  theme_minimal()

plotly::ggplotly(p)

closest_globe <- function(h) {
  idx <- which.min(abs(h_min_vector - h))  # index of closest value
  globe_vector[idx]  # return corresponding globe
}

# Example usage
closest_globe(0.012)


# Function to find globe for a given h via interpolation
globe_from_h <- function(h) {
  readed_data <- readRDS(
    file = here::here("data_files/h_min_max_vs_globe.RDS")
  )
  
  globe_vector <- readed_data$globe_vector
  h_min_vector <- readed_data$h_min_vector
  
  ord <- order(h_min_vector)
  h_sorted <- h_min_vector[ord]
  globe_sorted <- globe_vector[ord]
  
  return(round(approx(x = h_sorted, y = globe_sorted, xout = h)$y))
}

# Example usage
globe_from_h(h_min_vector[23])


# Function to get globe from h using spline
globe_from_h_spline <- function(h) {
  readed_data <- readRDS(
    file = here::here("data_files/h_min_max_vs_globe.RDS")
  )
  
  globe_vector <- readed_data$globe_vector
  h_min_vector <- readed_data$h_min_vector
  
  ord <- order(h_min_vector)
  h_sorted <- h_min_vector[ord]
  globe_sorted <- globe_vector[ord]
  
  # Create a smooth spline function: globe as function of h_max
  spline_func <- splinefun(x = h_sorted, 
                           y = globe_sorted, 
                           method = "monoH.FC")

  return(spline_func(h))
}
# Example usage
globe_from_h_spline(0.012)

hh <- 2^-seq(0,20, by = 1)
globe_hh <- globe_from_h_spline(hh)
globe_hh2 <- globe_from_h(hh)
plot(globe_vector, h_min_vector, col = "green", type = "l")
lines(globe_hh, hh, col = "red")
lines(globe_hh2, hh, col = "blue")

df_plot <- rbind(
  data.frame(globe = globe_vector, h = h_min_vector,
             curve = "h_min (green)"),
  data.frame(globe = globe_hh, h = hh,
             curve = "spline (red)"),
  data.frame(globe = globe_hh2, h = hh,
             curve = "direct (blue)")
)

# Plot
d <- ggplot(df_plot, aes(x = globe, y = h, color = curve)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1) +
  scale_color_manual(values = c(
    "h_min (green)" = "green",
    "spline (red)"  = "red",
    "direct (blue)" = "blue"
  )) +
  theme_minimal(base_size = 14) +
  labs(x = "globe", y = "h",
       color = "Curve")
plotly::ggplotly(d)
