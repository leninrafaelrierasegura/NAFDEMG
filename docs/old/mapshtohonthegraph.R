library(fmesher)
library(Matrix)
library(ggplot2)
library(tidyr)
library(dplyr)
library(MetricGraph)

gets.graph.tadpole <- function(flip_edge = FALSE){
  if(flip_edge) {
    edge1 <- rbind(c(0,0),c(1,0))[c(2,1),]
  } else {
    edge1 <- rbind(c(0,0),c(1,0))}
  theta <- seq(from=-pi,to=pi,length.out = 10000)
  edge2 <- cbind(1+1/pi+cos(theta)/pi,sin(theta)/pi)
  edges <- list(edge1, edge2)
  graph <- metric_graph$new(edges = edges, verbose = 0)
  graph$set_manual_edge_lengths(edge_lengths = c(1,2))
  #graph$build_mesh(h = h)
  return(graph)
}

graph <- gets.graph.tadpole(flip_edge = FALSE)


h_vector <-2^(-(10:1))
h_mih_vector <- numeric(length(h_vector))
h_max_vector <- numeric(length(h_vector))
for (i in 1:length(h_vector)) {
  h <- h_vector[i]
  graph$build_mesh(h = h)
  graph$compute_fem()
  h_mih_vector[i] <- mean(rowSums(graph$mesh$C))
  h_max_vector[i] <- mean(rowSums(graph$mesh$C))
  print(i)
}


# Create a data frame
df <- data.frame(
  index = seq_along(h_vector),
  h_input = h_vector,
  h_min = h_mih_vector,
  h_max = h_max_vector
)

# Convert to long format for ggplot
df_long <- df %>%
  pivot_longer(cols = c(h_input, h_min, h_max),
               names_to = "type",
               values_to = "h_value")

# Plot
p <- ggplot(df_long, aes(x = index, y = h_value, color = type)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("blue", "green", "red")) +
  labs(x = "Index", y = "h value", color = "Type",
       title = "Mesh spacing: h, h_min, h_max") +
  theme_minimal(base_size = 14)

plotly::ggplotly(p)
