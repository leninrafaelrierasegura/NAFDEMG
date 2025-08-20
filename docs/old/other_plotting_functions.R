# this plots the two meshes in 3D in a static way
graph.plotter.3d.two.meshes.static <- function(graph_finer, graph_coarser, f_on_graph_finer, f_on_graph_coarser){
  x_finer <- plotting.order(graph_finer$mesh$V[, 1], graph_finer)
  y_finer <- plotting.order(graph_finer$mesh$V[, 2], graph_finer)
  x_coarser <- plotting.order(graph_coarser$mesh$V[, 1], graph_coarser)
  y_coarser <- plotting.order(graph_coarser$mesh$V[, 2], graph_coarser)
  z_finer <- plotting.order(f_on_graph_finer, graph_finer)
  z_coarser <- plotting.order(f_on_graph_coarser, graph_coarser)
  x_range <- range(c(x_finer, x_coarser), na.rm = TRUE)
  y_range <- range(c(y_finer, y_coarser), na.rm = TRUE)
  z_range <- range(c(z_finer, z_coarser), na.rm = TRUE)
  p <- plot_ly() %>%
    add_trace(x = x_finer, y = y_finer, z = z_finer,
              type = "scatter3d", mode = "lines",
              line = list(color = "red", width = 4),
              name = "finer", showlegend = TRUE) %>%
    add_trace(x = x_coarser, y = y_coarser, z = z_coarser,
              type = "scatter3d", mode = "lines",
              line = list(color = "blue", width = 4),
              name = "coarser", showlegend = TRUE) %>%
    add_trace(x = x_finer, y = y_finer, z = rep(0, length(x_finer)),
              type = "scatter3d", mode = "lines",
              line = list(color = "black", width = 4),
              name = "Graph", showlegend = FALSE)
  p <- layout(p,
              scene = global.scene.setter(x_range, y_range, z_range),
              xaxis = list(visible = FALSE),
              yaxis = list(visible = FALSE))
  return(p)
}
graph.plotter.3d.two.meshes.static(overkill_graph, graph, projected_U_approx[,1], U_approx[,1])

# This function plots the two meshes in 3D with time animation
graph.plotter.3d.two.meshes.time <- function(graph_finer, graph_coarser, 
                                             time_seq, frame_val_to_display,
                                             f_on_graph_finer, f_on_graph_coarser) {
  # Spatial coordinates (ordered for plotting)
  x_finer <- plotting.order(graph_finer$mesh$V[, 1], graph_finer)
  y_finer <- plotting.order(graph_finer$mesh$V[, 2], graph_finer)
  x_coarser <- plotting.order(graph_coarser$mesh$V[, 1], graph_coarser)
  y_coarser <- plotting.order(graph_coarser$mesh$V[, 2], graph_coarser)
  
  # Apply plotting.order to each column (time step)
  z_finer <- apply(f_on_graph_finer, 2, plotting.order, graph = graph_finer)
  z_coarser <- apply(f_on_graph_coarser, 2, plotting.order, graph = graph_coarser)
  
  # Create data frames for plotting
  n_time <- ncol(z_finer)
  
  data_finer <- data.frame(
    x = rep(x_finer, times = n_time),
    y = rep(y_finer, times = n_time),
    z = as.vector(z_finer),
    frame = rep(time_seq, each = length(x_finer)),
    mesh = "finer"
  )
  
  data_coarser <- data.frame(
    x = rep(x_coarser, times = n_time),
    y = rep(y_coarser, times = n_time),
    z = as.vector(z_coarser),
    frame = rep(time_seq, each = length(x_coarser)),
    mesh = "coarser"
  )
  
  data_graph <- data.frame(
    x = rep(x_finer, times = n_time),  # use finer mesh for the baseline
    y = rep(y_finer, times = n_time),
    z = 0,
    frame = rep(time_seq, each = length(x_finer)),
    mesh = "baseline"
  )
  
  # Ranges
  x_range <- range(c(x_finer, x_coarser))
  y_range <- range(c(y_finer, y_coarser))
  z_range <- range(c(z_finer, z_coarser))
  
  # Plotly
  p <- plot_ly(frame = ~frame) %>%
    add_trace(data = data_finer,
              x = ~x, y = ~y, z = ~z,
              type = "scatter3d", mode = "lines",
              line = list(color = "red", width = 3),
              name = "finer") %>%
    add_trace(data = data_coarser,
              x = ~x, y = ~y, z = ~z,
              type = "scatter3d", mode = "lines",
              line = list(color = "blue", width = 3),
              name = "coarser") %>%
    add_trace(data = data_graph,
              x = ~x, y = ~y, z = ~z,
              type = "scatter3d", mode = "lines",
              line = list(color = "black", width = 2),
              name = "Graph", showlegend = FALSE)
  
  frame_name <- deparse(substitute(frame_val_to_display))
  
  p <- p %>%
    layout(
      scene = global.scene.setter(x_range, y_range, z_range),
      updatemenus = list(list(
        type = "buttons", showactive = FALSE,
        buttons = list(
          list(label = "Play", method = "animate",
               args = list(NULL, list(frame = list(duration = 2000 / length(time_seq), redraw = TRUE),
                                      fromcurrent = TRUE))),
          list(label = "Pause", method = "animate",
               args = list(NULL, list(mode = "immediate", 
                                      frame = list(duration = 0), redraw = FALSE)))
        )
      )),
      title = paste0(frame_name, ": ", formatC(frame_val_to_display[1], format = "f", digits = 4))
    ) %>%
    plotly_build()
  
  # Update frame titles
  for (i in seq_along(p$x$frames)) {
    p$x$frames[[i]]$layout <- list(
      title = paste0(frame_name, ": ", formatC(frame_val_to_display[i], format = "f", digits = 4))
    )
  }
  
  return(p)
}
graph.plotter.3d.two.meshes.time(overkill_graph, graph, time_seq, time_seq, projected_U_approx, U_approx)

# this is the same as before but adds vertical lines

graph.plotter.3d.two.meshes.time <- function(graph_finer, graph_coarser, 
                                             time_seq, frame_val_to_display,
                                             f_on_graph_finer, f_on_graph_coarser) {
  # Spatial coordinates (ordered for plotting)
  x_finer <- plotting.order(graph_finer$mesh$V[, 1], graph_finer)
  y_finer <- plotting.order(graph_finer$mesh$V[, 2], graph_finer)
  x_coarser <- plotting.order(graph_coarser$mesh$V[, 1], graph_coarser)
  y_coarser <- plotting.order(graph_coarser$mesh$V[, 2], graph_coarser)
  
  # Apply plotting.order to each column (time step)
  z_finer <- apply(f_on_graph_finer, 2, plotting.order, graph = graph_finer)
  z_coarser <- apply(f_on_graph_coarser, 2, plotting.order, graph = graph_coarser)
  
  # Create data frames for plotting
  n_time <- ncol(z_finer)
  
  data_finer <- data.frame(
    x = rep(x_finer, times = n_time),
    y = rep(y_finer, times = n_time),
    z = as.vector(z_finer),
    frame = rep(time_seq, each = length(x_finer)),
    mesh = "finer"
  )
  
  data_coarser <- data.frame(
    x = rep(x_coarser, times = n_time),
    y = rep(y_coarser, times = n_time),
    z = as.vector(z_coarser),
    frame = rep(time_seq, each = length(x_coarser)),
    mesh = "coarser"
  )
  
  data_graph <- data.frame(
    x = rep(x_finer, times = n_time),  # baseline on finer mesh
    y = rep(y_finer, times = n_time),
    z = 0,
    frame = rep(time_seq, each = length(x_finer)),
    mesh = "baseline"
  )
  
  # --------- Vertical lines ----------
  vertical_lines <- function(x, y, z, frame_vals) {
    do.call(rbind, lapply(seq_along(frame_vals), function(i) {
      idx <- ((i - 1) * length(x) + 1):(i * length(x))
      data.frame(
        x = rep(x, each = 3),
        y = rep(y, each = 3),
        z = as.vector(t(cbind(0, z[idx], NA))),
        frame = rep(frame_vals[i], each = length(x) * 3)
      )
    }))
  }
  
  vertical_finer   <- vertical_lines(data_finer$x[1:length(x_finer)],
                                     data_finer$y[1:length(y_finer)],
                                     data_finer$z,
                                     time_seq)
  
  vertical_coarser <- vertical_lines(data_coarser$x[1:length(x_coarser)],
                                     data_coarser$y[1:length(y_coarser)],
                                     data_coarser$z,
                                     time_seq)
  
  # Ranges
  x_range <- range(c(x_finer, x_coarser))
  y_range <- range(c(y_finer, y_coarser))
  z_range <- range(c(z_finer, z_coarser))
  
  # --------- Plotly object ----------
  p <- plot_ly(frame = ~frame) %>%
    add_trace(data = data_finer,
              x = ~x, y = ~y, z = ~z,
              type = "scatter3d", mode = "lines",
              line = list(color = "red", width = 3),
              name = "finer") %>%
    add_trace(data = data_coarser,
              x = ~x, y = ~y, z = ~z,
              type = "scatter3d", mode = "lines",
              line = list(color = "blue", width = 3),
              name = "coarser") %>%
    add_trace(data = data_graph,
              x = ~x, y = ~y, z = ~z,
              type = "scatter3d", mode = "lines",
              line = list(color = "black", width = 2),
              name = "Graph", showlegend = FALSE) %>%
    # Vertical lines finer
    add_trace(data = vertical_finer,
              x = ~x, y = ~y, z = ~z, frame = ~frame,
              type = "scatter3d", mode = "lines",
              line = list(color = "gray", width = 0.5),
              name = "Vertical finer", showlegend = FALSE) %>%
    # Vertical lines coarser
    add_trace(data = vertical_coarser,
              x = ~x, y = ~y, z = ~z, frame = ~frame,
              type = "scatter3d", mode = "lines",
              line = list(color = "darkgray", width = 0.5),
              name = "Vertical coarser", showlegend = FALSE)
  
  frame_name <- deparse(substitute(frame_val_to_display))
  
  p <- p %>%
    layout(
      scene = global.scene.setter(x_range, y_range, z_range),
      updatemenus = list(list(
        type = "buttons", showactive = FALSE,
        buttons = list(
          list(label = "Play", method = "animate",
               args = list(NULL, list(frame = list(duration = 2000 / length(time_seq), redraw = TRUE),
                                      fromcurrent = TRUE))),
          list(label = "Pause", method = "animate",
               args = list(NULL, list(mode = "immediate", 
                                      frame = list(duration = 0), redraw = FALSE)))
        )
      )),
      title = paste0(frame_name, ": ", formatC(frame_val_to_display[1], format = "f", digits = 4))
    ) %>%
    plotly_build()
  
  # Update frame titles
  for (i in seq_along(p$x$frames)) {
    p$x$frames[[i]]$layout <- list(
      title = paste0(frame_name, ": ", formatC(frame_val_to_display[i], format = "f", digits = 4))
    )
  }
  
  return(p)
}










