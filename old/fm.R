library(fmesher)

kappa <- 4
a <- 100
b <- 100

nx <- 10
ny <- 10

alpha <- 1.8
m = 4
beta <- alpha/2

T_final <- 0.5
time_step <- 0.001
time_seq <- seq(0, T_final, length.out = ((T_final - 0) / time_step + 1))
coeff <- 1
which_eig <- 9



FEM_matrices_coarse <- fem$fem_shifted_laplacian_sparse(a = a, b = b, nx = nx, ny = ny)
FEM_matrices_fine <-   fem$fem_shifted_laplacian_sparse(a = a, b = b, nx = 2*nx, ny = 2*ny)

mesh_coarse <- fm_mesh_2d(loc = FEM_matrices_coarse$mesh_loc, offset = 0, min.angle = 40)
mesh_fine <- fm_mesh_2d(loc = FEM_matrices_fine$mesh_loc, offset = 0, min.angle = 40)

true_eig <- compute_true_eigen_rectangle(
  a = a, 
  b = b, 
  loc = mesh_fine$loc,
  kappa = kappa,
  m_max = 5, 
  n_max = 5
)

eigvals <- true_eig$eigvals
eigfuncs <- true_eig$eigfuncs


U_0 <- coeff * eigfuncs[, which_eig]
U_true_matrix <- outer(U_0, exp(-(eigvals[which_eig]^(alpha/2)) * time_seq), FUN = "*")



plot_3d_slider_scatter(mesh_fine$loc, time_seq, U_true_matrix)


FEM_mat_coarse <- fm_fem(mesh_coarse, order = 2)
FEM_mat_fine <- fm_fem(mesh_fine, order = 2)

A_mat <- fm_basis(mesh_coarse, mesh_fine$loc)

Cl <- FEM_mat_coarse$c0
C <- FEM_mat_coarse$c1
G <- FEM_mat_coarse$g1
L <- kappa^2 * C + G

true_eig <- compute_true_eigen_rectangle(
  a = a, 
  b = b, 
  loc = mesh_coarse$loc,
  kappa = kappa,
  m_max = 5, 
  n_max = 5
)

eigvals <- true_eig$eigvals
eigfuncs <- true_eig$eigfuncs


U_0 <- coeff * eigfuncs[, which_eig]

my_op_frac <- my.fractional.operators.frac(L, beta, C, scale.factor = kappa^2, m = m, time_step)
U_approx_matrix <- solve_fractional_evolution(my_op_frac, time_step, time_seq, val_at_0 = U_0)

U_approx_fine_matrix <- A_mat %*% U_approx_matrix

diff <- U_true_matrix - U_approx_fine_matrix

plot_3d_slider_scatter(mesh_fine$loc, time_seq, diff)




# to plot on the sphere

f_vals <- eigfucs[,31]

# Compute normal vectors (for bumps)
normals <- cbind(x, y, z)
scale <- 0.1
x_bump <- x + scale * f_vals * normals[,1]
y_bump <- y + scale * f_vals * normals[,2]
z_bump <- z + scale * f_vals * normals[,3]

# Plotly
fig <- plot_ly() %>%
  # Original sphere
  add_mesh(
    x = x, y = y, z = z,
    i = tri[,1] - 1,
    j = tri[,2] - 1,
    k = tri[,3] - 1,
    color = I("lightblue"),
    opacity = 0.5,      # transparent
    name = "Sphere"
  ) %>%
  # Bumpy surface
  add_mesh(
    x = x_bump, y = y_bump, z = z_bump,
    i = tri[,1] - 1,
    j = tri[,2] - 1,
    k = tri[,3] - 1,
    intensity = f_vals,
    colorscale = "Viridis",
    flatshading = TRUE,
    opacity = 0.7,      # slightly transparent
    name = "Function surface"
  )

fig



plot_3d_slider_sphere_scatter <- function(mesh, eigvals, eigfuncs) {
  
  colorscale = "Viridis"
  
  x <- mesh$loc[,1]
  y <- mesh$loc[,2]
  z <- mesh$loc[,3]
  
  # Create frames for SCATTER
  frames <- lapply(seq_len(ncol(eigfuncs)), function(i) {
    list(
      name = as.character(i),
      data = list(list(
        x = x,
        y = y,
        z = z,
        type = "scatter3d",
        mode = "markers",
        marker = list(
          size = 5,
          color = eigfuncs[,i],
          colorscale = colorscale,
          showscale = TRUE
        )
      ))
    )
  })
  
  # Initial plot (frame 1)
  p <- plot_ly(
    x = x, y = y, z = z,
    type = "scatter3d",
    mode = "markers",
    marker = list(
      size = 5,
      color = eigfuncs[,1],
      colorscale = colorscale,
      showscale = TRUE
    ),
    frame = "1"
  )
  
  frame_name <- deparse(substitute(eigvals))
  
  p$x$frames <- frames
  
  # Layout and slider
  p <- p %>% layout(
    title = paste0(frame_name, ": ", eigvals[1]),
    sliders = list(
      list(
        active = 0,
        currentvalue = list(prefix = "Mode: "),
        pad = list(t = 50),
        steps = lapply(seq_len(ncol(eigfuncs)), function(i) {
          list(
            label = as.character(i),
            method = "animate",
            args = list(list(as.character(i)),
                        list(mode = "immediate",
                             frame = list(duration = 300, redraw = TRUE),
                             transition = list(duration = 0)))
          )
        })
      )
    ),
    updatemenus = list(
      list(
        type = "buttons",
        showactive = FALSE,
        y = 1,
        x = 1.15,
        xanchor = "right",
        yanchor = "top",
        buttons = list(
          list(label = "Play",
               method = "animate",
               args = list(NULL, list(frame = list(duration = 300, redraw = TRUE),
                                      fromcurrent = TRUE, mode = "immediate"))),
          list(label = "Pause",
               method = "animate",
               args = list(NULL, list(frame = list(duration = 0, redraw = FALSE),
                                      mode = "immediate")))
        )
      )
    )
  ) %>% plotly_build()
  
  # Update title in each frame
  for (i in seq_len(ncol(eigfuncs))) {
    p$x$frames[[i]]$layout <- list(
      title = paste0(frame_name, ": ", eigvals[i])
    )
  }
  
  return(p)
}



plot_3d_slider_sphere <- function(mesh, eigvals, eigfuncs) {
  
  colorscale = "Viridis"
  opacity = 1
  
  
  x <- mesh$loc[, 1]
  y <- mesh$loc[, 2]
  z <- mesh$loc[, 3]
  tri <- mesh$graph$tv
  
  
  # Create frames
  frames <- lapply(seq_len(ncol(eigfuncs)), function(i) {
    list(
      name = as.character(i),
      data = list(list(
        x = x,
        y = y,
        z = z,
        i = tri[,1] - 1,
        j = tri[,2] - 1,
        k = tri[,3] - 1,
        type = "mesh3d",
        intensity = eigfuncs[, i],
        colorscale = colorscale,
        opacity = opacity,
        flatshading = TRUE
      ))
    )
  })
  
  # Initial plot
  p <- plot_ly(
    x = x, y = y, z = z,
    i = tri[,1] - 1,
    j = tri[,2] - 1,
    k = tri[,3] - 1,
    type = "mesh3d",
    intensity = eigfuncs[,1],
    colorscale = colorscale,
    opacity = opacity,
    flatshading = TRUE,
    frame = "1"
  )
  
  frame_name <- deparse(substitute(eigvals))
  
  p$x$frames <- frames
  
  # Layout with slider
  p <- p %>% layout(
    title = paste0(frame_name, ": ", eigvals[1]),
    sliders = list(
      list(
        active = 0,
        currentvalue = list(prefix = "Frame: "),
        pad = list(t = 50),
        steps = lapply(seq_len(ncol(eigfuncs)), function(i) {
          list(
            label = as.character(i),
            method = "animate",
            args = list(list(as.character(i)),
                        list(mode = "immediate",
                             frame = list(duration = 300, redraw = TRUE),
                             transition = list(duration = 0)))
          )
        })
      )
    ),
    updatemenus = list(
      list(
        type = "buttons",
        showactive = FALSE,
        y = 1,
        x = 1.15,
        xanchor = "right",
        yanchor = "top",
        buttons = list(
          list(label = "Play",
               method = "animate",
               args = list(NULL, list(frame = list(duration = 300, redraw = TRUE),
                                      fromcurrent = TRUE, mode = "immediate"))),
          list(label = "Pause",
               method = "animate",
               args = list(NULL, list(frame = list(duration = 0, redraw = FALSE),
                                      mode = "immediate")))
        )
      )
    )
  ) %>% plotly_build()
  
  # Add titles for frames
  for (i in seq_len(ncol(eigfuncs))) {
    p$x$frames[[i]]$layout <- list(
      title = paste0(frame_name, ": ", eigvals[i])
    )
  }
  
  return(p)
}
