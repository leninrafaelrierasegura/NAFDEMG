library(fmesher)

kappa <- 1
a <- 20
b <- 10

nx <- 2
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

mesh_coarse <- fm_mesh_2d(loc = FEM_matrices_coarse$mesh_loc, offset = 0)
mesh_fine <- fm_mesh_2d(loc = FEM_matrices_fine$mesh_loc, offset = 0)

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



x <- mesh_fine$loc[,1]
y <- mesh_fine$loc[,2]
z <- eigfuncs[,9]
plotly::plot_ly(
  x = x,
  y = y,
  z = z,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 5, color = z, colorscale = "Viridis")
)



plot(mesh_coarse)
plot(mesh_fine)


FEM_mat_coarse <- fm_fem(mesh_coarse, order = 2)
FEM_mat_fine <- fm_fem(mesh_fine, order = 2)

A_mat <- fm_basis(mesh_coarse, mesh_fine$loc)

Cl <- FEM_mat_coarse$c0
C <- FEM_mat_coarse$c1
G <- FEM_mat_coarse$g1
L <- kappa^2 * C + G


