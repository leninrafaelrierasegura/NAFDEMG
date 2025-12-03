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


