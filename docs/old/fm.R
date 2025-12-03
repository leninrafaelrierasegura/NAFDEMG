library(fmesher)

a <- 100
b <- 100

nx <- 10
ny <- 10

FEM_matrices_coarse <- fem$fem_shifted_laplacian_sparse(
  a = a,
  b = b, 
  nx = nx, 
  ny = ny
)

FEM_matrices_fine <- fem$fem_shifted_laplacian_sparse(
  a = a,
  b = b, 
  nx = 2*nx, 
  ny = 2*ny
)


mesh_coarse <- fm_mesh_2d(loc = FEM_matrices_coarse$mesh_loc, offset = 0)
mesh_fine <- fm_mesh_2d(loc = FEM_matrices_fine$mesh_loc, offset = 0)



plot(mesh_coarse)
plot(mesh_fine)


FEM_mat_coarse <- fm_fem(mesh_coarse, order = 2)
FEM_mat_fine <- fm_fem(mesh_fine, order = 2)

A_mat <- fm_basis(mesh_coarse, mesh_fine$loc)

Cl <- FEM_mat_coarse$c0
C <- FEM_mat_coarse$c1
G <- FEM_mat_coarse$g1
L <- kappa^2 * C + G


