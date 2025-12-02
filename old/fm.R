library(fmesher)

mesh <- fm_mesh_2d(loc = FEM_matrices$coords, offset = 0)
plot(mesh)
FEM_mat <- fm_fem(mesh, order = 2)

C_lumped <- FEM_mat$c0
C <- FEM_mat$c1
G <- FEM_mat$g1
