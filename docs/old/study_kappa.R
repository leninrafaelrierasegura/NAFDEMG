# KAPPAS <- c(1,2,4,8,16,32)
# for (i in 1:length(KAPPAS)) {
#   kappa <- KAPPAS[i]
#   saveRDS(kappa, "old/kappa.RDS")
#   rmarkdown::render("conv_in_h.Rmd")
#   rmarkdown::render("conv_in_tau.Rmd")
#   rmarkdown::render("conv_in_m.Rmd")
# }


# KAPPAS <- c(1,2,4,8,16,32)
# DIVIDER <- c(1,2,4,6,8,10)
# for (i in 1:length(KAPPAS)) {
#   kappa <- KAPPAS[i]
#   saveRDS(kappa, "old/kappa.RDS")
#   for (j in 1:length(DIVIDER)) {
#     divider <- DIVIDER[j]
#     saveRDS(divider, "old/divider.RDS")
#     rmarkdown::render("control_conv_m.Rmd")
#   }
# }

rmarkdown::render("conv_in_m.Rmd")
KAPPAS <- c(1,2,4,8,16,32)
DIVIDER <- c(1,2,4,6,8,10)
for (i in 1:length(KAPPAS)) {
  kappa <- KAPPAS[i]
  saveRDS(kappa, "old/kappa.RDS")
  for (j in 1:length(DIVIDER)) {
    divider <- DIVIDER[j]
    saveRDS(divider, "old/divider.RDS")
    rmarkdown::render("control_conv_tau.Rmd")
  }
}

# for (i in 1:length(KAPPAS)) {
#   kappa <- KAPPAS[i]
#   saveRDS(kappa, "old/kappa.RDS")
#   for (j in 1:length(DIVIDER)) {
#     divider <- DIVIDER[j]
#     saveRDS(divider, "old/divider.RDS")
#     rmarkdown::render("control_conv_h.Rmd")
#   }
# }



rmarkdown::render("control_conv_h_with_iteration.Rmd")
rmarkdown::render("control_conv_m_with_iteration.Rmd")

KAPPAS <- c(4,8,16,32)
DIVIDER <- c(1,2,4,8,16)
for (i in 1:length(KAPPAS)) {
  kappa <- KAPPAS[i]
  saveRDS(kappa, "old/kappa.RDS")
  for (j in 1:length(DIVIDER)) {
    divider <- DIVIDER[j]
    saveRDS(divider, "old/divider.RDS")
    #rmarkdown::render("control_conv_h_with_iteration.Rmd")
    #rmarkdown::render("control_conv_tau_with_iteration.Rmd")
    rmarkdown::render("control_conv_m_with_iteration.Rmd")
  }
}


DIVIDER <- c(1,2,4,8,16)
for (j in 1:length(DIVIDER)) {
  divider <- DIVIDER[j]
  saveRDS(divider, "old/divider.RDS")
  rmarkdown::render("control_conv_m_with_iteration.Rmd")
}



M <- matrix(1:16,4,4)
M
matrix(pmax(4, pmin(10, M)), dim(M))


MS <- c(4,5,6,7,8)
for (j in 1:length(MS)) {
  m_cut <- MS[j]
  saveRDS(m_cut, "old/m_cut.RDS")
  rmarkdown::render("control_conv_h_with_iteration.Rmd")
  rmarkdown::render("control_conv_tau_with_iteration.Rmd")
}

rmarkdown::render("control_conv_tau_with_iteration.Rmd")
rmarkdown::render("control_conv_h_with_iteration.Rmd")




points <- 10^seq(log10(0.005), log10(0.1), length.out = 7)
points

v <- c(0.005, 0.0075, 0.01, 0.025, 0.05, 0.075, 0.1)
v <- seq(from = 0.005, to = 0.1, length.out = 7)
plot(log(points), rep(1, 7), xlab="Alpha", ylab="", yaxt="n", pch=19, cex=2)


