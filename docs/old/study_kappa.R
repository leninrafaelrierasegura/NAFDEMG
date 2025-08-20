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

for (i in 1:length(KAPPAS)) {
  kappa <- KAPPAS[i]
  saveRDS(kappa, "old/kappa.RDS")
  for (j in 1:length(DIVIDER)) {
    divider <- DIVIDER[j]
    saveRDS(divider, "old/divider.RDS")
    rmarkdown::render("control_conv_h.Rmd")
  }
}