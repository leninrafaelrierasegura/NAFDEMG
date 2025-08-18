KAPPAS <- c(1,2,4,8,16,32)
for (i in 1:length(KAPPAS)) {
  kappa <- KAPPAS[i]
  saveRDS(kappa, "old/kappa.RDS")
  rmarkdown::render("conv_in_h.Rmd")
  rmarkdown::render("conv_in_tau.Rmd")
  rmarkdown::render("conv_in_m.Rmd")
}