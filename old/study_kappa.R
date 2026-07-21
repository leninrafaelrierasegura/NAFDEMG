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


SIGMAS <- seq(from = 0.1, to = 4, by = 0.1)
RANGES <- seq(from = 0.1, to = 5, by = 0.2)
NUS <- seq(from = 0.1, to = 2.5, by = 0.1)
SIGMA.ES <- seq(from = 0.1, to = 5, by = 0.1)
N.REPS <- seq(from = 1, to = 10, by = 1)
for (i in 1:length(SIGMAS)) {
  sigma <- SIGMAS[i]
  saveRDS(sigma, "old/sigma.RDS")
  for (j in 1:length(RANGES)) {
    range <- RANGES[j]
    saveRDS(range, "old/range.RDS")
    for (k in 1:length(NUS)) {
      nu <- NUS[k]
      saveRDS(nu, "old/nu.RDS")
      for (l in 1:length(SIGMA.ES)) {
        sigma.e <- SIGMA.ES[l]
        saveRDS(sigma.e, "old/sigma_e.RDS")
        for (m in 1:length(N.REPS)) {
          n.rep <- N.REPS[m]
          saveRDS(n.rep, "old/n_rep.RDS")
          rmarkdown::render(here::here("old/Paper1_simulation.Rmd"))
        }
      }
    }
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
rmarkdown::render("control_conv_m_with_iteration.Rmd")
rmarkdown::render("control_conv_h_with_iteration.Rmd")



points <- 10^seq(log10(0.005), log10(0.1), length.out = 7)
points

v <- c(0.005, 0.0075, 0.01, 0.025, 0.05, 0.075, 0.1)
v <- seq(from = 0.005, to = 0.1, length.out = 7)
plot(log(points), rep(1, 7), xlab="Alpha", ylab="", yaxt="n", pch=19, cex=2)


rmarkdown::render("new_conv_in_h.Rmd")
rmarkdown::render("new_conv_in_tau.Rmd")
rmarkdown::render("new_conv_in_m.Rmd")


knitr::purl(here::here("conv_in_m_with_logtau_and_tauminustwo.Rmd"),
            output = here::here("old/conv_in_m_with_logtau_and_tauminustwo.R"))

