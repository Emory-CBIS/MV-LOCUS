library(MultiView.LOCUS)
library(testthat)



test_that("Multi View decomposition simulation runs", {
  skip_if_not_installed("MultiView.LOCUS")
  ICcorr<- readRDS(file.path(data_path, "ICcorr.rds"))
  S_x <- readRDS(file.path(data_path, "S_xreal.rds"))
  S_y <- readRDS(file.path(data_path, "S_yreal.rds"))

  sample1 <- sample(1:237, 100)
  mixing_cor_commonx <- t(ICcorr$M[sample1, c(5, 6)]) / apply(ICcorr$M[sample1, c(5, 6)], 2, sd)
  mixing_cor_commony <- mixing_cor_commonx + matrix(rnorm(200, 0, sqrt(1 / correlation_level^2 - 1)), nrow = 2)

  mixing_x <- t(rbind(mixing_cor_commonx, t(ICcorr$M[sample1, 8:9])))
  mixing_y <- t(rbind(mixing_cor_commony, t(ICcorr$M[sample1, 10:11])))

  noise <- noise_level
  Data_x <- mixing_x %*% S_x + matrix(rnorm(264 * 263 / 2 * 100, sd = noise), nrow = 100)
  Data_y <- mixing_y %*% S_y + matrix(rnorm(264 * 263 / 2 * 100, sd = noise), nrow = 100)

  res <- multi_view_decomposition(
    Y = list(Data_x, Data_y),
    q_common = 2,
    q = c(2, 2),
    V = 264,
    penalty = "SCAD",
    phi = 2,
    psi = 1,
    gamma = 2.1,
    rho = 0.95,
    espli1 = 1e-2,
    espli2 = 1e-2,
    MaxIteration = 500
  )
plot_conn(Ltrinv(res$S[[1]][4,],264,F))

})
