library(MultiView.LOCUS)
library(testthat)



test_that("Multi View decomposition simulation runs", {
  skip_if_not_installed("MultiView.LOCUS")
  system.file("data", package = "MultiView.LOCUS")
  set.seed(111)
  ICcorr<- readRDS(file.path(data_path, "ICcorr.rds"))
  S_x <- readRDS(file.path(data_path, "S_xreal.rds"))
  S_y <- readRDS(file.path(data_path, "S_yreal.rds"))
  noise_level = 0.005
  correlation_level = 0.4
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
    psi = 0.5,
    gamma = 2.1,
    rho = 0.95,
    espli1 = 1e-2,
    espli2 = 1e-2,
    MaxIteration = 500
  )

## Visualization of connectivity traits
library(gridExtra)  # for arranging plots

# Helper: plot a connectivity matrix from vector
plot_conn_wrapper <- function(vec, V) {
  mat <- Ltrinv(vec, V, d = FALSE)
  plot_conn(mat)
}

# Initialize a list to hold plots
plot_list <- list()

# Settings
V <- 264   # number of nodes
n_common <- 2   # number of common components
n_spe <- 2      # number of specific components per view

# View 1
for (i in 1:(n_common + n_spe)) {
  plot_list[[i]] <- plot_conn_wrapper(res$S_sparse[[1]][i, ], V)
}

# View 2
for (i in 1:(n_common + n_spe)) {
  plot_list[[i + 4]] <- plot_conn_wrapper(res$S_sparse[[2]][i, ], V)
}

# Arrange the plots: 2 rows × 4 columns
combined_plot = grid.arrange(grobs = plot_list, nrow = 2, ncol = 4)
#ggsave(
  filename = "multi_view_connectivity.png",
  plot = combined_plot,
  width = 4 * 4,   # 4 plots wide × 4 inches each
  height = 2 * 4,  # 2 rows × 4 inches each
  dpi = 300
)

})
