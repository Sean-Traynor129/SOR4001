library(spatstat)

# 1. Define observation window (10x10 = Area 100)
win <- square(10)

set.seed(42)

# 2. Random (Homogeneous Poisson Process)
# Expected points = lambda * Area = 3.0 * 100 = 300
pts_random <- rpoispp(lambda = 3.0, win = win)

# 3. Matérn Cluster Process
# Expected points = kappa * Area * mu = 0.2 * 100 * 15 = 300
pts_matern <- rMatClust(
  kappa = 0.17,   # Parent intensity (~17 parent clusters)
  scale = 0.5,   # Cluster radius
  mu = 15,       # Mean offspring per parent
  win = win
)

# 4. Simple Sequential Inhibition (SSI) Process
# Set r = 0.35 to allow 300 points to pack into a 10x10 space
pts_ssi <- rSSI(
  r = 0.35,      # Inhibition radius
  n = 300,       # Exact target point count
  win = win,
  giveup = 300  # Max iterations to attempt placement
)

# 5. Visual Comparison
par(mfrow = c(1, 3), mar = c(2, 2, 3, 1))

plot(pts_random, main = paste("Random (Poisson)\nN =", pts_random$n), cols = "darkblue",  pch = 16, cex = 0.6)
plot(pts_matern, main = paste("Clustered (Matérn)\nN =", pts_matern$n), cols = "darkgreen", pch = 16, cex = 0.6)
plot(pts_ssi,    main = paste("Inhibited (SSI)\nN =", pts_ssi$n),    cols = "darkred",   pch = 16, cex = 0.6)

# Reset graphics layout
par(mfrow = c(1, 1))

degraded_poisson <- rthin(pts_random, 0.99)
plot(degraded_poisson)

library(spatstat)

# 1. Structural Degradation (Cell Loss) using native rthin()
# P is the retention probability (1 - degradation rate)
degrade_pattern <- function(X, prop_remove) {
  if (prop_remove <= 0) return(X)
  rthin(X, P = 1 - prop_remove)[cite: 3, 16]
}

# 2. Positional Noise using native rjitter()
add_positional_noise <- function(X, sd_noise) {
  if (sd_noise <= 0) return(X)
  # rjitter applies random displacements to point coordinates
  rjitter(X, radius = sd_noise, distribution = "gaussian", retry = FALSE)[cite: 5, 16]
}

# 3. False Positives & Negatives using rthin() and superimpose()
add_fp_fn_noise <- function(X, fp_rate, fn_rate) {
  # Apply False Negatives (random removal via rthin)
  X_fn <- if (fn_rate > 0) rthin(X, P = 1 - fn_rate) else X[cite: 6, 16]
  
  # Apply False Positives (adding random points sampled from CSR)
  n_fp <- round(X$n * fp_rate)[cite: 6]
  if (n_fp > 0) {
    X_fp <- runifpoint(n = n_fp, win = Window(X))[cite: 6]
    return(superimpose(X_fn, X_fp))[cite: 6]
  }
  return(X_fn)
}