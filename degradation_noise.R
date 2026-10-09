library(spatstat)

set.seed(42)

window <- square(10)

degradation_levels <- seq(0.01, 0.20, by = 0.01)
noise_levels <- c(0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4)
fp_levels <- c(0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4)
fn_levels <- c(0,0.05,0.1,0.3)

random <- rpoispp(lambda = 3, 
                      win = window)

matern <- rMatClust(
  kappa = 0.2,   
  scale = 0.5,  
  mu = 15,
  win = window
)

ssi <- rSSI(
  r = 0.4,   
  n = 300,      
  win = window,
  giveup = 300
)

par(mfrow = c(1, 3), mar = c(2, 2, 3, 1))

plot(random, 
     main = paste("Random (Poisson)\nN =", random$n), 
     cols = "darkblue",  
     pch = 16, 
     cex = 0.6)

plot(matern, 
     main = paste("Clustered (Matérn)\nN =", matern$n), 
     cols = "darkgreen", 
     pch = 16, 
     cex = 0.6)

plot(ssi, 
     main = paste("Inhibited (SSI)\nN =", ssi$n),    
     cols = "darkred",   
     pch = 16,
     cex = 0.6)


for (deg in degradation_levels) {
  random_degrade <- rthin(random, 1-deg)
  plot(random_degrade, 
       main = paste("Random (Poisson)\nN =", random_degrade$n, "Degradation Level = ", deg), 
       cols = "darkblue",  
       pch = 16, 
       cex = 0.6)
}

for (n_level in noise_levels) {
  random_noise <- rjitter(random, n_level)
  plot(random_noise, 
       main = paste("Random (Poisson)\nN =", random_noise$n, "Noise Level = ", n_level),
       cols = "darkblue",  
       pch = 16, 
       cex = 0.6)
}

#I'm confused what the difference is between false negatives and degradation 
#- both seem to involve randomly removing points.

for (fn_level in fn_levels) {
  random_fn <- rthin(random, 1-fn_level)
  plot(random_fn, 
       main = paste("Random (Poisson)\nN = ", random_fn$n, "False Negative Level = ", fn_level),
       cols = "darkblue",  
       pch = 16, 
       cex = 0.6)
}

for (fp_level in fp_levels) {
  num_fp <- round(random$n*fp_level)
  fp <- runif(num_fp)
  FP <- rpoispp(lambda = num_fp/100, win=window)
  random_fp <- superimpose(random, FP)
  plot(random_fp, 
       main = paste("Random (Poisson)\nN = ", random_fp$n, "False Positive Level = ", fp_level),
       cols = "darkblue",  
       pch = 16, 
       cex = 0.6)
}

#How do you add/remove points from cluster, regular patterns while keeping their spatial distribution?
