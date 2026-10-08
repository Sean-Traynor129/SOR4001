library(spatstat)

set.seed(42)

window <- square(10)

pts_random <- rpoispp(lambda = 3, 
                      win = win)

pts_matern <- rMatClust(
  kappa = 0.2,   
  scale = 0.5,  
  mu = 15,
  win = window
)

pts_ssi <- rSSI(
  r = 0.4,   
  n = 300,      
  win = window,
  giveup = 300
)

plot(pts_random, 
     main = paste("Random (Poisson)\nN =", pts_random$n), 
     cols = "darkblue",  
     pch = 16, 
     cex = 0.6)

plot(pts_matern, 
     main = paste("Clustered (Matérn)\nN =", pts_matern$n), 
     cols = "darkgreen", 
     pch = 16, 
     cex = 0.6)

plot(pts_ssi, 
     main = paste("Inhibited (SSI)\nN =", pts_ssi$n),    
     cols = "darkred",   
     pch = 16,
     cex = 0.6)

help("rMatClust")
