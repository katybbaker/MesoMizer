library(multisensi)
library(sensitivity)
library(mizer)

setwd("C:/Users/Jackson Lab/Desktop/PhD")

#Load mizer details
species_params2 <- read.csv("mizer_params_new2.csv", header=TRUE, sep =",")
kappa <- 10
species_params2$R_max <- kappa * species_params2$w_inf^(-1.5)
species_params2$R_max <- species_params2$R_max * 100

#interaction matrix
inter <- read.csv("int_ls_new.csv", header = TRUE, sep = ",", row.names = NULL, stringsAsFactors = TRUE )
inter <- inter[-c(1)]
colnames(inter) <- c("Common Shiner", "Cisco", "Freshwater Drum", "Northern Pike", "Trout-perch", "Walleye", "Yellow Perch")
rownames(inter) <- colnames(inter)
inter <- as.matrix(inter)


m <- 10000

Xb <- data.frame(k_vb.6=runif(m,min=0.1818, max = 0.2222), w_inf.4 = runif(m,min=10023.99, max = 12251.55), k_vb.1 = runif(m, min = 0.5184, max = 0.6336),
                 k_vb.5 = runif(m, min = 0.39474, max = 0.48246), w_inf.1 = runif(m, min = 36.9, max = 45.1), R_max.5 = runif(m, min = 3.428205, max = 4.190028),
                 w_inf.6 = runif(m, min = 1927.233, max = 2355.507), w_inf.2 = runif(m, min = 636.5232, max = 777.9728), k_vb.4 = runif(m, min = 0.0918, max = 0.1122),
                 w_inf.5 = runif(m, min = 36.9, max = 45.1))

runsobol <- function(k_vb.6, w_inf.4, k_vb.1, k_vb.5, w_inf.1, R_max.5, w_inf.6, w_inf.2, k_vb.4, w_inf.5){
  species_params2$k_vb[6] <- k_vb.6
  species_params2$w_inf[4] <- w_inf.4
  species_params2$k_vb[1] <- k_vb.1
  species_params2$k_vb[5] <- k_vb.5
  species_params2$w_inf[1] <- w_inf.1
  species_params2$R_max[5] <- R_max.5
  species_params2$w_inf[6] <- w_inf.6
  species_params2$w_inf[2] <- w_inf.2
  species_params2$k_vb[4] <- k_vb.4
  species_params2$w_inf[5] <- w_inf.5
  
  params <- MizerParams(species_params = species_params2, interaction = inter, kappa = 10, w_pp_cutoff = 2)
  sim <- project(params, effort = 0, t_max=500, dt=0.1, t_save=1) 
  output <- getCommunitySlope(sim)[500,,]$slope
  return(output)
}


runsobol2 <- function(X){
  out <- matrix(nrow = nrow(X), ncol = 1, NA)
  for (i in 1:nrow(X)){
    out[i,] <- runsobol(X$k_vb.6[i], X$w_inf.4[i], X$k_vb.1[i], X$k_vb.5[i], X$w_inf.1[i], X$R_max.5[i], X$w_inf.6[i], X$w_inf.2[i], X$k_vb.4[i], X$w_inf.5[i])
  }
  out <- as.data.frame(out)
  return(out)
}

seq.sobol <- multisensi(design = sobol2007, model = runsobol2,
                                 reduction = NULL, analysis = analysis.sensitivity, center = TRUE,
                                 design.args = list(X1 = Xb[1:(m/2), ], X2 = Xb[(1 + m/2):m, ], nboot = 100),
                                 analysis.args = list(keep.outputs = FALSE))
