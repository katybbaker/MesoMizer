library(multisensi)
library(sensitivity)
library(mizer)

#Load mizer details
species_params2 <- read.csv("species_params.csv", header=TRUE, sep =",")
# kappa <- 10
# species_params2$R_max <- kappa * species_params2$w_inf^(-1.5)
# species_params2$R_max <- species_params2$R_max * 100

#interaction matrix
inter <- read.csv("int_ls_new.csv", header = TRUE, sep = ",", row.names = NULL, stringsAsFactors = TRUE )
inter <- inter[-c(1)]
colnames(inter) <- c("Common Shiner", "Cisco", "Freshwater Drum", "Northern Pike", "Trout-perch", "Walleye", "Yellow Perch")
rownames(inter) <- colnames(inter)
inter <- as.matrix(inter)


m <- 500

Xb <- data.frame(w_inf.4 = runif(m, min = 10023.99, max = 12251.55), w_inf.6 = runif(m, min= 1927.233 , max= 2355.507), k_vb.5 =runif(m, min= 0.39474, max = 0.48246),
                 k_vb.4 = runif(m, min = 0.0918, max = 0.1122), R_max.1 = runif(m, min=3.428205, max = 4.190028), k_vb.6 = runif(m, min= 0.1818, max=0.2222),
                 beta.6 = runif(m, min = 90, max =110), w_inf.3 = runif(m, min = 2742.03, max=3351.37), R_max.4 = runif(m, min = 0.000766, max= 0.000936),
                 k_vb.3 = runif(m, min =0.0603 , max = 0.0737))

runsobol <- function(w_inf.4, w_inf.6, k_vb.5, k_vb.4, R_max.1, k_vb.6, beta.6, w_inf.3, R_max.4, k_vb.3){
  species_params2$w_inf[4] <- w_inf.4
  species_params2$w_inf[6] <- w_inf.6
  species_params2$k_vb[5] <- k_vb.5
  species_params2$k_vb[4] <- k_vb.4
  species_params2$R_max[1] <- R_max.1
  species_params2$k_vb[6] <- k_vb.6
  species_params2$beta[6] <- beta.6
  species_params2$w_inf[3] <- w_inf.3
  species_params2$R_max[4] <- R_max.4
  species_params2$k_vb[3] <- k_vb.3
  
  params <- MizerParams(species_params = species_params2, interaction = inter, kappa = 10, w_pp_cutoff = 2)
  sim <- project(params, effort = 0, t_max=500, dt=0.1, t_save=1)
  output <- getCommunitySlope2(sim)[500,,]$mse
  return(output)
}

runsobol2 <- function(X){
  out <- matrix(nrow = nrow(X), ncol = 1, NA)
  for (i in 1:nrow(X)){
    out[i,] <- runsobol(X$w_inf.4[i], X$w_inf.6[i], X$k_vb.5[i], X$k_vb.4[i], X$R_max.1[i], X$k_vb.6[i], X$beta.6[i], X$w_inf.3[i], X$R_max.4[i], X$k_vb.3[i])
  }
  out <- as.data.frame(out)
  return(out)
}

seq.sobol <- multisensi(design = sobol2007, model = runsobol2,
                        reduction = NULL, analysis = analysis.sensitivity, center = TRUE,
                        design.args = list(X1 = Xb[1:(m/2), ], X2 = Xb[(1 + m/2):m, ], nboot = 100),
                        analysis.args = list(keep.outputs = FALSE))





############################################################################
#Load necessary functions
#load necessary functions to calculate measure of linearity (MSE)
library(dvmisc)

getCommunitySlope2 <- function(sim, species = 1:nrow(sim@params@species_params),
                               biomass = TRUE, ...) {
  check_species(sim, species)
  size_range <- get_size_range_array(sim@params, ...)
  # set entries for unwanted sizes to zero and sum over wanted species, giving
  # array (time x size)
  total_n <-
    apply(sweep(sim@n, c(2, 3), size_range, "*")[, species, , drop = FALSE],
          c(1, 3), sum)
  # numbers or biomass?
  if (biomass)
    total_n <- sweep(total_n, 2, sim@params@w, "*")
  # previously unwanted entries were set to zero, now set them to NA
  # so that they will be ignored when fitting the linear model
  total_n[total_n <= 0] <- NA
  # fit linear model at every time and put result in data frame
  slope <- plyr::adply(total_n, 1, function(x, w) {
    lm.fit <- lm(log(x) ~ log(w))
    summary_fit <- summary(lm.fit)
    out_df <- data.frame(
      slope = summary_fit$coefficients[2, 1],
      intercept = summary_fit$coefficients[1, 1],
      r2 = summary_fit$r.squared,
      mse = get_mse(lm.fit)
    )
  }, w = sim@params@w)
  dimnames(slope)[[1]] <- slope[, 1]
  slope <- slope[, -1]
  return(slope)
}


# internal
check_species <- function(object, species){
  if (!(is(species,"character") | is(species,"numeric")))
    stop("species argument must be either a numeric or character vector")
  if (is(species,"character"))
    check <- all(species %in% dimnames(object@n)$sp)  
  if (is(species,"numeric"))
    check <- all(species %in% 1:dim(object@n)[2])
  if (!check)
    stop("species argument not in the model species. species must be a character vector of names in the model, or a numeric vector referencing the species")
  return(check)
}
