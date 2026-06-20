######################################################*
# BIORESTAURACION 2016 - BIOMASS ANALYSIS
# BMAP - PERU LNG
#
# ESTIMATION OF PLANT BIOMASS BASED ON ALLOMETRIC 
# EQUATIONS DURING BIORESTORATION MONITORING 2016
#
# *Helper functions*
#
# Last version: 27-07-2020
# By: Hector Chuquillanqui / hchuquillanqui@gmail.
######################################################*

### 01. PLOT VALIDATION ####
plot_validation <- function(model, ...) {
  # old.par <- par(no.readonly = TRUE)
  par(mfrow = c(2,2))
  # https://stackoverflow.com/questions/14671172/how-to-convert-r-formula-to-text
  variables <- deparse(model$terms, width.cutoff = 500)
  model.name <- deparse(substitute(model))
  plot(model, which = c(1:4), ...)
  mtext(paste("[", model.name, "] : ", variables), 
        side = 3, line = -1.5, outer = TRUE, cex = 1, adj = 0.1)
  # par(old.par)
  par(mfrow = c(1,1))
}

### 02. ROOT MEAN SQUARED DEVIATION (RMSD) ####
# How much is difference between observed and predicted values
# (See Pineiro etal. 2008 about RMSD)
# Sometimes RMSD is refered as RMSE (Root mean squared error)
# https://en.wikipedia.org/wiki/Root-mean-square_deviation
# https://stats.stackexchange.com/questions/242787/how-to-interpret-root-mean-squared-error-rmse-vs-standard-deviation
# https://stats.stackexchange.com/questions/20741/mean-squared-error-vs-mean-squared-prediction-error
# https://stackoverflow.com/questions/43123462/how-to-obtain-rmse-out-of-lm-result
# https://stackoverflow.com/questions/40901445/function-to-calculate-r2-r-squared-in-r
get_rmsd <- function(mod) {
  agb_pre = predict(mod)
  rss = sum((agb_pre - mod$model[1])^2)
  mse = rss*1/(length(agb_pre)-1)
  sqrt(mse)
}

### 03. DISPLAY MODEL PERFORMANCE DESCRIPTORS ####
get_fit_descriptors <- function(mod) {
  require(MuMIn)
  require(broom)
  
  broom::glance(mod) %>% 
    select(r.squared, adj.r.squared, sigma, AIC) %>%
    ### AIC corrected for small samples (AICc)
    add_column(AICc = MuMIn::AICc(mod)) %>%
    ### Root mean square deviation 
    add_column(RMSD = get_rmsd(mod))
}

### 04. SET BIOMASS MODEL FROM SEVERAL EQUATION OPTIONS ####
eq_agb_test <- function(df,
                        equa.mod = NULL,
                        form = FALSE
) {
  # Use a particular equation according type.eq
  equation_model <- 
    switch(equa.mod,
           ii   = lm(log10(AGB) ~ log10(BA), data = df),
           iv   = lm(log10(AGB) ~ log10(HMAX), data = df),
           vi   = lm(log10(AGB) ~ log10(CA), data = df),
           viii = lm(log10(AGB) ~ log10(CA) + log10(HMAX), data = df),
           x    = lm(log10(AGB) ~ log10(BA) + log10(CA), data = df),
           xii  = lm(log10(AGB) ~ log10(BA) + log10(HMAX) + log10(CA), data = df)
    )
  
  if(form == TRUE) return(formula(equation_model))
  
  return(equation_model)
}


### 05. GET MODEL COEFFICIENT FROM OBS ~ PRED MODEL OR SELECTED EQUATION ####
get_coeffic <- function(data, indices, equa = NULL, set.form = NULL) {
  data <-  data[indices, ]
  
  if(!is.null(set.form)) {
    form <- formula(set.form)
  } else {
    form <- formula("agb_obs ~ agb_pre")
  }
  
  if(is.null(equa)) {
    lm.out <-  lm(form, data = data)
    return(lm.out$coefficients)
  }
  
  lm.out <- eq_agb_test(data, equa.mod = equa)
  return(lm.out$coefficients)
}

### 06. BOOTSTRAP CONFIDENCE INTERVAL OF COEFFICIENTS IN OBS ~ PRE MODEL ####
get_ci_obs_pre <- function(boot.obj) {
  ci_int <- boot::boot.ci(boot.obj, type = "bca", index = 1)
  ci_slp <- boot::boot.ci(boot.obj, type = "bca", index = 2)
  tibble(ci_low = c(ci_int$bca[4], ci_slp$bca[4]),
         ci_upp = c(ci_int$bca[5], ci_slp$bca[5])
  )
}


### 07. GET BOOTSTRAPING MEAN AND CONFIDENCE INTERVALS ####
# Make "boot" object
get_mean_boot <- 
  function(x) {
    Mboot <- boot::boot(x,
                        function(x, i)
                          mean(x[i]),
                        R = 10000)
  }

# Calculate CI from "boot" object
get_ci_boot <- function(xboot) {
  boot::boot.ci(xboot,
                conf = 0.95,
                type = c("bca")
  )$bca
}
