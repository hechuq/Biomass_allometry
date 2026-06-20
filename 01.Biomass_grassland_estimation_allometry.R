######################################################*
# BIORESTAURACION 2016 - BIOMASS ANALYSIS
# BMAP - PERU LNG
#
# ESTIMATION OF PLANT BIOMASS BASED ON ALLOMETRIC 
# EQUATIONS DURING BIORESTORATION MONITORING 2016
#
# *Updated version of "BMAP Biomass.R" script* 
#
# Last version: 27-07-2020
# By: Hector Chuquillanqui / hchuquillanqui@gmail.
######################################################*

library(tidyverse)
library(readxl)
library(stringr)
library(janitor) 
library(plotly)
library(broom)
library(MuMIn) # To calculate AIC corrected for small samples (AICc) 
library(caret) # To perform cross-validation

source("code/pairsPannelFunctions.R")
source("code/S01.helpers_functions.R")

##############################################################################*
### READ DATA AND DATA FORMAT                                             ####
##############################################################################*

# READ FROM EXCEL FILE
# files <- list.files(path = "data", full.names = TRUE)
biomass <- read_excel("data/Biomasa_2016.xlsx")

biomass %>% janitor::clean_names()

# RENAME COLUMNS
# Equivalent allometric variables measured in the field by Olivera etal. are
# indicated with comments
col_names <- c(
  "KP",   
  "transecto",
  "punto",
  "especie",
  "perimetro",
  "peso_fresco",
  "raiz_1",
  "peso_seco",
  "raiz_2",
  "basal_diam_1", # Longest basal diameter             = d1
  "basal_diam_2", # Basal diameter perpendicular to d1 = d2
  "folia_diam_1", # Longest tussock crown diameter     = dc1
  "folia_diam_2", # Diameter perpendicular to dc1      = dc2
  "altura_veget", # Height in the field                = h
                  # Note: hmax ("streched by hand") not measured
  "altura_inflo", 
  'altura_maxim', 
  "area_basal",   
  "area_folia",   
  "x19",
  "easting",
  "norting",
  "elevacion",
  "x23",
  "x24",
  "x25"
)

# FORMAT BIOMASS DATASET
biomass_dat <- biomass %>% 
  janitor::clean_names() %>%
  set_names(nm = col_names) %>%
  slice(1:56) %>% 
  select(-starts_with(c("x", "raiz"))) %>%
  mutate(peso_seco = as.numeric(peso_seco),
         id = row_number(),
         especie = str_replace(especie, ".*rigida", "Festuca rigida")) %>%
  separate(col = especie, into = c("genero", "epiteto"), remove = FALSE) %>%
  select(id, everything())

glimpse(biomass_dat)

# BASAL AND FOLIAR (CROWN) AREA CALCULATIONS
# Following Johnson etal. (1988), as cited by Olivera etal.
get_area_ellipse <- function(d_long, d_short) {
  pi * (d_long / 2) * (d_short / 2)
}

biomass_dat <-
  biomass_dat %>%
  mutate(basal_diam_2_new = ifelse(is.na(basal_diam_2), basal_diam_1, basal_diam_2),
         folia_diam_2_new = ifelse(is.na(folia_diam_2), folia_diam_1, folia_diam_2),
         area_basal_d11    = get_area_ellipse(basal_diam_1, basal_diam_1),
         area_basal_d12new = get_area_ellipse(basal_diam_1, basal_diam_2_new),
         area_folia_d12    = get_area_ellipse(folia_diam_1, folia_diam_2),
         area_folia_d12new = get_area_ellipse(folia_diam_1, folia_diam_2_new)
         )

##############################################################################*
### 1. EXPLORE DATA                                                       ####
##############################################################################*

### How above-ground biomass (AGB) data is distributed between species
biomass_dat %>%
  drop_na(peso_seco) %>%
  group_by(especie) %>%
  summarise(mean = mean(peso_seco),
            sample = length(peso_seco),
            total = sum(peso_seco)) %>%
  arrange(desc(mean))

### See relationships between AGB (peso_seco) and allometric measurements
biomass_dat %>%
  drop_na(peso_seco) %>%
  drop_na(altura_veget) %>%
  select(peso_seco, perimetro, 
         area_basal_d12new, area_folia_d12new,
         altura_veget, altura_maxim
         ) %>%
  mutate_all(.funs = ~log10(.)) %>%
  pairs(lower.panel = panel.smooth.2,
        diag.panel = panel.hist,
        upper.panel = panel.cor,
        lwd.smooth = 2)

### Visualize AGB relationship with selected allometric variable
sppxgen <- biomass_dat %>% distinct(genero, especie) %>% count(genero) 
spp <- biomass_dat %>% distinct(genero, especie) %>% arrange(genero)
spp.shp <- c(21, 22, 24)
gen.fil <- 
  c(
    "Calamagrostis" = "#e41a1c",
    "Festuca" = "#377eb8",
    "Jarava" = "#4daf4a",
    "Nasella" = "#984ea3",
    "Stipa" = "#ff7f00"
  )

spp.scales <- 
  purrr::map2(.x = sppxgen$n, .y = sppxgen$genero,
              .f = function(x, y) {
                tibble(
                  genero = y,
                  shp = spp.shp[1:x], 
                  fil = rep(gen.fil[y], x)
                )
              }
  ) %>%
  bind_rows() %>%
  add_column(especie = spp$especie)

spp.shp.manual <- spp.scales$shp %>% set_names(spp.scales$especie)
spp.fil.manual <- spp.scales$fil %>% set_names(spp.scales$especie)

plot_xy_relation <- 
  function(df, resp = "peso_seco", expl = NULL,
           make.plotly = FALSE) {

    p <-
      df %>%
      ggplot(aes_string(x = expl, y = resp)) +
      geom_smooth(method = "lm", color = "darkgrey") +
      geom_jitter(aes(shape = factor(especie), fill = genero),
                  size = 4, width = 0.02, height = 0.02, alpha = 0.7) +
      scale_x_log10() +
      scale_y_log10() +
      scale_shape_manual(values = spp.shp.manual) +
      scale_fill_manual(values = gen.fil) +
      guides(shape = guide_legend("Especies", 
                                  override.aes = list(shape = spp.shp.manual,
                                                      fill = spp.fil.manual)
      ),
      fill = guide_legend("Genero", 
                          override.aes = list(fill = gen.fil,
                                              shape = 21)
      )
      ) +
      labs(title = paste0("Relation of '", resp, "' with '", expl, "'"), 
           subtitle = "Response and explanatory were log10 transformed")

    if(make.plotly == TRUE) {
      p <- plotly::ggplotly(p)
    }

    return(p)
  }

# EXPLORE RELATIONS OF "peso_seco" WITH ALLOMETRIC VARIABLES IN LOG10 SCALE
plot_xy_relation(biomass_dat, resp = "peso_seco", expl = "area_basal_d12new")
plot_xy_relation(biomass_dat, resp = "peso_seco", expl = "area_folia_d12new")
plot_xy_relation(biomass_dat, resp = "peso_seco", expl = "altura_veget")
plot_xy_relation(biomass_dat, resp = "peso_seco", expl = "altura_maxim")


# PLOT INTERACTIVELY
# Caution: ggplotly show axes marks log10 transformed, but preserve relationship
plot_xy_relation(
  biomass_dat,
  resp = "peso_seco",
  expl = "area_basal_d12new",
  make.plotly = TRUE
)


##############################################################################*
### 2. MODELLING BIOMASS                                                 ####
###    Following Oliveras etal. 2013.
###
###
##############################################################################*

### 2.1 GRASS ALLOMETRY SENSU OLIVERAS  
# Oliveras etal. (2013) tested twelve models to estimate aboveground biomass (ABG)
# based on allometric variables: Basal area (BA), canopy area (CA) and plant
# maximun height (hmax) [Table 1]

# Equivalence between variables:
# Oliveras etal. ~ BMAP 2016
# AGB   ~  peso_seco
# BA    ~  area_basal_d12new
# CA    ~  area_folia_d12new
# hmax  ~  < no medido >
#          < alternativa > altura_maxim
# h     ~  altura_veget

# See description of variables in "data/Biomasa_2016.xlsx"

### 2.2 MODELS FOR TESTING
# First choices: Models IX, X, XI, XII (Table 1; Oliveras etal., 2013)
# [IX, XI] => Linear models
# [X, XII] => Non-linear models

# [IX]  => a + b(area_basal_d12new) + c(area_folia_d12new)
# [XI]  => a + b(area_basal_d12new) + c(altura_veget) + d(area_folia_d12new)
# [X]   => a(area_basal_d12new)^b * (area_folia_d12new)^c
# [XII] => a(area_basal_d12new)^b * (altura_veget)^c * (area_folia_d12new)^d

# Drop missing observations
bio_dat <- biomass_dat %>%
  drop_na(peso_seco) %>%   # 1 missing value
  drop_na(altura_veget)    # 2 missing value

# [IX, XI] => Linear models
mod_ix <- lm(peso_seco ~ area_basal_d12new + area_folia_d12new, 
             data = bio_dat)
mod_xi <- lm(peso_seco ~ area_basal_d12new + altura_veget + area_folia_d12new, 
             data = bio_dat)

summary(mod_ix)
anova(mod_ix)
plot_validation(mod_ix, col = 4)

summary(mod_xi)
anova(mod_xi)
plot_validation(mod_xi, col = 4)

# [X, XII] => Non-linear models
mod_x   <- lm(log10(peso_seco) ~ log10(area_basal_d12new) + log10(area_folia_d12new),
             data = bio_dat)
mod_xii <- lm(log10(peso_seco) ~ log10(area_basal_d12new) + log10(altura_veget) + log10(area_folia_d12new),
             data = bio_dat)

summary(mod_x)
anova(mod_x)
plot_validation(mod_x, col = 4)

summary(mod_xii)
anova(mod_xii)
plot_validation(mod_xii, col = 4)

### 2.3 MODEL COMPARISON
# Table of model fit descriptors
model_list <- list(mod_ix, mod_xi, mod_x, mod_xii)
names(model_list) <- c("mod_ix", "mod_xi", "mod_x", "mod_xii")
purrr::map_dfr(model_list, get_fit_descriptors, .id = "model")

# More parsimonious models
mod_x_simple <- update(mod_x, ~ . - log10(area_folia_d12new))
mod_xii_simple <- update(mod_xii, ~ . - log10(area_folia_d12new))

anova(mod_x_simple, mod_x)
anova(mod_xii_simple, mod_xii)
anova(mod_x_simple, mod_xii_simple)

model_simple_list <- list(mod_x, mod_xii, mod_x_simple, mod_xii_simple)
names(model_simple_list) <- c("mod_x", "mod_xii", "mod_x_simple", "mod_xii_simple") 
purrr::map_dfr(model_simple_list, get_fit_descriptors, .id = "model") 

### MODELS WITH ALTERNATIVE VARIABLES ("altura_maxim")
mod_xii_alt <- update(mod_xii, ~ . - log10(altura_veget) + log10(altura_maxim))
plot_validation(mod_xii_alt, col = 4)
summary(mod_xii_alt)
anova(mod_xii_alt)

mod_xii_alt_simple <- update(mod_xii_alt, ~ . - log10(area_folia_d12new))

# Compare with other models
model_simple_list_alt <- list(mod_x_simple, mod_xii_simple, mod_xii_alt_simple)
names(model_simple_list_alt) <- c("mod_x_simple", "mod_xii_simple", "mod_xii_alt_simple")
purrr::map_dfr(model_simple_list_alt, get_fit_descriptors, .id = "model")


##############################################################################*
### 3. TEST BIOMASS MODELS FROM PREVIOUS STUDIES                          ####
###    Use models proposed by Oliveras etal. 2013.
###
###
##############################################################################*

# Abreviated names for variables, following nomenclature from Oliveras etal.
bio_dat <-
  bio_dat %>%
  select(SPP = especie,
         GEN = genero,
         AGB = peso_seco, 
         BA = area_basal_d12new,
         CA = area_folia_d12new,
         HMAX = altura_veget) %>%
  mutate(agb_obs = log10(AGB))


### 3.1 FIT MODELS FROM OLIVERAS ETAL 2013 WITH BIORESTORATION DATA

# Best performance of non-linear equations in Oliveras etal.
eq_non_linear <- c("ii", "iv", "vi", "viii", "x", "xii")

# Store data and regression models for each equation in a data.frame
models_non_linear <-
  tibble(
    data = list(bio_dat),
    equation = eq_non_linear,
    lm.model = purrr::map(equation, .f = function(equa) eq_agb_test(bio_dat, equa))
  )

# Extract model descriptors, coeficients and predicted values
models_nl_details <-
  models_non_linear %>%
  mutate(mod_descriptors = purrr::map(lm.model, ~ get_fit_descriptors(.x))) %>%
  mutate(mod_coeficients = purrr::map(lm.model, ~ broom::tidy(.x, conf.int = TRUE))) %>%
  mutate(agb_pre = purrr::map(lm.model, ~ predict(.x)))

models_nl_details %>% unnest(mod_descriptors)
models_nl_details %>% unnest(mod_coeficients)


### 3.2 EXPLORE RELATION BETWEEN OBSERVED AND PREDICTED VALUES

# Observed and predicted data for each equation used
obs_pre_agb <-
  models_nl_details %>% unnest(c(data, agb_pre)) %>%
  mutate(agb_obs = log10(AGB)) %>%
  select(equation, agb_obs, agb_pre)

# Extract model descriptors, coefficients and bootstrap confidence interval
obs_pre_mods_details <-
  obs_pre_agb %>%
  group_by(equation) %>%
  nest() %>%
  mutate(obs_pre_mods = purrr::map(data, ~ lm(agb_obs ~ agb_pre, data = .x)),
         mods_descrip = purrr::map(obs_pre_mods, ~ get_fit_descriptors(.x)),
         tidy_coeffic = purrr::map(obs_pre_mods, ~ broom::tidy(.x)), 
         boot_coeffic = purrr::map(data, ~ boot::boot(.x, get_coeffic, R = 1000)),
         boot_ci_coef = purrr::map(boot_coeffic, ~ get_ci_obs_pre(.x)) 
         )

obs_pre_mods_details %>%
  unnest(mods_descrip) %>%
  select(-c(data, obs_pre_mods, boot_coeffic, tidy_coeffic, boot_ci_coef))
obs_pre_mods_details %>%
  unnest(c(tidy_coeffic, boot_ci_coef)) %>%
  select(-c(data, obs_pre_mods, mods_descrip, boot_coeffic))

# Update equation xii: simplify model
mod_xii_upd <- models_nl_details %>% filter(equation == "xii") %>% 
  select(lm.model) %>% pluck(1, 1)

mod_xii_upd_simple <- update(mod_xii_upd, . ~ . -log10(CA), data = bio_dat)
anova(mod_xii_upd, mod_xii_upd_simple)

bind_rows(models_non_linear,
          tibble(
            data = list(bio_dat),
            equation = "xii_upd",
            lm.model = list(mod_xii_upd_simple)
          )) %>%
  mutate(mod_descriptors = purrr::map(lm.model, ~ get_fit_descriptors(.x))) %>% 
  unnest(mod_descriptors)
  # --> Small constrast of 'r squared' and 'RMSD' between "xii" and "xii_upd"


### 3.3 VISUALIZE RELATION BETWEEN OBSERVED AND PREDICTED VALUES

obs_pre_agb %>%
  ggplot(aes(x = agb_pre, y = agb_obs, color = equation)) +
  geom_abline(slope = 1, linetype = 2, size = 1) +
  geom_point(size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, size = 1.5) +
  coord_fixed(xlim = c(0, 3.5), ylim = c(0,3.5)) +
  facet_wrap( . ~ equation) +
  scale_color_brewer(palette = "Set2") +
  ggtitle(label = bquote('Comparison of Oliveras equations: ' ~ bold('Observed') ~
                           'vs' ~ bold('Predicted')),
          subtitle = 'Coefficients estimated for Biorestoration data'
  ) +
  theme_bw()
  # --> NOTE 1: Dispersion in y-direction is related with values of RMSD, as a 
  #             measure of how much observed values depart from predicted ones.
  # --> NOTE 2: Regression line pass through 1:1 diagonal since this was 
  #             constructed over the same observed values that we used to fit
  #             the model (with their particular residuals).



##############################################################################*
### 4. EXPLORE PERFORMANCE OF PREVIOUS ALLOMETRIC MODELS                  ####
###    Calibrate and test models proposed by Oliveras etal. 2013
###    via cross validation with Biorestoration data
###
##############################################################################*

### 4.1 ITERATIVE PROCEDIMENT WITH CROSS-VALIDATION ----
# The following code calibrate a pre-defined set of equations by means of
# cross-validation method with caret::train, and repeat the proccess with
# a resampled 75% of Biorestoration data (train data), that changes at
# every iteration. After the calibration, the 25% left is used to assess
# model performance.

  # --> NOTE 1: Proccess takes long time (around 20 minutes) for 100 iterations
  #             and 6 equations with repeated cross-validation. Output object
  #             was stored in a .RDS file (biorestoration_biomass_cross-val.rds)
  #             which is called below in section 4.3.


# Working example with 10 iterations and 2 equations.
vec_iter <- 1:10
equations <- eq_non_linear[1:4]

# Full iterations (100) and all equations (6). Uncomment if necesary.
# vec_iter <- 1:100
# equations <- eq_non_linear
# After running the code below, save object with:
#  
# Remember: results can change a little because of random resampling

cross_validation_iter <-
  purrr::map(vec_iter, .f = function(i_var) {

    ### TRAINING AND TESTING DATA
    ### Split data in training (75%) and testing set (25%) via resampling
    inTraining <- caret::createDataPartition(bio_dat$AGB, p = .75, list = FALSE)
    training <- bio_dat[ inTraining, ]
    testing  <- bio_dat[-inTraining, ]

    testing_new <-
      testing %>% select(-c(SPP, GEN, agb_obs)) %>%
      set_names(paste0("log10(", names(.), ")")) %>%
      mutate_all(.funs = ~ log10(.))

    ### CUSTOMIZE REPEATED CROSS-VALIDATION
    train.control <-
      caret::trainControl(
        method = "repeatedcv",
        number = 5,    # Numbers of folds to split 75% data
        repeats = 20,  # Times to repeat CV
        savePredictions = "all",
        returnResamp = "all",
        returnData = TRUE,
        verboseIter = FALSE
      )

    mod_calibration <- 
      purrr::map(equations, .f = function(j_var) {

        ### TRAIN (CALIBRATE) MODEL 
        ### With 75% original data (4 folds for training, 1 fold for testing)
        equa.form <- eq_agb_test(df = bio_dat, equa.mod = j_var, form = TRUE)

        trained_model <-
          caret::train(
            form = equa.form,
            data = training,
            method = "lm",
            trControl = train.control 
          )

        cat("Model ", j_var, " -- ", "iteration : ", i_var, " DONE ", "\n")

        trained_output <- 
          trained_model[c("pred", "finalModel", "resample", "results")]

        ### PREDICT FROM TRAINED MODEL 
        ### With 25% original data
        predicted_agb <- 
          as_tibble(predict(trained_model$finalModel,
                            newdata = testing_new,
                            interval = "confidence"
          ))

        predicted_output <- 
          bind_cols(AGB_obs = testing_new$`log10(AGB)`, 
                    AGB_pre = predicted_agb$fit,
                    AGB_pre_low = predicted_agb$lwr,
                    AGB_pre_upp = predicted_agb$upr
          ) 

        ### CALIBRATION OUTPUTS
        return(list(trained = trained_output, 
                    predicted = predicted_output)
        )

      }) %>%
      set_names(paste0("equa_", equations))

    return(mod_calibration)

  }) %>%
  set_names(paste0("iteration_", vec_iter))

# INSPECT OUTPUT FROM ITERATIVE CROSS-VALIDATION OBJECT
str(cross_validation_iter, max.level = 2)

# Inspect output: 1 iteration out of 10
# Calibrated (train) output
str(cross_validation_iter$iteration_1$equa_ii$trained, max.level = 1)

# Predicted (test) output
str(cross_validation_iter$iteration_1$equa_ii$predicted)

# Calibrated (train) output objects
cross_validation_iter$iteration_1$equa_ii$trained$resample
cross_validation_iter$iteration_1$equa_ii$trained$pred
cross_validation_iter$iteration_1$equa_ii$trained$finalModel
cross_validation_iter$iteration_1$equa_ii$trained$results

### 4.2 EXTRACT DETAILED OUTPUTS FROM CV OBJECT ----
# FUNCTION: Extract cv outputs and present in tibble (dataframe) format
bind_cv_rows <- 
  function(cv.output,
           output = c("trained", "predicted"),
           trained.out = c("results", "finalModel", "pred", "resample")
           ) {

  output <- match.arg(output)
  trained.out <- match.arg(trained.out)

  if(output == "predicted") fun.formula <- formula(~ .x[[output]])
  if(output == "trained")   fun.formula <- formula(~ .x[[output]][[trained.out]])

  cv.output %>%
    purrr::map(
      .,
      .f = function(equa_lvl) {
        tmp <- equa_lvl %>% purrr::map(., fun.formula)
        if("lm" %in% class(tmp[[1]])) {
          tibble(equation = names(tmp), models = tmp)
        } else {
          tmp %>% bind_rows(., .id = "equation") %>% as_tibble
        }
      }
    ) %>%
    bind_rows(., .id = "iteration")

}

bind_cv_rows(cross_validation_iter, output = "trained", trained.out = "pred")
bind_cv_rows(cross_validation_iter, output = "trained", trained.out = "resample")
bind_cv_rows(cross_validation_iter, output = "trained", trained.out = "results")
bind_cv_rows(cross_validation_iter, output = "trained", trained.out = "finalModel")
bind_cv_rows(cross_validation_iter, output = "predicted")

bind_cv_rows(cross_validation_iter, output = "trained", trained.out = "finalModel")


### 4.3 EXPLORE COMPLETE CV OBJECT: 100 ITERATIONS WITH REPEATED CV ----
biorest_cross_val_iter <- readRDS("data/biorestoration_biomass_cross-val.rds")

### 4.3.1 TRAINED CROSS-VALIDATION RESULTS
# Based on 75% of original data

cv_biorest_train_results <-
  biorest_cross_val_iter %>%
  bind_cv_rows(output = "trained", trained.out = "results")

cv_biorest_results_mean <-
  cv_biorest_train_results %>%
  group_by(equation) %>%
  summarise_if(is.numeric, mean)

cv_biorest_results_boxplot <-
  cv_biorest_train_results %>%
  group_by(equation) %>%
  summarise(ymin.whisker  = median(RMSE) - 1.5*IQR(RMSE),
            y25 = quantile(RMSE, 0.25),
            y50 = mean(RMSE),
            y75 = quantile(RMSE, 0.75),
            ymax.whisker = median(RMSE) + 1.5*IQR(RMSE)
            )

cv_biorest_train_results %>%
  ggplot(aes(x = equation, y = RMSE)) +
  geom_violin(fill = "#addd8e") +
  geom_boxplot(
    data = cv_biorest_results_boxplot %>% left_join(cv_biorest_results_mean),
    aes(ymin = y25, lower = y25, 
        middle = y50,
        upper = y75, ymax = y75
        ),
    width = 0.05, fatten = NULL, fill = "NA",
    stat = "identity"
  ) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 3,
               fill = "#31a354") +
  scale_y_continuous(limits = c(0.2, 0.6), 
                     breaks = seq(0.2, 0.6, 0.05), 
                     expand = c(0,0)) +
  ggtitle(label = "Root-mean squared error (RMSE) from trained resamples",
          subtitle = paste0("Small values indicates better match between", "\n",
                            "*Observed* and *Predicted* values")) +
  theme_bw()

ggsave(device = "pdf", width = 6, filename = "figures/RMSE_iterations.pdf")


### 4.3.2 TESTED RESULTS FROM CROSS-VALIDATION 

# Based on 25% of original data
cv_biorest_test_predict <-
  biorest_cross_val_iter %>%
  bind_cv_rows(output = "predicted")

# Performance and coefficients for Observed ~ Predicted model
cv_biorest_test_perfor <-
  cv_biorest_test_predict %>%
  group_by(iteration, equation) %>%
  nest() %>%
  mutate(obs_pre_mods = purrr::map(data, ~ lm(AGB_obs ~ AGB_pre, data = .x)),
         mods_descrip = purrr::map(obs_pre_mods, ~ get_fit_descriptors(.x)),
         tidy_coeffic = purrr::map(obs_pre_mods, ~ broom::tidy(.x))
  )

# Mean of performance descriptors for each equation
cv_biorest_test_perfor_mean <-
  cv_biorest_test_perfor %>%
  unnest(mods_descrip) %>%
  group_by(equation) %>%
  summarise_if(is.numeric, mean)


### MODEL PERFORMANCE: RMSD
# Boxplot elements for RMSD (Root-mean square deviation) for each equation
cv_biorest_test_boxplot <-
  cv_biorest_test_perfor %>%
  unnest(mods_descrip) %>%
  group_by(equation) %>%
  summarise(ymin.whisker  = median(RMSD) - 1.5 * IQR(RMSD),
            y25 = quantile(RMSD, 0.25),
            y50 = mean(RMSD),
            y75 = quantile(RMSD, 0.75),
            ymax.whisker = median(RMSD) + 1.5 * IQR(RMSD)
  )

# Boostraping confidence interval (CI) for mean RMSD for each equation
cv_biorest_test_perfor_boot <-
  cv_biorest_test_perfor %>%
  unnest(mods_descrip) %>%
  group_by(equation) %>%
  summarise(mean = mean(get_mean_boot(RMSD)$t[,1]),
            ci_upp = get_ci_boot(get_mean_boot(RMSD))[5],
            ci_low = get_ci_boot(get_mean_boot(RMSD))[4]
  )

# Plot RMSD for all equations used
cv_biorest_test_perfor %>%
  unnest(mods_descrip) %>%
  select(-data, -obs_pre_mods) %>%
  ggplot(aes(x = equation, y = RMSD)) +
  geom_violin(fill = "#addd8e") +
  geom_boxplot(
    data = cv_biorest_test_boxplot %>% left_join(cv_biorest_test_perfor_mean),
    aes(ymin = y25, lower = y25, 
        middle = y50,
        upper = y75, ymax = y75
    ), 
    width = 0.2, fatten = NULL, fill = "NA",
    stat = "identity"
  ) +
  geom_errorbar(
    data = cv_biorest_test_perfor_boot %>% left_join(cv_biorest_test_perfor_mean),
    aes(ymin = ci_low, ymax = ci_upp),
    color = "grey50", width = 0, size = 1.5
  ) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 3,
               fill = "#31a354") +
  scale_y_continuous(limits = c(0, 0.7),
                     breaks = seq(0, 0.7, 0.05),
                     expand = c(0, 0)) +
  ggtitle(label = "Root-mean squared deviation (RMSD) from testing resamples",
          subtitle = paste0("Small values indicates better match between", "\n",
                            "*Observed* and *Predicted* values")
          ) +
  labs(caption = paste0("Green point = mean, ",
                        "black rectangle = Interquantile region, ",
                        "grey bar = bootstrap confidence intervals")
       ) +
  theme_bw()

ggsave(device = "pdf", width = 6, filename = "figures/RMSD_iterations_testing.pdf")


### OBSERVED ~ PREDICTED MODELS

# Plot regressions for all equations used
cv_biorest_test_predict %>%
  ggplot(aes(x = AGB_pre, y = AGB_obs, color = equation)) +
  geom_abline(slope = 1, linetype = 2, size = 0.5) +
  facet_wrap( ~ equation) +
  geom_smooth(aes(group = paste(iteration, equation)),
              method = "lm", se = FALSE, size = 0.1, alpha = 0.1, color = "grey") +
  geom_smooth(aes(group = equation),
              method = "lm", se = FALSE, size = 0.5, alpha = 1) +
  coord_fixed(xlim = c(0, 3), ylim = c(0, 3), expand = FALSE) +
  labs(x = "log10(AGB predicted)", y = "log10(ABG observed)") +
  ggtitle(label = "Observed vs Predicted AGB relationships",
          subtitle = paste0("Exact match indicated by 1:1 line", "\n",
                            "Regression lines in color represent all iterations")) +
  theme_bw()

ggsave(device = "pdf", width = 8, height = 5,
       filename = "figures/ObsPre_iterations_testing.pdf")


### ESTIMATED COEFFICIENTS

# Intercept and slope means and CI
cv_biorest_test_perfor_coef_boot <- 
  cv_biorest_test_perfor %>% 
  unnest(tidy_coeffic) %>%
  group_by(equation, term) %>%
  summarise(mean = mean(get_mean_boot(estimate)$t[,1]),
            ci_upp = get_ci_boot(get_mean_boot(estimate))[5],
            ci_low = get_ci_boot(get_mean_boot(estimate))[4]) %>%
  mutate(baseline = case_when(
    term == "(Intercept)" ~ 0,
    term == "AGB_pre" ~ 1
  ))

# Plot estimates with CI
cv_biorest_test_perfor_coef_boot %>%
  ggplot(aes(x = equation, y = mean)) +
  geom_hline(aes(yintercept = baseline), linetype = 2, color = "red") +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_upp),
    color = "darkgrey", width = 0.1, size = 1
  ) +
  geom_point(shape = 21, fill = "#31a354", size = 5) +
  facet_wrap(~ term, scales = "free_x") +
  coord_flip() +
  scale_y_continuous(breaks = scales::pretty_breaks(10), 
                     name = "Mean estimate (log10) by bootstraping") +
  ggtitle(label = "Estimated coefficients in obs ~ pre AGB relationship ",
          subtitle = paste0(
            "Mean and CI produced by bootstraping", "\n",
            "If CI intercept  not include 0, predictions are biased", "\n",
            "If CI AGB_pre not include 1, predictions are overestimated"
          )) +
  theme_bw()

ggsave(device = "pdf", width = 6, filename = "figures/COEFFICIENTS_iterations_testing.pdf")


### 4.4 TEST FOR DIFFERENCES IN PERFORMANCE AND COEFFICIENTS ----
# (Compare with plots produced in section 4.3.2)

# Test for differences in RMSD between equations
test_rmsd_equation <- 
  cv_biorest_test_perfor %>%
  unnest(mods_descrip) %>%
  ungroup() %>%
  mutate(equation = factor(equation)) %>%
  lm(RMSD ~ equation, data = .) %>%
  multcomp::glht(model = ., linfct = multcomp::mcp(equation = "Tukey")) %>%
  summary(., test = multcomp::adjusted("fdr"))

test_rmsd_equation %>% multcomp::cld(.) 
old.par <- par(mai=c(1,1,1.5,1), no.readonly=TRUE)
test_rmsd_equation %>% multcomp::cld(.) %>% plot()
par(old.par)

# Test for differences in slope and intercept between equations
test_term_equation <- 
  cv_biorest_test_perfor %>%
  unnest(tidy_coeffic) %>%
  ungroup() %>%
  mutate(equation = factor(equation)) %>%
  group_by(term) %>%
  nest() %>%
  mutate(model.coef = purrr::map(data, ~ lm(estimate ~ equation, data = .x))) %>%
  mutate(model.comp = purrr::map(model.coef,
                                 ~ multcomp::glht(
                                   model = .x, 
                                   linfct = multcomp::mcp(equation = "Tukey")))
         ) %>%
  mutate(model.summ = purrr::map(model.comp, 
                                ~ summary(.x,
                                          test = multcomp::adjusted("fdr"))))

test_term_equation$model.summ[[1]] %>% multcomp::cld()
test_term_equation$model.summ[[2]] %>% multcomp::cld()

plot(test_term_equation$model.summ[[1]] %>% multcomp::cld(), main = "intercept")
plot(test_term_equation$model.summ[[2]] %>% multcomp::cld(), main = "slope (AGB_pre)")


##############################################################################*
### 5. EXPLORE PERFORMANCE OF CALIBRATED EQUATIONS FROM OLIVERAS          ####
###     5.1. Test equations will full biorestoration data. 
###     5.2. Adjust coefficients from equation manually
###
##############################################################################*


### 5.1 TEST OLIVERAS EQUATIONS

# FUNCTION: Change between different equations
# Non-linear equations for multispecies pool (Oliveras etal. 2013)
eq_agb <- function(.BA, .HMAX = NULL, .CA = NULL,
                   equa.mod = c("ii", "viii", "x", "xii"),
                   festuca = FALSE,
                   modify.coeff = NULL, # Choices: "a", "b", "c", "d"
                   new.coeff = NULL     # New value to modify coeff
                   ) {

  # Coefficients for each equation
  if (equa.mod == "ii") {
    a = 0.921
    b = 0.756

    if (festuca == TRUE) {
      a = 0.158
      b = 0.659
    }
  }

  if(equa.mod == "viii") {
    a = -3.245
    b = 1.021
    c = 0.921

    if (festuca == TRUE) {
      a = -4.028
      b = 0.965
      c = 1.319
    }
  }

  if(equa.mod == "x") {
    a = 5.5*10^-3
    b = 0.478
    c = 0.675

    if (festuca == TRUE) {
      a = -1.669
      b = 0.492
      c = 0.463
    }
  }

  if(equa.mod == "xii") {
    a = 1.0*10^-3
    b = 0.480
    c = 0.935
    d = 0.373

    if (festuca == TRUE) {
      a = 8.1*10^-4
      b = 0.475
      c = 1.094
      d = 0.286
    }
  }

  # Change a selected coefficient, if desired
  if(!is.null(modify.coeff)) {
    assign(modify.coeff, value = new.coeff)
  }

  # Use a particular equation according type.eq
  switch(equa.mod,
         ii   = a + b*log10(.BA),
         viii = a + b*log10(.CA) + c*log10(.HMAX),
         x    = a + b*log10(.BA) + c*log10(.CA),
         xii  = a + b*log10(.BA) + c*log10(.HMAX) + d*log10(.CA)
         )

  }

# SET PREDICTIONS FROM "MULTIESPECIES" AND "FESTUCA" EQUATIONS
equation_use <- list("all_species_equation" = FALSE, "festuca_equation" = TRUE)

agb_pred_oliv_all <- 
  purrr::map(.x = equation_use,
             .f = function(use.equa) {
               bio_dat %>%
                 mutate(agb_obs = log10(AGB)) %>%
                 mutate(agb_xii  = eq_agb(.BA = BA, .HMAX = HMAX, .CA = CA, 
                                          festuca = use.equa, equa.mod = "xii"),
                        agb_x    = eq_agb(.BA = BA, .CA = CA, 
                                          festuca = use.equa, equa.mod = "x"),
                        agb_viii = eq_agb(.CA = CA, .HMAX = HMAX,
                                          festuca = use.equa, equa.mod = "viii"),
                        agb_ii   = eq_agb(.BA = BA, 
                                          festuca = use.equa, equa.mod = "ii")
                 ) %>%
                 select(contains("agb", ignore.case = FALSE)) %>%
                 gather(key = "equation", value = "agb_pre", -agb_obs)
             }) %>%
    bind_rows(.id = "equation_type")

# PLOT OBSERVED ~ PREDICTED RELATIONSHIP
agb_pred_oliv_all %>%
  ggplot(aes(x = agb_pre, y = agb_obs, color = equation)) +
  geom_abline(slope = 1, linetype = 2, size = 1) +
  geom_point(size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, size = 1.5) +
  coord_fixed(xlim = c(-1, 4), ylim = c(0,3.5)) +
  facet_grid(rows = vars(equation_type)) +
  scale_color_brewer(palette = "Set1") +
  ggtitle(label = bquote('Comparison of Oliveras equations: ' ~ bold('Observed') ~
                           'vs' ~ bold('Predicted')),
          subtitle = 'Using coefficients from Table 2 and Appendix S4 in Oliveras etal. (2013)'
  ) +
  theme_bw()

# GET PERFORMANCE OF DIFFERENT EQUATIONS
model_descriptors <- 
  agb_pred_oliv_all %>%
  group_by(equation_type, equation) %>%
  nest() %>%
  mutate(mods = purrr::map(data, .f = function(df) lm(agb_obs ~ agb_pre, data = df)),
         desc.mods = purrr::map(mods, get_fit_descriptors),
         rmsd = purrr::map(data, .f = function(df) {
           with(df,
                sqrt(sum((agb_pre - agb_obs)^2)*1/(length(agb_obs)-1))
                )
         })
         ) 
 
model_descriptors %>% unnest(c(desc.mods, rmsd))

# GET CONFIDENCE INTERVALS FOR COEFFICIENTES
model_estimates <- 
  model_descriptors %>% 
  mutate(tidy.coef = purrr::map(mods, .f = function(mods) broom::tidy(mods, conf.int = TRUE)),
         boot.coef = purrr::map(data, .f = function(df) boot::boot(df, get_coeffic, R = 1000)),
         ci_boot = purrr::map(boot.coef, .f = function(boot.obj) {
           ci_int <- boot::boot.ci(boot.obj, type = "bca", index = 1)
           ci_slp <- boot::boot.ci(boot.obj, type = "bca", index = 2)
           tibble(ci_low = c(ci_int$bca[4], ci_slp$bca[4]),
                  ci_upp = c(ci_int$bca[5], ci_slp$bca[5])
                  )
         })
         )

# GET MODELS WITH INTERCEPT AND SLOPE INCLUDING 0 OR 1
model_estimates %>% unnest(c(tidy.coef, ci_boot)) %>%
  select(-c(data, mods, desc.mods, rmsd, boot.coef)) %>%
  # https://stackoverflow.com/questions/40446004/select-grouped-rows-with-at-least-one-matching-criterion
  filter(any(ci_low < 0 & 0 < ci_upp | ci_low < 1 & 1 < ci_upp))


### 5.2 ADJUST COEFFICIENTS MANUALLY  

# FUNCTION: CHANGE MANUALLY COEFFICIENTS FROM A SELECTED EQUATION
plot_coef_adjust <- function(
  equation,
  coef_mod,
  coef_vec
) {
  
  
  equa_vec <- rep(equation_use, each = length(coef_vec)/2)
  
  agb_pred_oliv_all_coefmod <- 
    purrr::map2(.x = equa_vec,
                .y = coef_vec,
                .f = function(x, y) {
                  bio_dat %>%
                    mutate(agb_obs = log10(AGB)) %>%
                    mutate(agb_pre  = eq_agb(.BA = BA, .HMAX = HMAX, .CA = CA, 
                                             festuca = x, 
                                             equa.mod = equation,
                                             modify.coeff = coef_mod,
                                             new.coeff = y)
                    ) %>%
                    rename_at(vars(agb_pre), ~ paste0("agb_", equation)) %>%
                    mutate(coef_mod_val = y) %>%
                    mutate(coef_mod_lab = coef_mod) %>%
                    select(contains(c("agb", "coef"), ignore.case = FALSE)) %>%
                    gather(key = "equation", value = "agb_pre", 
                           -agb_obs, -coef_mod_val, -coef_mod_lab)
                }) %>%
    bind_rows(.id = "equation_type")
    
  agb_pred_oliv_all_coefmod %>%
    group_by(equation_type, equation, coef_mod_val) %>%
    nest() %>%
    mutate(mods = purrr::map(data, .f = function(df) lm(agb_obs ~ agb_pre, data = df)),
           # desc.mods = purrr::map(mods, get_fit_descriptors)
           rmsd = purrr::map(data, .f = function(df) {
             with(df,
                  sqrt(sum((agb_pre - agb_obs)^2)*1/(length(agb_obs)-1))
             )
           }),
           summ = purrr::map(mods, .f = function(mods) broom::tidy(mods))
           # conf = purrr::map(mods, .f = function(mods) as.data.frame(confint(mods)))
    ) %>%
    unnest(summ)
  
  labels_coeff_mod <-
    agb_pred_oliv_all_coefmod %>%
    split(f = list(agb_pred_oliv_all_coefmod$equation_type,
                   agb_pred_oliv_all_coefmod$coef_mod_val),
          sep = " " 
    ) %>%
    purrr::map(.f = function(df) {
      tibble(
        agb_obs = lm(agb_obs ~ agb_pre, data = df) %>%
          predict(newdata = data.frame(agb_pre = max(df$agb_pre))),
        agb_pre = max(df$agb_pre),
        coef_mod_lab = coef_mod
      )
    }) %>%
    bind_rows(.id = "eq_coef") %>%
    separate(col = eq_coef, into = c("equation_type", "coef_mod_val"), sep = " ")

  equa_cols <- list(
    "xii"  = c("#984ea3", "#dac2d9"),
    "x"    = c("#caf9c0", "#4daf4a"),
    "viii" = c("#ceffff", "#377eb8"),
    "ii"   = c("#ffdacf", "#e41a1c")
  )

  p <-
    agb_pred_oliv_all_coefmod %>%
    ggplot(aes(x = agb_pre, y = agb_obs, color = coef_mod_val)) +
    geom_abline(slope = 1, linetype = 2, size = 1) +
    geom_point(size = 2, alpha = 0.5) +
    geom_smooth(aes(color = coef_mod_val, group = coef_mod_val),
                method = "lm", se = FALSE, size = 1.5) +
    ggrepel::geom_label_repel(data = labels_coeff_mod,
                              aes(x= agb_pre, y = agb_obs,
                                  label = paste(coef_mod_lab, coef_mod_val, sep = ":")
                              ), 
                              direction = "x", seed = 123, nudge_y = 1.5, 
                              size = 2, label.size = 0.3, label.padding = 0.15,
                              min.segment.length = unit(0, "lines"),
                              segment.color = "darkgrey",
                              force = 1, color = "black") +
    coord_fixed(xlim = c(-1, 4), ylim = c(0,3.5)) +
    facet_grid(rows = vars(equation_type)) +
    scale_color_gradient(low = equa_cols[[equation]][1], 
                         high = equa_cols[[equation]][2],
                         name = "Coefficient") +
    # https://stackoverflow.com/questions/43696101/r-bquote-remove-the-space-before-approximately-equal-plotmath-symbol
    ggtitle(label = bquote('Change coefficient "' *bold(.(coef_mod))* '"' ~
                             'in equation' ~ bold(.(equation))),
            subtitle = 'Using coefficients from Table 2 and Appendix S4 in Oliveras etal. (2013)'
    ) +
    theme_bw()

  print(p)
}

# PLOT REGRESSION WITH MODIFIED COEFFICIENTS
plot_coef_adjust(equation = "x",
                 coef_mod = "a",
                 coef_vec = rep(seq(-1, 1.1, by = 0.2), 2) # Intercept trails for eq. "x"
                 )

plot_coef_adjust(equation = "xii",
                 coef_mod = "a",
                 coef_vec = rep(seq(-2, 0.1, by = 0.2), 2) # Intercept trails for eq. "xii"
)



##############################################################################*
### 6. COMPARE BIOMASS DATA FROM BIORESTORATION AND OLIVERAS              ####
###     A. Test equations will full biorestoration data. 
###     B. Adjust coefficients from equation manually
###    
##############################################################################*

### Biorestoration data: how much weight lost the sample after drying process 
biomass_dat %>%
  ggplot(aes(x = cut(peso_fresco, breaks = 15),
             y = (peso_fresco - peso_seco)/peso_fresco*100
  )) +
  geom_boxplot(color = "darkgrey") +
  geom_jitter(aes(fill = genero), size = 3, shape = 21) +
  coord_cartesian(ylim = c(0,100)) +
  scale_fill_manual(values = gen.fil) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle(paste0("AGB observations = ", nrow(biomass_dat)))


### Oliveras data extracted from figures and suplementary materials
oliveras_data <- read_excel("data/Default Dataset_fig2.xlsx")
oliveras_data %>%
  ggplot(aes(Observed_10exp_round)) +
  geom_histogram(bins = 25) +
  ggtitle(label = paste0("AGB observations = ", nrow(oliveras_data)))

oliveras_data_festuca <- read_excel("data/Default Dataset_sup002_apendixS07.xlsx")
oliveras_data_festuca %>%
  ggplot(aes(Observed_10exp_round)) +
  geom_histogram(bins = 25) +
  ggtitle(label = paste0("AGB observations = ", nrow(oliveras_data_festuca)))


### Compare datasets
biomass_full <-
  bio_dat %>%
  select(AGB) %>%
  set_names("biorest_2016") %>%
  bind_rows(oliveras_data %>%
              select(Observed_10exp) %>% 
              set_names("oliveras_multispp")) %>%
  bind_rows(oliveras_data_festuca %>%
              select(Observed_10exp) %>%
              set_names("oliveras_festuca")) %>%
  gather(key = "source", value = "peso_observado") %>%
  drop_na(peso_observado) %>%
  mutate(source = fct_rev(source))

biomass_full %>%
  ggplot(aes(peso_observado)) +
  geom_histogram(bins = 25, aes(fill = source), alpha = 0.5) +
  geom_text(data = biomass_full %>% count(source),
            aes(x = 750, y = 150, label = paste0("n = ", n))) +
  facet_grid(cols = vars(source)) + 
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  theme(legend.position = "bottom")

### FIN ###

