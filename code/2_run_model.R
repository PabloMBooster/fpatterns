# =====================================================================
# LENGTH-STRUCTURED MODEL
# North-Central anchoveta stock (Engraulis ringens)
# =====================================================================
library(dplyr)
source("code/0_internal_function.R")
source("code/1_model_functions.R")


# ---------------------------------------------------------------------
# 1. Inputs
# ---------------------------------------------------------------------

load("input/surveyFile.RData")  # Survey abundance by length class
load("input/catchFile.RData") # Weekly catch by length class
load("input/life_history_params_by_envir_scenario.RData")
load("input/species_data.RData")


# ---------------------------------------------------------------------
# 2. Setting
# ---------------------------------------------------------------------

stock                  <- "nc"
species_name           <- "anchoveta"
data_type              <- "boot"
survey_label           <- "survey"
environmental_scenario <- "neutral"

survey_unit <- 1e6  # Survey abundance is expressed in millions of individuals
catch_unit  <- 1    # Catch is expressed as number of individuals

length_marks <- seq(2, 20, by = 0.5)
length_labels <- as.character(length_marks)

# Length–weight relationship parameters: W = a × L^b
species_data$a <- 0.0072
species_data$b <- 2.9855

pars_wl <- list(a = species_data$a, b = species_data$b)

cat("Total catch (t): ", format(sum(catchFile$landing_t), big.mark = ","), "\n")
cat("Total abundance - survey: ", format(sum(survey), big.mark = ","), "\n")

# Catch at length
catch_weekly <- catchFile %>%
  filter(type == "Week")

available_length_columns <- intersect(length_labels, names(catch_weekly))

# Matrix: rows = lengths; columns = weeks
CatchData <- matrix(data = 0, nrow = length(length_labels), ncol = nrow(catch_weekly), 
                    dimnames = list(length_labels, paste0("Week_", catch_weekly$week)))

CatchData[available_length_columns, ] <- t(as.matrix(catch_weekly[, available_length_columns]))

# Abundance at length 
SurveyData <- matrix(data = 0, nrow = length(length_labels), ncol = ncol(survey),
                     dimnames = list(length_labels, colnames(survey)))

n_length_survey <- min(nrow(survey), nrow(SurveyData))

SurveyData[seq_len(n_length_survey), ] <- survey[seq_len(n_length_survey), , drop = FALSE]


# ---------------------------------------------------------------------
# 3.  projection
# ---------------------------------------------------------------------

week_labels <- colnames(CatchData)

output <- lengthBasedModel(survey_data  = SurveyData, 
                           data_type    = data_type, 
                           catch_data   = CatchData, 
                           scenario     = environmental_scenario, 
                           scenarios    = scenarios, 
                           specie       = species_data, 
                           week_labels  = week_labels, 
                           pars_wl      = pars_wl, 
                           survey_label = survey_label, 
                           survey_unit  = survey_unit, 
                           catch_unit   = catch_unit)


# ---------------------------------------------------------------------
# 4. Fishing mortality estimation
# ---------------------------------------------------------------------

M_l <- M_by_length(lengths   = length_marks, 
                   scenario  = environmental_scenario,
                   scenarios = scenarios)

F_model <- estimate_F_matrix(modelN     = output$modelN, 
                             CatchData  = CatchData, 
                             M_l        = M_l, 
                             unitSurvey = survey_unit, 
                             unitCatch  = catch_unit)

indicators <- F_indicators(F_model      = F_model, 
                           modelN       = output$modelN, 
                           adult_length = 12)


# ---------------------------------------------------------------------
# 5. outputs
# ---------------------------------------------------------------------

results_model <- list(CatchData   = CatchData, 
                      SurveyData  = SurveyData, 
                      model       = output, 
                      M_at_length = M_l, 
                      F_at_length = F_model, 
                      indicators  = indicators)

results_model$model$Biomass
results_model$indicators$F_stock
