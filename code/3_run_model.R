# =====================================================================
# LENGTH-STRUCTURED MODEL
# North-Central anchoveta stock (Engraulis ringens)
# =====================================================================
source("code/1_model_functions.R")
source("code/2_internal_function.R")
# ---------------------------------------------------------------------
# 1. Inputs
# ---------------------------------------------------------------------

load("input/abundance.RData")  # Abundancia del crucero por clase de talla
load("input/catch.RData")      # Captura semanal por clase de talla
load("input/life_history_parameters_by_environmental_scenario.RData")
load("input/species_data.RData")

library(dplyr)

# ---------------------------------------------------------------------
# 2. Setting
# ---------------------------------------------------------------------

stock <- "nc"
species_name <- "anchoveta"

data_type <- "boot"
survey_label <- "crucero"

environmental_scenario <- "neutro"

survey_unit <- 1e6  # Survey abundance is expressed in millions of individuals
catch_unit  <- 1    # Catch is expressed as number of individuals


length_marks <- seq(2, 20, by = 0.5)
length_labels <- as.character(length_marks)

# Length–weight relationship parameters: W = a × L^b
species_data$a <- 0.0072
species_data$b <- 2.9855

pars_wl <- list(
  a = species_data$a,
  b = species_data$b
)


cat("Captura total (t):", format(sum(catch$desembarque_t), big.mark = ","), "\n")
cat("Abundancia total del crucero:", format(sum(survey), big.mark = ","), "\n")

# Catch at length
catch_weekly <- catch %>%
  filter(tipo == "Semana")

available_length_columns <- intersect(length_labels, names(catch_weekly))

# Matriz: rows = lengths; columns = weeks
CatchData <- matrix(
  data = 0,
  nrow = length(length_labels),
  ncol = nrow(catch_weekly),
  dimnames = list(
    length_labels,
    paste0("Semana_", catch_weekly$semana)
  )
)

CatchData[available_length_columns, ] <- t(
  as.matrix(catch_weekly[, available_length_columns])
)


# Abundance at length 
SurveyData <- matrix(
  data = 0,
  nrow = length(length_labels),
  ncol = ncol(survey),
  dimnames = list(length_labels, colnames(survey))
)

n_length_survey <- min(nrow(survey), nrow(SurveyData))

SurveyData[seq_len(n_length_survey), ] <- survey[
  seq_len(n_length_survey), ,
  drop = FALSE
]

# ---------------------------------------------------------------------
# 3.  projection
# ---------------------------------------------------------------------

week_labels <- colnames(CatchData)

output <- lengthBasedModel(
  survey_data = SurveyData,
  data_type = data_type,
  catch_data = CatchData,
  scenario = environmental_scenario,
  scenarios = escenarios,
  specie = species_data,
  week_labels = week_labels,
  pars_wl = pars_wl,
  survey_label = survey_label,
  survey_unit = survey_unit,
  catch_unit = catch_unit
)

# ---------------------------------------------------------------------
# 4. Fishing mortality estimation
# ---------------------------------------------------------------------

M_l <- M_talla(
  tallas = length_marks,
  escenario = environmental_scenario,
  escenarios = escenarios
)

F_model <- estimar_F_matriz(
  modelN = output$modelN,
  CatchData = CatchData,
  M_l = M_l,
  unitSurvey = survey_unit,
  unitCatch = catch_unit
)

indicadores <- indicadores_F(
  F_model = F_model,
  modelN = output$modelN,
  talla_adulto = 12
)

# Tasa de explotación semanal basada en la biomasa total
indicadores$u <- catch_weekly$desembarque_t[-1] /
  output$Biomass[-1]

# ---------------------------------------------------------------------
# 5. outputs
# ---------------------------------------------------------------------

results_model <- list(
  CatchData = CatchData,
  SurveyData = SurveyData,
  model = output,
  M_at_length = M_l,
  F_at_length = F_model,
  indicators = indicadores
)

results_model$model$Biomass
results_model$indicators$F_stock
