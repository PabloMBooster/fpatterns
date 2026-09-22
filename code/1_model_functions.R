# -------------------------------------------------------------------
# R/model_projection.R
#
# Size-structured population projection.
# The implementation includes growth, natural mortality, maturity,
# biomass calculation, and bootstrap summaries.
# -------------------------------------------------------------------

# ===================================================================
# 1. model projection
# ===================================================================

projectmodel <- function(N, catch, a, b, k, Linf, sizeM, 
                         vectorM, freq = 52, sp = "anchoveta", 
                         Ts = 1, specie){
  
  # if(missing(species_data)){
  #   stop("You must provide 'species_data'.")
  # }
  
  # if(is.null(rownames(species_data)) ||
  #     !sp %in% rownames(species_data)){
  #   stop("The selected species was not found in 'species_data'.")
  # }
  
  # N should be a matrix: lengths x replicates.
  if(is.null(dim(N))){
    N <- matrix(N, ncol = 1)
  }
  
  N <- as.matrix(N)
  
  # if (is.null(rownames(N))){
  #   stop("N must contain the length classes as row names.")
  # }
  
  # specie <- getSpeciesInfo(
  #   sp = sp,
  #   species_data = species_data
  # )
  
  marks <- createMarks(specie)
  
  # if(nrow(N) != length(marks)){
  #   stop("The number of length classes in N does not match ",
  #        "the classes defined for the species.")
  # }
  #
  # if(length(catch) != nrow(N)){
  #   stop("catch must contain one observation per length class.")
  # }
  
  matrixN <- array(data = NA_real_, dim = c(Ts + 1, nrow(N), ncol(N)), 
                   dimnames = list(paste0("time_", 0:Ts), rownames(N), colnames(N)))
  
  matrixB <- matrix(data = NA_real_, nrow = Ts + 1, ncol = ncol(N), 
                    dimnames = list(paste0("time_", 0:Ts), colnames(N)))
  
  matrixBD <- matrix(data = NA_real_, nrow = Ts + 1, ncol = ncol(N), 
                     dimnames = list(paste0("time_", 0:Ts), colnames(N)))
  
  matrixBDR <- rep(NA_real_, ncol(N))
  
  for(i in seq_len(ncol(N))){
    simulation <- model(N0 = N[, i], catch = catch, a = a, b = b, k = k, 
                        Linf = Linf, sizeM = sizeM, vectorM = vectorM, 
                        freq = freq, specie = specie, Ts = Ts)
    
    matrixN[, , i] <- simulation$N
    matrixB[, i]   <- simulation$B
    matrixBD[, i]  <- simulation$BD
    matrixBDR[i]   <- simulation$BDR
  }
  
  N_median   <- apply(matrixN, MARGIN = c(1, 2), FUN = median, na.rm = TRUE)
  
  B_median   <- apply(matrixB, MARGIN = 1, FUN = median, na.rm = TRUE)
  
  BD_median  <- apply(matrixBD, MARGIN = 1, FUN = median, na.rm = TRUE)
  
  BDR_median <- median(matrixBDR, na.rm = TRUE)
  
  output <- list(N = N_median, B = B_median, BD = BD_median, BDR = BDR_median, 
                 raw = list(N   = matrixN, 
                            B   = matrixB,
                            BD  = matrixBD, 
                            BDR = matrixBDR))
  
  attr(output, "species")    <- sp
  attr(output, "frequency")  <- freq
  attr(output, "time_steps") <- Ts
  
  class(output) <- "surveyProj"
  
  return(output)
}


# ===================================================================
# 2. model's method
# ===================================================================

model <- function(N0, catch, a, b, k, Linf, sizeM, 
                  vectorM, freq, specie, Ts){
  
  A <- lengthProjMatrix(k = k, Linf = Linf, freq = freq, specie = specie)
  
  M <- naturalMortality(vectorM = vectorM, sizeM = sizeM, freq = freq, specie = specie)
  
  marks <- createMarks(specie)
  
  weight <- a * marks^b
  
  maturity <- maturity_ogive(specie)
  
  N <- matrix(data = NA_real_, nrow = Ts + 1, ncol = length(M), 
              dimnames = list(paste0("time_", 0:Ts), marks))
  
  N[1, ] <- N0
  
  for(t in seq_len(Ts)){
    
    # Natural mortality during the first half of the interval.
    survivors_half_step <- N[t, ]*exp(-M/2)
    
    # Model approximation: catch removal at the midpoint.
    remaining_abundance <- survivors_half_step - catch
    
    if(any(remaining_abundance < 0, na.rm = TRUE)){
      warning("Catch exceeds the available abundance in at least ", 
              "one length class. Negative values were truncated to zero.")
    }
    
    remaining_abundance <- pmax(remaining_abundance, 0)
    
    # Natural mortality during the second half of the interval and length transition.
    N[t + 1, ] <- as.numeric((remaining_abundance*exp(-M/2)) %*% A)
  }
  
  B   <- as.numeric(N %*% weight)
  BD  <- as.numeric(N %*% (maturity * weight))
  BDR <- tail(BD, n = 1)
  
  return(list(N = N, C = catch, B = B, BD = BD, BDR = BDR))
}


# ===================================================================
# 3. Weekly model population projection
# ===================================================================

lengthBasedModel <- function(survey_data, data_type, catch_data, 
                             scenario, specie, scenarios, survey_label, 
                             week_labels = NULL, pars_wl, 
                             survey_unit = 1e6, catch_unit = 1e3){
  
  # if(!scenario %in% names(scenarios)){
  #   stop("The selected environmental scenario was not found.")
  # }
  #
  # if(!all(rownames(survey_data) == rownames(catch_data))){
  #   stop("Survey and catch data must contain identical length classes.")
  # }
  
  scenario_parameters <- scenarios[[scenario]]
  
  n_lengths <- nrow(catch_data)
  n_times   <- ncol(catch_data) + 1
  n_boot    <- ncol(survey_data)
  
  # if(is.null(week_labels)){
  #   week_labels <- colnames(catch_data)
  # }
  
  model_array <- array(data = NA_real_, 
                       dim = c(n_lengths, n_times, n_boot), 
                       dimnames = list(rownames(catch_data), 
                                       c(survey_label, week_labels), 
                                       colnames(survey_data)))
  
  progress <- progress::progress_bar$new(format = " Projecting bootstrap samples [:bar] :percent eta: :eta", 
                                         total = n_boot, clear = FALSE, width = 60)
  
  for(bootstrap_id in seq_len(n_boot)){
    
    abundance_current <- matrix(survey_data[, bootstrap_id] * survey_unit, ncol = 1, 
                                dimnames = list(rownames(survey_data), survey_label))
    
    for(week_id in seq_len(ncol(catch_data))){
      
      model_projection <- projectmodel(N = abundance_current[, week_id, drop = FALSE], 
                                       catch = catch_data[, week_id] * catch_unit, 
                                       a = pars_wl$a, b = pars_wl$b, 
                                       k = scenario_parameters$k, 
                                       Linf = scenario_parameters$Linf, 
                                       sizeM = scenario_parameters$sizeM, 
                                       vectorM = scenario_parameters$vectorM,
                                       # Weekly temporal scale
                                       freq = 52, 
                                       specie = specie, 
                                       Ts = 1)
      
      abundance_next <- model_projection$N[2, ]
      
      abundance_current <- cbind(abundance_current, abundance_next)
    }
    
    colnames(abundance_current) <- c(survey_label, week_labels)
    
    model_array[, , bootstrap_id] <- abundance_current/survey_unit
    
    progress$tick()
  }
  
  # ---------------------------------------------------------------
  # Abundance summaries
  # ---------------------------------------------------------------
  
  # Dimensions: length x bootstrap replicate x time step
  model_array_N <- aperm(model_array, perm = c(1, 3, 2))
  
  model_median_N <- apply(model_array_N, MARGIN = c(1, 3), 
                          FUN = median, na.rm = TRUE)
  
  # ---------------------------------------------------------------
  # Biomass summaries
  # ---------------------------------------------------------------
  
  #species_info <- specie#modelProjection:::getSpeciesInfo(species)
  
  length_marks <- createMarks(specie)
  
  weight_at_length <- pars_wl$a * (length_marks ^ pars_wl$b)
  
  model_array_B <- apply(model_array_N, MARGIN = c(2, 3), 
                         FUN = function(x) x * weight_at_length)
  
  model_median_B <- apply(model_median_N, MARGIN = 2, 
                          FUN = function(x) x * weight_at_length)
  
  total_biomass <- colSums(model_median_B, na.rm = TRUE)
  
  # ---------------------------------------------------------------
  # Output
  # ---------------------------------------------------------------
  
  if(data_type == "boot"){
    output <- list(modelArrayN = model_array_N, modelN = model_median_N, 
                   modelArrayB = model_array_B, modelB = model_median_B, 
                   Biomass = total_biomass)
  }else{
    output <- list(modelN = model_median_N, modelB = model_median_B, 
                   Biomass = total_biomass)
  }
  
  attr(output, "freq") <- 52
  attr(output, "time_step") <- "weekly"
  
  message("Abundance unit: ", survey_unit, " individuals\n", 
          "Biomass unit: tonnes\n", 
          "Time step: weekly (52 steps per year)")
  
  return(output)
}
