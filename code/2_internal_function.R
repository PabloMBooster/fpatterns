
# ===================================================================
# Internal functions
# ===================================================================

# getSpeciesInfo <- function(sp,
#                            species_data) {
#   
#   required_columns <- c(
#     "Lmin",
#     "Lmax",
#     "bin",
#     "mat1",
#     "mat2"
#   )
#   
#   if (!all(required_columns %in% names(species_data))) {
#     stop(
#       "species_data debe incluir las columnas: ",
#       paste(required_columns, collapse = ", ")
#     )
#   }
#   
#   specie <- species_data[sp, required_columns, drop = FALSE]
#   
#   specie <- lapply(
#     specie,
#     function(x) as.numeric(as.character(x))
#   )
#   
#   return(specie)
# }


createMarks <- function(specie,
                         phi = FALSE) {
  
  marks <- seq(
    from = specie$Lmin,
    to = specie$Lmax,
    by = specie$bin
  )
  
  if (isTRUE(phi)) {
    
    marks_lower <- marks - 0.5 * specie$bin
    marks_upper <- marks + 0.5 * specie$bin
    
    marks <- sort(
      unique(
        c(marks_lower, marks_upper)
      )
    )
  }
  
  return(marks)
}


brody <- function(l,
                  Linf,
                  k) {
  
  Linf - (Linf - l) * exp(-k)
}


lengthProj <- function(l,
                        marks) {
  
  x <- sort(
    c(l, marks)
  )
  
  differences <- diff(x)
  
  position <- findInterval(
    l,
    marks
  ) + 1
  
  position <- seq(
    from = position[1],
    to = position[2]
  )
  
  proportions <- differences[position]
  proportions <- proportions / sum(proportions)
  
  n_lengths <- length(marks) - 1
  
  output <- numeric(
    n_lengths
  )
  
  position_output <- position - 1
  
  # Todo el crecimiento que excede la última clase
  # se acumula en la clase terminal.
  if (any(position_output > n_lengths)) {
    
    output[n_lengths] <- 1
    
  } else {
    
    output[position_output] <- proportions
  }
  
  return(output)
}


lengthProjMatrix <- function(
    k,
    Linf,
    freq,
    specie
) {
  
  dt <- 1 / freq
  
  marks <- createMarks(specie)
  marks_phi <- createMarks(specie, phi = TRUE)
  
  marks_lower <- marks - 0.5 * specie$bin
  marks_upper <- marks + 0.5 * specie$bin
  
  k_dt <- k * dt
  
  length_lower <- brody(
    l = marks_lower,
    Linf = Linf,
    k = k_dt
  )
  
  length_upper <- brody(
    l = marks_upper,
    Linf = Linf,
    k = k_dt
  )
  
  projected_bounds <- cbind(
    length_lower,
    length_upper
  )
  
  A <- t(
    apply(
      projected_bounds,
      MARGIN = 1,
      FUN = lengthProj,
      marks = marks_phi
    )
  )
  
  rownames(A) <- marks
  colnames(A) <- marks
  
  return(A)
}


naturalMortality <- function(
    vectorM,
    sizeM,
    freq,
    specie
) {
  
  dt <- 1 / freq
  
  marks <- createMarks(specie)
  
  mortality_group <- findInterval(
    marks,
    sizeM
  )
  
  mortality_group[mortality_group == 0] <- 1
  
  M <- vectorM[mortality_group] * dt
  
  names(M) <- marks
  
  return(M)
}


maturity_ogive <- function(specie) {
  
  marks <- createMarks(specie)
  
  maturity <- 1 / (
    1 + exp(specie$mat1 + specie$mat2 * marks)
  )
  
  return(maturity)
}


model <- function(
    N0,
    catch,
    a,
    b,
    k,
    Linf,
    sizeM,
    vectorM,
    freq,
    specie,
    Ts
) {
  
  A <- lengthProjMatrix(
    k = k,
    Linf = Linf,
    freq = freq,
    specie = specie
  )
  
  M <- naturalMortality(
    vectorM = vectorM,
    sizeM = sizeM,
    freq = freq,
    specie = specie
  )
  
  marks <- createMarks(specie)
  
  weight <- a * marks^b
  
  maturity <- maturity_ogive(specie)
  
  N <- matrix(
    data = NA_real_,
    nrow = Ts + 1,
    ncol = length(M),
    dimnames = list(
      paste0("time_", 0:Ts),
      marks
    )
  )
  
  N[1, ] <- N0
  
  for (t in seq_len(Ts)) {
    
    # Mortalidad natural durante la primera mitad del intervalo.
    survivors_half_step <- N[t, ] * exp(-M / 2)
    
    # Aproximación de model: extracción de la captura en el punto medio.
    remaining_abundance <- survivors_half_step - catch
    
    if (any(remaining_abundance < 0, na.rm = TRUE)) {
      warning(
        "La captura supera la abundancia disponible en al menos ",
        "una clase de talla. Los valores negativos se truncaron en cero."
      )
    }
    
    remaining_abundance <- pmax(
      remaining_abundance,
      0
    )
    
    # Mortalidad natural durante la segunda mitad y transición de talla.
    N[t + 1, ] <- as.numeric(
      (remaining_abundance * exp(-M / 2)) %*% A
    )
  }
  
  B <- as.numeric(
    N %*% weight
  )
  
  BD <- as.numeric(
    N %*% (maturity * weight)
  )
  
  BDR <- tail(
    BD,
    n = 1
  )
  
  return(
    list(
      N = N,
      C = catch,
      B = B,
      BD = BD,
      BDR = BDR
    )
  )
}



#-------------------------------------------------------
# Mortalidad natural por talla
#-------------------------------------------------------

M_talla <- function(tallas,
                    escenario = "neutro",
                    escenarios){
  
  # Verificar escenario
  if(!escenario %in% names(escenarios)){
    stop("Escenario no encontrado.")
  }
  
  esc <- escenarios[[escenario]]
  
  # Determinar en qué bloque cae cada talla
  bloque <- findInterval(tallas, esc$sizeM)
  
  # Asignar M correspondiente al bloque
  M <- esc$vectorM[bloque]
  
  names(M) <- tallas
  
  return(M)
}
#----------------

#-------------------------------------------------------
# Mortalidad por pesca por
#-------------------------------------------------------
#---------------------------------------------------------
# Estima F a partir de N, C y M usando Baranov
# N y C deben estar en las MISMAS unidades
#---------------------------------------------------------

estimar_F <- function(N, C, M,
                      dt = 1,
                      Fmax = 20){
  
  # Casos triviales
  if(is.na(N) || is.na(C) || is.na(M))
    return(NA_real_)
  
  if(N <= 0 || C <= 0)
    return(0)
  
  # La captura no puede superar la abundancia
  if(C >= N)
    return(NA_real_)
  
  fun <- function(F){
    
    N *
      (F/(F+M)) *
      (1-exp(-(F+M)*dt)) -
      C
    
  }
  
  # Verificar que exista una raíz
  f0   <- fun(0)
  fmax <- fun(Fmax)
  
  if(is.na(f0) || is.na(fmax))
    return(NA_real_)
  
  if(f0 * fmax > 0)
    return(NA_real_)
  
  uniroot(
    fun,
    interval = c(0, Fmax),
    tol = 1e-8
  )$root
  
}

#---------------------------------------------------------
# Mortalidad por pesca por talla y semana
#---------------------------------------------------------

estimar_F_matriz <- function(modelN,
                             CatchData,
                             M_l,
                             unitSurvey = 1e6,
                             unitCatch = 1,
                             dt = 1){
  
  # Ambas matrices quedan en las mismas unidades
  N <- modelN[, -1] * unitSurvey / unitCatch
  C <- CatchData
  
  F <- matrix(
    0,
    nrow = nrow(C),
    ncol = ncol(C),
    dimnames = dimnames(C)
  )
  
  for(i in seq_len(nrow(C))){
    
    for(j in seq_len(ncol(C))){
      
      F[i,j] <- estimar_F(
        N = N[i,j],
        C = C[i,j],
        M = M_l[i],
        dt = dt
      )
      
    }
    
  }
  
  return(F)
  
}


#----------------------------------------------------------
# Indicadores de mortalidad por pesca
#----------------------------------------------------------

indicadores_F <- function(F_model,
                          modelN,
                          talla_adulto = 12){
  
  tallas <- as.numeric(rownames(F_model))
  
  # Abundancia (elimina la columna del crucero)
  N <- modelN[, -1]
  
  adultos   <- tallas >= talla_adulto
  juveniles <- tallas < talla_adulto
  
  #-------------------------
  # F adultos (ponderado)
  #-------------------------
  
  F_adultos <- colSums(
    F_model[adultos, ] * N[adultos, ],
    na.rm = TRUE
  ) /
    colSums(
      N[adultos, ],
      na.rm = TRUE
    )
  
  #-------------------------
  # F juveniles (ponderado)
  #-------------------------
  
  F_juveniles <- colSums(
    F_model[juveniles, ] * N[juveniles, ],
    na.rm = TRUE
  ) /
    colSums(
      N[juveniles, ],
      na.rm = TRUE
    )
  
  #-------------------------
  # F stock (ponderado)
  #-------------------------
  
  F_stock <- colSums(
    F_model * N,
    na.rm = TRUE
  ) /
    colSums(
      N,
      na.rm = TRUE
    )
  
  #-------------------------
  # F máximo semanal
  #-------------------------
  
  F_max <- apply(
    F_model,
    2,
    max,
    na.rm = TRUE
  )
  
  #-------------------------
  # F instantaneo
  #-------------------------
  
  F_inst <- apply(
    F_model,
    2,
    sum,
    na.rm = TRUE
  )
  
  #-------------------------
  # Selectividad observada
  #-------------------------
  
  selectividad <- rowMeans(
    F_model,
    na.rm = TRUE
  )
  
  selectividad <- selectividad /
    max(selectividad, na.rm = TRUE)
  
  names(selectividad) <- tallas
  
  #-------------------------
  # F acumulado por talla
  #-------------------------
  
  F_acumulado <- rowSums(
    F_model,
    na.rm = TRUE
  )
  
  names(F_acumulado) <- tallas
  
  #-------------------------
  # Salida
  #-------------------------
  
  list(
    
    F_adultos     = F_adultos,
    F_juveniles   = F_juveniles,
    F_stock       = F_stock,
    F_max         = F_max,
    selectividad  = selectividad,
    F_acumulado   = F_acumulado
    
  )
  
}
