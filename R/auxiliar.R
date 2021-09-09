# @importFrom stats runif
# La función **simular_vector_convexo()** simula un vector de proporciones de suma uno
# (vector convexo) a partir de perturbar aleatoriamente un vector base de proporciones de
# suma uno. La perturbación, basada en una distribución uniforme, está modulada por el
# parámetro `d`.
simular_vector_convexo <- function(vector.convexo, d){
  output <- vector.convexo + stats::runif(length(vector.convexo), -d, d)
  output[output > 1L] <- 1L
  output[output < 0L] <- 0L
  if (sum(output) == 0L) output <- rep(1/length(vector.convexo), length(vector.convexo))
  output <- output/sum(output)
  return(output)
}

# La función **simular_matrices_fila_convexas()** simula un array de matrices de vectores
# fila convexos a partir de perturbar aleatoriamente una matriz base fila estandarizada.
# La perturbación, basada en una distribución uniforme, está modulada por el parámetro `d`.
simular_matrices_fila_convexas <- function(matriz, d, n.unidades){
  output <- array(NA, c(dim(matriz), n.unidades))
  for (i in 1L:n.unidades)
    output[, , i] <- t(apply(matriz, 1, simular_vector_convexo, d = d))
  return(output)
}

# La función **determinar_zeros_estructurales()** determina las coordenadas de las celdas
# que corresponden a ceros estructuales de una matriz de tranferencia resultado de aplicar
# la función **lphom()**.
determinar_zeros_estructurales <- function(lphom.object){
  output <- lphom.object$inputs$structural_zeros
  escenario <- lphom.object$inputs$new_and_exit_voters[1L]
  J.inic <- ncol(lphom.object$inputs$votes_election1)
  J.final <- ncol(lphom.object$origin)
  K.inic <- ncol(lphom.object$inputs$votes_election2)
  K.final <- ncol(lphom.object$destination)
  if ((escenario == "raw") & (J.inic < J.final) & (K.inic < K.final)){
    output[[length(output) + 1L]] <- c(J.final, K.final)
  }
  if ((escenario == "regular") & (K.inic < K.final) & (J.inic == J.final)){
    output[[length(output) + 1L]] <- c(J.final, K.final)
  }
  if ((escenario == "regular") & (K.inic < K.final) & (J.inic < J.final)){
    output[[length(output) + 1L]] <- c(J.final - 1L, K.final)
    output[[length(output) + 1L]] <- c(J.final, K.final)
  }
  if (escenario == "full"){
    output[[length(output) + 1L]] <- c(J.inic - 1L, K.inic)
    output[[length(output) + 1L]] <- c(J.inic, K.inic)
  }
  if (escenario == "gold"){
    output[[length(output) + 1L]] <- c(J.inic - 1L, K.inic)
    output[[length(output) + 1L]] <- c(J.inic - 1L, K.inic - 1L)
    output[[length(output) + 1L]] <- c(J.inic, K.inic)
    output[[length(output) + 1L]] <- c(J.inic, K.inic - 1L)
  }
  return(output)
}

# La función **simular_arrays_transferencia()** simula un array de matrices de transferencia
# perturbando aleatoriamente una matriz de transferencia calculada mediante **lphom()**.
# La perturbación, basada en una distribución uniforme, está modulada por el parámetro `d` y
# respeta las restricciones de ceros estructurales asociados a la matriz de transferencia
# obtenida por **lphom()**. Se simula una matriz de transferencia para cada unidad territorial.
simular_arrays_transferencia <- function(lphom.object, d){
  ceros <- determinar_zeros_estructurales(lphom.object)
  output <- simular_matrices_fila_convexas(matriz = lphom.object$VTM.complete,
                                           d = d,
                                           n.unidades = nrow(lphom.object$origin))
  if(!is.null(ceros)){
    for (ii in 1L:length(ceros)){
      output[ceros[[ii]][1L], ceros[[ii]][2L], ] <- 0L
    }
    for (ii in 1L:nrow(lphom.object$origin)){
      output[, , ii] <- output[, , ii]/rowSums(output[, , ii])
    }
  }
  return(output)
}

# La función **simular_arrays_votos()** simula un array de matrices de votos aplicando
# sobre los vectores de votos de la elección de origen en cada unidad territorial un array
# de proporciones obtenido perturbando aleatoriamente una matriz de transferencia calculada
# mediante **lphom()** utilizando la función **simular_arrays_transferencia()**.
simular_arrays_votos <- function(lphom.object, d){
  output <- simular_arrays_transferencia(lphom.object = lphom.object, d = d)
  for (ii in 1L:dim(output)[3L]){
    for (jj in 1L:dim(output)[1L]){
      output[jj, , ii] <- lphom.object$origin[ii, jj] * output[jj, , ii]
    }
  }
  return(output)
}

# Dado un array de matrices origen-destino de votos correspondientes a un conjunto de unidades
# territoriales, la función **HET_MT.votos_MT.prop_Y()** calcula: (i) el índice de heterogeneidad
# como está definido en la ecuación (10) de Romero et al. (2020), (ii) la matriz de origen-destino
# global que corresponde al array, (iii) la matriz de transferencia global que corresponde al
# array y (iv) la matriz de destino de votos en cada una de las unidades territoriales.
# The heterogeneity index accounts for the percentage of voters which should be shifted to
# match perfectly the homogeneity hypothesis. Es decir, lo distintas que son las matrices de
# transferencia en cada unidad territorial.
HET_MT.votos_MT.prop_Y <- function(array.votos){
  mt.votos <- apply(array.votos, c(1,2), sum)
  mt.prop <- mt.votos/rowSums(mt.votos)
  # En algunas remuestras boosting podrían haber opciones sin votos
  mt.prop[is.nan(mt.prop)] <- 1/ncol(mt.prop)
  origen <- t(apply(array.votos, c(1,3), sum))
  array.homogeneo <- array(NA, dim(array.votos))
  for (ii in 1L:dim(array.homogeneo)[3L]){
    for (jj in 1L:dim(array.homogeneo)[1L]){
      array.homogeneo[jj, , ii] <- origen[ii, jj] * mt.prop[jj, ]
    }
  }
  HET <- 100*(0.5*sum(abs(array.homogeneo - array.votos)))/sum(array.votos)
  EHet <- t(apply(array.homogeneo - array.votos, c(2,3), sum))
  Error.H <- apply(abs(array.homogeneo - array.votos), 3, sum)
  output <- list("HET" = HET, "MT.votos" = mt.votos, "MT.prop" = mt.prop, "EHet" = EHet,
                 "Error.H" = Error.H, "Y" = t(apply(array.votos, c(2,3), sum)))
  return(output)
}

# La función **EI_index()** calcula el índice de discrepancia o de error, siguiendo la
# definición de la ecuación (12) de Romero et al. (2020), entre la matriz origen-destino
# de votos real (que es conocida en cada escenario simulado), y una matriz origen-destino
# de votos estimada, por ejmplo, mediante **lphom()**.
EI_index <- function(matriz.real, matriz.estimada){
  EI <- 50*sum(abs(matriz.real - matriz.estimada))/sum(matriz.real)
  return(EI)
}

# La función **simular_y_resumir()** genera un escenario simulado y resume el escenario
# simulado con cuatro estadísticos: (i) el estadístico HET calculado utilizando la función
# **HET_MT.votos_MT.prop_Y()**, (ii) el estadístico HETe calculado a partir de los datos
# del escenario simulado calculado utilizando **lphom()**, (iii) el estadístico IE calculado
# utilizando la función **IE_index()** y (iv) el valor de perturbación d utilizado en el escenario.
simular_y_resumir <- function(lphom.object, d){
  ceros <- determinar_zeros_estructurales(lphom.object)
  array.votos <- simular_arrays_votos(lphom.object, d)
  resumenes <- HET_MT.votos_MT.prop_Y(array.votos)
  # Ajuste de decimales
  resumenes$Y[,1L] <- resumenes$Y[,1L] + (rowSums(lphom.object$origin)-rowSums(resumenes$Y))
  negativos <- resumenes$Y[, 1L] < 0
  origen <- lphom.object$origin
  origen[negativos, 1L] <- origen[negativos, 1L] - resumenes$Y[negativos, 1L]
  resumenes$Y[negativos, 1L] <- 0
  # Estimacion
  estimacion <- suppressMessages(lphom(votes_election1 = origen,
                                       votes_election2 = resumenes$Y,
                                       new_and_exit_voters = "raw",
                                       structural_zeros = ceros,
                                       counts = lphom.object$inputs$counts,
                                       verbose = FALSE,
                                       lphom.object$inputs$solver))
  estimacion$VTM.complete <- estimacion$VTM.complete[1L:nrow(resumenes$MT.prop),
                                                     1L:ncol(resumenes$MT.prop)]
  estimacion$VTM.complete <- estimacion$VTM.complete/rowSums(estimacion$VTM.complete)
  matriz.estimada <- estimacion$VTM.complete*colSums(lphom.object$origin)
  EI <- EI_index(resumenes$MT.votos, matriz.estimada)
  resumen <- c(resumenes$HET, estimacion$HETe, EI, d)
  names(resumen) <- c("HET", "HETe", "EI", "d")
  output <- list("resumen" = resumen, "MT.real" = resumenes$MT.prop,
                 "MT.estimada" = estimacion$VTM.complete)
  return(output)
}

# La función **simular_escenarios()** simula, para cada elemento de un mallado de perturbaciones,
# **M.d**, **B** escenarios para cada valor de **d** y salva los estadísticos de resumen asociados
# a cada escenario simulado y las matrices de transferencia reales simuladas y estimadas.
simular_escenarios <- function(lphom.object, M.d, B = 30){
  estadisticos <- NULL
  MT.reales <- MT.estimadas <- array(NA, c(dim(lphom.object$VTM.complete), length(M.d)*B))
  for (ii in 1L:length(M.d)){
    for(bb in 1L:B){
      escenario <- simular_y_resumir(lphom.object, M.d[ii])
      estadisticos <- rbind( estadisticos, escenario$resumen )
      MT.reales[, , (ii-1L)*B + bb] <- escenario$MT.real
      MT.estimadas[, , (ii-1L)*B + bb] <- escenario$MT.estimada
    }
  }
  output <- list("estadisticos" = estadisticos, "MT.reales" = MT.reales, "MT.estimadas" = MT.estimadas)
  return(output)
}

# Función **hallar_d0()**. Para poder obtener una aproximación de la incertidumbre (error)
# asociada a la estimación de una matriz de transferencia obtenida mediante **lphom()**,
# medida a través del coeficiente EI, se utiliza un procedimiento basado en simulación.
# En concreto, (i) se estudia la relación entre la distancia (el índice de error EI)
# entre las matrices simuladas (próximas a la matriz obtenida por **lphom()**) y
# las estimadas por **lphom()** de tales matrices simuladas y el índice de heterogeneidad
# estimado por **lphom()**, HETe, para esas mismas matrices simuladas, y (ii) se utiliza
# la relación obtenida para lograr una estimación de EI a partir del coeficiente HETe obtenido
# en nuestra elección de estudio (HETe0). La relación debe estudiarse en un entorno del
# coeficiente HETe0. La función **hallard0()** busca, para nuestro caso de estudio, el valor d0
# (para las simulaciones) que corresponde aproximadamente a HETe0 y el incremento id0 de la
# perturbación **d** que corresponde aproximadamente a un incremento de 1 punto en el HETe resultante.
hallar_d0 <- function(lphom.object){
  M.inic <- c(.01, .03, .05, .08, .1, .15, .2, .25, .3)
  estadisticos <- simular_escenarios(lphom.object = lphom.object,
                                     M.d = M.inic, B = 3)$estadisticos[, c(2,4)]
  iter <- 1L
  while(lphom.object$HETe > max(estadisticos[,1L]) & (iter < 5)){
    M.inic <- 2*M.inic
    estadisticos <- simular_escenarios(lphom.object = lphom.object,
                                       M.d = M.inic, B = 3)$estadisticos[, c(2,4)]
    iter <- iter + 1L
  }
  ajuste <- stats::lm(estadisticos[,1L] ~ estadisticos[, 2L])
  d0 <- (lphom.object$HETe - ajuste$coef[1L])/ajuste$coef[2L]
  id0 <- 1/ajuste$coef[2L]
  output <- c(d0, id0)
  names(output) <- c("d0", "id0")
  return(output)
}


# La función **simular_errores_estimacion_pjk()** está basada en la función **simular_escenarios()**,
# guarda los índices de heterogeneidad estimados, los verdaderos valores simulados, los valores estimados
# utilizando **lphom()** y los errores de estimación asociados a un coeficiente pjk específico
# (coordenadas j y k). Toma una salida de la función **simular_escenarios()** y genera una matriz
# de coeficientes de heterogenidad, simulados, estimados y errores.

simular_errores_estimacion_pjk <- function(sim_esc.object, j, k){
  HETe <- sim_esc.object$estadisticos[,2L]
  pjk.simulado <- sim_esc.object$MT.reales[j, k, ]
  pjk.estimado <- sim_esc.object$MT.estimadas[j, k, ]
  sesgojk <- pjk.estimado - pjk.simulado
  output <- cbind(HETe, pjk.simulado, pjk.estimado, sesgojk)
  return(output)
}

# La función **modelos_ajuste_estimacion_pjk()** toma un conjunto de escenarios simulados y
# ajusta las regresiones que sirven de base para construir los intervalos de confianza de
# los coeficientes pjk estimados.
modelos_ajuste_estimacion_pjk <- function(matriz.errores, HETe0){
  he <- matriz.errores[, 1L]
  he2 <- he^2
  error <- matriz.errores[, 4L]
  ajuste.lineal <- stats::lm(error ~ he + he2 + 0L)
  eom <- ajuste.lineal$coef[1L]*HETe0 + ajuste.lineal$coef[2L]*HETe0^2
  resid <- error - ajuste.lineal$coef[1L]*he + ajuste.lineal$coef[2L]*he2
  df <- data.frame(resid, he, he2)
  df <- df[stats::complete.cases(df),]
  ajuste.varianza <- stats::lm(resid ~ he + he2 + 0L, data = df)
  seo <- sqrt(max(ajuste.varianza$coef[1L]*HETe0 + ajuste.varianza$coef[2L]*HETe0^2, 0L))
  output <- c(eom, seo)
  return(output)
}

#' @importFrom Rsymphony Rsymphony_solve_LP
#' @importFrom lpSolve lp
# La función **calculo_MT_unidad()** calcula la matriz de transferencia que corresponde a una
# unidad de votación. Para ello toma como base una matriz de transferencia global y tratando de
# mover esta solución lo menos posible busca otra que verique las restricciones asociadas a la unidad
# de votación. Es decir, que aplicada sobre los resultados de origen de la unidad genere exactamente
# los resultados de destino de la unidad. Esta solucion puede ser indeterminada.
calculo_MT_unidad <- function(lphom.object, iii, solver){
  xt <- lphom.object$origin[iii, ]
  yt <- lphom.object$destination[iii,]
  pg <- lphom.object$VTM.complete/rowSums(lphom.object$VTM.complete)
  ceros <- determinar_zeros_estructurales(lphom.object)
  # Parameters
  nj <- length(xt)
  nk <- length(yt)
  njk <- nj * nk

  # restricciones sum(pjk)=1
  a1 <- kronecker(diag(nj), t(rep(1L, nk)))
  b1 <- rep(1L, nj)

  # restricciones cuadre total de votos para los nk partidos
  at <- t(kronecker(xt, diag(nk)))
  bt <- yt

  # restricciones para definir los ejk
  ajk <- cbind(kronecker(diag(xt), diag(nk)), t(kronecker(diag(njk), c(1L,-1L))))
  bjk <- as.vector(t(xt * pg))

  # Sintesis
  a <- rbind(cbind(rbind(a1, at), matrix(0L, nj+nk,2L*njk)), ajk)
  b <- c(b1, bt, bjk)

  # restricciones ceros estructurales
  if (length(ceros) > 0){
    ast = matrix(0L, length(ceros), ncol(a))
    bst = rep(0L, length(ceros))
    for (i in 1L:length(ceros)){
      ast[i, nk*(ceros[[i]][1L]-1L) + ceros[[i]][2L]] <- 1L
    }
    a <- rbind(a, ast)
    b <- c(b, bst)
  }
  # calculo bounds
#  lb <- matrix(0, 3*njk ,1)
#  ub <- matrix(Inf, 3*njk, 1)
  # Objective function, to minimize
  fun.obj <- c(rep(0L, njk), rep(1L, 2L*njk))
  # Solution
  if (solver == "lp_solve"){
    sol <- suppressWarnings(lpSolve::lp('min', fun.obj, a, rep('=', length(b)), b))
  } else {
    sol <- Rsymphony::Rsymphony_solve_LP(obj = fun.obj,
                                       mat = a,
                                       dir = rep('==', length(b)),
                                       rhs = b)
  }
#  z <- sol$solution
  # Output programa lineal 1
  sol1 <- matrix(sol$solution[1:njk], nj, nk, TRUE,
                 dimnames=dimnames(lphom.object$VTM.complete))
  return(sol1)
}


#' @importFrom Rsymphony Rsymphony_solve_LP
#' @importFrom lpSolve lp
# La función **calculo_MT_unidad_abs()** calcula la matriz de transferencia que corresponde a una
# unidad de votación. Para ello toma como base una matriz de transferencia global y tratando de
# mover esta solución lo menos posible busca otra que verique las restricciones asociadas a la unidad
# de votación. Es decir, que aplicada sobre los resultados de origen de la unidad genere exactamente
# los resultados de destino de la unidad. Esta solucion, que puede ser indeterminada, es refinada en una
# segunda etapa eligiendo entre todas las matrices que cumplen las restricciones,
# la matriz con menor distancia L1 con la matriz original.
calculo_MT_unidad_abs <- function(lphom.object, iii, solver){
  xt <- lphom.object$origin[iii, ]
  yt <- lphom.object$destination[iii,]
  pg <- lphom.object$VTM.complete/rowSums(lphom.object$VTM.complete)
  ceros <- determinar_zeros_estructurales(lphom.object)
  # Parameters
  nj <- length(xt)
  nk <- length(yt)
  njk <- nj * nk

  # restricciones sum(pjk)=1
  a1 <- kronecker(diag(nj), t(rep(1L, nk)))
  b1 <- rep(1L, nj)

  # restricciones cuadre total de votos para los nk partidos
  at <- t(kronecker(xt, diag(nk)))
  bt <- yt

  # restricciones para definir los ejk
  ajk <- cbind(kronecker(diag(xt), diag(nk)), t(kronecker(diag(njk), c(1L,-1L))))
  bjk <- as.vector(t(xt * pg))

  # Sintesis
  a <- rbind(cbind(rbind(a1, at), matrix(0, nj+nk,2*njk)), ajk)
  b <- c(b1, bt, bjk)

  # restricciones ceros estructurales
  if (length(ceros) > 0){
    ast = matrix(0L, length(ceros), ncol(a))
    bst = rep(0L, length(ceros))
    for (i in 1L:length(ceros)){
      ast[i, nk*(ceros[[i]][1L] - 1L) + ceros[[i]][2L]] <- 1L
    }
    a <- rbind(a, ast)
    b <- c(b, bst)
  }
  # calculo bounds
#  lb <- matrix(0, 3*njk ,1)
#  ub <- matrix(Inf, 3*njk, 1)
  # Objective function, to minimize
  fun.obj <- c(rep(0L, njk), rep(1L, 2L*njk))
  # Solution
  if (solver == "lp_solve"){
    sol <- suppressWarnings(lpSolve::lp('min', fun.obj, a, rep('=', length(b)), b))
  } else {
    sol <- Rsymphony::Rsymphony_solve_LP(obj = fun.obj,
                                       mat = a,
                                       dir = rep('==', length(b)),
                                       rhs = b)
  }
  #z <- sol$solution
  # Output programa lineal 1
  sol1 <- matrix(sol$solution[1L:njk], nj, nk, TRUE,
                 dimnames=dimnames(lphom.object$VTM.complete))
  # Segundo programa lineal
  a <- rbind(a, fun.obj)
  b <- c(b, sol$objval)
  aa <- cbind(diag(njk),
               matrix(0L, njk, 2L*njk),
               a[(nj+nk+1L):(nj+nk+njk),(njk+1L):(3L*njk)])
  na <- cbind(a, matrix(0L, nrow(a), 2L*njk))
  na <- rbind(na, aa)
  pgp <- t(pg)
  nb <- c(b, as.vector(pgp))
  nf <- c(rep(0L, 3L*njk), rep(1L, 2L*njk))
  if (solver == "lp_solve"){
    nsol <- suppressWarnings(lpSolve::lp('min', nf, na, rep('=', length(nb)), nb))
  }else {
    nsol <- Rsymphony::Rsymphony_solve_LP(obj = nf,
                                       mat = na,
                                       dir = rep('==', length(nb)),
                                       rhs = nb)
  }
  #  nz <- nsol$solution
  # Output programa lineal 2
  sol2 <- matrix(nsol$solution[1L:njk], nj, nk, TRUE,
                 dimnames=dimnames(lphom.object$VTM.complete))
  return(sol2)
}

#' @importFrom Rsymphony Rsymphony_solve_LP
#' @importFrom lpSolve lp
# La función **calculo_MT_unidad_max()** calcula la matriz de transferencia que corresponde a una
# unidad de votación. Para ello toma como base una matriz de transferencia global y tratando de
# mover esta solución lo menos posible busca otra que verique las restricciones asociadas a la unidad
# de votación. Es decir, que aplicada sobre los resultados de origen de la unidad genere exactamente
# los resultados de destino de la unidad. Esta solucion, que puede ser indeterminada, es refinada en una
# segunda etapa eligiendo entre todas las matrices que cumplen las restricciones,
# la matriz con menor distancia L_Inf con la matriz original.
calculo_MT_unidad_max <- function(lphom.object, iii, solver){
  xt <- lphom.object$origin[iii, ]
  yt <- lphom.object$destination[iii,]
  pg <- lphom.object$VTM.complete/rowSums(lphom.object$VTM.complete)
  ceros <- determinar_zeros_estructurales(lphom.object)
  # Parameters
  nj <- length(xt)
  nk <- length(yt)
  njk <- nj * nk

  # restricciones sum(pjk)=1
  a1 <- kronecker(diag(nj), t(rep(1L, nk)))
  b1 <- rep(1L, nj)

  # restricciones cuadre total de votos para los nk partidos
  at <- t(kronecker(xt, diag(nk)))
  bt <- yt

  # restricciones para definir los ejk
  ajk <- cbind(kronecker(diag(xt), diag(nk)), t(kronecker(diag(njk), c(1L,-1L))))
  bjk <- as.vector(t(xt * pg))

  # Sintesis
  a <- rbind(cbind(rbind(a1, at), matrix(0L, nj+nk, 2L*njk)), ajk)
  b <- c(b1, bt, bjk)

  # restricciones ceros estructurales
  if (length(ceros) > 0){
    ast = matrix(0L, length(ceros), ncol(a))
    bst = rep(0L, length(ceros))
    for (i in 1L:length(ceros)){
      ast[i, nk*(ceros[[i]][1L]-1L) + ceros[[i]][2L]] <- 1L
    }
    a <- rbind(a, ast)
    b <- c(b, bst)
  }
  # Objective function, to minimize
  fun.obj <- c(rep(0L, njk), rep(1L, 2L*njk))
  # Solution
  if (solver == "lp_solve"){
    sol <- suppressWarnings(lpSolve::lp('min', fun.obj, a, rep('=', length(b)), b))
  } else {
    sol <- Rsymphony::Rsymphony_solve_LP(obj = fun.obj,
                                       mat = a,
                                       dir = rep('==', length(b)),
                                       rhs = b)
  }
  #z <- sol$solution
  # Output programa lineal 1
  sol1 <- matrix(sol$solution[1L:njk], nj, nk, TRUE,
                 dimnames=dimnames(lphom.object$VTM.complete))

  # Segundo programa lineal
  a <- rbind(a, fun.obj)
  b <- c(b, sol$objval)
  a <- cbind(a, rep(0L, nrow(a)))
  ad <- cbind(diag(njk),
              matrix(0L,njk,2L*njk),
              -rep(1L, njk))
  bd <- rep(0L, njk)
  nf <- c(rep(0L, 3L*njk), 1L)
  na <- rbind(a, ad)
  nb <- c(b, bd)
  if (solver == "lp_solve"){
    nsol <- suppressWarnings(lpSolve::lp('min', nf, na,
                                       c(rep('=', length(b)), rep('<=', length(bd))), nb))
  } else {
    nsol <- Rsymphony::Rsymphony_solve_LP(obj = nf,
                                       mat = na,
                                       dir = c(rep('==', length(b)), rep('<=', length(bd))),
                                       rhs = nb)
  }
  # nz <- nsol$solution
  # Output programa lineal 2
  sol2 <- matrix(nsol$solution[1L:njk], nj, nk, TRUE,
                 dimnames=dimnames(lphom.object$VTM.complete))
  return(sol2)
}

#' @importFrom Rsymphony Rsymphony_solve_LP
#' @importFrom lpSolve lp
# La función **lphom_local()** calcula la matriz de transferencia que corresponde a una unidad de
# votación. Para ello toma como base la matriz de transferencia global obtenida por **lphom()** y
# tratando de mover esta solución lo menos posible busca otra que verique las restricciones asociadas
# a la unidad de votación. Es decir, que aplicada sobre los resultados de origen de la unidad genere
# exactamente los resultados de destino de la unidad. Esta solucion puede ser indeterminada.
lphom_local <- function(lphom.object, iii, solver){
  xt <- lphom.object$origin[iii, ]
  yt <- lphom.object$destination[iii,]
  pg <- lphom.object$VTM.complete/rowSums(lphom.object$VTM.complete)
  ceros <- determinar_zeros_estructurales(lphom.object)
  # Parameters
  nj <- length(xt)
  nk <- length(yt)
  njk <- nj * nk

  # restricciones sum(pjk)=1
  a1 <- kronecker(diag(nj), t(rep(1L, nk)))
  b1 <- rep(1L, nj)

  # restricciones cuadre total de votos para los nk partidos
  at <- t(kronecker(xt, diag(nk)))
  bt <- yt

  # restricciones para definir los ejk
  ajk <- cbind(kronecker(diag(xt), diag(nk)), t(kronecker(diag(njk), c(1L,-1L))))
  bjk <- as.vector(t(xt * pg))

  # Sintesis
  a <- rbind(cbind(rbind(a1, at), matrix(0L, nj+nk, 2L*njk)), ajk)
  b <- c(b1, bt, bjk)

  # Restricciones de proporciones constantes para salidas
  # Escenarios
  escenario <- lphom.object$inputs$new_and_exit_voters[1L]
  J.inic <- ncol(lphom.object$inputs$votes_election1)
  K.inic <- ncol(lphom.object$inputs$votes_election2)
  # Raw scenario. Constraints pjk(j,K) constant,
  if ((escenario == "raw") & (J.inic < nj) & (K.inic < nk)){
    pb <- yt[nk]/sum(xt[1L:(nj - 1L)])
    ab <- matrix(0L, nj,3L*njk)
    bb <- rep(0L, nj)
    for (j in 1L:nj){
      ab[j,j*nk] <- 1L
      bb[j] <- pb*(j < nj)
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }
  if ((escenario == "raw") & (J.inic == nj) & (K.inic < nk)){
    pb <- yt[nk]/sum(xt[1L:nj])
    ab <- matrix(0L, nj,3L*njk)
    bb <- rep(0L, nj)
    for (j in 1L:nj){
      ab[j,j*nk] <- 1L
      bb[j] <- pb
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }

  # Regular scenario. Constraints pjk(j,K) constant.
  if ((escenario == "regular") & (K.inic < nk) & (J.inic < nj)){
    pb <- yt[nk]/sum(xt[1L:(nj-2L)])
    ab <- matrix(0L, nj, 3L*njk)
    bb <- rep(0L, nj)
    for (j in 1L:nj){
      ab[j,j*nk] <- 1
      bb[j] <- pb*(j < (nj-1L))
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }
  if ((escenario == "regular") & (J.inic == nj) & (K.inic < nk)){
    pb <- yt[nk]/sum(xt[1L:(nj-1L)])
    ab <- matrix(0L, nj,3L*njk)
    bb <- rep(0L, nj)
    for (j in 1L:nj){
      ab[j,j*nk] = 1L
      bb[j]=pb*(j < nj)
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }

  # Full scenario. Constraints pjk(j,K) constant.
  if (escenario == "full"){
    pb <- yt[nk]/sum(xt[1L:(nj - 2L)])
    ab <- matrix(0L, nj, 3L*njk)
    bb <- rep(0L, nj)
    for (j in 1L:nj){
      ab[j,j*nk] <- 1L
      bb[j] <- pb*(j < (nj-1L))
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }

  # Gold scenario. Constraints pjk(j,K) constant.
  if (escenario == "gold"){
    # Column K-1
    pb <- yt[nk-1L]/sum(xt[1L:(nj-2L)])
    ab <- matrix(0L, nj, 3L*njk)
    bb <- rep(0L, nj)
    for (j in 1L:nj){
      ab[j,(nk-1L) + (j-1L)*nk] <- 1L
      bb[j] <- pb*(j < (nj-1L))
    }
    a = rbind(a,ab)
    b = c(b,bb)
    # Column K
    pb <- yt[nk]/sum(xt[1L:(nj-2L)])
    ab <- matrix(0L, nj, 3L*njk)
    bb <- rep(0L, nj)
    for (j in 1L:nj){
      ab[j, j*nk] <- 1L
      bb[j] <- pb*(j < (nj-1L))
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }

  # Restricciones de ceros estructurales introducidos por el usuario
  if (length(ceros) > 0){
    ast = matrix(0L, length(ceros), ncol(a))
    bst = rep(0L, length(ceros))
    for (i in 1L:length(ceros)){
      ast[i, nk*(ceros[[i]][1L]-1L) + ceros[[i]][2L]] <- 1L
    }
    a <- rbind(a, ast)
    b <- c(b, bst)
  }

  # calculo bounds
#  lb <- matrix(0, 3*njk ,1)
#  ub <- matrix(Inf, 3*njk, 1)
  # Objective function, to minimize
  fun.obj <- c(rep(0L, njk), rep(1L, 2L*njk))
  # Solution
  if (solver == "lp_solve"){
    sol <- suppressWarnings(lpSolve::lp('min', fun.obj, a, rep('=', length(b)), b))
  #z <- sol$solution
  } else {
    sol <- Rsymphony::Rsymphony_solve_LP(obj = fun.obj,
                                       mat = a,
                                       dir = rep('==', length(b)),
                                       rhs = b)
  }
  # Output programa lineal 1
  sol1 <- matrix(sol$solution[1:njk], nj, nk, TRUE,
                 dimnames=dimnames(lphom.object$VTM.complete))
  return(sol1)
}



#' @importFrom Rsymphony Rsymphony_solve_LP
#' @importFrom lpSolve lp
# La función **lphom_local_abs()** calcula la matriz de transferencia que corresponde a una unidad de
# votación. Para ello toma como base la matriz de transferencia global obtenida por **lphom()** y
# tratando de mover esta solución lo menos posible busca otra que verique las restricciones asociadas
# a la unidad de votación. Es decir, que aplicada sobre los resultados de origen de la unidad genere
# exactamente los resultados de destino de la unidad. Esta solucion, que puede ser indeterminada, es refinada
# en una segunda etapa eligiendo entre todas las matrices que cumplen las restricciones,
# la matriz con menor distancia L1 con la matriz original.
lphom_local_abs <- function(lphom.object, iii, solver){
  xt <- lphom.object$origin[iii, ]
  yt <- lphom.object$destination[iii,]
  pg <- lphom.object$VTM.complete/rowSums(lphom.object$VTM.complete)
  ceros <- determinar_zeros_estructurales(lphom.object)

  # Parameters
  nj <- length(xt)
  nk <- length(yt)
  njk <- nj * nk

  # restricciones sum(pjk)=1
  a1 <- kronecker(diag(nj), t(rep(1L, nk)))
  b1 <- rep(1L, nj)

  # restricciones cuadre total de votos para los nk partidos
  at <- t(kronecker(xt, diag(nk)))
  bt <- yt

  # restricciones para definir los ejk
  ajk <- cbind(kronecker(diag(xt), diag(nk)), t(kronecker(diag(njk), c(1L,-1L))))
  bjk <- as.vector(t(xt * pg))

  # Sintesis
  a <- rbind(cbind(rbind(a1, at), matrix(0L, nj+nk, 2L*njk)), ajk)
  b <- c(b1, bt, bjk)

  # Restricciones de proporciones constantes para salidas
  # Escenarios
  escenario <- lphom.object$inputs$new_and_exit_voters[1L]
  J.inic <- ncol(lphom.object$inputs$votes_election1)
  K.inic <- ncol(lphom.object$inputs$votes_election2)
  # Raw scenario. Constraints pjk(j,K) constant,
  if ((escenario == "raw") & (J.inic < nj) & (K.inic < nk)){
    pb <- yt[nk]/sum(xt[1L:(nj-1L)])
    ab <- matrix(0L, nj, 3L*njk)
    bb <- rep(0L, nj)
    for (j in 1L:nj){
      ab[j,j*nk] <- 1L
      bb[j] <- pb*(j < nj)
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }
  if ((escenario == "raw") & (J.inic == nj) & (K.inic < nk)){
    pb <- yt[nk]/sum(xt[1L:nj])
    ab <- matrix(0L, nj, 3L*njk)
    bb <- rep(0L, nj)
    for (j in 1L:(nj)){
      ab[j,j*nk] <- 1L
      bb[j] <- pb
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }

  # Regular scenario. Constraints pjk(j,K) constant.
  if ((escenario == "regular") & (K.inic < nk) & (J.inic < nj)){
    pb <- yt[nk]/sum(xt[1L:(nj-2L)])
    ab <- matrix(0L, nj, 3L*njk)
    bb <- rep(0L, nj)
    for (j in 1L:nj){
      ab[j,j*nk] <- 1L
      bb[j] <- pb*(j < (nj-1L))
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }
  if ((escenario == "regular") & (J.inic == nj) & (K.inic < nk)){
    pb <- yt[nk]/sum(xt[1L:(nj-1L)])
    ab <- matrix(0L, nj, 3L*njk)
    bb <- rep(0L, nj)
    for (j in 1L:nj){
      ab[j,j*nk] = 1L
      bb[j]=pb*(j < nj)
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }

  # Full scenario. Constraints pjk(j,K) constant.
  if (escenario == "full"){
    pb <- yt[nk]/sum(xt[1L:(nj-2L)])
    ab <- matrix(0, nj, 3L*njk)
    bb <- rep(0L, nj)
    for (j in 1L:nj){
      ab[j,j*nk] <- 1L
      bb[j] <- pb*(j < (nj-1L))
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }

  # Gold scenario. Constraints pjk(j,K) constant.
  if (escenario == "gold"){
    # Column K-1
    pb <- yt[nk-1L]/sum(xt[1L:(nj-2L)])
    ab <- matrix(0L, nj, 3L*njk)
    bb <- rep(0L, nj)
    for (j in 1L:nj){
      ab[j,(nk-1L) + (j-1L)*nk] <- 1L
      bb[j] <- pb*(j < (nj-1L))
    }
    a = rbind(a,ab)
    b = c(b,bb)
    # Column K
    pb <- yt[nk]/sum(xt[1L:(nj-2L)])
    ab <- matrix(0L, nj, 3L*njk)
    bb <- rep(0L, nj)
    for (j in 1L:nj){
      ab[j, j*nk] <- 1L
      bb[j] <- pb*(j < (nj-1L))
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }

  # Restricciones de ceros estructurales introducidos por el usuario
  if (length(ceros) > 0){
    ast = matrix(0L, length(ceros), ncol(a))
    bst = rep(0L, length(ceros))
    for (i in 1L:length(ceros)){
      ast[i, nk*(ceros[[i]][1L]-1L) + ceros[[i]][2L]] <- 1L
    }
    a <- rbind(a, ast)
    b <- c(b, bst)
  }

  # calculo bounds
#  lb <- matrix(0, 3*njk ,1)
#  ub <- matrix(Inf, 3*njk, 1)
  # Objective function, to minimize
  fun.obj <- c(rep(0L, njk), rep(1L, 2L*njk))
  # Solution
  if (solver == "lp_solve"){
    sol <- suppressWarnings(lpSolve::lp('min', fun.obj, a, rep('=', length(b)), b))
  } else {
    sol <- Rsymphony::Rsymphony_solve_LP(obj = fun.obj,
                                       mat = a,
                                       dir = rep('==', length(b)),
                                       rhs = b)
  }
  # z <- sol$solution
  # Output programa lineal 1
  sol1 <- matrix(sol$solution[1:njk], nj, nk, TRUE,
                 dimnames=dimnames(lphom.object$VTM.complete))

  # Segundo programa lineal
  a <- rbind(a, fun.obj)
  b <- c(b, sol$objval)
  aa <- cbind(diag(njk),
              matrix(0L, njk, 2L*njk),
              a[(nj+nk+1L):(nj+nk+njk),(njk+1L):(3L*njk)])
  na <- cbind(a, matrix(0L, nrow(a), 2L*njk))
  na <- rbind(na, aa)
  pgp <- t(pg)
  nb <- c(b, as.vector(pgp))
  nf <- c(rep(0L, 3L*njk), rep(1L, 2L*njk))
  if (solver == "lp_solve"){
    nsol <- suppressWarnings(lpSolve::lp('min', nf, na, rep('=', length(nb)), nb))
  } else {
    nsol <- Rsymphony::Rsymphony_solve_LP(obj = nf,
                                       mat = na,
                                       dir = rep('==', length(nb)),
                                       rhs = nb)
  }
  #nz <- nsol$solution
  # Output programa lineal 2
  sol2 <- matrix(nsol$solution[1:njk], nj, nk, TRUE,
                 dimnames=dimnames(lphom.object$VTM.complete))
  return(sol2)
}


#' @importFrom Rsymphony Rsymphony_solve_LP
#' @importFrom lpSolve lp
# La función **lphom_local_max()** calcula la matriz de transferencia que corresponde a una unidad de
# votación. Para ello toma como base la matriz de transferencia global obtenida por **lphom()** y
# tratando de mover esta solución lo menos posible busca otra que verique las restricciones asociadas
# a la unidad de votación. Es decir, que aplicada sobre los resultados de origen de la unidad genere
# exactamente los resultados de destino de la unidad. Esta solucion, que puede ser indeterminada, es refinada
# en una segunda etapa eligiendo entre todas las matrices que cumplen las restricciones,
# la matriz con menor distancia L_Inf con la matriz original.
lphom_local_max <- function(lphom.object, iii, solver){
  xt <- lphom.object$origin[iii, ]
  yt <- lphom.object$destination[iii,]
  pg <- lphom.object$VTM.complete/rowSums(lphom.object$VTM.complete)
  ceros <- determinar_zeros_estructurales(lphom.object)

  # Parameters
  nj <- length(xt)
  nk <- length(yt)
  njk <- nj * nk

  # restricciones sum(pjk)=1
  a1 <- kronecker(diag(nj), t(rep(1L, nk)))
  b1 <- rep(1L, nj)

  # restricciones cuadre total de votos para los nk partidos
  at <- t(kronecker(xt, diag(nk)))
  bt <- yt

  # restricciones para definir los ejk
  ajk <- cbind(kronecker(diag(xt), diag(nk)), t(kronecker(diag(njk), c(1L,-1L))))
  bjk <- as.vector(t(xt * pg))

  # Sintesis
  a <- rbind(cbind(rbind(a1, at), matrix(0L, nj+nk, 2L*njk)), ajk)
  b <- c(b1, bt, bjk)

  # Restricciones de proporciones constantes para salidas
  # Escenarios
  escenario <- lphom.object$inputs$new_and_exit_voters[1]
  J.inic <- ncol(lphom.object$inputs$votes_election1)
  K.inic <- ncol(lphom.object$inputs$votes_election2)
  # Raw scenario. Constraints pjk(j,K) constant,
  if ((escenario == "raw") & (J.inic < nj) & (K.inic < nk)){
    pb <- yt[nk]/sum(xt[1L:(nj-1L)])
    ab <- matrix(0L, nj, 3L*njk)
    bb <- rep(0L, nj)
    for (j in 1L:nj){
      ab[j,j*nk] <- 1L
      bb[j] <- pb*(j < nj)
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }
  if ((escenario == "raw") & (J.inic == nj) & (K.inic < nk)){
    pb <- yt[nk]/sum(xt[1L:nj])
    ab <- matrix(0L, nj, 3L*njk)
    bb <- rep(0L, nj)
    for (j in 1L:(nj)){
      ab[j,j*nk] <- 1L
      bb[j] <- pb
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }

  # Regular scenario. Constraints pjk(j,K) constant.
  if ((escenario == "regular") & (K.inic < nk) & (J.inic < nj)){
    pb <- yt[nk]/sum(xt[1L:(nj-2L)])
    ab <- matrix(0L, nj, 3L*njk)
    bb <- rep(0L, nj)
    for (j in 1L:nj){
      ab[j,j*nk] <- 1L
      bb[j] <- pb*(j < (nj-1L))
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }
  if ((escenario == "regular") & (J.inic == nj) & (K.inic < nk)){
    pb <- yt[nk]/sum(xt[1L:(nj-1L)])
    ab <- matrix(0L, nj, 3L*njk)
    bb <- rep(0L, nj)
    for (j in 1L:nj){
      ab[j,j*nk] = 1L
      bb[j]=pb*(j < nj)
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }

  # Full scenario. Constraints pjk(j,K) constant.
  if (escenario == "full"){
    pb <- yt[nk]/sum(xt[1L:(nj-2L)])
    ab <- matrix(0L, nj, 3L*njk)
    bb <- rep(0L, nj)
    for (j in 1L:nj){
      ab[j,j*nk] <- 1L
      bb[j] <- pb*(j < (nj - 1L))
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }

  # Gold scenario. Constraints pjk(j,K) constant.
  if (escenario == "gold"){
    # Column K-1
    pb <- yt[nk-1L]/sum(xt[1L:(nj-2L)])
    ab <- matrix(0L, nj, 3L*njk)
    bb <- rep(0L, nj)
    for (j in 1L:nj){
      ab[j,(nk-1L) + (j-1L)*nk] <- 1L
      bb[j] <- pb*(j < (nj-1L))
    }
    a = rbind(a,ab)
    b = c(b,bb)
    # Column K
    pb <- yt[nk]/sum(xt[1L:(nj-2L)])
    ab <- matrix(0L, nj, 3L*njk)
    bb <- rep(0L, nj)
    for (j in 1L:nj){
      ab[j, j*nk] <- 1L
      bb[j] <- pb*(j < (nj-1L))
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }

  # Restricciones de ceros estructurales introducidos por el usuario
  if (length(ceros) > 0){
    ast = matrix(0L, length(ceros), ncol(a))
    bst = rep(0L, length(ceros))
    for (i in 1L:length(ceros)){
      ast[i, nk*(ceros[[i]][1L]-1L) + ceros[[i]][2L]] <- 1L
    }
    a <- rbind(a, ast)
    b <- c(b, bst)
  }

  # calculo bounds
#  lb <- matrix(0, 3*njk ,1)
#  ub <- matrix(Inf, 3*njk, 1)
  # Objective function, to minimize
  fun.obj <- c(rep(0L, njk), rep(1L, 2L*njk))
  # Solution
  if (solver == "lp_solve"){
    sol <- suppressWarnings(lpSolve::lp('min', fun.obj, a, rep('=', length(b)), b))
  } else {
    sol <- Rsymphony::Rsymphony_solve_LP(obj = fun.obj,
                                       mat = a,
                                       dir = rep('==', length(b)),
                                       rhs = b)
  }
  # Output programa lineal 1
  sol1 <- matrix(sol$solution[1L:njk], nj, nk, TRUE,
                 dimnames=dimnames(lphom.object$VTM.complete))

  # Segundo programa lineal
  a <- rbind(a, fun.obj)
  b <- c(b, sol$objval)
  a <- cbind(a, rep(0L, nrow(a)))
  ad <- cbind(diag(njk),
              matrix(0L, njk, 2L*njk),
              -rep(1L, njk))
  bd <- rep(0L, njk)
  nf <- c(rep(0L, 3L*njk), 1L)
  na <- rbind(a, ad)
  nb <- c(b, bd)
  if (solver == "lp_solve"){
    nsol <- suppressWarnings(lpSolve::lp('min', nf, na,
                                       c(rep('=', length(b)), rep('<=', length(bd))), nb))
  } else {
    nsol <- Rsymphony::Rsymphony_solve_LP(obj = nf,
                                       mat = na,
                                       dir = c(rep('==', length(b)), rep('<=', length(bd))),
                                       rhs = nb)
  }
  # Output programa lineal 2
  sol2 <- matrix(nsol$solution[1L:njk], nj, nk, TRUE,
                 dimnames=dimnames(lphom.object$VTM.complete))
  return(sol2)
}


# La función **lp_solver_local()** determmina la función a emplear como función de ajuste local
# dadas las especificaciones introducidas por el usuario
lp_solver_local <- function(uniform, distance.local){
  if (uniform) {
    if (distance.local[1L] == "abs") {
      lphom_unit <- lphom_local_abs
    } else if (distance.local[1L] == "max") {
      lphom_unit <- lphom_local_max
    } else if (distance.local[1L] == "none") {
      lphom_unit <- lphom_local
    }
  } else {
    if (distance.local[1L] == "abs") {
      lphom_unit <- calculo_MT_unidad_abs
    } else if (distance.local[1L] == "max") {
      lphom_unit <- calculo_MT_unidad_max
    } else if (distance.local[1L] == "none") {
      lphom_unit <- calculo_MT_unidad
    }
  }
  return(lphom_unit)
}


# La funcion dec2counts encuentra la matriz con entradas enteras más próxima de una
# matriz de transferencia con entradas decimales.
dec2counts<-function(matriz, vector.fila, vector.columna){
  #
  # aprox.entera(matriz,vector.columna,vector.fila)
  #
  # Funcion para buscar la solucion entera mas proxima de una matriz con entradas
  # no enteras debe cumplir restricciones de agregacion enteras por filas y columnas.
  #
  # INPUT:
  #       matriz: matriz mxn inicial (no necesariamente cuadrada) con entradas no
  #               necesariamente enteras que se busca mover lo menos posible
  #               para una matriz entera verificando las agregaciones por filas
  #               y columnas dadas por vector.fila y vector.columna.
  #       vector.fila: vector de m componentes con lo que deben sumar las filas de
  #                    matriz despues de la aproximacion entera.
  #       vector.columna: vector de n componentes con lo que deben sumar las columnas
  #                       de matriz despues de la aproximacion entera.
  # OUTPUT:
  #       Una matriz origen-destino de numeros enteros cumpliendo las restricciones
  #       de suma de filas y columnas.
  #

  # funcion objetivo
  objetivo<-rep(c(1L, 1L, 0L),length(vector.fila)*length(vector.columna))

  # Restricciones de que los coeficientes de la matriz mas el valor positivo
  # menos el valor negativo debe ser igual al entero mas proximo
  R1 <- t(kronecker(diag(length(vector.fila)*length(vector.columna)),
                    c(1L, -1L, 1L)))
  c1 <- as.vector(t(matriz))
  # Restricciones de suma de filas
  R2 <- t(kronecker(diag(length(vector.fila)),rep(c(0L,0L,1L),length(vector.columna))))
  c2 <- vector.fila
  # Restricciones de suma de columnas
  R3<- kronecker(t(rep(1L,length(vector.fila))),
                 kronecker(diag(length(vector.columna)),t(c(0L,0L,1L))))
  c3 <- vector.columna
  # Conjunto de todas las restricciones
  R <- rbind(R1,R2,R3)
  c0 <- c(c1,c2,c3)
  # Tipo de restricciones
  direc <- rep("==",length(c0))
  # Indices de las variables que han de ser enteras
  #    indices <- which(rep(c(0,0,1),length(vector.fila)*length(vector.columna))==1)
  tipos <- rep(c("C","C","I"),length(vector.fila)*length(vector.columna))
  indices <- which(tipos=="I")
  # Matriz de transferencia con valores enteros
  output <-matrix(Rsymphony::Rsymphony_solve_LP(obj = objetivo,
                                                mat = R,
                                                dir = direc,
                                                rhs= c0,
                                                types = tipos)$solution[indices],
                  length(vector.fila),
                  length(vector.columna), TRUE)
  return(output)
}


# La funcion model_LP calcula todos los elementos (matrices y vectores) que definen
# el programa lineal basico asociado a lphom, bajo supuesto de elecciones simultaneas.
model_LP <- function(votes_election1, votes_election2) {
  x <- as.matrix(votes_election1)
  y <- as.matrix(votes_election2)
  if (nrow(x) != nrow(y)) {
    stop("The number of spatial units is different in origin and destination")
  }
  if (min(x, y) < 0) {
    stop("Negative values for voters (electors) are not allowed")
  }
  d <- rowSums(y) - rowSums(x)
  if (any(d != 0L)) {
    stop("At least in a unit the sums (by row) of both elections differ")
  }

  # Parameters
  I <- nrow(x)
  J <- ncol(x)
  K <- ncol(y)
  JK <- J * K
  IK <- I * K

  # Constraints sum(pjk)=1
  a1 <- cbind(kronecker(diag(J), t(rep(1L, K))), matrix(0L, J, 2L*IK))
  b1 <- rep(1L, J) # Segundos miembros asociados

  # Constraints total match for K parties
  xt <- colSums(x)
  at <- cbind(t(kronecker(xt, diag(K))), matrix(0L, K, 2L*IK))
  bt <- colSums(y) # Segundos miembros asociados

  # Constraints to match votes in the I spatial units and K parties
  bp <- as.vector(t(y))
  ap <- cbind(kronecker(x, diag(K)), t(kronecker(diag(IK), c(1L, -1L))))

  # Joining the three sets of constraints
  a <- rbind(a1, at, ap)
  b <- c(b1, bt, bp)

  # Bounds
  # lb <- matrix(0, JK + 2 * IK, 1)
  # ub <- matrix(Inf, JK + 2 * IK, 1)

  # Objective function, to minimize
  fun.obj <- c(rep(0L, JK), rep(1L, 2L * IK))

  output <- list("a" = a, "b" = b,
                 #"lb" = lb, "ub" = ub,
                 "fun.obj" = fun.obj)
  return(output)
}


# La funcion model_local_LP calcula todos los elementos (matrices y vectores) que definen
# el programa lineal basico asociado a lphom_local, bajo supuesto de elecciones simultaneas.
model_local_LP <- function(MT, marginal_fila, marginal_columna){

  xt <- as.vector(as.matrix(marginal_fila))
  yt <- as.vector(as.matrix(marginal_columna))
  pg <- MT/rowSums(MT)

  # Parameters
  nj <- length(xt)
  nk <- length(yt)
  njk <- nj * nk

  # restricciones sum(pjk)=1
  a1 <- kronecker(diag(nj), t(rep(1, nk)))
  b1 <- rep(1, nj)

  # restricciones cuadre total de votos para los nk partidos
  at <- t(kronecker(xt, diag(nk)))
  bt <- yt

  # restricciones para definir los ejk
  ajk <- cbind(kronecker(diag(xt), diag(nk)), t(kronecker(diag(njk), c(1L,-1L))))
  bjk <- as.vector(t(xt * pg))

  # Sintesis
  a <- rbind(cbind(rbind(a1, at), matrix(0, nj+nk,2*njk)), ajk)
  b <- c(b1, bt, bjk)

  fun.obj <- c(rep(0L, njk), rep(1L, 2L*njk))

  output <- list("a" = a, "b" = b,
                 #"lb" = lb, "ub" = ub,
                 "fun.obj" = fun.obj)
  return(output)
}


# La funcion MT_joint_local calcula la matriz de transferencia de votos que corresponde
# a una unidad considerando conjuntamente las dos posibles ordenaciones de filas y columnas,
# siendo la solucion congruente y no dependiendo con solucion unica como se proyecte
# a filas o columnas cada clasificacion
MT_joint_local_0 <- function(MT, marginal_fila, marginal_columna, solver){

  xt <- as.vector(as.matrix(marginal_fila))
  yt <- as.vector(as.matrix(marginal_columna))

  # Parameters
  J <- nrow(MT)
  K <- ncol(MT)
  JK <- J*K

  # Etapa 1
  m.local.12 <- model_local_LP(MT, marginal_fila, marginal_columna)
  m.local.21 <- model_local_LP(t(MT), marginal_columna, marginal_fila)

  a <- rbind(cbind(m.local.12$a, matrix(0L, nrow(m.local.12$a), ncol(m.local.21$a))),
             cbind(matrix(0L, nrow(m.local.21$a), ncol(m.local.12$a)), m.local.21$a))
  b <- c(m.local.12$b, m.local.21$b)
  fun.obj <- c(m.local.12$fun.obj, m.local.21$fun.obj)

  # We add the las J*K constraints of congruence
  ap <- kronecker(diag(xt), diag(K))
  aq <- matrix(0L, JK, JK)
  for (k in 1L:K){
    aq <- aq + kronecker(t(diag(K)[,k]),
                         kronecker(-diag(J), yt[k]*diag(K)[,k]))
  }

  a <- rbind(a,
             cbind(ap, matrix(0L, JK, 2L*JK),
                   aq, matrix(0L, JK, 2L*JK)))
  b <- c(b, rep(0L, JK))

  # names1 <- colnames(marginal_fila)
  # names2 <- colnames(marginal_columna)

  # Solution Etapa 1
  if (solver == "lp_solve"){
    sol <- suppressWarnings(lpSolve::lp('min', fun.obj, a, rep('=', length(b)), b))
  } else if (solver == "symphony") {
    sol <- Rsymphony::Rsymphony_solve_LP(obj = fun.obj,
                                         mat = a,
                                         dir = rep('==', length(b)),
                                         rhs = b)
  }
  z <- sol$solution

  # Etapa 2
  a <- rbind(a, fun.obj)
  b <- c(b, sol$objval)
  # matriz pjk
  aa <- cbind(diag(JK),
              matrix(0L, JK, 2L*JK), matrix(0L, JK, 3L*JK),
              a[(J+K+1L):(J+K+JK),(JK+1L):(3L*JK)], matrix(0L,JK,2L*JK))
  na <- cbind(a, matrix(0L, nrow(a),4L*JK))
  na <- rbind(na, aa)
  pgp <- t(MT/rowSums(MT))
  nb <- c(b, as.vector(pgp))
  # matriz qkj
  aa2 <- cbind(matrix(0L, JK, 3L*JK), diag(JK),
               matrix(0L, JK, 2L*JK), matrix(0L, JK, 2L*JK),
               a[(2L*(J+K)+JK+1L):(2L*(J+K)+2L*JK),(JK+3L*JK+1L):(6L*JK)])
  na2 <- rbind(na, aa2)
  pgp2 <- t(t(MT)/colSums(MT))
  nb2 <- c(nb, as.vector(pgp2))
  nf <- c(rep(0L, 6L*JK), rep(1L, 4L*JK))
  if (solver == "lp_solve"){
    nsol <- suppressWarnings(lpSolve::lp('min', nf, na2, rep('=', length(nb2)), nb2))
  } else {
    nsol <- Rsymphony::Rsymphony_solve_LP(obj = nf,
                                          mat = na2,
                                          dir = rep('==', length(nb2)),
                                          rhs = nb2)
  }
  # Output programa lineal 2
  VTM.votos <- matrix(nsol$solution[1L:JK], J, K,
                      TRUE, dimnames = dimnames(MT)) * xt

  return(VTM.votos)
}


# La funcion MT_joint_local garantiza la congruencia de la solucion en las unidades
# donde existe mas de una solucion
MT_joint_local <- function(MT, marginal_fila, marginal_columna, solver){
  M12 <- MT_joint_local_0(round(MT, 4L),
                          marginal_fila,
                          marginal_columna,
                          solver)
  M21 <- MT_joint_local_0(round(t(MT), 4L),
                          marginal_columna,
                          marginal_fila,
                          solver)
  output <- (M12+t(M21))/2
  return(output)
}



# La funcion HET_joint calcula el indice de heterogeniedad
# asociado a la ecuacion (33) de Symmetry estimating RxC tables by ecological inference
HET_joint <- function(array.votos){
  mt.votos <- apply(array.votos, c(1,2), sum)
  mt.prop1 <- mt.votos/rowSums(mt.votos)
  mt.prop2 <- t(mt.votos)/rowSums(t(mt.votos))
  filas <- t(apply(array.votos, c(1,3), sum))
  columnas <- t(apply(array.votos, c(2,3), sum))
  array.homogeneo1 <- array(NA, dim(array.votos))
  array.homogeneo2 <- array(NA, dim(array.votos))

  for (ii in 1L:dim(array.homogeneo1)[3L]){
    for (jj in 1L:dim(array.homogeneo1)[1L]){
      array.homogeneo1[jj, , ii] <- filas[ii, jj] * mt.prop1[jj, ]
    }
    for (jj in 1L:dim(array.homogeneo2)[2L]){
      array.homogeneo2[, jj, ii] <- columnas[ii, jj] * mt.prop2[jj, ]
    }
  }

  HETe1 <- 100*(0.5*sum(abs(array.homogeneo1 - array.votos)))/sum(array.votos)
  HETe2 <- 100*(0.5*sum(abs(array.homogeneo2 - array.votos)))/sum(array.votos)
  HETe <- (HETe1 + HETe2)/2

  output <- list("HETe" = HETe, "HETe1" = HETe1, "HETe2" = HETe2)
  return(output)
}

