# @importFrom stats runif
# La función **simular_vector_convexo()** simula un vector de proporciones de suma uno
# (vector convexo) a partir de perturbar aleatoriamente un vector base de proporciones de
# suma uno. La perturbación, basada en una distribución uniforme, está modulada por el
# parámetro `d`.
simular_vector_convexo <- function(vector.convexo, d){
  output <- vector.convexo + stats::runif(length(vector.convexo), -d, d)
  output[output > 1] <- 1
  output[output < 0] <- 0
  if (sum(output) == 0) output <- rep(1/length(vector.convexo), length(vector.convexo))
  output <- output/sum(output)
  return(output)
}

# La función **simular_matrices_fila_convexas()** simula un array de matrices de vectores
# fila convexos a partir de perturbar aleatoriamente una matriz base fila estandarizada.
# La perturbación, basada en una distribución uniforme, está modulada por el parámetro `d`.
simular_matrices_fila_convexas <- function(matriz, d, n.unidades){
  output <- array(NA, c(dim(matriz), n.unidades))
  for (i in 1:n.unidades)
    output[, , i] <- t(apply(matriz, 1, simular_vector_convexo, d = d))
  return(output)
}

# La función **determinar_zeros_estructurales()** determina las coordenadas de las celdas
# que corresponden a ceros estructuales de una matriz de tranferencia resultado de aplicar
# la función **lphom()**.
determinar_zeros_estructurales <- function(lphom.object){
  output <- lphom.object$inputs$structural_zeros
  escenario <- lphom.object$inputs$new_and_exit_voters[1]
  J.inic <- ncol(lphom.object$inputs$votes_election1)
  J.final <- ncol(lphom.object$origin)
  K.inic <- ncol(lphom.object$inputs$votes_election2)
  K.final <- ncol(lphom.object$destination)
  if ((escenario == "raw") & (J.inic < J.final) & (K.inic < K.final)){
    output[[length(output) + 1]] <- c(J.final, K.final)
  }
  if ((escenario == "regular") & (K.inic < K.final) & (J.inic == J.final)){
    output[[length(output) + 1]] <- c(J.final, K.final)
  }
  if ((escenario == "regular") & (K.inic < K.final) & (J.inic < J.final)){
    output[[length(output) + 1]] <- c(J.final - 1, K.final)
    output[[length(output) + 1]] <- c(J.final, K.final)
  }
  if (escenario == "full"){
    output[[length(output) + 1]] <- c(J.inic - 1, K.inic)
    output[[length(output) + 1]] <- c(J.inic, K.inic)
  }
  if (escenario == "gold"){
    output[[length(output) + 1]] <- c(J.inic - 1, K.inic)
    output[[length(output) + 1]] <- c(J.inic - 1, K.inic - 1)
    output[[length(output) + 1]] <- c(J.inic, K.inic)
    output[[length(output) + 1]] <- c(J.inic, K.inic - 1)
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
    for (ii in 1:length(ceros)){
      output[ceros[[ii]][1], ceros[[ii]][2], ] <- 0
    }
    for (ii in 1:nrow(lphom.object$origin)){
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
  for (ii in 1:dim(output)[3]){
    for (jj in 1:dim(output)[1]){
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
  for (ii in 1:dim(array.homogeneo)[3]){
    for (jj in 1:dim(array.homogeneo)[1]){
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
  resumenes$Y[,1] <- resumenes$Y[,1] + (rowSums(lphom.object$origin)-rowSums(resumenes$Y))
  negativos <- resumenes$Y[, 1] < 0
  origen <- lphom.object$origin
  origen[negativos, 1] <- origen[negativos, 1] - resumenes$Y[negativos, 1]
  resumenes$Y[negativos, 1] <- 0
  # Estimacion
  estimacion <- suppressMessages(lphom(origen, resumenes$Y, "raw", ceros, FALSE))
  estimacion$VTM.complete <- estimacion$VTM.complete[1:nrow(resumenes$MT.prop),
                                                     1:ncol(resumenes$MT.prop)]
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
  for (ii in 1:length(M.d)){
    for(bb in 1:B){
      escenario <- simular_y_resumir(lphom.object, M.d[ii])
      estadisticos <- rbind( estadisticos, escenario$resumen )
      MT.reales[, , (ii-1)*B + bb] <- escenario$MT.real
      MT.estimadas[, , (ii-1)*B + bb] <- escenario$MT.estimada
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
  iter <- 1
  while(lphom.object$HETe > max(estadisticos[,1]) & (iter < 5)){
    M.inic <- 2*M.inic
    estadisticos <- simular_escenarios(lphom.object = lphom.object,
                                       M.d = M.inic, B = 3)$estadisticos[, c(2,4)]
    iter <- iter + 1
  }
  ajuste <- stats::lm(estadisticos[,1] ~ estadisticos[, 2])
  d0 <- (lphom.object$HETe - ajuste$coef[1])/ajuste$coef[2]
  id0 <- 1/ajuste$coef[2]
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
  HETe <- sim_esc.object$estadisticos[,2]
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
  he <- matriz.errores[, 1]
  he2 <- he^2
  error <- matriz.errores[, 4]
  ajuste.lineal <- stats::lm(error ~ he + he2 + 0)
  eom <- ajuste.lineal$coef[1]*HETe0 + ajuste.lineal$coef[2]*HETe0^2
  ajuste.varianza <- stats::lm(ajuste.lineal$residuals^2 ~ he + he2 + 0)
  seo <- sqrt(max(ajuste.varianza$coef[1]*HETe0 + ajuste.varianza$coef[2]*HETe0^2, 0))
  output <- c(eom, seo)
  return(output)
}

#' @importFrom lpSolve lp
# La función **calculo_MT_unidad()** calcula la matriz de transferencia que corresponde a una
# unidad de votación. Para ello toma como base una matriz de transferencia global y tratando de
# mover esta solución lo menos posible busca otra que verique las restricciones asociadas a la unidad
# de votación. Es decir, que aplicada sobre los resultados de origen de la unidad genere exactamente
# los resultados de destino de la unidad. Esta solucion puede ser indeterminada.
calculo_MT_unidad <- function(lphom.object, iii){
  xt <- lphom.object$origin[iii, ]
  yt <- lphom.object$destination[iii,]
  pg <- lphom.object$VTM.complete/rowSums(lphom.object$VTM.complete)
  ceros <- determinar_zeros_estructurales(lphom.object)
  # Parameters
  nj <- length(xt)
  nk <- length(yt)
  njk <- nj * nk
  # restricciones sum(pjk)=1
  a1 <- matrix(0, nj, njk)
  for (j in 1:nj){
    a1[j,((j-1)*nk+1):(j*nk)] <- 1
  }
  b1 <- rep(1, nj)
  # restricciones cuadre total de votos para los nk partidos
  at <- NULL
  for (j in 1:nj){
    at <- cbind(at, xt[j]*diag(nk))
  }
  bt <- yt
  a <- rbind(a1, at)
  b <- c(b1, bt)
  a <- cbind(a, matrix(0, nj+nk,2*njk))
  # restricciones para definir los ejk
  ajk <- matrix(0, njk,njk)
  ejk <- matrix(0, njk,2*njk)
  bjk <- rep(0, njk)
  for (j in 1:nj){
    for (k in 1:nk){
      fi <- (j-1)*nk + k
      ajk[fi,fi] <- xt[j]
      co <- 2*(fi-1)+ c(1, 2)
      ejk[fi,co] <- c(1, -1)
      bjk[fi] <- xt[j] * pg[j,k]
    }
  }
  ajk <- cbind(ajk, ejk)
  a <- rbind(a, ajk)
  b <- c(b, bjk)
  # restricciones ceros estructurales
  if (length(ceros) > 0){
    ast = matrix(0, length(ceros), ncol(a))
    bst = rep(0, length(ceros))
    for (i in 1:length(ceros)){
      ast[i, nk*(ceros[[i]][1]-1) + ceros[[i]][2]] <- 1
    }
    a <- rbind(a, ast)
    b <- c(b, bst)
  }
  # calculo bounds
  lb <- matrix(0, 3*njk ,1)
  ub <- matrix(Inf, 3*njk, 1)
  # Objective function, to minimize
  fun.obj <- c(rep(0,njk), rep(1,2*njk))
  # Solution
  sol <- suppressWarnings(lpSolve::lp('min', fun.obj, a, rep('=', length(b)), b))
  z <- sol$solution
  # Output programa lineal 1
  sol1 <- matrix(z[1:njk], nj, nk, TRUE, dimnames=dimnames(lphom.object$VTM.complete))
  return(sol1)
}


#' @importFrom lpSolve lp
# La función **calculo_MT_unidad_abs()** calcula la matriz de transferencia que corresponde a una
# unidad de votación. Para ello toma como base una matriz de transferencia global y tratando de
# mover esta solución lo menos posible busca otra que verique las restricciones asociadas a la unidad
# de votación. Es decir, que aplicada sobre los resultados de origen de la unidad genere exactamente
# los resultados de destino de la unidad. Esta solucion, que puede ser indeterminada, es refinada en una
# segunda etapa eligiendo entre todas las matrices que cumplen las restricciones,
# la matriz con menor distancia L1 con la matriz original.
calculo_MT_unidad_abs <- function(lphom.object, iii){
  xt <- lphom.object$origin[iii, ]
  yt <- lphom.object$destination[iii,]
  pg <- lphom.object$VTM.complete/rowSums(lphom.object$VTM.complete)
  ceros <- determinar_zeros_estructurales(lphom.object)
  # Parameters
  nj <- length(xt)
  nk <- length(yt)
  njk <- nj * nk
  # restricciones sum(pjk)=1
  a1 <- matrix(0, nj, njk)
  for (j in 1:nj){
    a1[j,((j-1)*nk+1):(j*nk)] <- 1
  }
  b1 <- rep(1, nj)
  # restricciones cuadre total de votos para los nk partidos
  at <- NULL
  for (j in 1:nj){
    at <- cbind(at, xt[j]*diag(nk))
  }
  bt <- yt
  a <- rbind(a1, at)
  b <- c(b1, bt)
  a <- cbind(a, matrix(0, nj+nk,2*njk))
  # restricciones para definir los ejk
  ajk <- matrix(0, njk,njk)
  ejk <- matrix(0, njk,2*njk)
  bjk <- rep(0, njk)
  for (j in 1:nj){
    for (k in 1:nk){
      fi <- (j-1)*nk + k
      ajk[fi,fi] <- xt[j]
      co <- 2*(fi-1)+ c(1, 2)
      ejk[fi,co] <- c(1, -1)
      bjk[fi] <- xt[j] * pg[j,k]
    }
  }
  ajk <- cbind(ajk, ejk)
  a <- rbind(a, ajk)
  b <- c(b, bjk)
  # restricciones ceros estructurales
  if (length(ceros) > 0){
    ast = matrix(0, length(ceros), ncol(a))
    bst = rep(0, length(ceros))
    for (i in 1:length(ceros)){
      ast[i, nk*(ceros[[i]][1]-1) + ceros[[i]][2]] <- 1
    }
    a <- rbind(a, ast)
    b <- c(b, bst)
  }
  # calculo bounds
  lb <- matrix(0, 3*njk ,1)
  ub <- matrix(Inf, 3*njk, 1)
  # Objective function, to minimize
  fun.obj  <- c(rep(0,njk), rep(1,2*njk))
  # Solution
  sol <- suppressWarnings(lpSolve::lp('min', fun.obj, a, rep('=', length(b)), b))
  z <- sol$solution
  # Output programa lineal 1
  sol1 <- matrix(z[1:njk], nj, nk, TRUE, dimnames=dimnames(lphom.object$VTM.complete))
  # Segundo programa lineal
  a <- rbind(a, fun.obj)
  b <- c(b, sol$objval)
  aa <- cbind(diag(njk),
               matrix(0,njk,2*njk),
               a[(nj+nk+1):(nj+nk+njk),(njk+1):(3*njk)])
  na <- cbind(a, matrix(0, 1+nj+nk+njk,2*njk))
  na <- rbind(na, aa)
  pgp <- t(pg)
  nb <- c(b, as.vector(pgp))
  nf <- c(rep(0, 3*njk), rep(1,2*njk))
  nsol <- suppressWarnings(lpSolve::lp('min', nf, na, rep('=', length(nb)), nb))
  nz <- nsol$solution
  # Output programa lineal 2
  sol2 <- matrix(nz[1:njk], nj, nk, TRUE, dimnames=dimnames(lphom.object$VTM.complete))
  return(sol2)
}

#' @importFrom lpSolve lp
# La función **calculo_MT_unidad_max()** calcula la matriz de transferencia que corresponde a una
# unidad de votación. Para ello toma como base una matriz de transferencia global y tratando de
# mover esta solución lo menos posible busca otra que verique las restricciones asociadas a la unidad
# de votación. Es decir, que aplicada sobre los resultados de origen de la unidad genere exactamente
# los resultados de destino de la unidad. Esta solucion, que puede ser indeterminada, es refinada en una
# segunda etapa eligiendo entre todas las matrices que cumplen las restricciones,
# la matriz con menor distancia L_Inf con la matriz original.
calculo_MT_unidad_max <- function(lphom.object, iii){
  xt <- lphom.object$origin[iii, ]
  yt <- lphom.object$destination[iii,]
  pg <- lphom.object$VTM.complete/rowSums(lphom.object$VTM.complete)
  ceros <- determinar_zeros_estructurales(lphom.object)
  # Parameters
  nj <- length(xt)
  nk <- length(yt)
  njk <- nj * nk
  # restricciones sum(pjk)=1
  a1 <- matrix(0, nj, njk)
  for (j in 1:nj){
    a1[j,((j-1)*nk+1):(j*nk)] <- 1
  }
  b1 <- rep(1, nj)
  # restricciones cuadre total de votos para los nk partidos
  at <- NULL
  for (j in 1:nj){
    at <- cbind(at, xt[j]*diag(nk))
  }
  bt <- yt
  a <- rbind(a1, at)
  b <- c(b1, bt)
  a <- cbind(a, matrix(0, nj+nk,2*njk))
  # restricciones para definir los ejk
  ajk <- matrix(0, njk,njk)
  ejk <- matrix(0, njk,2*njk)
  bjk <- rep(0, njk)
  for (j in 1:nj){
    for (k in 1:nk){
      fi <- (j-1)*nk + k
      ajk[fi,fi] <- xt[j]
      co <- 2*(fi-1)+ c(1, 2)
      ejk[fi,co] <- c(1, -1)
      bjk[fi] <- xt[j] * pg[j,k]
    }
  }
  ajk <- cbind(ajk, ejk)
  a <- rbind(a, ajk)
  b <- c(b, bjk)
  # restricciones ceros estructurales
  if (length(ceros) > 0){
    ast = matrix(0, length(ceros), ncol(a))
    bst = rep(0, length(ceros))
    for (i in 1:length(ceros)){
      ast[i, nk*(ceros[[i]][1]-1) + ceros[[i]][2]] <- 1
    }
    a <- rbind(a, ast)
    b <- c(b, bst)
  }
  # calculo bounds
  lb <- matrix(0, 3*njk ,1)
  ub <- matrix(Inf, 3*njk, 1)
  # Objective function, to minimize
  fun.obj <- c(rep(0,njk), rep(1,2*njk))
  # Solution
  sol <- suppressWarnings(lpSolve::lp('min', fun.obj, a, rep('=', length(b)), b))
  z <- sol$solution
  # Output programa lineal 1
  sol1 <- matrix(z[1:njk], nj, nk, TRUE, dimnames=dimnames(lphom.object$VTM.complete))

  # Segundo programa lineal
  a <- rbind(a, fun.obj)
  b <- c(b, sol$objval)
  a <- cbind(a, rep(0, nrow(a)))
  ad <- cbind(diag(njk),
              matrix(0,njk,2*njk),
              -rep(1, njk))
  bd <- rep(0, njk)
  nf <- c(rep(0, 3*njk), 1)
  na <- rbind(a, ad)
  nb <- c(b, bd)
  nsol <- suppressWarnings(lpSolve::lp('min', nf, na,
                                       c(rep('=', length(b)), rep('<=', length(bd))), nb))
  nz <- nsol$solution
  # Output programa lineal 2
  sol2 <- matrix(nz[1:njk], nj, nk, TRUE, dimnames=dimnames(lphom.object$VTM.complete))
  return(sol2)
}

#' @importFrom lpSolve lp
# La función **lphom_local()** calcula la matriz de transferencia que corresponde a una unidad de
# votación. Para ello toma como base la matriz de transferencia global obtenida por **lphom()** y
# tratando de mover esta solución lo menos posible busca otra que verique las restricciones asociadas
# a la unidad de votación. Es decir, que aplicada sobre los resultados de origen de la unidad genere
# exactamente los resultados de destino de la unidad. Esta solucion puede ser indeterminada.
lphom_local <- function(lphom.object, iii){
  xt <- lphom.object$origin[iii, ]
  yt <- lphom.object$destination[iii,]
  pg <- lphom.object$VTM.complete/rowSums(lphom.object$VTM.complete)
  ceros <- determinar_zeros_estructurales(lphom.object)
  # Parameters
  nj <- length(xt)
  nk <- length(yt)
  njk <- nj * nk
  # restricciones sum(pjk)=1
  a1 <- matrix(0, nj, njk)
  for (j in 1:nj){
    a1[j,((j-1)*nk+1):(j*nk)] <- 1
  }
  b1 <- rep(1, nj)
  # restricciones cuadre total de votos para los nk partidos
  at <- NULL
  for (j in 1:nj){
    at <- cbind(at, xt[j]*diag(nk))
  }
  bt <- yt
  a <- rbind(a1, at)
  b <- c(b1, bt)
  a <- cbind(a, matrix(0, nj+nk,2*njk))
  # restricciones para definir los ejk
  ajk <- matrix(0, njk,njk)
  ejk <- matrix(0, njk,2*njk)
  bjk <- rep(0, njk)
  for (j in 1:nj){
    for (k in 1:nk){
      fi <- (j-1)*nk + k
      ajk[fi,fi] <- xt[j]
      co <- 2*(fi-1)+ c(1, 2)
      ejk[fi,co] <- c(1, -1)
      bjk[fi] <- xt[j] * pg[j,k]
    }
  }
  ajk <- cbind(ajk, ejk)
  a <- rbind(a, ajk)
  b <- c(b, bjk)

  # Restricciones de proporciones constantes para salidas
  # Escenarios
  escenario <- lphom.object$inputs$new_and_exit_voters[1]
  J.inic <- ncol(lphom.object$inputs$votes_election1)
  K.inic <- ncol(lphom.object$inputs$votes_election2)
  # Raw scenario. Constraints pjk(j,K) constant,
  if ((escenario == "raw") & (J.inic < nj) & (K.inic < nk)){
    pb <- yt[nk]/sum(xt[1:(nj-1)])
    ab <- matrix(0, nj,3*njk)
    bb <- rep(0,nj)
    for (j in 1:nj){
      ab[j,j*nk] <- 1
      bb[j] <- pb*(j < nj)
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }
  if ((escenario == "raw") & (J.inic == nj) & (K.inic < nk)){
    pb <- yt[nk]/sum(xt[1:nj])
    ab <- matrix(0, nj,3*njk)
    bb <- rep(0,nj)
    for (j in 1:(nj)){
      ab[j,j*nk] <- 1
      bb[j] <- pb
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }

  # Regular scenario. Constraints pjk(j,K) constant.
  if ((escenario == "regular") & (K.inic < nk) & (J.inic < nj)){
    pb <- yt[nk]/sum(xt[1:(nj-2)])
    ab <- matrix(0, nj,3*njk)
    bb <- rep(0,nj)
    for (j in 1:nj){
      ab[j,j*nk] <- 1
      bb[j] <- pb*(j < (nj-1))
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }
  if ((escenario == "regular") & (J.inic == nj) & (K.inic < nk)){
    pb <- yt[nk]/sum(xt[1:(nj-1)])
    ab <- matrix(0, nj,3*njk)
    bb <- rep(0,nj)
    for (j in 1:nj){
      ab[j,j*nk] = 1
      bb[j]=pb*(j < nj)
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }

  # Full scenario. Constraints pjk(j,K) constant.
  if (escenario == "full"){
    pb <- yt[nk]/sum(xt[1:(nj-2)])
    ab <- matrix(0, nj,3*njk)
    bb <- rep(0,nj)
    for (j in 1:nj){
      ab[j,j*nk] <- 1
      bb[j] <- pb*(j < (nj-1))
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }

  # Gold scenario. Constraints pjk(j,K) constant.
  if (escenario == "gold"){
    # Column K-1
    pb <- yt[nk-1]/sum(xt[1:(nj-2)])
    ab <- matrix(0, nj,3*njk)
    bb <- rep(0,nj)
    for (j in 1:nj){
      ab[j,(nk-1) + (j-1)*nk] <- 1
      bb[j] <- pb*(j < (nj-1))
    }
    a = rbind(a,ab)
    b = c(b,bb)
    # Column K
    pb <- yt[nk]/sum(xt[1:(nj-2)])
    ab <- matrix(0, nj,3*njk)
    bb <- rep(0,nj)
    for (j in 1:nj){
      ab[j, j*nk] <- 1
      bb[j] <- pb*(j < (nj-1))
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }

  # Restricciones de ceros estructurales introducidos por el usuario
  if (length(ceros) > 0){
    ast = matrix(0, length(ceros), ncol(a))
    bst = rep(0, length(ceros))
    for (i in 1:length(ceros)){
      ast[i, nk*(ceros[[i]][1]-1) + ceros[[i]][2]] <- 1
    }
    a <- rbind(a, ast)
    b <- c(b, bst)
  }

  # calculo bounds
  lb <- matrix(0, 3*njk ,1)
  ub <- matrix(Inf, 3*njk, 1)
  # Objective function, to minimize
  fun.obj <- c(rep(0,njk), rep(1,2*njk))
  # Solution
  sol <- suppressWarnings(lpSolve::lp('min', fun.obj, a, rep('=', length(b)), b))
  z <- sol$solution

  # Output programa lineal 1
  sol1 <- matrix(z[1:njk], nj, nk, TRUE, dimnames=dimnames(lphom.object$VTM.complete))
  return(sol1)
}



#' @importFrom lpSolve lp
# La función **lphom_local_abs()** calcula la matriz de transferencia que corresponde a una unidad de
# votación. Para ello toma como base la matriz de transferencia global obtenida por **lphom()** y
# tratando de mover esta solución lo menos posible busca otra que verique las restricciones asociadas
# a la unidad de votación. Es decir, que aplicada sobre los resultados de origen de la unidad genere
# exactamente los resultados de destino de la unidad. Esta solucion, que puede ser indeterminada, es refinada
# en una segunda etapa eligiendo entre todas las matrices que cumplen las restricciones,
# la matriz con menor distancia L1 con la matriz original.
lphom_local_abs <- function(lphom.object, iii){
  xt <- lphom.object$origin[iii, ]
  yt <- lphom.object$destination[iii,]
  pg <- lphom.object$VTM.complete/rowSums(lphom.object$VTM.complete)
  ceros <- determinar_zeros_estructurales(lphom.object)

  # Parameters
  nj <- length(xt)
  nk <- length(yt)
  njk <- nj * nk
  # restricciones sum(pjk)=1
  a1 <- matrix(0, nj, njk)
  for (j in 1:nj){
    a1[j,((j-1)*nk+1):(j*nk)] <- 1
  }
  b1 <- rep(1, nj)
  # restricciones cuadre total de votos para los nk partidos
  at <- NULL
  for (j in 1:nj){
    at <- cbind(at, xt[j]*diag(nk))
  }
  bt <- yt
  a <- rbind(a1, at)
  b <- c(b1, bt)
  a <- cbind(a, matrix(0, nj+nk,2*njk))
  # restricciones para definir los ejk
  ajk <- matrix(0, njk,njk)
  ejk <- matrix(0, njk,2*njk)
  bjk <- rep(0, njk)
  for (j in 1:nj){
    for (k in 1:nk){
      fi <- (j-1)*nk + k
      ajk[fi,fi] <- xt[j]
      co <- 2*(fi-1)+ c(1, 2)
      ejk[fi,co] <- c(1, -1)
      bjk[fi] <- xt[j] * pg[j,k]
    }
  }
  ajk <- cbind(ajk, ejk)
  a <- rbind(a, ajk)
  b <- c(b, bjk)

  # Restricciones de proporciones constantes para salidas
  # Escenarios
  escenario <- lphom.object$inputs$new_and_exit_voters[1]
  J.inic <- ncol(lphom.object$inputs$votes_election1)
  K.inic <- ncol(lphom.object$inputs$votes_election2)
  # Raw scenario. Constraints pjk(j,K) constant,
  if ((escenario == "raw") & (J.inic < nj) & (K.inic < nk)){
    pb <- yt[nk]/sum(xt[1:(nj-1)])
    ab <- matrix(0, nj,3*njk)
    bb <- rep(0,nj)
    for (j in 1:nj){
      ab[j,j*nk] <- 1
      bb[j] <- pb*(j < nj)
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }
  if ((escenario == "raw") & (J.inic == nj) & (K.inic < nk)){
    pb <- yt[nk]/sum(xt[1:nj])
    ab <- matrix(0, nj,3*njk)
    bb <- rep(0,nj)
    for (j in 1:(nj)){
      ab[j,j*nk] <- 1
      bb[j] <- pb
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }

  # Regular scenario. Constraints pjk(j,K) constant.
  if ((escenario == "regular") & (K.inic < nk) & (J.inic < nj)){
    pb <- yt[nk]/sum(xt[1:(nj-2)])
    ab <- matrix(0, nj,3*njk)
    bb <- rep(0,nj)
    for (j in 1:nj){
      ab[j,j*nk] <- 1
      bb[j] <- pb*(j < (nj-1))
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }
  if ((escenario == "regular") & (J.inic == nj) & (K.inic < nk)){
    pb <- yt[nk]/sum(xt[1:(nj-1)])
    ab <- matrix(0, nj,3*njk)
    bb <- rep(0,nj)
    for (j in 1:nj){
      ab[j,j*nk] = 1
      bb[j]=pb*(j < nj)
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }

  # Full scenario. Constraints pjk(j,K) constant.
  if (escenario == "full"){
    pb <- yt[nk]/sum(xt[1:(nj-2)])
    ab <- matrix(0, nj,3*njk)
    bb <- rep(0,nj)
    for (j in 1:nj){
      ab[j,j*nk] <- 1
      bb[j] <- pb*(j < (nj-1))
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }

  # Gold scenario. Constraints pjk(j,K) constant.
  if (escenario == "gold"){
    # Column K-1
    pb <- yt[nk-1]/sum(xt[1:(nj-2)])
    ab <- matrix(0, nj,3*njk)
    bb <- rep(0,nj)
    for (j in 1:nj){
      ab[j,(nk-1) + (j-1)*nk] <- 1
      bb[j] <- pb*(j < (nj-1))
    }
    a = rbind(a,ab)
    b = c(b,bb)
    # Column K
    pb <- yt[nk]/sum(xt[1:(nj-2)])
    ab <- matrix(0, nj,3*njk)
    bb <- rep(0,nj)
    for (j in 1:nj){
      ab[j, j*nk] <- 1
      bb[j] <- pb*(j < (nj-1))
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }

  # Restricciones de ceros estructurales introducidos por el usuario
  if (length(ceros) > 0){
    ast = matrix(0, length(ceros), ncol(a))
    bst = rep(0, length(ceros))
    for (i in 1:length(ceros)){
      ast[i, nk*(ceros[[i]][1]-1) + ceros[[i]][2]] <- 1
    }
    a <- rbind(a, ast)
    b <- c(b, bst)
  }

  # calculo bounds
  lb <- matrix(0, 3*njk ,1)
  ub <- matrix(Inf, 3*njk, 1)
  # Objective function, to minimize
  fun.obj <- c(rep(0,njk), rep(1,2*njk))
  # Solution
  sol <- suppressWarnings(lpSolve::lp('min', fun.obj, a, rep('=', length(b)), b))
  z <- sol$solution

  # Output programa lineal 1
  sol1 <- matrix(z[1:njk], nj, nk, TRUE, dimnames=dimnames(lphom.object$VTM.complete))

  # Segundo programa lineal
  a <- rbind(a, fun.obj)
  b <- c(b, sol$objval)
  aa <- cbind(diag(njk),
              matrix(0,njk,2*njk),
              a[(nj+nk+1):(nj+nk+njk),(njk+1):(3*njk)])
  na <- cbind(a, matrix(0, nrow(a),2*njk))
  na <- rbind(na, aa)
  pgp <- t(pg)
  nb <- c(b, as.vector(pgp))
  nf <- c(rep(0, 3*njk), rep(1,2*njk))
  nsol <- suppressWarnings(lpSolve::lp('min', nf, na, rep('=', length(nb)), nb))
  nz <- nsol$solution
  # Output programa lineal 2
  sol2 <- matrix(nz[1:njk], nj, nk, TRUE, dimnames=dimnames(lphom.object$VTM.complete))
  return(sol2)
}


#' @importFrom lpSolve lp
# La función **lphom_local_max()** calcula la matriz de transferencia que corresponde a una unidad de
# votación. Para ello toma como base la matriz de transferencia global obtenida por **lphom()** y
# tratando de mover esta solución lo menos posible busca otra que verique las restricciones asociadas
# a la unidad de votación. Es decir, que aplicada sobre los resultados de origen de la unidad genere
# exactamente los resultados de destino de la unidad. Esta solucion, que puede ser indeterminada, es refinada
# en una segunda etapa eligiendo entre todas las matrices que cumplen las restricciones,
# la matriz con menor distancia L_Inf con la matriz original.
lphom_local_max <- function(lphom.object, iii){
  xt <- lphom.object$origin[iii, ]
  yt <- lphom.object$destination[iii,]
  pg <- lphom.object$VTM.complete/rowSums(lphom.object$VTM.complete)
  ceros <- determinar_zeros_estructurales(lphom.object)

  # Parameters
  nj <- length(xt)
  nk <- length(yt)
  njk <- nj * nk
  # restricciones sum(pjk)=1
  a1 <- matrix(0, nj, njk)
  for (j in 1:nj){
    a1[j,((j-1)*nk+1):(j*nk)] <- 1
  }
  b1 <- rep(1, nj)
  # restricciones cuadre total de votos para los nk partidos
  at <- NULL
  for (j in 1:nj){
    at <- cbind(at, xt[j]*diag(nk))
  }
  bt <- yt
  a <- rbind(a1, at)
  b <- c(b1, bt)
  a <- cbind(a, matrix(0, nj+nk,2*njk))
  # restricciones para definir los ejk
  ajk <- matrix(0, njk,njk)
  ejk <- matrix(0, njk,2*njk)
  bjk <- rep(0, njk)
  for (j in 1:nj){
    for (k in 1:nk){
      fi <- (j-1)*nk + k
      ajk[fi,fi] <- xt[j]
      co <- 2*(fi-1)+ c(1, 2)
      ejk[fi,co] <- c(1, -1)
      bjk[fi] <- xt[j] * pg[j,k]
    }
  }
  ajk <- cbind(ajk, ejk)
  a <- rbind(a, ajk)
  b <- c(b, bjk)

  # Restricciones de proporciones constantes para salidas
  # Escenarios
  escenario <- lphom.object$inputs$new_and_exit_voters[1]
  J.inic <- ncol(lphom.object$inputs$votes_election1)
  K.inic <- ncol(lphom.object$inputs$votes_election2)
  # Raw scenario. Constraints pjk(j,K) constant,
  if ((escenario == "raw") & (J.inic < nj) & (K.inic < nk)){
    pb <- yt[nk]/sum(xt[1:(nj-1)])
    ab <- matrix(0, nj,3*njk)
    bb <- rep(0,nj)
    for (j in 1:nj){
      ab[j,j*nk] <- 1
      bb[j] <- pb*(j < nj)
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }
  if ((escenario == "raw") & (J.inic == nj) & (K.inic < nk)){
    pb <- yt[nk]/sum(xt[1:nj])
    ab <- matrix(0, nj,3*njk)
    bb <- rep(0,nj)
    for (j in 1:(nj)){
      ab[j,j*nk] <- 1
      bb[j] <- pb
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }

  # Regular scenario. Constraints pjk(j,K) constant.
  if ((escenario == "regular") & (K.inic < nk) & (J.inic < nj)){
    pb <- yt[nk]/sum(xt[1:(nj-2)])
    ab <- matrix(0, nj,3*njk)
    bb <- rep(0,nj)
    for (j in 1:nj){
      ab[j,j*nk] <- 1
      bb[j] <- pb*(j < (nj-1))
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }
  if ((escenario == "regular") & (J.inic == nj) & (K.inic < nk)){
    pb <- yt[nk]/sum(xt[1:(nj-1)])
    ab <- matrix(0, nj,3*njk)
    bb <- rep(0,nj)
    for (j in 1:nj){
      ab[j,j*nk] = 1
      bb[j]=pb*(j < nj)
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }

  # Full scenario. Constraints pjk(j,K) constant.
  if (escenario == "full"){
    pb <- yt[nk]/sum(xt[1:(nj-2)])
    ab <- matrix(0, nj,3*njk)
    bb <- rep(0,nj)
    for (j in 1:nj){
      ab[j,j*nk] <- 1
      bb[j] <- pb*(j < (nj-1))
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }

  # Gold scenario. Constraints pjk(j,K) constant.
  if (escenario == "gold"){
    # Column K-1
    pb <- yt[nk-1]/sum(xt[1:(nj-2)])
    ab <- matrix(0, nj,3*njk)
    bb <- rep(0,nj)
    for (j in 1:nj){
      ab[j,(nk-1) + (j-1)*nk] <- 1
      bb[j] <- pb*(j < (nj-1))
    }
    a = rbind(a,ab)
    b = c(b,bb)
    # Column K
    pb <- yt[nk]/sum(xt[1:(nj-2)])
    ab <- matrix(0, nj,3*njk)
    bb <- rep(0,nj)
    for (j in 1:nj){
      ab[j, j*nk] <- 1
      bb[j] <- pb*(j < (nj-1))
    }
    a = rbind(a,ab)
    b = c(b,bb)
  }

  # Restricciones de ceros estructurales introducidos por el usuario
  if (length(ceros) > 0){
    ast = matrix(0, length(ceros), ncol(a))
    bst = rep(0, length(ceros))
    for (i in 1:length(ceros)){
      ast[i, nk*(ceros[[i]][1]-1) + ceros[[i]][2]] <- 1
    }
    a <- rbind(a, ast)
    b <- c(b, bst)
  }

  # calculo bounds
  lb <- matrix(0, 3*njk ,1)
  ub <- matrix(Inf, 3*njk, 1)
  # Objective function, to minimize
  fun.obj <- c(rep(0,njk), rep(1,2*njk))
  # Solution
  sol <- suppressWarnings(lpSolve::lp('min', fun.obj, a, rep('=', length(b)), b))
  z <- sol$solution

  # Output programa lineal 1
  sol1 <- matrix(z[1:njk], nj, nk, TRUE, dimnames=dimnames(lphom.object$VTM.complete))

  # Segundo programa lineal
  a <- rbind(a, fun.obj)
  b <- c(b, sol$objval)
  a <- cbind(a, rep(0, nrow(a)))
  ad <- cbind(diag(njk),
              matrix(0,njk,2*njk),
              -rep(1, njk))
  bd <- rep(0, njk)
  nf <- c(rep(0, 3*njk), 1)
  na <- rbind(a, ad)
  nb <- c(b, bd)
  nsol <- suppressWarnings(lpSolve::lp('min', nf, na,
                                       c(rep('=', length(b)), rep('<=', length(bd))), nb))
  nz <- nsol$solution
  # Output programa lineal 2
  sol2 <- matrix(nz[1:njk], nj, nk, TRUE, dimnames=dimnames(lphom.object$VTM.complete))
  return(sol2)
}


