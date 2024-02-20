## Auxiliary functions of lp_apriori

# Function to test that all inputs of lp_apriori are correct
test_inputs_lp_apriori <- function(argg){
  
  # Test integers  
  argg.1 <- as.matrix(argg$votes_election1)
  argg.2 <- as.matrix(argg$votes_election2)
  condicion <- max(abs(argg.1 - round(argg.1))) + max(abs(argg.2 - round(argg.2)))
  if(condicion > 0L & argg$integers)
    stop("Integer solutions cannot be computed. At least a value included in 'votes_election1' or 'votes_election2' is decimal.")
  
  # Test names 1
  if(is.character(argg$weights)){
    if (!(argg$weights[1L] %in% c("x", "xy", "constant", "expected", "counts", "sqrt", "sd")))
      stop('Not allowed string for argument "weights". The only allowed strings for "weights" are "constant", "x", "xy", "expected", "counts", "sqrt" and "sd".')
  }
  
  # Test names 2
  if (!(argg$new_and_exit_voters[1L] %in% c("raw", "regular", "ordinary", "semifull", "enriched",
                                           "simultaneous", "full", "fullreverse", "gold", "adjust1", "adjust2")))
    stop('Not allowed string for argument "new_and_exit_voters". The only allowed strings for "new_and_exit_voters" are "raw", "simultaneous", "regular", "ordinary", "enriched", , "adjust1", "adjust2", "semifull", "full", "fullreverse" and "gold".')
  
  # Test names 3
  if (!(argg$solver[1L] %in% c("lp_solve", "symphony")))
    stop('Not allowed string for argument "solver". The only allowed strings for "solver" are "lp_solve" and "symphony".')

  # Test names 3
  if (!(argg$integers.solver[1L] %in% c("lp_solve", "symphony")))
    stop('Not allowed string for argument "integers.solver". The only allowed strings for "integers.solver" are "lp_solve" and "symphony".')
  
  # Tests numeric data
  x <- as.matrix(argg$votes_election1)
  y <- as.matrix(argg$votes_election2)
  if (ncol(x) == 1L) x <- t(x)
  if (ncol(y) == 1L) y <- t(y)
  if (nrow(x) != nrow(y))
    stop('The number of spatial units is different in origin and destination.')
  if (argg$new_and_exit_voters[1L] %in% c("semifull", "simultaneous", "full", "fullreverse", "gold")){
    if (!identical(round(rowSums(x)), round(rowSums(y)))){
      texto <- paste0('The number of voters (electors) in Election 1 and ',
                      'Election 2 differ (in at least a territorial unit). This is not ',
                      'allowed in a \"', argg$new_and_exit_voters[1L], '\" scenario.')
      stop(texto)
    } 
  }
  if (min(x,y) < 0L) stop('Negative values for voters (electors) are not allowed')
  if (min(argg$apriori, na.rm = T) < 0L) 
    stop('Negative values for a priori proportions (counts) are not allowed')
  
  return(list("x" = x, "y" = y))
}


# Function to compute net entries and exits for lp_apriori
lp_apriori_net <- function(x0, y0){
  NET_ENTRIES <- NET_EXITS <- rep(0L, nrow(x0))
  x <- x0
  y <- y0
  
  # Estimation of net entries/exits in the election census
  d <- rowSums(y) - rowSums(x)
  if (any(d != 0L)) {
    NET_ENTRIES[d > 0L] <- d[d > 0L]
    NET_EXITS[d < 0L] <- -d[d < 0L]
  }
  
  # Net entries and exits
  if (sum(NET_ENTRIES) > 0L){
    x <- cbind(x, NET_ENTRIES)
  }
  if (sum(NET_EXITS) > 0L){
    y <- cbind(y0, NET_EXITS)
  }
  
  x0 <- colSums(x0)
  y0 <- colSums(y0)
  x <- colSums(x)
  y <- colSums(y)
  return(list("x0" = x0, "y0" = y0, "x" = x, "y" = y))
}

test_apriori_weights <- function(x0, y0, P0, pesos0, scenario){

  P0 <- as.matrix(P0)
  if(is.character(pesos0)){
    pesos <- pesos0[1]
  } else {
    pesos0 <- as.matrix(pesos0)
    if(!isTRUE(all.equal(dim(P0), dim(pesos0)))) 
      stop("The apriori and weights matrices have different orders.")
  }  
  
  if(!is.null(dim(x0))){
    x0 <- colSums(x0)
    y0 <- colSums(y0)
  }
  
  J1 <- length(x0)
  K2 <- length(y0)
  J0 <- nrow(P0)
  K0 <- ncol(P0)
  
  # Se completa con NAs en caso hay insuficiente apriori info
  if(J0 < J1){
    nf <- J1 - J0
    P0 <- rbind(P0, matrix(NA, nf, ncol(P0)))
  }
  
  if(K0 < K2){
    nf <- K2 - K0
    P0 <- cbind(P0, matrix(NA, nrow(P0), nf))
  }
  
  # Test cuando cuadran
  if(scenario %in% c("simultaneous", "semifull", "full", "fullreverse", "gold")){
    if (K2 < ncol(P0)){
      warning("A priori probabilities have been provided for more destination options than available 
              in 'votes_election2'. The excess of columns in 'apriori' has been omitted.")
      P0 <- P0[, 1L:K2]
    }
    
    if(J1 < nrow(P0)){
      warning("A priori probabilities have been provided for more origin options than available 
              in 'votes_election1'. The excess of columns in of 'apriori' has been omitted.")
      P0 <- P0[1L:J1, ]
    }
  }
  
  # Test cuando pueden no cuadrar
  if (scenario %in% c("raw", "regular", "ordinary", "enriched")){
    if((K2 + 1L) < ncol(P0)){
      warning("A priori probabilities have been provided for more destination options than posible 
              in 'votes_election2'. The excess of columns in 'apriori' has been omitted.")
      P0 <- P0[, 1L:(K2 + 1L)]
    }
    
    if((J1 + 1L) < nrow(P0)){
      warning("A priori probabilities have been provided for more origin options than possible 
              in 'votes_election1'. The excess of columns in of 'apriori' has been omitted.")
      P0 <- P0[1L:(J1 + 1L), ]
    }
  }
  
  if(!is.character(pesos0)){
    pesos0 <- as.matrix(pesos0)
    pesos <- matrix(NA, nrow(P0), ncol(P0))
    nf <- min(nrow(P0), nrow(pesos0))
    nc <- min(ncol(P0), ncol(pesos0))
    pesos[1L:nf, 1L:nc] <- pesos0[1L:nf, 1L:nc]
  }
  
  return(list("apriori" = P0, "weights" = pesos))
}


## Funtion to compute the weights of the pjk0
calcular_weights <- function(x, y, P0, weights){
  # Test 1
  if (is.character(weights)){
    if( !(weights %in% c("x", "xy", "constant", "expected", "counts", "sqrt", "sd"))){
      stop('Only "x", "xy", "constant", "expected", "counts", "sqrt" or "sd" are strings allowed for the argument "weights"') 
    } else {
      output <- calcular_weights_string(x = x, y = y, P0 = P0, weights = weights)
    }
  } else {
    output <- as.matrix(weights)
#    if (min(weights, na.rm = T) < 0L)
#      stop("Negative weights are not allowed")
#    if(!isTRUE(all.equal(dim(P0), dim(weights)))){ 
#      stop("P0 and weights have different orders")
#    } else {
#      output <- matrix(NA, length(x), length(y))
#      nf <- min(length(x),nrow(weights))
#      nc <- min(length(y), ncol(weights))
#      output[1:nrow(weights), 1:ncol(weights)] <- weights
      # output[is.na(output)] <- Inf
#    }
  }
  return(output)
}


calcular_weights_string <- function(x, y, P0, weights){
  if(weights == "x"){
    output <- matrix(rep(x, length(y)), length(x), length(y))
    if(min(output[output != 0], na.rm = TRUE) > 1L) 
      output[output != 0] <- output[output != 0]/min(output[output != 0], na.rm = TRUE)
  }  
  
  if (weights == "xy"){
    output <- kronecker(x, t(y))
    if(min(output[output != 0], na.rm = TRUE) > 1L) 
      output[output != 0] <- output[output != 0]/min(output[output != 0], na.rm = TRUE)
  }
  
  if (weights == "constant")
    output <- matrix(1L, length(x), length(y))
  
  if (weights == "expected"){
    if (min(P0, na.rm = TRUE) < 0L)
      stop("Negative values in 'apriori' are not allowed")
    min.v <- min(nrow(P0), length(x))
    output <- P0[1:min.v, ]*x[1:min.v]
    # output[output == 0] <- Inf
    if(min(output[output != 0], na.rm = TRUE) > 1L) 
      output[output != 0] <- output[output != 0]/min(output[output != 0], na.rm = TRUE)
  }
  
  if (weights == "counts"){
    if (min(P0, na.rm = TRUE) < 0L)
      stop("Negative values in 'apriori' are not allowed")
    output <- P0 + .5
    if(min(output[output != 0], na.rm = TRUE) > 1L) 
      output[output != 0] <- output[output != 0]/min(output[output != 0], na.rm = TRUE)
  }
  
  if (weights == "sqrt"){
    if (min(P0, na.rm = TRUE) < 0L)
      stop("Negative values in 'apriori' are not allowed")
    output <- sqrt(P0 + .5)
    if(min(output[output != 0], na.rm = TRUE) > 1L) 
      output[output != 0] <- output[output != 0]/min(output[output != 0], na.rm = TRUE)
  }
  
  if (weights == "sd"){
    if (min(P0, na.rm = TRUE) < 0L)
      stop("Negative values in 'apriori' are not allowed")
    P0t <- P0/rowSums(P0, na.rm = TRUE)
    output <- sqrt((P0 + .5)/(P0t*(1L-P0t)))
    output[P0t > 0.999 & !is.na(P0t)] <- sqrt((P0[P0t > 0.999 & !is.na(P0t)] + .5)/(0.995*(1L - 0.995)))
    output[is.infinite(output)] <- sqrt(2)               
    if(min(output[output != 0], na.rm = TRUE) > 1L) 
      output[output != 0] <- output[output != 0]/min(output[output != 0], na.rm = TRUE)
  }
  output0 <- matrix(NA, length(x), length(y))
  output0[1:nrow(output), 1:ncol(output)] <- output
  return(output0)
}


## Function for simultaneous or type I scenarios 
lp_apriori_simultaneous <- function(x, y, P0, weights, solver = "lp_solve", ...){
  #argg <- c(as.list(environment()), list(...))
  pesos.u <- calcular_weights(x = x, y = y, P0 = P0, weights = weights)
  sumo.filas <- rowSums(P0, na.rm = TRUE)
  for (i in 1L:length(sumo.filas)){
    if(sumo.filas[i] > 1L){
      P0[i, ] <- P0[i, ]/sumo.filas[i]
    }
  }

  J <- length(x)
  K <- length(y)
  JK <- J*K
  
  A2 <- cbind(kronecker(diag(J), t(rep(1L, K))), matrix(0L, J, 2*JK))
  A3 <- cbind(kronecker(t(x), diag(K)), matrix(0L, K, 2*JK))
  A4 <- cbind(diag(JK), kronecker(t(c(-1L, 1L)), diag(JK)))
  
  b1 <- c(rep(1L, J), y)
  b2 <- as.vector(t(P0))
  
  f1 <- c(rep(0L, JK))
  f2 <- (as.vector(t(pesos.u)))
  
  # Dealing with missing values
  falta <- which(is.na(b2))
  A4[falta, ] <- 0L
  b2[falta] <- 0L
  f2[falta] <- 0L
  f2[is.na(f2)] <- 0L
  
  A <- rbind(A2, A3, A4)
  b <- c(b1, b2)
  f <- c(f1, f2, f2)
  
  if (solver == "lp_solve")
    solucion <- lpSolve::lp('min', f, A, rep('=', length(b)), b)
  if (solver == "symphony")
    solucion <- Rsymphony::Rsymphony_solve_LP(obj = f,
                                              mat = A,
                                              dir = rep('==', length(b)),
                                              rhs = b)
  mt <- matrix(solucion$solution[1L:JK], J, K, TRUE)
  
  for (j in 1L:J){
    if(x[j] == 0L) mt[j, ] <- 0L
  }
  for (k in 1L:K){
    if(y[k] == 0L) mt[, k] <- 0L
  }
  
  VTM <- round(mt*100, 2)
  
  output <- list("VTM" = VTM, "VTM.complete" = mt,  "error" = solucion$objval,
                 "weights" = pesos.u)
  return(output)
}



## Function for raw1 or type II scenarios
lp_apriori_raw1 <- function(x, y, P0, weights, solver = "lp_solve", ...){
  pesos.u <- calcular_weights(x = x, y = y, P0 = P0, weights = weights)
  sumo.filas <- rowSums(P0, na.rm = TRUE)
  for (i in 1L:length(sumo.filas)){
    if(sumo.filas[i] > 1L){
      P0[i, ] <- P0[i, ]/sumo.filas[i]
    }
  }
  
  J <- length(x)
  K <- length(y)
  JK <- J*K
  
  A2 <- cbind(kronecker(diag(J), t(rep(1L, K))), matrix(0L, J, 2*JK))
  A3 <- cbind(kronecker(t(x), diag(K)), matrix(0L, K, 2*JK))
  A4 <- cbind(diag(JK), kronecker(t(c(-1L,1L)), diag(JK)))
  A5 <- cbind(kronecker(cbind(rep(1L, J - 1L), -diag(J - 1L)), t(c(rep(0L, K - 1L), 1L))), 
              matrix(0L, J - 1L, 2*JK))  
  
  b1 <- c(rep(1L, J), y)
  b2 <- as.vector(t(P0))
  b3 <- rep(0L, J - 1L)
  
  f1 <- c(rep(0L, JK))
  f2 <- (as.vector(t(pesos.u)))
  
  # Dealing with missing values
  falta <- which(is.na(b2))
  A4[falta, ] <- 0L
  b2[falta] <- 0L
  f2[falta] <- 0L
  f2[is.na(f2)] <- 0L
  
  A <- rbind(A2, A3, A4, A5)
  b <- c(b1, b2, b3)
  f <- c(f1, f2, f2)
  
  if (solver == "lp_solve")
    solucion <- lpSolve::lp('min', f, A, rep('=', length(b)), b)
  if (solver == "symphony")
    solucion <- Rsymphony::Rsymphony_solve_LP(obj = f,
                                              mat = A,
                                              dir = rep('==', length(b)),
                                              rhs = b)
  mt <- matrix(solucion$solution[1L:JK], J, K, TRUE)
  
  for (j in 1L:J){
    if(x[j] == 0L) mt[j, ] <- 0L
  }
  for (k in 1L:K){
    if(y[k] == 0L) mt[, k] <- 0L
  }
  return(mt)   
}


## Function for lp_apriori_raw2() for raw2 or type III (when uniform = T) 
##  or type IV (when uniform = F) scenarios
lp_apriori_raw2 <- function(x, y, P0, weights, uniform = TRUE, solver = "lp_solve"){
  pesos.u <- calcular_weights(x = x, y = y, P0 = P0, weights = weights)
  sumo.filas <- rowSums(P0, na.rm = TRUE)
  for (i in 1L:length(sumo.filas)){
    if(sumo.filas[i] > 1L){
      P0[i, ] <- P0[i, ]/sumo.filas[i]
    }
  }
  
  J <- length(x)
  K <- length(y)
  JK <- J*K
  
  A2 <- cbind(kronecker(diag(J), t(rep(1L, K))), matrix(0L, J, 2*JK))
  A3 <- cbind(kronecker(t(x), diag(K)), matrix(0L, K, 2*JK))
  A4 <- cbind(diag(JK), kronecker(t(c(-1L, 1L)), diag(JK)))
  A5 <- cbind(kronecker(cbind(rep(1L, J - 2), -diag(J - 2), rep(0L, J - 2)), t(c(rep(0L, K - 1L), 1L))), 
              matrix(0L, J - 2, 2*JK))  
  A6 <- t(c(rep(0L, JK - 1L), 1L, rep(0L, 2*JK)))
  
  b1 <- c(rep(1L, J), y)
  b2 <- as.vector(t(P0))
  b3 <- rep(0L, J - 1L)
  
  f1 <- c(rep(0L, JK))
  f2 <- (as.vector(t(pesos.u)))
  
  # Dealing with missing values
  falta <- which(is.na(b2))
  A4[falta, ] <- 0L
  b2[falta] <- 0L
  f2[falta] <- 0L
  f2[is.na(f2)] <- 0L
  
  if (uniform){
    A <- rbind(A2, A3, A4, A5, A6)
    b <- c(b1, b2, b3)
    f <- c(f1, f2, f2)
  } else {
    b3 <- 0L
    A <- rbind(A2, A3, A4, A6)
    b <- c(b1, b2, b3)
    f <- c(f1, f2, f2)
  } 
  
  if (solver == "lp_solve")
    solucion <- lpSolve::lp('min', f, A, rep('=', length(b)), b)
  if (solver == "symphony")
    solucion <- Rsymphony::Rsymphony_solve_LP(obj = f,
                                              mat = A,
                                              dir = rep('==', length(b)),
                                              rhs = b)
  mt <- matrix(solucion$solution[1L:JK], J, K, TRUE)
  
  for (j in 1L:J){
    if(x[j] == 0L) mt[j, ] <- 0L
  }
  for (k in 1L:K){
    if(y[k] == 0L) mt[, k] <- 0L
  }
  return(mt)   
}


### Function for raw scenarios: lp_apriori_raw()
lp_apriori_raw <- function(x0, y0, x, y, P0, weights, uniform = TRUE, solver = "lp_solve", ...){
  J0 <- length(x0)
  J <- length(x)
  K0 <- length(y0)
  K <- length(y)
  JK <- J*K
  
  if(J == J0 & J < nrow(P0)){
    warning("A priori probabilities have been provided for new entries to the census, 
            but no net new entries are derived from the provided election results.
            The last row of 'apriori' is therefore omitted.")
    P0 <- P0[-nrow(P0), ]
  }
  
  if(K == K0 & K < ncol(P0)){
    warning("A priori probabilities have been provided for new exits to the census, 
            but no net new exits are derived from the provided election results.
            The last row of 'apriori' is therefore omitted.")
    P0 <- P0[, -ncol(P0)]
  }
    
  if(K > ncol(P0)){
    P0 <- cbind(P0, rep(NA, nrow(P0)))
  }
  if(J > nrow(P0)){
    P0 <- rbind(P0, rep(NA, ncol(P0)))
  }
  
  if (is.matrix(weights)){
    if(nrow(P0) < nrow(weights)) weights <- weights[-nrow(weights), ]
    if(ncol(P0) < ncol(weights)) weights <- weights[, -ncol(weights)]
    if(nrow(P0) > nrow(weights)) weights <- rbind(weights, rep(NA, ncol(weights)))
    if(ncol(P0) > ncol(weights)) weights <- cbind(weights, rep(NA, nrow(weights)))
  }
  
  if (J == J0 & K == K0){
    mt <- lp_apriori_simultaneous(x = x, y = y, P0 = P0, weights = weights, 
                                  solver = solver)$VTM.complete
  }
  if (J == J0 & K > K0){
    if (!uniform){
      mt <- lp_apriori_simultaneous(x = x, y = y, P0 = P0, weights = weights, 
                                    solver = solver)$VTM.complete
    } else {
      mt <- lp_apriori_raw1(x = x, y = y, P0 = P0, weights = weights, 
                            solver = solver)
    }
  }
  
  if (J > J0 & K == K0){
    mt <- lp_apriori_raw2(x = x, y = y, P0 = P0, weights = weights, 
                          uniform = FALSE, solver = solver)
  }
  
  if (J > J0 & K > K0){
    mt <- lp_apriori_raw2(x = x, y = y, P0 = P0, weights = weights, 
                          uniform = uniform, solver = solver)
  }
  
  for (j in 1L:J){
    if(x[j] == 0L) mt[j, ] <- 0L
  }
  for (k in 1L:K){
    if(y[k] == 0L) mt[, k] <- 0L
  }
  
  VTM <- round(mt*100, 2)
  
  output <- list("VTM" = VTM, "VTM.complete" = mt, "weights" = weights)
  
  return(output)
}



### Function for lp_apriori_regular2(), regular2 or type V (when uniform = T)
### or type VI (when uniform = F) scenarios
lp_apriori_regular2 <- function(x, y, P0, weights, uniform, solver = solver){
  pesos.u <- calcular_weights(x = x, y = y, P0 = P0, weights = weights)
  sumo.filas <- rowSums(P0, na.rm = TRUE)
  for (i in 1L:length(sumo.filas)){
    if(sumo.filas[i] > 1L){
      P0[i, ] <- P0[i, ]/sumo.filas[i]
    }
  }
  
  J <- length(x)
  K <- length(y)
  JK <- J*K
  
  A2 <- cbind(kronecker(diag(J), t(rep(1L, K))), matrix(0L, J, 2*JK))
  A3 <- cbind(kronecker(t(x), diag(K)), matrix(0L, K, 2*JK))
  A4 <- cbind(diag(JK), kronecker(t(c(-1L, 1L)), diag(JK)))
  A5 <- cbind(kronecker(cbind(rep(1L, J - 3), -diag(J - 3), matrix(0L, J - 3, 2)), t(c(rep(0L, K - 1L), 1L))), 
              matrix(0L, J - 3, 2*JK))
  A6 <- t(c(rep(0L, (J - 1L)*K - 1L), 1L, rep(0L, 2*JK + K)))
  A7 <- t(c(rep(0L, JK - 1L), 1L, rep(0L, 2*JK)))
  b1 <- c(rep(1L, J), y)
  b2 <- as.vector(t(P0))
  b3 <- rep(0L, J - 1L)
  
  f1 <- c(rep(0L, JK))
  f2 <- (as.vector(t(pesos.u)))
  
  # Dealing with missing values
  falta <- which(is.na(b2))
  A4[falta, ] <- 0L
  b2[falta] <- 0L
  f2[falta] <- 0L
  f2[is.na(f2)] <- 0L
  
  if (uniform){
    A <- rbind(A2, A3, A4, A5, A6, A7)
    b <- c(b1, b2, b3)
    f <- c(f1, f2, f2)
  } else {
    b3 <- rep(0L, 2)
    A <- rbind(A2, A3, A4, A6, A7)
    b <- c(b1, b2, b3)
    f <- c(f1, f2, f2)
  }
  
  if (solver == "lp_solve")
    solucion <- lpSolve::lp('min', f, A, rep('=', length(b)), b)
  if (solver == "symphony")
    solucion <- Rsymphony::Rsymphony_solve_LP(obj = f,
                                              mat = A,
                                              dir = rep('==', length(b)),
                                              rhs = b)
  mt <- matrix(solucion$solution[1L:JK], J, K, TRUE)
  
  for (j in 1L:J){
    if(x[j] == 0L) mt[j, ] <- 0L
  }
  for (k in 1L:K){
    if(y[k] == 0L) mt[, k] <- 0L
  }
  return(mt)   
}


### Function for regular scenarios: lp_apriori_regular()
lp_apriori_regular <- function(x0, y0, x, y, P0, weights, uniform = TRUE, solver = "lp_solve", ...){
  J0 <- length(x0)
  J <- length(x)
  K0 <- length(y0)
  K <- length(y)
  JK <- J*K
  
  if(J == J0 & J < nrow(P0)){
    warning("A priori probabilities for new entries to the census (different from new voters by age) 
            have been provided, but no net new entries are derived from the provided
            election results. The last row of 'apriori' is therefore omitted.")
    P0 <- P0[-nrow(P0), ]
  }
  if(K == K0 & K < ncol(P0)){
    warning("A priori probabilities have been provided for new exits to the census, 
            but no net new exits are derived from the provided election results.
            The last row of 'apriori' is therefore omitted.")
    P0 <- P0[, -ncol(P0)]
  }
  
  if(K > ncol(P0)){
    P0 <- cbind(P0, rep(NA, nrow(P0)))
  }
  if(J > nrow(P0)){
    P0 <- rbind(P0, rep(NA, ncol(P0)))
  }
  
  if (is.matrix(weights)){
    if(nrow(P0) < nrow(weights)) weights <- weights[-nrow(weights), ]
    if(ncol(P0) < ncol(weights)) weights <- weights[, -ncol(weights)]
    if(nrow(P0) > nrow(weights)) weights <- rbind(weights, rep(NA, ncol(weights)))
    if(ncol(P0) > ncol(weights)) weights <- cbind(weights, rep(NA, nrow(weights)))
  }
  
  if (K == K0){
    if (J == J0){
      mt <- lp_apriori_raw2(x = x, y = y, P0 = P0, weights = weights, 
                            uniform = FALSE, solver = solver)
    } else{
      mt <- lp_apriori_regular2(x = x, y = y, P0 = P0, weights = weights, 
                                uniform = FALSE, solver = solver)
    }
  } else {
    if (J == J0){
      mt <- lp_apriori_raw2(x = x, y = y, P0 = P0, weights = weights, 
                            uniform = uniform, solver = solver)
    } else{
      mt <- lp_apriori_regular2(x = x, y = y, P0 = P0, weights = weights, 
                                uniform = uniform, solver = solver)
    }
  }
  
  for (j in 1L:J){
    if(x[j] == 0L) mt[j, ] <- 0L
  }
  for (k in 1L:K){
    if(y[k] == 0L) mt[, k] <- 0L
  }
  
  VTM <- round(mt*100, 2)
  
  output <- list("VTM" = VTM, "VTM.complete" = mt, "weights" = weights)
  
  return(output)
}


### Function lp_apriori_ordinary1(), for ordinary1 or type VII scenarios
lp_apriori_ordinary1 <- function(x, y, P0, weights, solver = "lp_solve", ...){
  pesos.u <- calcular_weights(x = x, y = y, P0 = P0, weights = weights)
  sumo.filas <- rowSums(P0, na.rm = TRUE)
  for (i in 1L:length(sumo.filas)){
    if(sumo.filas[i] > 1L){
      P0[i, ] <- P0[i, ]/sumo.filas[i]
    }
  }
  
  J <- length(x)
  K <- length(y)
  JK <- J*K
  
  A2 <- cbind(kronecker(diag(J), t(rep(1L, K))), matrix(0L, J, 2*JK))
  A3 <- cbind(kronecker(t(x), diag(K)), matrix(0L, K, 2*JK))
  A4 <- cbind(diag(JK), kronecker(t(c(-1L,1L)), diag(JK)))
  A5.1 <- cbind(kronecker(cbind(rep(1L, J - 1L), -diag(J - 1L)), t(c(rep(0L, K - 2), 1L, 0L))), 
                matrix(0L, J - 1L, 2*JK))
  A5.2 <- cbind(kronecker(cbind(rep(1L, J - 1L), -diag(J - 1L)), t(c(rep(0L, K - 1L), 1L))), 
                matrix(0L, J - 1L, 2*JK))  
  
  b1 <- c(rep(1L, J), y)
  b2 <- as.vector(t(P0))
  b3.1 <- rep(0L, J - 1L)
  b3.2 <- rep(0L, J - 1L)
  
  f1 <- c(rep(0L, JK))
  f2 <- (as.vector(t(pesos.u)))
  
  # Dealing with missing values
  falta <- which(is.na(b2))
  A4[falta, ] <- 0L
  b2[falta] <- 0L
  f2[falta] <- 0L
  f2[is.na(f2)] <- 0L
  
  A <- rbind(A2, A3, A4, A5.1, A5.2)
  b <- c(b1, b2, b3.1, b3.2)
  f <- c(f1, f2, f2)
  
  if (solver == "lp_solve")
    solucion <- lpSolve::lp('min', f, A, rep('=', length(b)), b)
  if (solver == "symphony")
    solucion <- Rsymphony::Rsymphony_solve_LP(obj = f,
                                              mat = A,
                                              dir = rep('==', length(b)),
                                              rhs = b)
  mt <- matrix(solucion$solution[1L:JK], J, K, TRUE)
  
  for (j in 1L:J){
    if(x[j] == 0L) mt[j, ] <- 0L
  }
  for (k in 1L:K){
    if(y[k] == 0L) mt[, k] <- 0L
  }
  return(mt)   
}


### Function lp_apriori_ordinary2(), for ordinary2 or type VIII (when uniform = T)
### or type IX (when uniform = F) scenarios
lp_apriori_ordinary2 <- function(x, y, P0, weights, uniform, solver = solver){
  pesos.u <- calcular_weights(x = x, y = y, P0 = P0, weights = weights)
  sumo.filas <- rowSums(P0, na.rm = TRUE)
  for (i in 1L:length(sumo.filas)){
    if(sumo.filas[i] > 1L){
      P0[i, ] <- P0[i, ]/sumo.filas[i]
    }
  }
  
  J <- length(x)
  K <- length(y)
  JK <- J*K
  
  A2 <- cbind(kronecker(diag(J), t(rep(1L, K))), matrix(0L, J, 2*JK))
  A3 <- cbind(kronecker(t(x), diag(K)), matrix(0L, K, 2*JK))
  A4 <- cbind(diag(JK), kronecker(t(c(-1L, 1L)), diag(JK)))
  A5.1 <- cbind(kronecker(cbind(rep(1L, J - 2), -diag(J - 2), rep(0L, J-2)), t(c(rep(0L, K - 2), 1L, 0L))), 
                matrix(0L, J - 2, 2*JK))
  A5.2 <- cbind(kronecker(cbind(rep(1L, J - 2), -diag(J - 2), rep(0L, J-2)), t(c(rep(0L, K - 1L), 1L))), 
                matrix(0L, J - 2, 2*JK))  
  A6.1 <- t(c(rep(0L, JK - 2), 1L, rep(0L, 2*JK + 1L)))
  A6.2 <- t(c(rep(0L, JK - 1L), 1L, rep(0L, 2*JK)))
  
  b1 <- c(rep(1L, J), y)
  b2 <- as.vector(t(P0))
  b3.1 <- rep(0L, J - 1L)
  b3.2 <- rep(0L, J - 1L)
  
  f1 <- c(rep(0L, JK))
  f2 <- (as.vector(t(pesos.u)))
  
  # Dealing with missing values
  falta <- which(is.na(b2))
  A4[falta, ] <- 0L
  b2[falta] <- 0L
  f2[falta] <- 0L
  f2[is.na(f2)] <- 0L
  
  if (uniform){
    A <- rbind(A2, A3, A4, A5.1, A5.2, A6.1, A6.2)
    b <- c(b1, b2, b3.1, b3.2)
    f <- c(f1, f2, f2)
  } else {
    b3 <- rep(0L, 2)
    A <- rbind(A2, A3, A4, A6.1, A6.2)
    b <- c(b1, b2, b3)
    f <- c(f1, f2, f2)
  }
  
  if (solver == "lp_solve")
    solucion <- lpSolve::lp('min', f, A, rep('=', length(b)), b)
  if (solver == "symphony")
    solucion <- Rsymphony::Rsymphony_solve_LP(obj = f,
                                              mat = A,
                                              dir = rep('==', length(b)),
                                              rhs = b)
  mt <- matrix(solucion$solution[1L:JK], J, K, TRUE)
  
  for (j in 1L:J){
    if(x[j] == 0L) mt[j, ] <- 0L
  }
  for (k in 1L:K){
    if(y[k] == 0L) mt[, k] <- 0L
  }
  return(mt)   
}


### Function lp_apriori_ordinary(),  for ordinary scenarios
lp_apriori_ordinary <- function(x0, y0, x, y, P0, weights, uniform = TRUE, solver = "lp_solve", ...){
  J0 <- length(x0)
  J <- length(x)
  K0 <- length(y0)
  K <- length(y)
  JK <- J*K
  
  if(K == K0 & K < ncol(P0)){
    warning("A priori probabilities for new exists to the census (different from deaths)
            have been provided, but no net new exists are derived from the provided
            election results. The last row of 'apriori' is therefore omitted.")
    P0 <- P0[, -ncol(P0)]
  }
  if(J == J0 & J < nrow(P0)){
    warning("A priori probabilities have been provided for new entries to the census, 
            but no net new entries are derived from the provided election results.
            The last row of 'apriori' is therefore omitted.")
    P0 <- P0[-nrow(P0), ]
  }
  
  if(K > ncol(P0)){
    P0 <- cbind(P0, rep(NA, nrow(P0)))
  }
  if(J > nrow(P0)){
    P0 <- rbind(P0, rep(NA, ncol(P0)))
  }
  
  if (is.matrix(weights)){
    if(nrow(P0) < nrow(weights)) weights <- weights[-nrow(weights), ]
    if(ncol(P0) < ncol(weights)) weights <- weights[, -ncol(weights)]
    if(nrow(P0) > nrow(weights)) weights <- rbind(weights, rep(NA, ncol(weights)))
    if(ncol(P0) > ncol(weights)) weights <- cbind(weights, rep(NA, nrow(weights)))
  }
  
  if (K == K0 & J == J0){
    if (uniform){
      mt <- lp_apriori_raw1(x = x, y = y, P0 = P0, weights = weights, solver = solver)
    } else {
      mt <- lp_apriori_simultaneous(x = x, y = y, P0 = P0, 
                                    weights = weights, solver = solver)$VTM.complete  
    }
  }
  
  if (K == K0 & J > J0){  
    mt <- lp_apriori_raw2(x = x, y = y, P0 = P0, weights = weights, 
                          uniform = uniform, solver = solver)
  }
  
  if (K > K0 & J == J0){
    if (uniform){
      mt <- lp_apriori_ordinary1(x = x, y = y, P0 = P0, weights = weights, solver = solver)
    } else{
      mt <- lp_apriori_simultaneous(x = x, y = y, P0 = P0, 
                                    weights = weights, solver = solver)$VTM.complete  
    }
  }
  
  if (K > K0 & J > J0){
    mt <- lp_apriori_ordinary2(x = x, y = y, P0 = P0, weights = weights, 
                               uniform = uniform, solver = solver)
  }
  
  for (j in 1L:J){
    if(x[j] == 0L) mt[j, ] <- 0L
  }
  for (k in 1L:K){
    if(y[k] == 0L) mt[, k] <- 0L
  }
  
  VTM <- round(mt*100, 2)
  
  output <- list("VTM" = VTM, "VTM.complete" = mt, "weights" = weights)
  
  return(output)
}


### Función lp_apriori_enriched1(), for enriched or type X (when uniform = T)
### or type XI (when uniform = F) scenarios
lp_apriori_enriched1 <- function(x, y, P0, weights, uniform, solver = solver){
  pesos.u <- calcular_weights(x = x, y = y, P0 = P0, weights = weights)
  sumo.filas <- rowSums(P0, na.rm = TRUE)
  for (i in 1L:length(sumo.filas)){
    if(sumo.filas[i] > 1L){
      P0[i, ] <- P0[i, ]/sumo.filas[i]
    }
  }
  
  J <- length(x)
  K <- length(y)
  JK <- J*K
  
  A2 <- cbind(kronecker(diag(J), t(rep(1L, K))), matrix(0L, J, 2*JK))
  A3 <- cbind(kronecker(t(x), diag(K)), matrix(0L, K, 2*JK))
  A4 <- cbind(diag(JK), kronecker(t(c(-1L, 1L)), diag(JK)))
  A5.1 <- cbind(kronecker(cbind(rep(1L, J - 3), -diag(J - 3), matrix(0L, J - 3, 2)),
                          t(c(rep(0L, K - 2), 1L, 0L))), 
                matrix(0L, J - 3, 2*JK))
  A5.2 <- cbind(kronecker(cbind(rep(1L, J - 3), -diag(J - 3), matrix(0L, J - 3, 2)),
                          t(c(rep(0L, K - 1L), 1L))), 
                matrix(0L, J - 3, 2*JK))
  
  A6.1 <- t(c(rep(0L, (J - 1L)*K - 2), 1L, rep(0L, 2*JK + K + 1L)))
  A6.2 <- t(c(rep(0L, (J - 1L)*K - 1L), 1L, rep(0L, 2*JK + K)))
  A7.1 <- t(c(rep(0L, JK - 2), 1L, rep(0L, 2*JK + 1L)))
  A7.2 <- t(c(rep(0L, JK - 1L), 1L, rep(0L, 2*JK)))
  
  b1 <- c(rep(1L, J), y)
  b2 <- as.vector(t(P0))
  b3.1 <- rep(0L, J - 1L)
  b3.2 <- rep(0L, J - 1L)
  
  f1 <- c(rep(0L, JK))
  f2 <- (as.vector(t(pesos.u)))
  
  # Dealing with missing values
  falta <- which(is.na(b2))
  A4[falta, ] <- 0L
  b2[falta] <- 0L
  f2[falta] <- 0L
  f2[is.na(f2)] <- 0L
  
  if (uniform){
    A <- rbind(A2, A3, A4, A5.1, A5.2, A6.1, A6.2, A7.1, A7.2)
    b <- c(b1, b2, b3.1, b3.2)
    f <- c(f1, f2, f2)
  } else {
    b3 <- rep(0L, 4)
    A <- rbind(A2, A3, A4, A6.1, A6.2, A7.1, A7.2)
    b <- c(b1, b2, b3)
    f <- c(f1, f2, f2)
  }
  
  if (solver == "lp_solve")
    solucion <- lpSolve::lp('min', f, A, rep('=', length(b)), b)
  if (solver == "symphony")
    solucion <- Rsymphony::Rsymphony_solve_LP(obj = f,
                                              mat = A,
                                              dir = rep('==', length(b)),
                                              rhs = b)
  mt <- matrix(solucion$solution[1L:JK], J, K, TRUE)
  
  for (j in 1L:J){
    if(x[j] == 0L) mt[j, ] <- 0L
  }
  for (k in 1L:K){
    if(y[k] == 0L) mt[, k] <- 0L
  }
  return(mt)   
}


### Function lp_apriori_enriched(), for enriched scenarios
lp_apriori_enriched <- function(x0, y0, x, y, P0, weights, uniform = TRUE, solver = "lp_solve", ...){
  J0 <- length(x0)
  J <- length(x)
  K0 <- length(y0)
  K <- length(y)
  JK <- J*K
  
  if(K == K0 & K < ncol(P0)){
    warning("A priori probabilities for new exists to the census (different from deaths) 
            have been provided, but no net new exists are derived from the provided
            election results. The last row of 'apriori' is therefore omitted.")
    P0 <- P0[, -ncol(P0)]
  }
  
  if(J == J0 & J < nrow(P0)){
    warning("A priori probabilities for new entries to the census (different from new voters by age) 
            have been provided, but no net new entries are derived from the provided
            election results. The last row of 'apriori' is therefore omitted.")
    P0 <- P0[-nrow(P0), ]
  }
  
  if(K > ncol(P0)){
    P0 <- cbind(P0, rep(NA, nrow(P0)))
  }
  if(J > nrow(P0)){
    P0 <- rbind(P0, rep(NA, ncol(P0)))
  }
  
  if (is.matrix(weights)){
    if(nrow(P0) < nrow(weights)) weights <- weights[-nrow(weights), ]
    if(ncol(P0) < ncol(weights)) weights <- weights[, -ncol(weights)]
    if(nrow(P0) > nrow(weights)) weights <- rbind(weights, rep(NA, ncol(weights)))
    if(ncol(P0) > ncol(weights)) weights <- cbind(weights, rep(NA, nrow(weights)))
  }
  
  if (K == K0 & J == J0){
    mt <- lp_apriori_raw2(x = x, y = y, P0 = P0, weights = weights, 
                          uniform = uniform, solver = solver)  
  }
  
  if (K == K0 & J > J0){  
    mt <- lp_apriori_regular2(x = x, y = y, P0 = P0, weights = weights, 
                              uniform = uniform, solver = solver)
  }
  
  if (K > K0 & J == J0){
    mt <- lp_apriori_ordinary2(x = x, y = y, P0 = P0, weights = weights, 
                               uniform = uniform, solver = solver)
  }
  
  if (K > K0 & J > J0){
    mt <- lp_apriori_enriched1(x = x, y = y, P0 = P0, weights = weights, 
                               uniform = uniform, solver = solver)
  }
  
  for (j in 1L:J){
    if(x[j] == 0L) mt[j, ] <- 0L
  }
  for (k in 1L:K){
    if(y[k] == 0L) mt[, k] <- 0L
  }
  
  VTM <- round(mt*100, 2)
  
  output <- list("VTM" = VTM, "VTM.complete" = mt, "weights" = weights)
  
  return(output)
}


## Function lp_apriori_semifull, for semifull scenarios
lp_apriori_semifull <- function(x, y, P0, weights, uniform = TRUE, solver = "lp_solve", ...){
  J <- length(x)
  K <- length(y)
  JK <- J*K
  
  if(K > ncol(P0)){
    P0 <- cbind(P0, rep(NA, nrow(P0)))
  }
  if(J > nrow(P0)){
    P0 <- rbind(P0, rep(NA, ncol(P0)))
  }
  
  if (is.matrix(weights)){
    if(nrow(P0) < nrow(weights)) weights <- weights[-nrow(weights), ]
    if(ncol(P0) < ncol(weights)) weights <- weights[, -ncol(weights)]
    if(nrow(P0) > nrow(weights)) weights <- rbind(weights, rep(NA, ncol(weights)))
    if(ncol(P0) > ncol(weights)) weights <- cbind(weights, rep(NA, nrow(weights)))
  }
  
  mt <- lp_apriori_raw2(x = x, y = y, P0 = P0, weights = weights, 
                            uniform = uniform, solver = solver)
  
  for (j in 1L:J){
    if(x[j] == 0L) mt[j, ] <- 0L
  }
  for (k in 1L:K){
    if(y[k] == 0L) mt[, k] <- 0L
  }
  
  VTM <- round(mt*100, 2)
  
  output <- list("VTM" = VTM, "VTM.complete" = mt, "weights" = weights)
  
  return(output)
}


## Function lp_apriori_fullreverse, for fullreverse scenarios 
lp_apriori_fullreverse <- function(x, y, P0, weights, uniform = TRUE, solver = "lp_solve", ...){
  J <- length(x)
  K <- length(y)
  JK <- J*K
  
  if(K > ncol(P0)){
    P0 <- cbind(P0, rep(NA, nrow(P0)))
  }
  if(J > nrow(P0)){
    P0 <- rbind(P0, rep(NA, ncol(P0)))
  }
  
  if (is.matrix(weights)){
    if(nrow(P0) < nrow(weights)) weights <- weights[-nrow(weights), ]
    if(ncol(P0) < ncol(weights)) weights <- weights[, -ncol(weights)]
    if(nrow(P0) > nrow(weights)) weights <- rbind(weights, rep(NA, ncol(weights)))
    if(ncol(P0) > ncol(weights)) weights <- cbind(weights, rep(NA, nrow(weights)))
  }
  
  mt <- lp_apriori_ordinary2(x = x, y = y, P0 = P0, weights = weights, 
                             uniform = uniform, solver = solver)
  
  for (j in 1L:J){
    if(x[j] == 0L) mt[j, ] <- 0L
  }
  for (k in 1L:K){
    if(y[k] == 0L) mt[, k] <- 0L
  }
  
  VTM <- round(mt*100, 2)
  
  output <- list("VTM" = VTM, "VTM.complete" = mt, "weights" = weights)
  
  return(output)
}


## Function lp_apriori_full, for full scenarios 
lp_apriori_full <- function(x, y, P0, weights, uniform = TRUE, solver = "lp_solve", ...){
  J <- length(x)
  K <- length(y)
  JK <- J*K
  
  if(K > ncol(P0)){
    P0 <- cbind(P0, rep(NA, nrow(P0)))
  }
  if(J > nrow(P0)){
    P0 <- rbind(P0, rep(NA, ncol(P0)))
  }
  
  if (is.matrix(weights)){
    if(nrow(P0) < nrow(weights)) weights <- weights[-nrow(weights), ]
    if(ncol(P0) < ncol(weights)) weights <- weights[, -ncol(weights)]
    if(nrow(P0) > nrow(weights)) weights <- rbind(weights, rep(NA, ncol(weights)))
    if(ncol(P0) > ncol(weights)) weights <- cbind(weights, rep(NA, nrow(weights)))
  }
  
  mt <- lp_apriori_regular2(x = x, y = y, P0 = P0, weights = weights, 
                            uniform = uniform, solver = solver)
  
  for (j in 1L:J){
    if(x[j] == 0L) mt[j, ] <- 0L
  }
  for (k in 1L:K){
    if(y[k] == 0L) mt[, k] <- 0L
  }
  
  VTM <- round(mt*100, 2)
  
  output <- list("VTM" = VTM, "VTM.complete" = mt, "weights" = weights)
  
  return(output)
}

## Function lp_apriori_gold, for gold scenarios
lp_apriori_gold <- function(x, y, P0, weights, uniform = TRUE, solver = "lp_solve", ...){
  J <- length(x)
  K <- length(y)
  JK <- J*K
  
  if(K > ncol(P0)){
    P0 <- cbind(P0, rep(NA, nrow(P0)))
  }
  if(J > nrow(P0)){
    P0 <- rbind(P0, rep(NA, ncol(P0)))
  }
  
  if (is.matrix(weights)){
    if(nrow(P0) < nrow(weights)) weights <- weights[-nrow(weights), ]
    if(ncol(P0) < ncol(weights)) weights <- weights[, -ncol(weights)]
    if(nrow(P0) > nrow(weights)) weights <- rbind(weights, rep(NA, ncol(weights)))
    if(ncol(P0) > ncol(weights)) weights <- cbind(weights, rep(NA, nrow(weights)))
  }
  
  mt <- lp_apriori_enriched1(x = x, y = y, P0 = P0, weights = weights, 
                             uniform = uniform, solver = solver)
  
  for (j in 1:J){
    if(x[j] == 0L) mt[j, ] <- 0L
  }
  for (k in 1:K){
    if(y[k] == 0L) mt[, k] <- 0L
  }
  
  VTM <- round(mt*100, 2)
  
  output <- list("VTM" = VTM, "VTM.complete" = mt, "weights" = weights)
  
  return(output)
}


