#' Integer-adjusts outputs of the lphom-family functions
#'
#' @description Takes as input an object generated with an algorithm of the lphom-family 
#' (lphom, tslphom, nslphom, tslphom_dual, nslphom_joint, ....) and returns 
#' as output an object of the same class as the input object with all their relevant estimated (local and global) transfer matrices 
#' of counts updated to their closest integer matrices. The rest of main components of the object are also accordingly updated.
#'
#' @author Jose M. Pavia, \email{pavia@@uv.es}
#' @references ...
#'
#' @param x An object output of a lphom family algorithm
#' @param ... Other arguments passed on the method. Not currently used. 
#'                        
#' @details
#' The updating of the matrices is performed using integer linear programming after imposing all the row- and column-constraints. 
#'                                                 
#' @return
#' An object of the same class and components as `x` with its components properly updated after adjusting the estimated count matrices in `x` using integer linear programming
#' 
#' @export
#' 
#' @examples
#' mt.ns <- nslphom(France2017P[, 1:8] , France2017P[, 9:12], new_and_exit_voters= "raw")
#' mt.ns <- adjust2integers(mt.ns)
#

adjust2integers <- function(x, ...){
  UseMethod("adjust2integers")
}
  
#' @export
adjust2integers.lphom <- function(x, ...) 
{ 
  
  if (inherits(x, "ei_lp")){
    ajustar_con <- adjust2integers_ei_lp
  } else if (inherits(x, "ei_dual")){
    ajustar_con <- adjust2integers_ei_dual 
  } else if (inherits(x, "ei_joint")){
    ajustar_con <- adjust2integers_ei_joint
  } else {
    ajustar_con <- adjust2integers_default
  }
  
  ajustar_con(x = x, ...)
}


adjust2integers_default <- function(x, ...){
  warning(paste0("adjust2integers does not know how to handle an object of class ", 
                class(x), 
                ".\n  It has been devised to deal with outputs of functions of the lphom-family."))
  
}

adjust2integers_lphom <- function(x, ...){
  y <- x
  y$VTM.complete.votes <- dec2counts(x$VTM.complete.votes, 
                                     colSums(x$origin),
                                     colSums(x$destination))
  dimnames(y$VTM.complete.votes) <- dimnames(x$VTM.complete.votes)
  y$VTM.votes <- y$VTM.complete.votes[1L:nrow(x$VTM.votes), 1L:ncol(x$VTM.votes)]
  y$VTM.complete <- y$VTM.complete.votes/rowSums(y$VTM.complete.votes)
  y$VTM <- round(100*y$VTM.complete[1L:nrow(x$VTM), 1L:ncol(x$VTM)], 2)
  y$OTM <- t(y$VTM.complete.votes)/colSums(y$VTM.complete.votes)
  y$OTM <- round(100*y$OTM[1L:nrow(x$OTM), 1L:ncol(x$OTM)], 2)
  y$EHet <- x$destination - x$origin %*% y$VTM.complete
  y$HETe <- 100*sum(abs(y$EHet))/sum(y$VTM.complete.votes)#
  
  filas0 <- which(rowSums(y$VTM.votes) == 0)
  colum0 <- which(colSums(y$VTM.votes) == 0)
  y$VTM[filas0, ] <- 0
  y$VTM.complete[filas0, ] <- 0
  y$OTM[colum0, ] <- 0
  
  return(y)
}

adjust2integers_ei_lp <- function(x, ...){
  if (class(x)[1] == "lphom"){
    adjust2integers_lphom(x = x, ...)
  } else {
    y <- x
    for (i in 1L:nrow(x$origin)){
      y$VTM.votes.units[, , i] <- dec2counts(x$VTM.votes.units[, , i],
                                             x$origin[i,], x$destination[i,])
      y$VTM.prop.units[, , i] <- y$VTM.votes.units[, , i]/rowSums(y$VTM.votes.units[, , i])
    }
    y$VTM.prop.units[is.na(y$VTM.prop.units)] <- 0L
    y$VTM.complete.votes <- apply(y$VTM.votes.units, c(1, 2), sum)
    y$VTM.votes <- y$VTM.complete.votes[1L:nrow(x$VTM.votes), 1L:ncol(x$VTM.votes)]
    y$VTM.complete <- y$VTM.complete.votes/rowSums(y$VTM.complete.votes)
    y$VTM <- round(100*y$VTM.complete[1L:nrow(x$VTM), 1L:ncol(x$VTM)], 2)
    y$OTM <- t(y$VTM.complete.votes)/colSums(y$VTM.complete.votes)
    y$OTM <- round(100*y$OTM[1L:nrow(x$OTM), 1L:ncol(x$OTM)], 2)
    y$EHet <- x$destination - x$origin %*% y$VTM.complete
    y$HETe <- HET_MT.votos_MT.prop_Y(y$VTM.votes.units)$HET

    filas0 <- which(rowSums(y$VTM.complete.votes) == 0)
    colum0 <- which(colSums(y$VTM.complete.votes) == 0)
    y$VTM[filas0, ] <- 0
    y$VTM.complete[filas0, ] <- 0
    y$OTM[colum0, ] <- 0
    return(y)
  }
}

adjust2integers_lphom_dual <- function(x, ...){
  y <- x
  y$VTM.votes.w <- dec2counts(x$VTM.votes.w, 
                              colSums(x$lphom.object.12$origin),
                              colSums(x$lphom.object.12$destination))
  y$VTM.votes.a <- dec2counts(x$VTM.votes.a, 
                              colSums(x$lphom.object.12$origin),
                              colSums(x$lphom.object.12$destination))
  dimnames(y$VTM.votes.w) <- dimnames(y$VTM.votes.a) <- dimnames(x$VTM.votes.w)
  y$VTM12.w <- y$VTM.votes.w/rowSums(y$VTM.votes.w)
  y$VTM21.w <- t(y$VTM.votes.w)/colSums(y$VTM.votes.w)
  y$VTM12.a <- y$VTM.votes.a/rowSums(y$VTM.votes.a)
  y$VTM21.a <- t(y$VTM.votes.a)/colSums(y$VTM.votes.a)
  
  HETe1 <- sum(abs(as.matrix(x$lphom.object.12$origin) %*% y$VTM12.a - x$lphom.object.12$destination))
  HETe2 <- sum(abs(as.matrix(x$lphom.object.12$destination) %*% y$VTM21.a - x$lphom.object.12$origin))
  y$HETe.a <- 50*(HETe1 + HETe2)/sum(y$VTM.votes.a)
  HETe1 <- sum(abs(as.matrix(x$lphom.object.12$origin) %*% y$VTM12.w - x$lphom.object.12$destination))
  HETe2 <- sum(abs(as.matrix(x$lphom.object.12$destination) %*% y$VTM21.w - x$lphom.object.12$origin))
  y$HETe.w <- 50*(HETe1 + HETe2)/sum(y$VTM.votes.w)
  
  return(y)
}

adjust2integers_ei_dual <- function(x, ...){
  if (class(x)[1] == "lphom_dual"){
    adjust2integers_lphom_dual(x = x, ...)
  } else {
    y <- x
    for (i in 1L:dim(x$VTM.votes.units.w)[3]){
      y$VTM.votes.units.w[, , i] <- dec2counts(x$VTM.votes.units.w[, , i],
                                               round(rowSums(x$VTM.votes.units.w[, , i])), 
                                               round(colSums(x$VTM.votes.units.w[, , i])))
  
      y$VTM.votes.units.a[, , i] <- dec2counts(x$VTM.votes.units.a[, , i],
                                               round(rowSums(x$VTM.votes.units.a[, , i])), 
                                               round(colSums(x$VTM.votes.units.a[, , i])))
    }
    dimnames(y$VTM.votes.units.w) <- dimnames(y$VTM.votes.units.a) <- c(dimnames(x$VTM.votes.w),
                                                                        list(rownames(x$tslphom.object.12$origin)))
    
    y$VTM.votes.w <- apply(y$VTM.votes.units.w, c(1, 2), sum)
    y$VTM.votes.a <- apply(y$VTM.votes.units.a, c(1, 2), sum)
   
    y$VTM12.w <- y$VTM.votes.w/rowSums(y$VTM.votes.w)
    y$VTM21.w <- t(y$VTM.votes.w)/colSums(y$VTM.votes.w)
    y$VTM12.a <- y$VTM.votes.a/rowSums(y$VTM.votes.a)
    y$VTM21.a <- t(y$VTM.votes.a)/colSums(y$VTM.votes.a)
    
    y$HETe.a <- HET_MT.votos_MT.prop_Y(y$VTM.votes.units.a)$HET
    y$HETe.w <- HET_MT.votos_MT.prop_Y(y$VTM.votes.units.w)$HET
    
    return(y)
  }
}


adjust2integers_lphom_joint <- function(x, ...){
  y <- x
  y$VTM.votes <- dec2counts(x$VTM.votes, 
                            colSums(x$inputs$votes_election1),
                            colSums(x$inputs$votes_election2))
  dimnames(y$VTM.votes) <- dimnames(x$VTM.votes)
  y$VTM12 <- y$VTM.votes/rowSums(y$VTM.votes)
  y$VTM21 <- t(y$VTM.votes)/colSums(y$VTM.votes)
  
  y$EHet.12 <- as.matrix(x$inputs$votes_election2) - as.matrix(x$inputs$votes_election1) %*% y$VTM12
  y$EHet.21 <- as.matrix(x$inputs$votes_election1) - as.matrix(x$inputs$votes_election2) %*% y$VTM21

  y$HETe <- 50*(sum(abs(y$EHet.12)) + sum(abs(y$EHet.21)))/sum(y$VTM.votes)
  
  return(y)
}

adjust2integers_ei_joint <- function(x, ...){
  if (class(x)[1] == "lphom_joint"){
    adjust2integers_lphom_joint(x = x, ...)
  } else {
    # invisible(validObject(x))
    y <- x
    for (i in 1L:nrow(x$inputs$votes_election1)){
      y$VTM.votes.units[, , i] <- dec2counts(x$VTM.votes.units[, , i],
                                             x$inputs$votes_election1[i,], 
                                             x$inputs$votes_election2[i,])
    }
    
    y$VTM.votes <- apply(y$VTM.votes.units, c(1, 2), sum)
    
    y$VTM12 <- y$VTM.votes/rowSums(y$VTM.votes)
    y$VTM21 <- t(y$VTM.votes)/colSums(y$VTM.votes)
    y$EHet.12 <- as.matrix(x$inputs$votes_election2) - as.matrix(x$inputs$votes_election1) %*% y$VTM12
    y$EHet.21 <- as.matrix(x$inputs$votes_election1) - as.matrix(x$inputs$votes_election2) %*% y$VTM21
    
    y$HETe <- HET_joint(y$VTM.votes.units)$HETe
    
    return(y)
  }
}
