#' 	Summarize a lphom-family object
#'
#' @description Summary method for objects obtained with an algorithm of the lphom-family (lphom, tslphom, nslphom, tslphom_dual, nslphom_joint, ....).
#'
#' @author Jose M. Pavia, \email{pavia@@uv.es}
#'
#' @param object An object output of a **lphom** family algorithm.
#' @param ... Other arguments passed on to methods. Not currently used.
#' 
#' @return
#' An object of class `"summary_lphom"`. 
#' 
#'
#' @export
#' @method summary lphom
#' @examples
#' mt.ns <- nslphom(France2017P[, 1:8] , France2017P[, 9:12], new_and_exit_voters= "raw")
#' summary(mt.ns)
#'

summary.lphom <- function(object, ...){
  
  if (inherits(object, "ei_lp")){
    resumir_con <- summary_ei_lp
  } else if (inherits(object, "ei_dual")){
    resumir_con <- summary_ei_dual 
  } else if (inherits(object, "ei_joint")){
    resumir_con <-summary_ei_joint
  }
  
  resumir_con(object = object, ...)
}



summary_ei_lp  <- function(object, ...) 
{
  
  output <- NULL
  output$prop.matrix <- object$VTM
  output$counts.matrix <- object$VTM.votes

  output$row.margins <- round(rowSums(object$VTM.complete.votes))
  output$col.margins <- round(colSums(object$VTM.complete.votes))
  
  output$HETe <- object$HETe
  
  class(output) <- "summary.lphom"
  output
    
}

summary_ei_dual  <- function(object, ...)
{
  
  output <- NULL
  output$prop.matrix.w <- object$VTM12.w*100
  output$counts.matrix.w <- object$VTM.votes.w

  output$prop.matrix.a <- object$VTM12.a*100
  output$counts.matrix.a <- object$VTM.votes.a
  
  output$row.margins <- round(rowSums(object$VTM.votes.w))
  output$col.margins <- round(colSums(object$VTM.votes.w))
  
  output$HETe.w <- object$HETe.w
  output$HETe.a <- object$HETe.a
  
  class(output) <- "summary.lphom"
  output
}

summary_ei_joint  <- function(object, ...)
{
  
  output <- NULL
  output$prop.matrix <- object$VTM12*100
  output$counts.matrix <- object$VTM.votes
  
  output$row.margins <- round(rowSums(object$VTM.votes))
  output$col.margins <- round(colSums(object$VTM.votes))
  
  output$HETe <- object$HETe
  
  class(output) <- "summary.lphom"
  output
}
