#' Implements lphom_dual algorithm
#'
#' @description Estimates RxC vote transfer matrices (ecological contingency tables) with lphom_dual
#'
#' @author Jose M. Pavia, \email{pavia@@uv.es}
#' @author Rafael Romero \email{rromero@@eio.upv.es}
#' @references Pavia, JM and Romero, R (2024). Symmetry estimating RxC vote transfer matrices from aggregate data. *Journal of the Royal Statistical Society, Series A – Statistics in Society*, forthcoming.  \doi{10.1093/jrsssa/qnae013}
#'
#' @param votes_election1 data.frame (or matrix) of order IxJ with the counts to be initially
#'                        mapped to rows. When estimating vote transfer matrices, the votes gained by 
#'                        the *J* political options competing on election 1 (or origin) in the *I*
#'                        territorial units considered.  The sum by rows of `votes_election1` and
#'                        `votes_election2` must coincide.
#'
#' @param votes_election2 data.frame (or matrix) of order IxK with the counts to be initially mapped 
#'                        to columns. When estimating vote transfer matrices, the votes gained by
#'                        the *K* political options competing on election 2 (or destination) in the *I* 
#'                        territorial units considered. The sum by rows of `votes_election1` and
#'                        `votes_election2` must coincide.
#'
#' @param integers A TRUE/FALSE value that indicates whether the LP solution of counts (votes) must be approximate
#'                 to the closest integer solution using ILP. Default, FALSE.
#'
#' @param solver A character string indicating the linear programming solver to be used, only
#'               `lp_solve` and `symphony` are allowed. By default, `lp_solve`. The package `Rsymphony`
#'               needs to be installed for the option `symphony` to be used.
#'
#' @param integers.solver A character string indicating the linear programming solver to be used to approximate
#'                        to the closest integer solution, only `symphony` and `lp_solve` are allowed.
#'                        By default, `symphony`. The package `Rsymphony` needs to be installed for the option `symphony` 
#'                        to be used. Only used when `integers = TRUE`. 
#'                        
#' @param ... Other arguments to be passed to the function. Not currently used.
#'  

#' @return
#' A list with the following components
#'    \item{VTM.votes.w}{ The matrix of order JxK with the estimated cross-distribution of votes of elections 1 and 2, attained weighting the two dual solutions using as weights the corresponding HTEe estimates.}
#'    \item{VTM.votes.a}{ The matrix of order JxK with the estimated cross-distribution of votes of elections 1 and 2, attained simple averaging the two dual solutions.}
#'    \item{HETe.w}{ Estimated heterogeneity index associated to the `VTM.votes.w` solution.}
#'    \item{HETe.a}{ Estimated heterogeneity index associated to the `VTM.votes.a` solution.}
#'    \item{VTM12.w}{ The matrix of order JxK with the estimated row-standardized proportions of vote transitions from election 1 to election 2 associated to the `VTM.votes.w` solution.}
#'    \item{VTM21.w}{ The matrix of order KxJ with the estimated row-standardized proportions of vote transitions from election 2 to election 1 associated to the `VTM.votes.w` solution.}
#'    \item{VTM12.a}{ The matrix of order JxK with the estimated row-standardized proportions of vote transitions from election 1 to election 2 associated to the `VTM.votes.a` solution.}
#'    \item{VTM21.a}{ The matrix of order KxJ with the estimated row-standardized proportions of vote transitions from election 2 to election 1 associated to the `VTM.votes.a` solution.}
#'    \item{lphom.object.12}{ The output of the \code{\link{lphom}} function attained solving the problem X --> Y. That is, mapping `votes_election1` to rows and `votes_election2` to columns.}
#'    \item{lphom.object.21}{ The output of the \code{\link{lphom}} function attained solving the problem Y --> X. That is, mapping `votes_election2` to rows and `votes_election1` to columns.}
#'    \item{inputs}{ A list containing all the objects with the values used as arguments by the function.}
#' @export
#'
#' @family linear programing ecological inference functions
#'
#' @seealso \code{\link{lphom}} \code{\link{tslphom_dual}} \code{\link{nslphom_dual}} \code{\link{lphom_joint}} \code{\link{tslphom_joint}} \code{\link{nslphom_joint}}
#'
#' @examples
#' x <- France2017P[, 1:8]
#' y <- France2017P[, 9:12]
#' y[,1] <- y[,1]  - (rowSums(y) - rowSums(x))
#' mt <- lphom_dual(x, y)
#' mt$VTM.votes.w
#' mt$HETe.w
#
#' @importFrom lpSolve lp
#


lphom_dual <- function(votes_election1,
                       votes_election2,
                       integers = FALSE,
                       solver = "lp_solve",
                       integers.solver = "symphony",
                       ...){

  inputs <- c(as.list(environment()), list(...))
  integers <- inputs$integers <- test_integers(argg = inputs)
  
  if (integers.solver == "lp_solve"){
    dec2counts <- dec2counts_lp
  } else {
    dec2counts <- dec2counts_symphony
  }
  
  # library(lphom)
  lphom.object.12 <- lphom(votes_election1, votes_election2, "simultaneous",
                           integers = integers, solver = solver,
                           integers.solver = integers.solver)
  lphom.object.21 <- lphom(votes_election2, votes_election1, "simultaneous",
                           integers = integers, solver = solver,
                           integers.solver = integers.solver)
  votos12 <- lphom.object.12$VTM.complete.votes
  votos21 <- t(lphom.object.21$VTM.complete.votes)
  VTM.votos <- (votos12 + votos21)/2
  VTM.votos.weigthed <- (votos12*lphom.object.12$HETe^-1 + votos21*lphom.object.21$HETe^-1)/
    (lphom.object.12$HETe^-1 + lphom.object.21$HETe^-1)
  VTM1 <- VTM.votos/rowSums(VTM.votos)
  VTM2 <- t(VTM.votos)/colSums(VTM.votos)
  VTM1.weighted <- VTM.votos.weigthed/rowSums(VTM.votos.weigthed)
  VTM2.weighted <- t(VTM.votos.weigthed)/colSums(VTM.votos.weigthed)

  if (integers){
    VTM.votos <- dec2counts(VTM.votos, rowSums(VTM.votos), colSums(VTM.votos))
    VTM.votos.weigthed <- dec2counts(VTM.votos.weigthed, rowSums(VTM.votos.weigthed),
                                     colSums(VTM.votos.weigthed))
    VTM1 <- VTM.votos/rowSums(VTM.votos)
    VTM2 <- t(VTM.votos)/colSums(VTM.votos)
    VTM1.weighted <- VTM.votos.weigthed/rowSums(VTM.votos.weigthed)
    VTM2.weighted <- t(VTM.votos.weigthed)/colSums(VTM.votos.weigthed)
  }

  HETe1 <- sum(abs(as.matrix(votes_election1) %*% VTM1 - votes_election2))
  HETe2 <- sum(abs(as.matrix(votes_election2) %*% VTM2 - votes_election1))
  HETe.a <- 50*(HETe1 + HETe2)/sum(VTM.votos)
  HETe1 <- sum(abs(as.matrix(votes_election1) %*% VTM1.weighted - votes_election2))
  HETe2 <- sum(abs(as.matrix(votes_election2) %*% VTM2.weighted - votes_election1))
  HETe.w <- 50*(HETe1 + HETe2)/sum(VTM.votos.weigthed)

  dimnames(VTM.votos.weigthed) <- dimnames(VTM.votos) <-
    dimnames(VTM1.weighted) <- dimnames(VTM1) <- dimnames(lphom.object.12$VTM)
  dimnames(VTM2.weighted) <- dimnames(VTM2) <- dimnames(lphom.object.12$OTM)

  output <- list("VTM.votes.w" = VTM.votos.weigthed, "VTM.votes.a" = VTM.votos, "HETe.w" = HETe.w,
              "HETe.a" = HETe.a, "VTM12.w" = VTM1.weighted, "VTM21.w" = VTM2.weighted,
              "VTM12.a" = VTM1, "VTM21.a" = VTM2, "lphom.object.12" = lphom.object.12,
              "lphom.object.21" = lphom.object.21, "inputs" = inputs)
  class(output) <- c("lphom_dual", "ei_dual", "lphom")
  return(output)
}
