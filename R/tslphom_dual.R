#' Implements the tslphom_dual algorithm
#'
#' @description Estimates RxC vote transfer matrices (ecological contingency tables) with tslphom_dual
#'
#' @author Jose M. Pavia, \email{pavia@@uv.es}
#' @author Rafael Romero \email{rromero@@eio.upv.es}
#' @references Pavia, JM and Romero, R (2021). Symmetry estimating R×C vote transfer matrices from aggregate data, mimeo.
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
#'                        territorial units considered. In general, The sum by rows of `votes_election1` and
#'                        `votes_election2` must coincide.
#'
#' @param integers A TRUE/FALSE value that indicates whether the problem is solved in integer values
#'                 in both iterations: zero (lphom) and final (including unit) solutions. If TRUE, the LP matrices
#'                 are approximated to the closest integer solution solving the corresponding Integer Linear Program.
#'                 Default, FALSE.
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
#'  \item{VTM.votes.w}{ The matrix of order JxK with the estimated cross-distribution of votes of elections 1 and 2,
#'                      attained weighting the two dual solutions using as weights the corresponding HTEe estimates.}
#'  \item{VTM.votes.units.w}{ The array of order JxKxI with the local estimated cross-distributions of votes of elections 1 and 2 by unit,
#'                      attained weighting the two dual solutions using as weights the corresponding HTEe estimates.}
#'  \item{VTM.votes.a}{ The matrix of order JxK with the estimated cross-distribution of votes of elections 1 and 2,
#'                      attained simple averaging the two dual solutions.}
#'  \item{VTM.votes.units.a}{ The matrix of order JxKxI with the estimated cross-distributions of votes of elections 1 and 2 by unit,
#'                      attained weighting the two dual solutions using as weights the corresponding HTEe estimates.}
#'  \item{HETe.w}{ Estimated heterogeneity index associated to the `VTM.votes.w` solution.}
#'  \item{HETe.a}{ Estimated heterogeneity index associated to the `VTM.votes.a` solution.}
#'  \item{VTM12.w}{ The matrix of order JxK with the estimated row-standardized proportions of vote transitions from election 1
#'                  to election 2 associated to the `VTM.votes.w` solution.}
#'  \item{VTM21.w}{ The matrix of order KxJ with the estimated row-standardized proportions of vote transitions from election 2
#'                  to election 1 associated to the `VTM.votes.w` solution.}
#'  \item{VTM12.a}{ The matrix of order JxK with the estimated row-standardized proportions of vote transitions from election 1
#'                  to election 2 associated to the `VTM.votes.a` solution.}
#'  \item{VTM21.a}{ The matrix of order KxJ with the estimated row-standardized proportions of vote transitions from election 2
#'                  to election 1 associated to the `VTM.votes.a` solution.}
#'  \item{tslphom.object.12}{ The output of the \code{\link{tslphom}} function attained solving the problem X --> Y,
#'                          that is, mapping `votes_election1` to rows and `votes_election2` to columns.}
#'  \item{tslphom.object.21}{ The output of the \code{\link{tslphom}} function attained solving the problem Y --> X,
#'                          that is, mapping `votes_election2` to rows and `votes_election1` to columns.}
#'  \item{inputs}{ A list containing all the objects with the values used as arguments by the function.}
#' @export
#'
#' @family linear programing ecological inference functions
#' @seealso \code{\link{tslphom}} \code{\link{lphom_dual}} \code{\link{nslphom_dual}} \code{\link{lphom_joint}} \code{\link{tslphom_joint}} \code{\link{nslphom_joint}}
#'
#' @examples
#' x <- France2017P[, 1:8]
#' y <- France2017P[, 9:12]
#' y[,1] <- y[,1]  - (rowSums(y) - rowSums(x))
#' mt <- tslphom_dual(x, y)
#' mt$VTM.votes.w
#' mt$HETe.w
#'

#' @importFrom lpSolve lp
#


tslphom_dual <- function(votes_election1,
                         votes_election2,
                         integers = FALSE,
                         solver = "lp_solve",
                         integers.solver = "symphony",
                         ...){
  
  argg <- c(as.list(environment()), list(...))
  integers <- test_integers(argg)
  
  if (integers.solver == "lp_solve"){
    dec2counts <- dec2counts_lp
  } else {
    dec2counts <- dec2counts_symphony
  }
  
  inputs <- list("votes_election1" = votes_election1, "votes_election2" = votes_election2,
                 "integers" = integers, "solver" = solver, 
                 "integers.solver" = integers.solver)

  lphom.object.12 <- tslphom(votes_election1, votes_election2, "simultaneous",
                           integers = integers, solver = solver,
                           integers.solver = integers.solver)
  lphom.object.21 <- tslphom(votes_election2, votes_election1, "simultaneous",
                           integers = integers, solver = solver,
                           integers.solver = integers.solver)

  votos.units.a <- (lphom.object.12$VTM.votes.units +
                      aperm(lphom.object.21$VTM.votes.units, c(2, 1, 3)))/2
  votos.units.w <- (lphom.object.12$VTM.votes.units*lphom.object.12$HETe^-1 +
                      aperm(lphom.object.21$VTM.votes.units, c(2, 1, 3))*lphom.object.21$HETe^-1)/
    (lphom.object.12$HETe^-1 + lphom.object.21$HETe^-1)

  if (integers){
    for (i in 1L:nrow(lphom.object.12$origin)){
      votos.units.a[, , i] <- dec2counts(votos.units.a[, , i],
                                         lphom.object.12$origin[i,],
                                         lphom.object.12$destination[i,])
      votos.units.w[, , i] <- dec2counts(votos.units.w[, , i],
                                         lphom.object.12$origin[i,],
                                         lphom.object.12$destination[i,])
    }
  }

  VTM.votos <- apply(votos.units.a, c(1,2), sum)
  VTM.votos.weigthed <- apply(votos.units.w, c(1,2), sum)

  VTM1 <- VTM.votos/rowSums(VTM.votos)
  VTM2 <- t(VTM.votos)/colSums(VTM.votos)
  VTM1.weighted <- VTM.votos.weigthed/rowSums(VTM.votos.weigthed)
  VTM2.weighted <- t(VTM.votos.weigthed)/colSums(VTM.votos.weigthed)

  HETe.a <- HET_MT.votos_MT.prop_Y(votos.units.a)$HET
  HETe.w <- HET_MT.votos_MT.prop_Y(votos.units.w)$HET

  dimnames(VTM.votos.weigthed) <- dimnames(VTM.votos) <-
    dimnames(VTM1.weighted) <- dimnames(VTM1) <- dimnames(lphom.object.12$VTM)
  dimnames(VTM2.weighted) <- dimnames(VTM2) <- dimnames(lphom.object.12$OTM)

  dimnames(votos.units.a) <- dimnames(votos.units.w) <- dimnames(lphom.object.12$VTM_votes)

  output <- list("VTM.votes.w" = VTM.votos.weigthed, "VTM.votes.units.w" = votos.units.w,
              "VTM.votes.a" = VTM.votos, "VTM.votes.units.a" = votos.units.a,
              "HETe.w" = HETe.w, "HETe.a" = HETe.a, "VTM12.w" = VTM1.weighted,
              "VTM21.w" = VTM2.weighted, "VTM12.a" = VTM1, "VTM21.a" = VTM2,
              "tslphom.object.12" = lphom.object.12, "tslphom.object.21" = lphom.object.21,
              "inputs" = inputs)
  class(output) <- c("tslphom_dual", "ei_dual", "lphom")
  return(output)
}
