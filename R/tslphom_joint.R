#' Implements the tslphom_joint algorithm
#'
#' @description Estimates RxC vote transfer matrices (ecological contingency tables) with tslphom_joint
#'
#' @author Jose M. Pavia, \email{pavia@@uv.es}
#' @author Rafael Romero \email{rromero@@eio.upv.es}
#' @references Pavia, JM and Romero, R (2021). Symmetry estimating RxC vote transfer matrices from aggregate data, mimeo.
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
#'               `lp_solve` and `symphony` are allowed. By default, `lp_solve`.
#'
#' @param ... Other arguments to be passed to the function. Not currently used.
#' 

#' @return
#' A list with the following components
#'    \item{VTM.votes}{ A matrix of order JxK with the estimated cross-distribution of votes of elections 1 and 2.}
#'    \item{HETe}{ The estimated heterogeneity index associated to the `VTM.votes` solution.}
#'    \item{VTM12}{ The matrix of order JxK with the estimated row-standardized proportions of vote transitions from election 1 to election 2 associated to the `VTM.votes` solution.}
#'    \item{VTM21}{ The matrix of order KxJ with the estimated row-standardized proportions of vote transitions from election 2 to election 1 associated to the `VTM.votes` solution.}
#'    \item{VTM.votes.units}{ An array of order JxKxI with the estimated matrix of cross-distributions of votes of elections 1 and 2 attained for each unit after congruently adjusting the **lphom_joint()** initial estimate.}
#'    \item{EHet12}{ A matrix of order IxK measuring in each unit a distance to the homogeneity hypothesis. That is, the differences under the homogeneity hypothesis between the actual recorded results and the expected results in each territorial unit for each option of election two. The matrix Eik.}
#'    \item{EHet21}{ A matrix of order IxJ measuring in each unit a distance to the homogeneity hypothesis. That is, the differences under the homogeneity hypothesis between the actual recorded results and the expected results in each territorial unit for each option of election one. The matrix Eij.}
#'    \item{inputs}{ A list containing all the objects with the values used as arguments by the function.}
#'    \item{solution_init}{ A list with the main outputs produced by **lphom_joint()**.}
#'
#' @export
#'
#' @family linear programing ecological inference functions
#'
#' @seealso \code{\link{tslphom}} \code{\link{lphom_dual}} \code{\link{tslphom_dual}} \code{\link{nslphom_dual}} \code{\link{lphom_joint}} \code{\link{nslphom_joint}}
#'
#' @examples
#' x <- France2017P[, 1:8]
#' y <- France2017P[, 9:12]
#' y[,1] <- y[,1]  - (rowSums(y) - rowSums(x))
#' mt <- tslphom_joint(x, y)
#' mt$VTM.votes
#' mt$HETe
#
#' @importFrom Rsymphony Rsymphony_solve_LP
#' @importFrom lpSolve lp
#

tslphom_joint <- function(votes_election1,
                        votes_election2,
                        integers = FALSE,
                        solver = "lp_solve",
                        ...){
  
  argg <- c(as.list(environment()), list(...))
  integers <- test_integers(argg)
  
  # Calculo de la solucion inicial con lphom_joint
  lphom_inic <- lphom_joint(votes_election1 = votes_election1, votes_election2 = votes_election2,
                            integers = integers, solver = solver)

  # Calculo de las soluciones locales
  VTM.votes_units <- array(NA, c(dim(lphom_inic$VTM.votes), nrow(votes_election1)))
  for (i in 1L:nrow(votes_election1)){
    VTM.votes_units[, , i] <- MT_joint_local(MT = lphom_inic$VTM.votes,
                                             marginal_fila = votes_election1[i, ],
                                             marginal_columna = votes_election2[i, ],
                                             solver = solver)
  }
  if (integers){
    for (i in 1L:nrow(votes_election1)){
      VTM.votes_units[, , i] <- dec2counts(VTM.votes_units[, , i],
                                           votes_election1[i, ],
                                           votes_election2[i, ])
    }
  }

  VTM.votes <- apply(VTM.votes_units, c(1,2), sum)
  HETe <- HET_joint(VTM.votes_units)$HETe
  VTM1 <- VTM.votes/rowSums(VTM.votes)
  VTM2 <- t(VTM.votes)/rowSums(t(VTM.votes))

  dimnames(VTM.votes) <- dimnames(VTM1) <- dimnames(lphom_inic$VTM.votes)
  dimnames(VTM2) <- dimnames(lphom_inic$VTM2)
  dimnames(VTM.votes_units) <- c(dimnames(lphom_inic$VTM.votes), list(rownames(votes_election1)))
  
  eik <- votes_election2 - as.matrix(votes_election1) %*% VTM1
  eij <- votes_election1 - as.matrix(votes_election2) %*% VTM2

  inic <- lphom_inic[1L:6L]
  names(inic) <- paste0(names(inic), "_init")

  output <- list("VTM.votes" = VTM.votes, "HETe" = HETe, "VTM12" = VTM1,
              "VTM21" = VTM2, "VTM.votes.units" = VTM.votes_units, "EHet.12" = eik, 
              "EHet.21" = eij, "inputs" = lphom_inic$inputs, "solution_init" = inic)
  class(output) <- c("tslphom_joint", "ei_joint", "lphom")
  return(output)
}
