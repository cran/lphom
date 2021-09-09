#' Implements the nslphom_joint algorithm
#'
#' @description Estimates RxC vote transfer matrices (ecological contingency tables) with nslphom_joint
#'
#' @author Jose M. Pavia, \email{pavia@@uv.es}
#' @author Rafael Romero \email{rromero@@eio.upv.es}
#' @references Pavia, JM and Romero, R (2021). Symmetry estimating R×C vote transfer matrices from aggregate data, mimeo.
#'
#' @param votes_election1 data.frame (or matrix) of order IxJ with the votes gained by the *J*
#'                        political options competing on election 1 (or origin) in the *I*
#'                        territorial units considered. In general, the counts to be initially
#'                        mapped to columns. The sum by rows of `votes_election1` and
#'                        `votes_election2` must coincide.
#'
#' @param votes_election2 data.frame (or matrix) of order IxK with the votes gained by
#'                        the *K* political options competing on election 2 (or destination)
#'                        in the *I* territorial units considered. In general, the counts to be
#'                        initially mapped to columns. The sum by rows of `votes_election1` and
#'                        `votes_election2` must coincide.
#'
#' @param iter.max Maximum number of iterations to be performed. The process ends independently when either
#'                 the number of iterations reaches iter.max or when the maximum variation between two
#'                 consecutive estimates of both ways probability transfer matrices are less than `tol`.
#'                 By default, 10.
#'
#' @param min.first A TRUE/FALSE value. If FALSE, the matrix associated with the minimum `HETe` after
#'                  performing `iter.max` iterations is taken as solution.
#'                  If TRUE, the associated matrix to the instant in which the first decrease of `HETe` occurs
#'                  is taken as solution. The process stops at that moment. In this last scenario
#'                  (when `min.first = TRUE`), `iter.max` is is forced to be at least 100. Default, FALSE.
#'
#' @param counts A TRUE/FALSE value that indicates whether the problem is solved in integer values (counts) in
#'               each iteration: zero (lphom) and intermediate and final (including unit) solutions.
#'               If TRUE, the initial LP matrices are approximated in each iteration to the closest integer solution
#'               solving the corresponding Integer Linear Program. Default, FALSE.
#'
#' @param solver A character string indicating the linear programming solver to be used, only
#'               `lp_solve` and `symphony` are allowed. By default, `lp_solve`.
#'
#' @param tol Maximum deviation allowed between two consecutive iterations. The process ends when the maximum
#'            variation between the estimated cross-distributions of votes between two consecutive
#'            iterations is less than `tol` or the maximum number of iterations has been reached. By default, 0.001.
#'
#' @return
#' A list with the following components
#'    \item{VTM.votes}{ A matrix of order JxK with the estimated cross-distribution of votes of elections 1 and 2.}
#'    \item{HETe}{ The estimated heterogeneity index associated to the `VTM.votes` solution.}
#'    \item{VTM12}{ The matrix of order JxK with the estimated row-standardized proportions of vote transitions from election 1 to election 2 associated to the `VTM.votes` solution.}
#'    \item{VTM21}{ The matrix of order KxJ with the estimated row-standardized proportions of vote transitions from election 2 to election 1 associated to the `VTM.votes` solution.}
#'    \item{VTM.votes.units}{ An array of order JxKxI with the estimated matrix of cross-distributions of votes of elections 1 and 2 attained for each unit in iteration of the solution.}
#'    \item{iter}{ The real final number of iterations performed before ending the process.}
#'    \item{iter.min}{ Number of the iteration associated to the selected `VTM.votes` solution.}
#'    \item{inputs}{ A list containing all the objects with the values used as arguments by the function.}
#'    \item{solution_init}{ A list with the main outputs produced by **lphom_joint()**.}
#'
#' @export
#'
#' @family linear programing ecological inference functions
#'
#' @seealso \code{\link{nslphom}} \code{\link{lphom_dual}} \code{\link{tslphom_dual}} \code{\link{nslphom_dual}} \code{\link{lphom_joint}} \code{\link{tslphom_joint}}
#'
#' @examples
#' x <- France2017P[, 1:8]
#' y <- France2017P[, 9:12]
#' y[,1] <- y[,1]  - (rowSums(y) - rowSums(x))
#' mt <- nslphom_joint(x, y)
#' mt$VTM.votes
#' mt$HETe
#
#' @importFrom Rsymphony Rsymphony_solve_LP
#' @importFrom lpSolve lp
#

nslphom_joint <- function(votes_election1,
                        votes_election2,
                        iter.max = 10,
                        min.first = FALSE,
                        counts = FALSE,
                        solver = "lp_solve",
                        tol = 0.001){

  # Calculo de la solucion inicial con lphom_joint
  lphom_inic <- lphom_joint(votes_election1 = votes_election1, votes_election2 = votes_election2,
                            counts = counts, solver = solver)

  if (min.first) iter.max <- max(100L, iter.max)

  # Inicio proceso iterativo
  VTM.votes.sequence <- array(NA, c(dim(lphom_inic$VTM.votes), iter.max + 1L))
  VTM.votes_units.iter <- array(NA, c(dim(lphom_inic$VTM.votes),
                                       nrow(votes_election1)))

  iter <- 0L
  HETe.iter <- Inf
  dif.max <- Inf
  VTM.iter <- lphom_inic$VTM.votes
  VTM.votes.sequence[, , iter + 1L] <- VTM.iter
  HETe.sequence <- lphom_inic$HETe

  # Proceso iterativo
  while (iter < iter.max & dif.max > tol){
    # Calculo de las soluciones locales
    for (i in 1L:nrow(votes_election1)){
      VTM.votes_units.iter[, , i] <- MT_joint_local(MT = VTM.iter,
                                                    marginal_fila = votes_election1[i, ],
                                                    marginal_columna = votes_election2[i, ],
                                                    solver = solver)
    }
    if (counts){
      for (i in 1L:nrow(votes_election1)){
        VTM.votes_units.iter[, , i] <- dec2counts(VTM.votes_units.iter[, , i],
                                                  votes_election1[i, ],
                                                  votes_election2[i, ])
      }
    }

    HETe <- HET_joint(VTM.votes_units.iter)$HETe
    iter <- iter + 1L
    VTM.votos <- apply(VTM.votes_units.iter, c(1,2), sum)
    dif.max <- max(abs(VTM.votos - VTM.iter))
    VTM.iter <- VTM.votos
    VTM.votes.sequence[, , iter + 1L] <- VTM.iter
    HETe.sequence <- c(HETe.sequence, HETe)

    if (HETe < HETe.iter){
      VTM.votes <- VTM.votos
      VTM.votes_units <- VTM.votes_units.iter
      VTM1 <- VTM.votes/rowSums(VTM.votes)
      VTM2 <- t(VTM.votes)/rowSums(t(VTM.votes))
      HETe.iter <- HETe
      iter.min <- iter
    }
    if (min.first & (HETe.sequence[iter + 1L] > HETe.sequence[iter])) dif.max <- -Inf
  } # End while

  # Eliminar sobrantes en arrays si converga antes de alcanzar el máximo de iteraciones
  VTM.votes.sequence <- VTM.votes.sequence[, , 1L:(iter + 1L)]

  # Mejorando la salida
  dimnames(VTM.votes) <- dimnames(VTM1) <- dimnames(lphom_inic$VTM.votes)
  dimnames(VTM2) <- dimnames(lphom_inic$VTM2)
  dimnames(VTM.votes_units) <- c(dimnames(lphom_inic$VTM.votes), list(rownames(votes_election1)))
  dimnames(VTM.votes.sequence) <- c(dimnames(VTM.votes), list(c(0:iter)))

  inic <- lphom_inic[1:6]
  names(inic) <- paste0(names(inic), "_init")

  inputs <- lphom_inic$inputs
  inputs$iter.max <- iter.max
  inputs$min.first <- min.first
  inputs$tol <- tol

  return(list("VTM.votes" = VTM.votes, "HETe" = HETe.iter, "VTM12" =VTM1, "VTM21" = VTM2,
              "HETe.sequence" = HETe.sequence, "VTM.votes.sequence" = VTM.votes.sequence,
              "VTM.votes.units" = VTM.votes_units, "iter" = iter, "iter.min" = iter.min,
              "inputs" = inputs, "solution_init" = inic))
}
