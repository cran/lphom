#' Implements the lphom_joint algorithm
#'
#' @description Estimates RxC vote transfer matrices (ecological contingency tables) with lphom_joint
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
#'    \item{VTM.votes}{ A matrix of order JxK with the estimated cross-distribution of votes of elections 1 and 2.}
#'    \item{HETe}{ The estimated heterogeneity index associated to the `VTM.votes` solution.}
#'    \item{VTM12}{ The matrix of order JxK with the estimated row-standardized proportions of vote transitions from election 1 to election 2 associated to the `VTM.votes` solution.}
#'    \item{VTM21}{ The matrix of order KxJ with the estimated row-standardized proportions of vote transitions from election 2 to election 1 associated to the `VTM.votes` solution.}
#'    \item{EHet12}{ A matrix of order IxK measuring in each unit a distance to the homogeneity hypothesis. That is, the differences under the homogeneity hypothesis between the actual recorded results and the expected results in each territorial unit for each option of election two. The matrix Eik.}
#'    \item{EHet21}{ A matrix of order IxJ measuring in each unit a distance to the homogeneity hypothesis. That is, the differences under the homogeneity hypothesis between the actual recorded results and the expected results in each territorial unit for each option of election one. The matrix Eij.}
#'    \item{deterministic.bounds}{ A list of two matrices of order JxK containing for each vote transition the lower and upper proportions allowed given the observed aggregates.}
#'    \item{inputs}{ A list containing all the objects with the values used as arguments by the function.}
#' @export
#'
#' @family linear programing ecological inference functions
#' @seealso \code{\link{lphom}} \code{\link{lphom_dual}} \code{\link{tslphom_dual}} \code{\link{nslphom_dual}} \code{\link{tslphom_joint}} \code{\link{nslphom_joint}}
#'
#' @examples
#' x <- France2017P[, 1:8]
#' y <- France2017P[, 9:12]
#' y[,1] <- y[,1]  - (rowSums(y) - rowSums(x))
#' mt <- lphom_joint(x, y)
#' mt$VTM.votes
#' mt$HETe
#
#' @importFrom lpSolve lp
#

lphom_joint <- function(votes_election1,
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
  
  # Parameters
  I <- nrow(votes_election1);
  J <- ncol(votes_election1);
  K <- ncol(votes_election2)
  JK <- J*K; IK <- I*K; IJ <- I*J

  modelo12 <- model_LP(votes_election1, votes_election2)
  modelo21 <- model_LP(votes_election2, votes_election1)

  a <- rbind(cbind(modelo12$a, matrix(0L, nrow(modelo12$a), ncol(modelo21$a))),
             cbind(matrix(0L, nrow(modelo21$a), ncol(modelo12$a)), modelo21$a))
  b <- c(modelo12$b, modelo21$b)
  f <- c(modelo12$f, modelo21$f)

  # We add the las J*K constraints of congruence
  yt <- colSums(votes_election2)
  ap <- kronecker(diag(colSums(votes_election1)), diag(K))
  aq <- matrix(0L, JK, JK)
  for (k in 1L:K){
    aq <- aq + kronecker(t(diag(K)[,k]),
                         kronecker(diag(J), -yt[k]*diag(K)[,k]))
  }
  a <- rbind(a,
             cbind(ap, matrix(0L, JK, 2L*IK),
                   aq, matrix(0L, JK, 2L*IJ)))
  b <- c(b, rep(0L, JK))

  names1 <- colnames(votes_election1)
  names2 <- colnames(votes_election2)

  # Solution
  if (solver == "lp_solve"){
    sol <- suppressWarnings(lpSolve::lp('min', f, a, rep('=', length(b)), b))
  } else if (solver == "symphony") {
    sol = Rsymphony::Rsymphony_solve_LP(obj = f,
                                        mat = a,
                                        dir = rep('==', length(b)),
                                        rhs = b)
  }
  z <- sol$solution

  VTM1 <- matrix(z[1L:JK], J, K, TRUE, dimnames = list(names1, names2))
  VTM2 <- matrix(z[ncol(modelo12$a) + (1:JK)], K, J, TRUE, dimnames = list(names2, names1))
  VTM.votos <- VTM1 * colSums(votes_election1)

  if (integers){
    VTM.votos <- dec2counts(VTM.votos, rowSums(VTM.votos), colSums(VTM.votos))
    VTM1 <- VTM.votos/rowSums(VTM.votos)
    VTM2 <- t(VTM.votos)/colSums(VTM.votos)
    dimnames(VTM.votos) <- dimnames(VTM1) <- list(names1, names2)
    dimnames(VTM2) <- list(names2, names1)
  }

  eik <- votes_election2 - as.matrix(votes_election1) %*% VTM1
  eij <- votes_election1 - as.matrix(votes_election2) %*% VTM2

 # e <- z[(JK+1):(JK+2*IK)]
 # e <- matrix(e,2,IK)
 # e <- e[1,] - e[2,]
 # eik <- t(matrix(e,K,I))
 # colnames(eik) <- names2
 # rownames(eik) <- rownames(votes_election1)

 # e <- z[ncol(modelo12$a) + ((JK+1):(JK+2*IJ))]
 # e <- matrix(e,2,IJ)
 # e <- e[1,] - e[2,]
 # eij <- t(matrix(e, J, I))
 # colnames(eij) <- names1
 # rownames(eij) <- rownames(votes_election1)

  HETe <- 50*(sum(abs(eik)) + sum(abs(eij)))/sum(VTM.votos)
  
  det.bounds <- bounds_compound(origin = votes_election1, 
                                destination = votes_election2, zeros = NULL)[c(1,2)]

  output <- list("VTM.votes" = VTM.votos, "HETe" = HETe, "VTM12" = VTM1, "VTM21" = VTM2,
                 "EHet.12" = eik, "EHet.21" = eij, "deterministic.bounds" = det.bounds, 
                 "inputs" = inputs)
  class(output) <- c("lphom_joint", "ei_joint", "lphom")
  return(output)
}
