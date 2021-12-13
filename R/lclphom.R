#' Implements lclphom algorithm
#'
#' @description  Estimates RxC vote transfer matrices (ecological contingency tables) with lclphom
#'
#' @author Jose M. Pavia, \email{pavia@@uv.es}
#' @author Rafael Romero \email{rromero@@eio.upv.es}
#'
#' @param votes_election1 data.frame (or matrix) of order IxJ (likely of final order IxJ-1
#'                        in `regular` and `raw` scenarios) with the votes gained by the *J*
#'                        political options competing on election 1 (or origin) in the *I*
#'                        territorial units considered. In general, the row marginals 
#'                        of the *I* tables.
#'
#' @param votes_election2 data.frame (or matrix) of order IxK (likely of final order IxK-1
#'                        in `regular` and `raw` scenarios) with the votes gained by
#'                        the *K* political options competing on election 2 (or destination)
#'                        in the *I* territorial units considered. In general, the column marginals 
#'                        of the *I* tables.
#'
#' @param new_and_exit_voters A character string indicating the level of information available
#'                            regarding new entries and exits of the election censuses between the
#'                            two elections. This argument captures the different options discussed
#'                            on Section 3 of Romero et al. (2020). This argument admits five values:
#'                            `raw`, `regular`, `simultaneous`, `full` and `gold`. Default, `raw`.
#'                            The argument `simultaneous` should be used in a typical ecological inference 
#'                            problem.
#'
#' @param structural_zeros Default, NULL. A list of vectors of length two, indicating the election options
#'                         for which no transfer of votes are allowed between election 1 and election 2.
#'                         For instance, when new_and_exit_voters is set to `"regular"`,
#'                         lphom implicitly `states structural_zeros = list(c(J, K))` in case exits and/or
#'                         entries are computed because the sum by rows of `votes_election1` and
#'                         `votes_election2` does not coincide.
#'
#' @param iter.max Maximum number of iterations to be performed. The process ends when either the
#'                 number of iterations reaches `iter.max` or when there is no error reduction in any
#'                 local unit between two consecutive iterations. By default, 1000.
#'
#' @param uniform A TRUE/FALSE value that indicates if census exits affects all the electoral options in a
#'                (relatively) similar fashion in each voting unit: equation (13) of Pavia and Romero (2021a).
#'                Default, TRUE.
#'
#' @param distance.local A string argument that indicates whether the second step of the lphom_local algorithm
#'                       should be performed to solve potential indeterminacies of local solutions.
#'                       Default, `"abs"`.
#'                       If `distance.local = "abs"` lphom_local selects in its second step the matrix
#'                       closer to the temporary global solution under L_1 norm, among the first step compatible matrices.
#'                       If `distance.local = "max"` lphom_local selects in its second step the matrix
#'                       closer to the temporary global solution under L_Inf norm, among the first step compatible matrices.
#'                       If `distance.local = "none"`, the second step of lphom_local is not performed.
#'
#' @param counts A TRUE/FALSE value that indicates whether the problem is solved in integer values (counts) in
#'               each iteration, including zero (lphom) and local solutions. Initial LP matrices are
#'               approximate to the closest integer solution solving the corresponding Integer Linear Program.
#'               Default, FALSE.
#'
#' @param solver A character string indicating the linear programming solver to be used, only
#'               `symphony` and `lp_solve` are allowed. By default, `lp_solve`.
#'
#' @param verbose A TRUE/FALSE value that indicates if the main outputs of the function should be
#'                printed on the screen. Default, FALSE.
#'
#' @details Description of the `new_and_exit_voters` argument in more detail.
#' \itemize{
#'  \item{`raw`: }{The default value. This argument accounts for the most plausible scenario when
#'                 estimating vote transfer matrices: A scenario with two elections elapsed at least some
#'                 months where only the raw election data recorded in the *I* territorial units, 
#'                 in which the area under study is divided, are available. 
#'                 In this scenario, net exits (basically deaths) and net entries (basically 
#'                 new young voters) are estimated according to equation (7) of Romero et al. (2020). 
#'                 Constraints defined by equations (8) and (9) of Romero et al. (2020) and (13) of Pavia
#'                 and Romero (2021a) are imposed. 
#'                 In this scenario, when net exits and/or net entries are negligible (such as between 
#'                 the first- and second-round of French Presidential elections), they are omitted in 
#'                 the outputs.}
#'  \item{`regular`: }{ For estimating vote transfer matrices, this value accounts for a scenario with 
#'                 two elections elapsed at least some months where (i) the column *J* of `votes_election1` 
#'                 corresponds to new young electors who have the right to vote for the first time and (ii)
#'                 net exits (basically a consequence of mortality), and eventually net entries,
#'                 are computed according equation (7) of Romero et al. (2020), and (iii) within each unit 
#'                 it is assummed that net exits affect equally all the first *J-1* options 
#'                 of election 1. Hence when `uniform = TRUE` equation (13) of Pavia and 
#'                 Romero (2021a) is applied. Constraints (8) and (9) of Romero et al. (2020)
#'                 are imposed to start the process.}
#'  \item{`simultaneous`: }{ This is the value to be used in a classical ecological inference problems, 
#'                such as for racial voting, and in a scenario with two simultaneous elections. 
#'                In this scenario, the sum by rows of `votes_election1` and `votes_election2` must coincide. 
#'                Constraints defined by equations (8) and (9) of Romero et al. (2020) and (13) of Pavia and 
#'                Romero (2021a) are not included in the model.}
#'  \item{`full`: }{This value accounts for a scenario with two elections elapsed at least some
#'                months, where: (i) the column *J-1* of votes_election1 totals new young
#'                electors that have the right to vote for the first time; (ii) the column *J*
#'                of votes_election1 measures new immigrants that have the right to vote; and
#'                (iii) the column *K* of votes_election2 corresponds to total exits of the census
#'                lists (due to death or emigration). In this scenario, the sum by rows of
#'                `votes_election1` and `votes_election2` must agree and constraints (8)
#'                and (9) of Romero et al. (2020) are imposed.}
#'  \item{`gold`: }{This value accounts for a scenario similar to full, where total exits are
#'               separated out between exits due to emigration (column *K-1* of `votes_election2`)
#'               and death (column *K* of `votes_election2`). In this scenario, the sum by rows
#'               of `votes_election1` and `votes_election2` must agree. The same restrictions
#'               as in the above scenario apply but for both columns *K-1* and *K* of the vote
#'               transition probability matrix}
#' }
#'
#' @return
#' A list with the following components
#'  \item{VTM}{ A matrix of order JxK with the estimated percentages of row-standardized vote transitions from election 1 to election 2.}
#'  \item{VTM.votes}{ A matrix of order JxK with the estimated vote transitions from election 1 to election 2.}
#'  \item{OTM}{ A matrix of order KxJ with the estimated percentages of the origin of the votes obtained for the different options of election 2.}
#'  \item{HETe}{ The estimated heterogeneity index as defined in equation (15) of Pavia and Romero (2021a).}
#'  \item{VTM.complete}{ A matrix of order J'xK' with the estimated proportions of row-standardized vote transitions from election 1 to election 2, including in `regular` and `raw` scenarios the row and the column corresponding to net_entries and net_exits even when they are really small, less than 1% in all units.}
#'  \item{VTM.complete.votes}{ A matrix of order J'xK' with the estimated vote transitions from election 1 to election 2, including in `regular` and `raw` scenarios the row and the column corresponding to net_entries and net_exits even when they are really small, less than 1% in all units.}
#'  \item{VTM.prop.units}{ An array of order J'xK'xI with the estimated proportions of vote transitions from election 1 to election 2 attained for each unit in the solution.}
#'  \item{VTM.votes.units}{ An array of order J'xK'xI with the estimated matrix of vote transitions from election 1 to election 2 attained for for each unit in the solution.}
#'  \item{VTM.complete.last.iter}{ A matrix of order J'xK' with the estimated proportions of vote transitions from election 1 to election 2, including in `regular` and `raw` scenarios the row and the column corresponding to net_entries and net_exits even when they are really small, less than 1% in all units, corresponding to the final iteration.}
#'  \item{VTM.sequence}{ Array of order J'xK'x(iter+1) (where `iter` is the efective number of iterations performed) of the intermediate estimated matrices corresponding to each iteration.}
#'  \item{HETe.sequence}{ Numeric vector of length `iter+1` with the `HETe` coefficients corresponding to the matrices in `VTM.sequence`.}
#'  \item{VTM.prop.units.last.iter}{ An array of order J'xK'xI with the estimated proportions of vote transitions from election 1 to election 2 attained for each unit in the final iteration.}
#'  \item{VTM.votes.units.last.iter}{ An array of order J'xK'xI with the estimated matrix of vote transitions from election 1 to election 2 attained for each unit in the final iteration.}
#'  \item{zeros}{ A list of vectors of length two, indicating the election options for which no transfer of votes are allowed between election 1 and election 2.}
#'  \item{iter}{ The real final number of iterations performed before ending the process.}
#'  \item{iter.units}{ A matrix of order Ix(iter+1) with the number of iteration corresponding to the solution selected for each unit
#'                    in each iteration.}
#'  \item{errors}{ A vector of length I with the minimal error observed in the sequence for each unit. It corresponds to
#'                the unit-error associated with the solution linked with either `VTM.prop.units` or `VTM.votes.units`.}
#'  \item{inputs}{ A list containing all the objects with the values used as arguments by the function.}
#'  \item{origin}{ A matrix with the final data used as votes of the origin election after taking into account the level of information available regarding to new entries and exits of the election censuses between the two elections.}
#'  \item{destination}{ A matrix with the final data used as votes of the origin election after taking into account the level of information available regarding to new entries and exits of the election censuses between the two elections.}
#'  \item{EHet}{ A matrix of order IxK measuring in each spatial unit a distance to the homogeneity hypothesis, that is, the differences under the homogeneity hypothesis between the actual recorded results and the expected results with the solution in each territorial unit for each option of election 2.}
#'  \item{solution_init}{ A list with the main outputs produced by **lphom()**.}
#'  \itemize{
#'  \item{`VTM_init`:}{ A matrix of order JxK with the estimated percentages of vote transitions from election 1 to election 2 initially obtained by **lphom()**.}
#'  \item{`VTM.votes_init`:}{ A matrix of order JxK with the estimated vote transitions from election 1 to election 2 initially obtained by **lphom()**.}
#'  \item{`OTM_init`:}{ A matrix of order KxJ with the estimated percentages of the origin of the votes obtained for the different options of election 2 initially obtained by **lphom()**.}
#'  \item{`HETe_init`:}{ The estimated heterogeneity index defined in equation (10) of Romero et al. (2020). }
#'  \item{`EHet_init`:}{ A matrix of order IxK measuring in each spatial unit the distance to the homogeneity hypothesis, that is, the differences under the homogeneity hypothesis between the actual recorded results and the expected results, using the **lphom()** solution, in each territorial unit for each option of election 2.}
#'  \item{`VTM.complete_init`:}{ A matrix of order J'xK' with the estimated proportions of vote transitions from election 1 to election 2 initially obtained by **lphom()**, including in `regular` and `raw` scenarios the row and the column corresponding to net_entries and net_exits even when they are really small, less than 1% in all units.}
#'  \item{`VTM.complete.votes_init`:}{ A matrix of order J'xK' with the estimated vote transitions from election 1 to election 2 initially obtained by **lphom()**, including in `regular` and `raw` scenarios the row and the column corresponding to net_entries and net_exits even when they are really small, less than 1% in all units.}
#' }
#'
#' @export
#'
#' @family linear programing ecological inference functions
#' @seealso \code{\link{lphom}} \code{\link{tslphom}} \code{\link{nslphom}}
#'
#' @examples
#' mt.lc <- lclphom(France2017P[, 1:8] , France2017P[, 9:12], new_and_exit_voters= "raw")
#' mt.lc$VTM
#' mt.lc$HETe
#' mt.lc$solution_init$HETe_init
#'
#' @importFrom Rsymphony Rsymphony_solve_LP
#' @importFrom lpSolve lp
#'

lclphom <- function(votes_election1,
                    votes_election2,
                    new_and_exit_voters = c("raw", "regular", "simultaneous", "full", "gold"),
                    structural_zeros = NULL,
                    iter.max = 1000,
                    uniform = TRUE,
                    distance.local = c("abs", "max", "none"),
                    counts = FALSE,
                    solver = "lp_solve",
                    verbose = FALSE ){

  if (iter.max < 0 | iter.max%%1 > 0)
    stop('iter.max must be a positive integer')
  if (!(distance.local[1L] %in% c("abs", "max", "none")))
    stop('Not allowed string for argument "distance.local".
         The only allowed strings for "distance.local" are "abs", "max" and "none".')

  # Funcion local a aplicar
  lphom_unit <- lp_solver_local(uniform = uniform,
                                distance.local = distance.local[1])

  # Calculo de la solucion inicial
  lphom_inic <- lphom(votes_election1 = votes_election1, votes_election2 = votes_election2,
                      new_and_exit_voters = new_and_exit_voters, structural_zeros,
                      counts = counts, verbose = FALSE, solver = solver)
  lphom0 <- lphom_inic

  zeros <- determinar_zeros_estructurales(lphom_inic)

  # Inicio proceso iterativo
  VTM.sequence <- array(NA, c(dim(lphom_inic$VTM.complete), iter.max + 1L))
  iter <- 0L
  n.changes <- Inf
  VTM.iter <- lphom0$VTM.complete <- lphom_inic$VTM.complete
  VTM.sequence[, , iter+1L] <- VTM.iter
  HETe.sequence <- lphom_inic$HETe
  iter.units <- matrix(0L,  nrow(lphom_inic$origin), iter.max + 1L)
  errors <- rep(Inf, nrow(lphom_inic$origin))
  VTM_units.iter <- votos_units.iter <- array(NA, c(dim(lphom_inic$VTM.complete), nrow(lphom_inic$origin)))

  while (iter < iter.max & n.changes > 0){
    # Calculo de las soluciones locales
    VTM_units.t <- votos_units.t <- array(NA, c(dim(lphom_inic$VTM.complete), nrow(lphom_inic$origin)))
    for (i in 1L:nrow(lphom_inic$origin)){
      VTM_units.t[, , i] <- lphom_unit(lphom.object = lphom0, iii = i, solver = solver)
      if (counts){
        votos_units.t[, , i] <- dec2counts(VTM_units.t[, , i]*lphom_inic$origin[i, ],
                                          lphom_inic$origin[i,], lphom_inic$destination[i,])
        VTM_units.t[, , i] <- votos_units.t[, , i]/rowSums(votos_units.t[, , i])
      } else {
        votos_units.t[, , i] <- VTM_units.t[, , i]/rowSums(VTM_units.t[, , i])*lphom_inic$origin[i, ]
      }
      VTM_units.t[lphom_inic$origin[i, ] == 0L, , i] <- 0L
    }
    votos_units.t[is.na(votos_units.t)] <- 0L
    Errors.t <- HET_MT.votos_MT.prop_Y(votos_units.t)$Error.H
    iter <- iter + 1L
    # Unidades para las que mejora el error asociado a la solución local
    n.changes <- Errors.t < errors
    iter.units[, iter + 1L] <- iter.units[, iter]
    iter.units[n.changes, iter + 1L] <- iter
    errors[n.changes] <- Errors.t[n.changes]
    VTM_units.iter[, , n.changes] <- VTM_units.t[, , n.changes]
    votos_units.iter[, , n.changes] <- votos_units.t[, , n.changes]
    VTM_votos_homogeneos <- HET_MT.votos_MT.prop_Y(votos_units.iter)
    VTM.iter <- lphom0$VTM.complete <- VTM_votos_homogeneos$MT.pro
    VTM.sequence[, , iter + 1L] <- VTM.iter
    HETe.sequence <- c(HETe.sequence, VTM_votos_homogeneos$HET)
    n.changes <- sum(n.changes)
  } # End while
  if (iter != iter.max) iter <- iter - 1L
  iter.units <- iter.units[, 1L:(iter + 1L)]
  VTM.sequence <- VTM.sequence[, , 1L:(iter + 1L)]
  dimnames(iter.units) <- list(rownames(lphom_inic$origin),
                            paste0("iter = ", 0L:iter))
  dimnames(VTM.sequence) <- c(dimnames(lphom_inic$VTM.complete),
                              list(paste0("iter = ", 0L:iter)))
  HETe.sequence <- HETe.sequence[1L:(iter + 1L)]
  VTM.complete.iter <- VTM.iter
  dimnames(VTM_units.iter) <- dimnames(votos_units.iter) <- c(dimnames(lphom_inic$VTM.complete),
                                                              list(rownames(lphom_inic$origin)))

  # Obtaining the final solution
  VTM_units <- votos_units <- array(NA, c(dim(lphom_inic$VTM.complete), nrow(lphom_inic$origin)))
  for (i in 1L:nrow(lphom_inic$origin)){
    VTM_units[, , i] <- lphom_unit(lphom.object = lphom0, iii = i, solver = solver)
    votos_units[, , i] <- VTM_units.t[, , i]/rowSums(VTM_units.t[, , i])*lphom_inic$origin[i, ]
    VTM_units[lphom_inic$origin[i, ] == 0L, , i] <- 0L
  }
  votos_units[is.na(votos_units)] <- 0L
  VTM_votos_homogeneos <- HET_MT.votos_MT.prop_Y(votos_units)
  VTM_votos <- VTM_votos_homogeneos$MT.votos
  OTM <- round(t(VTM_votos)/colSums(VTM_votos)*100, 2)
  OTM <- OTM[c(1L:nrow(lphom_inic$OTM)), c(1:ncol(lphom_inic$OTM))]
  VTM.complete <- VTM_votos_homogeneos$MT.pro

  # Improving the solution
  dimnames(OTM) <- dimnames(lphom_inic$OTM)
  HETe <- HETe.sequence[iter + 1L]
  dimnames(VTM_units) <- dimnames(votos_units) <- c(dimnames(lphom_inic$VTM.complete),
                                                    list(rownames(lphom_inic$origin)))
  EHet <- VTM_votos_homogeneos$EHet
  dimnames(EHet) <- dimnames(lphom_inic$EHet)
  dimnames(VTM.complete) <- dimnames(VTM.complete.iter) <- dimnames(VTM_votos) <- dimnames(lphom_inic$VTM.complete)
  VTM <- round(VTM.complete[c(1L:nrow(lphom_inic$VTM)), c(1L:ncol(lphom_inic$VTM))]*100, 2)
  VTM.votes <- VTM_votos[c(1L:nrow(lphom_inic$VTM)), c(1L:ncol(lphom_inic$VTM))]

  if (verbose){
    cat("\n\nEstimated Heterogeneity Index HETe:",round(HETe, 2),"%\n")
    cat("\n\n Matrix of vote transitions (in %) from Election 1 to Election 2\n\n")
    print(VTM)
    cat("\n\n Origin (%) of the votes obtained in Election 2\n\n")
    print(OTM)
  }
  lphom_inic$inputs$verbose <- verbose
  inputs <- c(lphom_inic$inputs, "iter.max" = iter.max, "uniform" = uniform,
              "distance.local" = distance.local)
  inic <- lphom_inic[c(1L:6L, 10L)]
  names(inic) <- paste0(names(inic), "_init")
  
  # Caso de filas o columnas con cero votos
  filas0 <- which(rowSums(VTM_votos) == 0)
  colum0 <- which(colSums(VTM_votos) == 0)
  VTM[filas0, ] <- 0
  VTM.complete[filas0, ] <- 0
  OTM[colum0, ] <- 0
  
  output <- list("VTM" = VTM, "VTM.votes" = VTM.votes, "OTM" = OTM, "HETe" = HETe,
              "VTM.complete" = VTM.complete, "VTM.complete.votes" = VTM_votos,
              "VTM.prop.units" = VTM_units, "VTM.votes.units" = votos_units,
              "VTM.complete.last.iter" = VTM.complete.iter, "VTM.sequence" = VTM.sequence,
              "HETe.sequence" = HETe.sequence, "VTM.prop.units.last.iter" = VTM_units.iter,
              "VTM.votes.units.last.iter" = votos_units.iter, "zeros" = zeros,
              "iter" = iter, "iter.units" = iter.units, "errors" = errors, "EHet" = EHet,
              "inputs" = inputs, "origin" = lphom_inic$origin,
              "destination" = lphom_inic$destination, "solution_init" = inic)
  class(output) <- c("lclphom", "ei_lp", "lphom")
  return(output)
}
