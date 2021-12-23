#' Implements tslphom algorithm
#'
#' @description  Estimates RxC vote transfer matrices (ecological contingency tables) with tslphom
#'
#' @author Jose M. Pavia, \email{pavia@@uv.es}
#' @author Rafael Romero \email{rromero@@eio.upv.es}
#' @references Pavia, JM, and Romero, R (2021). Improving estimates accuracy of voter transitions. Two new algorithms for ecological inference based on linear programming. \doi{10.31124/advance.14716638.v1}.
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
#' @param structural_zeros Default NULL. A list of vectors of length two, indicating the election options
#'                         for which no transfer of votes are allowed between election 1 and election 2.
#'                         For instance, when new_and_exit_voters is set to `"regular"`,
#'                         lphom implicitly `states structural_zeros = list(c(J, K))` in case exits and/or
#'                         entries are computed because the sum by rows of `votes_election1` and
#'                         `votes_election2` does not coincide.
#'
#' @param uniform A TRUE/FALSE value that indicates if census exits affects all the electoral options in a
#'                (relatively) similar fashion in each voting unit: equation (13) of Pavia and Romero (2021).
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
#' @param integers A TRUE/FALSE value that indicates whether the problem is solved in integer values
#'                 in both iterations, including iteration zero (lphom) and final (including unit) solutions. If TRUE, the LP matrices
#'                 are approximated to the closest integer solution solving the corresponding Integer Linear Program.
#'                 Default, FALSE.
#'
#' @param verbose A TRUE/FALSE value that indicates if the main outputs of the function should be
#'                printed on the screen. Default, FALSE.
#'
#' @param solver A character string indicating the linear programming solver to be used, only
#'               `lp_solve` and `symphony` are allowed. By default, `lp_solve`.
#'
#' @param ... Other arguments to be passed to the function. Not currently used.
#'  

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
#'  \item{HETe}{ The estimated heterogeneity index as defined in equation (15) of Pavia and Romero (2021).}
#'  \item{VTM.complete}{ A matrix of order J'xK' with the estimated proportions of row-standardized vote transitions from election 1 to election 2, including in `regular` and `raw` scenarios the row and the column corresponding to net_entries and net_exits even when they are really small, less than 1% in all units.}
#'  \item{VTM.complete.votes}{ A matrix of order J'xK' with the estimated vote transitions from election 1 to election 2, including in `regular` and `raw` scenarios the row and the column corresponding to net_entries and net_exits even when they are really small, less than 1% in all units.}
#'  \item{VTM.prop.units}{ An array of order J'xK'xI with the estimated proportions of vote transitions from election 1 to election 2 attained for each unit after adjusting the **lphom()** initial estimate.}
#'  \item{VTM.votes.units}{ An array of order J'xK'xI with the estimated matrix of vote transitions from election 1 to election 2 attained for each unit after adjusting the **lphom()** initial estimate.}
#'  \item{zeros}{ A list of vectors of length two, indicating the election options for which no transfer of votes are allowed between election 1 and election 2.}
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
#'
#' @family linear programing ecological inference functions
#' @seealso \code{\link{lphom}} \code{\link{nslphom}} \code{\link{lclphom}}
#'
#' @examples
#' mt.ts <- tslphom(France2017P[, 1:8] , France2017P[, 9:12], new_and_exit_voters= "raw")
#' mt.ts$VTM
#' mt.ts$HETe
#' mt.ts$solution_init$HETe_init
#'
#' @importFrom Rsymphony Rsymphony_solve_LP
#' @importFrom lpSolve lp
#'

#
tslphom <- function(votes_election1,
                    votes_election2,
                    new_and_exit_voters = c("raw", "regular", "simultaneous", "full", "gold"),
                    structural_zeros = NULL,
                    uniform = TRUE,
                    distance.local = c("abs", "max", "none"),
                    integers = FALSE,
                    verbose = FALSE,
                    solver = "lp_solve",
                    ...){

  argg <- c(as.list(environment()), list(...))
  integers <- test_integers(argg)
  
  if (!(distance.local[1L] %in% c("abs", "max", "none")))
    stop('Not allowed string for argument "distance.local".
         The only allowed strings for "distance.local" are "abs", "max" and "none".')

  # Calculo de la solucion inicial
  lphom_inic <- lphom(votes_election1 = votes_election1, votes_election2 = votes_election2,
                      new_and_exit_voters = new_and_exit_voters, structural_zeros,
                      integers = integers, verbose = FALSE, solver = solver)

  # Funcion local a aplicar
  lphom_unit <- lp_solver_local(uniform = uniform,
                                distance.local = distance.local)

  # Calculo de las soluciones locales
  zeros <- determinar_zeros_estructurales(lphom_inic)
  VTM_units <- votos_units <- array(NA, c(dim(lphom_inic$VTM.complete), nrow(lphom_inic$origin)))

  for (i in 1L:nrow(lphom_inic$origin)){
      VTM_units[, , i] <- lphom_unit(lphom.object = lphom_inic, iii = i, solver = solver)
      if (integers){
        votos_units[, , i] <-dec2counts(VTM_units[, , i]*lphom_inic$origin[i, ],
                                        lphom_inic$origin[i,], lphom_inic$destination[i,])
        VTM_units[, , i] <- votos_units[, , i]/rowSums(votos_units[, , i])
      } else {
        votos_units[, , i] <- VTM_units[, , i]/rowSums(VTM_units[, , i])*lphom_inic$origin[i, ]
      }
      VTM_units[lphom_inic$origin[i, ] == 0L, , i] <- 0L
  }
  votos_units[is.na(votos_units)] <- 0L
  VTM_votos_homogeneos <- HET_MT.votos_MT.prop_Y(votos_units)
  VTM_votos <- VTM_votos_homogeneos$MT.votos
  VTM.complete <- VTM_votos_homogeneos$MT.pro
  EHet <- VTM_votos_homogeneos$EHet
  dimnames(VTM.complete) <- dimnames(lphom_inic$VTM.complete)
  dimnames(VTM_votos) <- dimnames(lphom_inic$VTM.complete)
  OTM <- round(t(VTM_votos)/colSums(VTM_votos)*100, 2)
  OTM <- OTM[c(1L:nrow(lphom_inic$OTM)), c(1L:ncol(lphom_inic$OTM))]
  dimnames(OTM) <- dimnames(lphom_inic$OTM)
  HETe <- VTM_votos_homogeneos$HET
  dimnames(VTM_units) <- dimnames(votos_units) <- c(dimnames(VTM.complete),
                                                    list(rownames(lphom_inic$origin)))
  dimnames(EHet) <- dimnames(lphom_inic$EHet)
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
  inputs <- c(lphom_inic$inputs, "uniform" = uniform, "distance.local" = distance.local)
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
                 "VTM.prop.units" = VTM_units, "VTM.votes.units" = votos_units, "zeros" = zeros,
                 "inputs" = inputs, "origin" = lphom_inic$origin,
                 "destination" = lphom_inic$destination, "EHet" = EHet, "solution_init" = inic)
  class(output) <- c("tslphom", "ei_lp", "lphom")
  return(output)
}
