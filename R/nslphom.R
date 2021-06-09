#' Implements nslphom algorithm
#'
#' @description  Estimates RxC vote transfer matrices (ecological contingency tables) with nslphom
#'
#' @author Jose M. Pavia, \email{pavia@@uv.es}
#' @author Rafael Romero \email{rromero@@eio.upv.es}
#' @references Pavia, JM, and Romero, R (2021). Improving estimates accuracy of voter transitions. Two new algorithms for ecological inference based on linear programming, mimeo.
#'
#' @param votes_election1 data.frame (or matrix) of order IxJ (likely of final order IxJ-1
#'                        in `regular` and `raw` scenarios) with the votes gained by the *J*
#'                        political options competing on election 1 (or origin) in the *I*
#'                        territorial units considered.
#'
#' @param votes_election2 data.frame (or matrix) of order IxK (likely of final order IxK-1
#'                        in `regular` and `raw` scenarios) with the votes gained by
#'                        the *K* political options competing on election 2 (or destination)
#'                        in the *I* territorial units considered.
#'
#' @param new_and_exit_voters A character string indicating the level of information available
#'                            regarding new entries and exits of the election censuses between the
#'                            two elections. This argument captures the different options discussed
#'                            on Section 3 of Romero et al. (2020). This argument admits five values:
#'                            `regular`, `raw`, `simultaneous`, `full` and `gold`. Default, `regular`.
#'
#' @param structural_zeros Default, NULL. A list of vectors of length two, indicating the election options
#'                         for which no transfer of votes are allowed between election 1 and election 2.
#'                         For instance, when new_and_exit_voters is set to `"regular"`,
#'                         lphom implicitly `states structural_zeros = list(c(J, K))` in case exits and/or
#'                         entries are computed because the sum by rows of `votes_election1` and
#'                         `votes_election2` does not coincide.
#'
#' @param iter.max Maximum number of iterations to be performed. The process ends when either the
#'                 number of iterations reaches iter.max or when the maximum variation between two consecutive
#'                 estimates of the probability transfer matrix is less than `tol`. By default, 10.
#'
#' @param min.first A TRUE/FALSE value. If FALSE, the matrix associated with the minimum `HETe` after
#'                  performing `iter.max` iterations is taken as solution.
#'                  If TRUE, the associated matrix to the instant in which the first decrease of `HETe` occurs
#'                  is taken as solution. The process stops at that moment. In this last scenario
#'                  (when `min.first = TRUE`), `burnin = 0` is forced and `iter.max` is at least 100. Default, FALSE.
#'
#' @param uniform A TRUE/FALSE value that indicates if census exits affects all the electoral options in a
#'                (relatively) similar fashion in each voting unit: equation (13) of Pavia and Romero (2021).
#'                Default, TRUE.
#'
#' @param distance.local A string argument that indicates whether the second step of the lphom_local algorithm
#'                       should be performed in order to resolve potential indeterminacies of local solutions.
#'                       Default, `"abs"`.
#'                       If `distance.local = "abs"` lphom_local selects in its the second step the matrix
#'                       closer to the temporary global solution under L_1 norm, among the first step compatible matrices.
#'                       If `distance.local = "abs"` lphom_local selects in its the second step the matrix
#'                       closer to the temporary global solution under L_Inf norm, among the first step compatible matrices.
#'                       If `distance.local = "none"`, the second step of lphom_local is not performed.
#'
#' @param burnin Number of initial solutions to be discarded before determining the final solution. By default, 0.
#'
#' @param verbose A TRUE/FALSE value that indicates if the main outputs of the function should be
#'                printed on the screen. Default, FALSE.
#'
#' @param tol Maximum deviation allowed between two consecutive iterations. The process ends when the maximum
#'            variation between two proportions for the estimation of the transfer matrix between two consecutive
#'            iterations is less than `tol` or the maximum number of iterations has been reached. By default, 0.00001.
#'
#' @details Description of the `new_and_exit_voters` argument in more detail.
#' \itemize{
#'  \item{`regular`: }{The default value. This argument accounts for the most plausible scenario.
#'                    A scenario with two elections elapsed at least some months.
#'                    In this scenario, (i) the column *J* of  `votes_election1` corresponds to
#'                    new young electors who have the right to vote for the first time and (ii)
#'                    net exits (basically a consequence of mortality), and eventually net entries,
#'                    are computed according equation (7) of Romero et al. (2020), and (iii) we
#'                    assume net exits affect equally all the first *J-1* options of election 1,
#'                     hence (8) and (9) constraints of Romero et al. (2020) are imposed.}
#'  \item{`raw`: }{This value accounts for a scenario with two elections where only the raw
#'                 election data recorded in the *I* territorial units, in which the area
#'                 under study is divided, are available. In this scenario, net exits
#'                 (basically deaths) and net entries (basically new young voters) are estimated
#'                 according to equation (7) of Romero et al. (2020). Constraints defined by
#'                 equations (8) and (9) of Romero et al. (2020) are imposed. In this scenario,
#'                 when net exits and/or net entries are negligible (such as between the first- and
#'                 second-round of French Presidential elections), they are omitted in the outputs.}
#'  \item{`simultaneous`: }{This value accounts for either a scenario with two simultaneous elections
#'                 or a classical ecological inference problem. In this scenario, the sum by rows
#'                 of `votes_election1` and `votes_election2` must coincide. Constraints
#'                 defined by equations (8) and (9) of Romero et al. (2020) are not included
#'                 in the model.}
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
#'  \item{VTM}{ A matrix of order JxK with the estimated percentages of vote transitions from election 1 to election 2.}
#'  \item{OTM}{ A matrix of order KxJ with the estimated percentages of the origin of the votes obtained for the different options of election 2.}
#'  \item{HETe}{ The estimated heterogeneity index as defined in equation (15) of Pavia and Romero (2021).}
#'  \item{VTM.complete}{ A matrix of order J'xK' with the estimated proportions of vote transitions from election 1 to election 2, including in `regular` and `raw` scenarios the row and the column corresponding to net_entries and net_exits even when they are really small, less than 1% in all units.}
#'  \item{VTM.sequence}{ Array of order J'xK'x(iter+1) (where `iter` is the efective number of iterations performed) of the estimated matrices corresponding to each iteration.}
#'  \item{HETe.sequence}{ Numeric vector of length `iter+1` with the `HETe` coefficients corresponding to the matrices in `VTM.sequence`.}
#'  \item{VTM_units}{ An array of order J'xK'xI with the estimated proportions of vote transitions from election 1 to election 2 attained for each unit in the selected iteration.}
#'  \item{VTM_votes}{ An array of order J'xK'xI with the estimated matrix of vote transitions from election 1 to election 2 attained for for each unit in the selected iteration.}
#'  \item{zeros}{ A list of vectors of length two, indicating the election options for which no transfer of votes are allowed between election 1 and election 2.}
#'  \item{iter}{ The real final number of iterations performed before ending the process.}
#'  \item{iter.min}{ Number of the iteration associated to the selected `VTM` solution.}
#'  \item{inputs}{ A list containing all the objects with the values used as arguments by the function.}
#'  \item{origin}{ A matrix with the final data used as votes of the origin election after taking into account the level of information available regarding to new entries and exits of the election censuses between the two elections.}
#'  \item{destination}{ A matrix with the final data used as votes of the origin election after taking into account the level of information available regarding to new entries and exits of the election censuses between the two elections.}
#'  \item{EHet}{ A matrix of order IxK measuring in each spatial unit a distance to the homogeneity hypothesis, that is, the differences under the homogeneity hypothesis between the actual recorded results and the expected results with the solution in each territorial unit for each option of election two.}
#'  \item{solution_init}{ A list with the main outputs produced by **lphom()**.}
#'  \itemize{
#'  \item{`VTM_init`:}{ A matrix of order JxK with the estimated percentages of vote transitions from election 1 to election 2 initially obtained by **lphom()**.}
#'  \item{`OTM_init`:}{ A matrix of order KxJ with the estimated percentages of the origin of the votes obtained for the different options of election 2 initially obtained by **lphom()**.}
#'  \item{`HETe_init`:}{ The estimated heterogeneity index defined in equation (10) of Romero et al. (2020). }
#'  \item{`EHet_init`:}{ A matrix of order IxK measuring in each spatial unit the distance to the homogeneity hypothesis, that is, the differences under the homogeneity hypothesis between the actual recorded results and the expected results, using the **lphom()** solution, in each territorial unit for each option of election two.}
#'  \item{`VTM.complete_init`:}{ A matrix of order J'xK' with the estimated proportions of vote transitions from election 1 to election 2 initially obtained by **lphom()**, including in `regular` and `raw` scenarios the row and the column corresponding to net_entries and net_exits even when they are really small, less than 1% in all units.}
#' }
#'
#' @export
#'
#'
#' @family linear programing ecological inference functions
#' @seealso \code{\link{lphom}} \code{\link{tslphom}}
#'
#' @examples
#' mt.ns <- nslphom(France2017P[, 1:8] , France2017P[, 9:12],
#'                  new_and_exit_voters= "raw")
#' mt.ns$VTM
#' mt.ns$HETe
#' mt.ns$solution_init$HETe_init
#
nslphom <- function(votes_election1, votes_election2,
                    new_and_exit_voters = c("regular", "raw", "simultaneous", "full", "gold"),
                    structural_zeros = NULL, iter.max = 10, min.first = FALSE,
                    uniform = TRUE, distance.local = c("abs", "max", "none"),
                    burnin = 0, verbose = FALSE, tol = 10^-5 ){

  # Calculo de la solucion inicial
  lphom_inic <- lphom(votes_election1 = votes_election1, votes_election2 = votes_election2,
                      new_and_exit_voters = new_and_exit_voters, structural_zeros, verbose = FALSE)
  lphom0 <- lphom_inic

  if (!(distance.local[1] %in% c("abs", "max", "none")))
    stop('Not allowed string for argument "distance.local".
         The only allowed strings for "distance.local" are "abs", "max" and "none".')

  # Funcion local a aplicar
  if (uniform) {
    if (distance.local[1] == "abs") {
      lphom_unit <- lphom_local_abs
    } else if (distance.local[1] == "max") {
      lphom_unit <- lphom_local_max
    } else if (distance.local[1] == "none") {
      lphom_unit <- lphom_local
    }
  } else {
    if (distance.local[1] == "abs") {
      lphom_unit <- calculo_MT_unidad_abs
    } else if (distance.local[1] == "max") {
      lphom_unit <- calculo_MT_unidad_max
    } else if (distance.local[1] == "none") {
      lphom_unit <- calculo_MT_unidad
    }
  }

  if (min.first){
    burnin <- 0
    iter.max <- max(iter.max, 100)
  }
  zeros <- determinar_zeros_estructurales(lphom_inic)

  # Inicio proceso iterativo
  VTM.sequence <- array(NA, c(dim(lphom_inic$VTM.complete), iter.max + 1))
  VTM_votos.sequence <- array(NA, c(dim(lphom_inic$VTM.complete), iter.max + 1))
  VTM_units.sequence <- array(NA, c(dim(lphom_inic$VTM.complete),
                                    nrow(lphom_inic$origin), iter.max + 1))
  votos_units.sequence <- array(NA, c(dim(lphom_inic$VTM.complete),
                                      nrow(lphom_inic$origin), iter.max + 1))
  EHet.sequence <- array(NA, c(dim(lphom_inic$EHet), iter.max + 1))
  iter <- 0
  dif.max <- Inf
  VTM.iter <- lphom0$VTM.complete <- lphom_inic$VTM.complete
  VTM.sequence[, , iter+1] <- VTM.iter
  HETe.sequence <- lphom_inic$HETe
  EHet.sequence[, , iter+1] <- lphom_inic$HETe

  while (iter < iter.max & dif.max > tol){
    # Calculo de las soluciones locales
    VTM_units <- votos_units <- array(NA, c(dim(lphom_inic$VTM.complete), nrow(lphom_inic$origin)))
    for (i in 1:nrow(lphom_inic$origin)){
      VTM_units[, , i] <- lphom_unit(lphom.object = lphom0, iii = i)
      votos_units[, , i] <- VTM_units[, , i]/rowSums(VTM_units[, , i])*lphom_inic$origin[i, ]
      VTM_units[lphom_inic$origin[i, ] == 0, , i] <- 0
    }
    VTM_votos_homogeneos <- HET_MT.votos_MT.prop_Y(votos_units)
    VTM_votos <- VTM_votos_homogeneos$MT.votos
    VTM.complete <- VTM_votos_homogeneos$MT.pro
    iter <- iter + 1
    dif.max <- max(abs(VTM.complete - VTM.iter))
    #      dif0 <- max(VTM.complete - abs(lphom_inic$VTM.complete))
    #      print(iter)
    #      print(dif0)
    #      print(dif.max)
    VTM.iter <- lphom0$VTM.complete <- VTM.complete
    VTM.sequence[, , iter+1] <- VTM.iter
    HETe.sequence <- c(HETe.sequence, VTM_votos_homogeneos$HET)
    VTM_votos.sequence[, , iter+1] <- VTM_votos
    VTM_units.sequence[, , , iter+1] <- VTM_units
    votos_units.sequence[, , , iter+1] <- votos_units
    EHet.sequence[, , iter+1] <- VTM_votos_homogeneos$EHet
    if (min.first & (HETe.sequence[iter+1] > HETe.sequence[iter])) dif.max <- -Inf
  } # End while
  dimnames(VTM.sequence) <- c(dimnames(lphom_inic$VTM.complete),
                              list(paste0("iter = ", 0:iter.max)))
  #  names(HETe.sequence) <- paste0("iter = ", 0:iter.max)

  # En caso de que converga antes de alcanzar el máximo de iteraciones
  VTM.sequence <- VTM.sequence[, , 1:(iter+1)]
  VTM_votos.sequence <- VTM_votos.sequence[, , 1:(iter+1)]
  VTM_units.sequence <- VTM_units.sequence[, , , 1:(iter+1)]
  votos_units.sequence <- votos_units.sequence[, , , 1:(iter+1)]

  # Solution
  if (iter < burnin) burnin <- iter - 1
  iter.select <- which.min(HETe.sequence[(burnin+2):(iter+1)])
  VTM.complete <- VTM.sequence[, , burnin+1+iter.select]
  VTM_votos <- VTM_votos.sequence[, , burnin+1+iter.select]
  VTM_units <- VTM_units.sequence[, , , burnin+1+iter.select]
  votos_units <- votos_units.sequence[, , , burnin+1+iter.select]
  OTM <- round(t(VTM_votos)/colSums(VTM_votos)*100, 2)
  OTM <- OTM[c(1:nrow(lphom_inic$OTM)), c(1:ncol(lphom_inic$OTM))]
  EHet <- EHet.sequence[, , burnin+1+iter.select]

  # Improving the solution
  dimnames(OTM) <- dimnames(lphom_inic$OTM)
  HETe <- HETe.sequence[burnin+1+iter.select]
  dimnames(VTM_units) <- dimnames(votos_units) <- c(dimnames(VTM.complete),
                                                    list(rownames(lphom_inic$origin)))
  dimnames(EHet) <- dimnames(lphom_inic$EHet)

  VTM <- round(VTM.complete[c(1:nrow(lphom_inic$VTM)), c(1:ncol(lphom_inic$VTM))]*100, 2)

  if (verbose){
    cat("\n\nEstimated Heterogeneity Index HETe:",round(HETe, 2),"%\n")
    cat("\n\n Matrix of vote transitions (in %) from Election 1 to Election 2\n\n")
    print(VTM)
    cat("\n\n Origin (%) of the votes obtained in Election 2\n\n")
    print(OTM)
  }
  inputs <- lphom_inic$inputs
  inputs$iter.max <- iter.max
  inputs$min.first <- min.first
  inputs$uniform <- uniform
  inputs$distance.local <- distance.local
  inputs$burnin <- burnin
  inputs$tol <- tol
  inic <- lphom_inic[c(1:4,8)]
  names(inic) <- paste0(names(inic), "_init")
  return(list("VTM" = VTM, "OTM" = OTM, "HETe" = HETe, "VTM.complete" = VTM.complete,
              "VTM.sequence" = VTM.sequence, "HETe.sequence" = HETe.sequence,
              "VTM_units" = VTM_units, "VTM_votes" = votos_units, "zeros" = zeros,
              "iter" = iter, "iter.min" = burnin + iter.select, "EHet" = EHet,
              "inputs" = inputs, "origin" = lphom_inic$origin,
              "destination" = lphom_inic$destination, "solution_init" = inic))
}
