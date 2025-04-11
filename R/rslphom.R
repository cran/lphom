#' Implements rslphom algorithm
#'
#' @description  Estimates RxC (JxK) vote transfer matrices (ecological contingency tables) with rslphom
#'
#' @author Jose M. Pavia, \email{pavia@@uv.es}
#' @references Pavia, JM, and Romero, R (2024). Improving estimates accuracy of voter transitions. Two new algorithms for ecological inference based on linear programming, *Sociological Methods & Research*, 53(4), 1491–1533. \doi{10.1177/00491241221092725}.
#' @references Pavia, JM, and Penades, A (2024). A bottom-up approach for ecological inference.
#'
#' @param votes_election1 data.frame (or matrix) of order IxJ1 with the votes gained by 
#'                        (or the counts corresponding to) the J1 political options competing
#'                        (available) on election 1 (or origin) in the I units considered.
#'                        In general, the row marginals of the I tables corresponding
#'                        to the units.
#'                        
#' @param votes_election2 data.frame (or matrix) of order IxK2
#'                        with the votes gained by (or the counts corresponding to) the K2
#'                        political options competing (available) on election 2 (or destination) 
#'                        in the I (territorial) units considered. In general, the column marginals 
#'                        of the I tables corresponding to the units.
#'
#' @param emphasis A numerical vector of values between 0 and 1 informing of the weights/emphasis to be used
#'                 to promote each unit when estimating its transfer matrix. Default, 0.995. When the length
#'                 of `emphasis` is one, only a weight (a level of emphasis) is analyzed. When the length of
#'                 `emphasis` is higher than one, as many as different weights/emphasis as the length of
#'                 emphasis are tried in the estimation of the transfer matrix of each unit. In each unit,
#'                 the local solution selected corresponds to the transfer matrix with lower expected error.
#'
#' @param new_and_exit_voters A character string indicating the level of information available
#'                            in `votes_election1` and `votes_election2` regarding new entries 
#'                            and exits of the election censuses between the two elections. 
#'                            This argument allows, in addition to the options discussed in Pavia 
#'                            (2023), three more options. This argument admits eleven different values: 
#'                            `raw`, `regular`, `ordinary`, `enriched`, `adjust1`, `adjust2`, 
#'                            `simultaneous`, `semifull`, `full`, `fullreverse` and `gold`. 
#'                            Default, `raw`.
#'                            
#' @param apriori data.frame (or matrix) of order J0xK0 with an initial estimate of the
#'                (row-standarized) global voter transition proportions/fractions, pjk0, between
#'                the first J0 (election) options of election 1 and the first K0 (election) options
#'                of election 2. This matrix can contain some missing values. When no a priori
#'                information is available `apriori` is a null object. Default, `NULL`.
#'
#' @param lambda A number between 0 and 1 informing the relative weight the user assigns to the
#'               `apriori` information. Setting `lambda = 0` is equivalent to not having a priori
#'               information (i.e., `apriori = NULL`). Default, `0.5`.
#'
#' @param uniform A `TRUE/FALSE` value that informs whether census exits impact all the electoral options
#'                in a (relatively) similar fashion in all iterations, including iteration 0 and 
#'                when deriving units tables. If `uniform = TRUE` typically at least one of the equations  
#'                among equations (6) to (11) of Pavia (2023) is included in the underlying model.
#'                This parameter has no effect in `simultaneous` scenarios. It also has not impact 
#'                in `raw` and `regular` scenarios when no net exits are estimated by the function
#'                from the provided information. Default, `TRUE`.
#'
#' @param structural_zeros Default `NULL`. A list of vectors of length two, indicating the election options
#'                         for which no transfer of votes are allowed between election 1 and election 2.
#'                         For instance, when new_and_exit_voters is set to `"semifull"`,
#'                         lphom implicitly states `structural_zeros = list(c(J1, K2))`.
#'
#' @param integers A `TRUE/FALSE` value that indicates whether the problem is solved in integer values
#'                 in all the steps, including lphom intermediate solutions and unit solutions.
#'                 If `integers = TRUE`, the LP matrices are approximated to the closest integer solution solving
#'                 the corresponding Integer Linear Program. Default, `FALSE`.
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
#' @param save.local.by.emphasis A `TRUE/FALSE` value that indicates if the estimated matrices obtained in each unit
#'                               for each value of emphasis should be saved. Default, `FALSE`.
#'
#' @param verbose A `TRUE/FALSE` value that indicates if a summary of the results of the computations performed
#'                to estimate net entries and exits should be printed on the screen. Default, `TRUE`.
#'
#' @param solver A character string indicating the linear programming solver to be used, only
#'               `lp_solve` and `symphony` are allowed. By default, `lp_solve`. The package `Rsymphony`
#'               needs to be installed for the option `symphony` to be used.
#'
#' @param integers.solver A character string indicating the linear programming solver to be used for
#'                        approximating the LP solution to the closest integer solution.
#'                        Only `symphony` and `lp_solve` are allowed. By default, `symphony`.
#'                        The package `Rsymphony` needs to be installed for the
#'                        option `symphony` to be used. Only used when `integers = TRUE`.
#'
#' @param ... Other arguments to be passed to the function. Not currently used.
#'

#'
#' @details Description of the `new_and_exit_voters` argument in more detail.
#' \itemize{
#'   \item `raw`: The default value. This argument accounts for the most plausible scenario when
#'                 estimating vote transfer matrices. A scenario with two elections elapsed at least
#'                 some months where only the raw election data recorded in the I (territorial) units,
#'                 in which the electoral space under study is divided, are available.
#'                 In this scenario, net exits and net entries are estimated according to
#'                 equation (7) of Romero et al. (2020). When both net entries and exits are no
#'                 null, constraint (15) of Pavia (2023) applies. If there are net exits and `uniform = TRUE`
#'                 either constraints (6) or (8) and (15) of Pavia (2023) are imposed. In this scenario,
#'                 J could be equal to J1 or J1 + 1 and K equal to K2 or K2 + 1.
#'  \item `regular`: This value accounts for a scenario with
#'                 two elections elapsed at least some months where (i) the column J1
#'                 of `votes_election1` corresponds to new young electors who have the right
#'                 to vote for the first time, (ii) net exits and maybe other additional
#'                 net entries are computed according to equation (7) of Romero et al. (2020), and
#'                 (iii) we can (or not) assume that net exits impact equally all the first J1 - 1
#'                 options of election 1. When both net entries and exits are no null, constraints
#'                 (13) and (15) of Pavia (2023) apply. If `uniform = TRUE` and there are net exits either
#'                 constraints (8) or (11) of Pavia (2023), depending on whether there are or not net
#'                 entries, are also imposed. In this scenario, J could be equal to J1 or J1 + 1 and
#'                 K equal to K2 or K2 + 1. Note that this scenario could be also used if
#'                 column J1 of `votes_election1` would correspond to immigrants instead of 
#'                 new young electors.
#'  \item `ordinary`: This value accounts for a scenario
#'                 with two elections elapsed at least some months where (i) the column K1
#'                 of `votes_election2` corresponds to electors who died in the period between
#'                 elections, (ii) net entries and maybe other additional net exits are
#'                 computed according to equation (7) of Romero et al. (2020), and (iii) we can
#'                 assume (or not) that exits impact equally all the J1 options of election 1.
#'                 When both net entries and exits are no null, constraints (14) and
#'                 (15) of Pavia (2023) apply and if `uniform = TRUE` either constraints
#'                 (8) and (9) or, without net entries, (6) and (7) of Pavia (2023) are also imposed.
#'                 In this scenario, J could be equal to J1 or J1 + 1 and K equal to K2 or K2 + 1.
#'                 Note that this scenario could be also used if column K1 of 
#'                 `votes_election2` would correspond to emigrants instead of deaths.
#'  \item `enriched`: This value accounts for a scenario that somehow combine `regular` and
#'                 `ordinary` scenarios. We consider two elections elapsed at least some months where
#'                 (i) the column J1 of `votes_election1` corresponds to new young electors
#'                  who have the right to vote for the first time, (ii) the column K2 of
#'                 `votes_election2` corresponds to electors who died in the interperiod
#'                 election, (iii) other (net) entries and (net) exits are computed according
#'                 to equation (7) of Romero et al. (2020), and (iv) we can assume
#'                 (or not) that exits impact equally all the J1 - 1 options of election 1.
#'                 When both net entries and exits are no null, constraints (12) to
#'                 (15) of Pavia (2023) apply and if `uniform = TRUE` constraints
#'                 (10) and (11) of Pavia (2023) are also imposed. In this scenario, J could be equal
#'                 to J1 or J1 + 1 and K equal to K2 or K2 + 1. Note that this scenario could be also used if
#'                 the column J1 of `votes_election1` would correspond to immigrants instead of 
#'                 new young electors and/or if column K1 of `votes_election2` would correspond
#'                 to emigrants instead of deaths.
#'  \item `adjust1`: This value accounts for a scenario
#'                 with two elections elapsed at least some months where the census in 
#'                 each of the I polling units of the first election (the row-sums of `votes_election1`) are 
#'                 proportionally adjusted to match the corresponding census of the polling units in the 
#'                 second election (the row-sums of `votes_election2`).
#'                 If `integers = TRUE`, each row in `votes_election1` is proportionally adjusted to the closest integer
#'                 vector whose sum is equal to the sum of the corresponding row in `votes_election2`.
#'  \item `adjust2`: This value accounts for a scenario
#'                 with two elections elapsed at least some months where the census in 
#'                 each of the I polling units of the second election (the row-sums of `votes_election2`) 
#'                 are proportionally adjusted to match the corresponding census of the polling units 
#'                 in the first election (the row-sums of `votes_election1`).
#'                 If `integers = TRUE`, each row in `votes_election2` is adjusted to the closest integer
#'                 vector whose sum is equal to the sum of the corresponding row in `votes_election1`.
#'  \item `simultaneous`: This is the value to be used in classical ecological inference problems,
#'                such as in ecological studies of racial voting, and in scenarios with two simultaneous elections.
#'                In this scenario, the sum by rows of `votes_election1` and `votes_election2` must coincide.
#'                Constraints defined by equations (8) and (9) of Romero et al. (2020) are not included in
#'                the model. In this case, the lphom function just implements the basic model defined,
#'                for instance, by equations (1) to (5) of Pavia (2024).
#'  \item `semifull`: This value accounts for a scenario with two elections elapsed at least some
#'                months, where: (i) the column J1 = J of `votes_election1` totals new
#'                electors (young and immigrants) that have the right to vote for the first time and
#'                (ii) the column K2 = K of `votes_election2` corresponds to total exits of the census
#'                lists (due to death or emigration). In this scenario, the sum by rows of
#'                `votes_election1` and `votes_election2` must agree and constraint (15)
#'                of Pavia (2023) apply. Additionally, if `uniform = TRUE` constraints
#'                (8) of Pavia (2023) are also imposed.
#'  \item `full`: This value accounts for a scenario with two elections elapsed at least some
#'                months, where (i) the column J - 1 of `votes_election1` totals new young
#'                electors that have the right to vote for the first time, (ii) the column J (=J1)
#'                of `votes_election1` measures new immigrants that have the right to vote and
#'                (iii) the column K (=K2) of `votes_election2` corresponds to total exits of the census
#'                lists (due to death or emigration). In this scenario, the sum by rows of
#'                `votes_election1` and `votes_election2` must agree and constraints (13)
#'                and (15) of Pavia (2023) apply.  Additionally, if `uniform = TRUE` constraints
#'                (11) of Pavia (2023) are also imposed.
#'  \item `fullreverse`: This value is somehow the mirror version of `full`. 
#'                It accounts for a scenario with two elections elapsed at least some
#'                months, where (i) the column J1 = J of `votes_election1` totals new
#'                electors (young and immigrants) that have the right to vote for the first time and
#'                (ii) where total exits are separated out between exits due to emigration 
#'                (column K - 1 of `votes_election2`) and death (column K of `votes_election2`). 
#'               In this scenario, the sum by rows of `votes_election1` and `votes_election2` must 
#'               agree and constraints (14) and (15) of Pavia (2023) apply. 
#'                Additionally, if `uniform = TRUE` constraints (8) and (9) of Pavia (2023) are also imposed.
#'  \item `gold`: This value accounts for a scenario similar to `full`, where total exits are
#'               separated out between exits due to emigration (column K - 1 of `votes_election2`)
#'               and death (column K of `votes_election2`). In this scenario, the sum by rows
#'               of `votes_election1` and `votes_election2` must agree. Constraints (12) to
#'               (15) of Pavia (2023) apply and if `uniform = TRUE` constraints (10) and (11)
#'               of Pavia (2023) are also imposed.
#' }
#'
#' @return
#' A list with the following components
#'  \item{VTM}{ A matrix of order J'xK' (where J'=J-1 or J and K'= K-1 or K) with the estimated percentages of row-standardized vote transitions from election 1 to election 2.
#'             In `raw`, `regular`, `ordinary` and `enriched` scenarios when the percentage of net entries is small, less than 1% of the census in all units,
#'             net entries are omitted (i.e., the number of rows of `VTM` is equal to J1) even when estimates for net entries different from zero are obtained. Likewise, in the same scenarios when the percentage of net exits is small, less than 1%
#'             of the census in all units, net exits are omitted (i.e., the number of rows of `VTM` is equal to K2) even when estimates for net exits different from zero are obtained.}
#'  \item{VTM.votes}{ A matrix of order J'xK' (where J'=J-1 or J and K'= K-1 or K) with the estimated vote transitions from election 1 to election 2.
#'             In `raw`, `regular`, `ordinary` and `enriched` scenarios when the percentage of net entries is small, less than 1% of the census,
#'             net entries are omitted (i.e., J = J1) even when estimates for net entries different from zero are obtained. Likewise, in the same scenarios when the percentage of net exits is small, less than 1%
#'             of the census, net exits are omitted (i.e., K = K2) even when estimates for net exits different from zero are obtained.}
#'  \item{OTM}{ A matrix of order KxJ with the estimated percentages of the origin of the votes obtained for the different options of election 2.}
#'  \item{HETe}{ The estimated heterogeneity index as defined in equation (15) of Pavia and Romero (2022).}
#'  \item{VTM.complete}{ A matrix of order JxK with the estimated proportions of row-standardized vote transitions from election 1 to election 2, including in `raw`, `regular`, `ordinary` and `enriched` scenarios the row and the column corresponding to net_entries and net_exits even when they are really small, less than 1% in all units.}
#'  \item{VTM.complete.votes}{ A matrix of order JxK with the estimated vote transitions from election 1 to election 2, including in `raw`, `regular`, `ordinary` and `enriched` scenarios the row and the column corresponding to net_entries and net_exits even when they are really small, less than 1% in all units.}
#'  \item{VTM.prop.units}{ An array of order JxKxI with the estimated proportions of vote transitions from election 1 to election 2 attained for each unit after adjusting the **lphom()** initial estimate.}
#'  \item{VTM.votes.units}{ An array of order JxKxI with the estimated matrix of vote transitions from election 1 to election 2 attained for each unit after adjusting the **lphom()** initial estimate.}
#'  \item{VTM.sequence}{ Array of order JxKxlength(emphasis) with the global estimated matrices corresponding to each weight.}
#'  \item{zeros}{ A list of vectors of length two, indicating the election options for which no transfer of votes are allowed between election 1 and election 2.}
#'  \item{errors}{ A matrix of order Ixlength(emphasis) with the expected errors for each unit and weight. The
#'                solution determined by `VTM.prop.units` or `VTM.votes.units` is the one obtained combining the unit solutions
#'                corresponding to the minimum observed errors.}
#'  \item{VTM.prop.units.by.emphasis}{ An array of order JxKxIxlength(emphasis) with the estimated proportions of vote transitions from election 1 to election 2 attained in each unit for each weight.
#'                                    This is a `NULL` array if `save.local.by.emphasis = FALSE`.}
#'  \item{deterministic.bounds}{ A list of two matrices of order JxK and two arrays of order JxKxI containing for each vote transition the lower and upper allowed proportions given the observed aggregates.}
#'  \item{inputs}{ A list containing all the objects with the values used as arguments by the function.}
#'  \item{origin}{ A matrix with the final data used as votes of the origin election after taking into account the level of information available regarding to new entries and exits of the election censuses between the two elections.}
#'  \item{destination}{ A matrix with the final data used as votes of the origin election after taking into account the level of information available regarding to new entries and exits of the election censuses between the two elections.}
#'  \item{EHet}{ A matrix of order IxK measuring in each spatial unit a distance to the homogeneity hypothesis, that is, the differences under the homogeneity hypothesis between the actual recorded results and the expected results with the solution in each territorial unit for each option of election 2.}
#'  \item{solution_init}{ A list with the main outputs produced by **lphom()**.}
#'  \itemize{
#'  \item `VTM_init`: A matrix of order J'xK' with the estimated percentages of vote transitions from election 1 to election 2 initially obtained by **lphom()** with the raw data, without promoting any unit.
#'  \item `VTM.votes_init`: A matrix of order J'xK' with the estimated vote transitions from election 1 to election 2 initially obtained by **lphom()**  with the raw data, without promoting any unit.
#'  \item `OTM_init`: A matrix of order KxJ with the estimated percentages of the origin of the votes obtained for the different options of election 2 initially obtained by **lphom()**  with the raw data, without promoting any unit.
#'  \item `HETe_init`: The estimated heterogeneity index defined in equation (10) of Romero et al. (2020).
#'  \item `EHet_init`: A matrix of order IxK measuring in each spatial unit the distance to the homogeneity hypothesis, that is, the differences under the homogeneity hypothesis between the actual recorded results and the expected results, using the **lphom()** solution  with the raw data, without promoting any unit, in each territorial unit for each option of election 2.
#'  \item `VTM.complete_init`: A matrix of order JxK with the estimated proportions of vote transitions from election 1 to election 2 initially obtained by **lphom()**, including in `raw`, `regular`, `ordinary` and `enriched` scenarios the row and the column corresponding to net_entries and net_exits even when they are really small, less than 1% in all units.
#'  \item `VTM.complete.votes_init`: A matrix of order JxK with the estimated vote transitions from election 1 to election 2 initially obtained by **lphom()**, including in `raw`, `regular`, `ordinary` and `enriched` scenarios the row and the column corresponding to net_entries and net_exits even when they are really small, less than 1% in all units.
#' }

#'
#' @export
#'
#'
#' @family linear programing ecological inference functions
#' @seealso \code{\link{lphom}} \code{\link{tslphom}} \code{\link{nslphom}} \code{\link{lclphom}}
#'
#' @examples
#' mt.rs <- rslphom(France2017P[, 1:8] , France2017P[, 9:12], emphasis = 0.5)
#' mt.rs$VTM
#'
#' @importFrom lpSolve lp
#'
#
rslphom <- function(votes_election1,
                    votes_election2,
                    emphasis = 0.995,
                    new_and_exit_voters = c("raw", "regular", "ordinary", "enriched", 
                                            "adjust1", "adjust2", "simultaneous", 
                                            "semifull", "full", "fullreverse", "gold"),
                    apriori = NULL,
                    lambda = 0.5,
                    uniform = TRUE,
                    structural_zeros = NULL,
                    integers = FALSE,
                    distance.local = c("abs", "max", "none"),
                    save.local.by.emphasis = FALSE,
                    verbose = TRUE,
                    solver = "lp_solve",
                    integers.solver = "symphony",
                    ...){
  
  argg <- c(as.list(environment()), list(...))
  integers <- test_integers(argg)
  weights <- emphasis
  
  if (!(distance.local[1L] %in% c("abs", "max", "none")))
    stop('Not allowed string for argument "distance.local".
         The only allowed strings for "distance.local" are "abs", "max" and "none".')
  
  if (max(weights) >= 1 | min(weights) < 0)
    stop('Invalid "emphasis" values. The values of "emphasis" must be smaller than one and non-negative.')
  
  if (integers.solver == "lp_solve"){
    dec2counts <- dec2counts_lp
  } else {
    dec2counts <- dec2counts_symphony
  }
  
  # Calculo de la solucion inicial
  lphom_inic <- lphom(votes_election1 = votes_election1, votes_election2 = votes_election2,
                      new_and_exit_voters = new_and_exit_voters,
                      apriori = apriori, lambda = lambda, uniform = uniform,
                      structural_zeros = structural_zeros, integers = integers,
                      verbose = verbose, solver = solver, integers.solver = integers.solver)
  
  # Funcion local a aplicar
  lphom_unit <- lp_solver_local(uniform = uniform,
                                distance.local = distance.local)
  
  # Ceros estructurales
  zeros <- determinar_zeros_estructurales(lphom_inic)
  
  # Objects for saving main outputs
  n.weights <- length(weights)
  VTM_units <- votos_units <- array(NA, c(dim(lphom_inic$VTM.complete), nrow(lphom_inic$origin)))
  VTM_units.sol <- votos_units.sol <- VTM_units
  errors <- matrix(NA, nrow(lphom_inic$origin), n.weights)
  VTM.sequence <- array(NA, c(dim(lphom_inic$VTM.complete), n.weights))
  VTM.prop.units.by.emphasis <- NULL
  if (save.local.by.emphasis)
    VTM.prop.units.by.emphasis <- array(NA, c(dim(lphom_inic$VTM.complete), nrow(lphom_inic$origin), n.weights))
  
  # Local solutions
  errors.min <- rep(Inf, nrow(lphom_inic$origin))
  for (w in 1L:n.weights){
    peso <- weights[w]
    for (i in 1L:nrow(lphom_inic$origin)){
      # Solution lphom tras reescalado de observaciones individuales
      lphom.temp <- rescaled(lphom.object = lphom_inic, w = peso, ii = i, ceros = zeros)
      VTM_units[, , i] <- lphom_unit(lphom.object = lphom.temp, iii = i, solver = solver)
      if (integers){
        votos_units[, , i] <-dec2counts(VTM_units[, , i]*lphom_inic$origin[i, ],
                                        lphom_inic$origin[i,], lphom_inic$destination[i,])
        VTM_units[, , i] <- votos_units[, , i]/rowSums(votos_units[, , i])
      } else {
        votos_units[, , i] <- VTM_units[, , i]/rowSums(VTM_units[, , i])*lphom_inic$origin[i, ]
      }
      VTM_units[lphom_inic$origin[i, ] == 0L, , i] <- 0L
      
      errors[i, w] <- sum(abs(lphom_inic$destination[i, ] -
                                lphom_inic$origin[i, ] %*% lphom.temp$VTM.complete))
      
      if(errors[i, w] < errors.min[i]){
        errors.min[i] <- errors[i, w]
        VTM_units.sol[, , i] <- VTM_units[, , i]
        votos_units.sol[, , i] <- votos_units[, , i]
      }
    } # End loop units
    if(save.local.by.emphasis){
      VTM.prop.units.by.emphasis[, , , w] <- VTM_units
    }
    VTM.sequence[, , w] <- HET_MT.votos_MT.prop_Y(votos_units)$MT.pro
  } # End loop weights
  
  votos_units.sol[is.na(votos_units.sol)] <- 0L
  VTM_votos_homogeneos <- HET_MT.votos_MT.prop_Y(votos_units.sol)
  VTM_votos <- VTM_votos_homogeneos$MT.votos
  VTM.complete <- VTM_votos_homogeneos$MT.pro
  EHet <- VTM_votos_homogeneos$EHet
  dimnames(VTM.complete) <- dimnames(lphom_inic$VTM.complete)
  dimnames(VTM_votos) <- dimnames(lphom_inic$VTM.complete)
  OTM <- round(t(VTM_votos)/colSums(VTM_votos)*100, 2)
  OTM <- OTM[c(1L:nrow(lphom_inic$OTM)), c(1L:ncol(lphom_inic$OTM))]
  dimnames(OTM) <- dimnames(lphom_inic$OTM)
  HETe <- VTM_votos_homogeneos$HET
  dimnames(VTM_units.sol) <- dimnames(votos_units.sol) <- c(dimnames(VTM.complete),
                                                            list(rownames(lphom_inic$origin)))
  dimnames(EHet) <- dimnames(lphom_inic$EHet)
  VTM <- round(VTM.complete[c(1L:nrow(lphom_inic$VTM)), c(1L:ncol(lphom_inic$VTM))]*100, 2)
  VTM.votes <- VTM_votos[c(1L:nrow(lphom_inic$VTM)), c(1L:ncol(lphom_inic$VTM))]
  
  dimnames(errors) <- c(list(rownames(lphom_inic$origin)),
                        list(paste("emphasis = ", weights)))
  
  dimnames(VTM.sequence) <- c(dimnames(VTM.complete),
                              list(paste("emphasis = ", weights)))
  
  inputs <- c(lphom_inic$inputs, "emphasis" = weights, "save.local.by.emphasis" = save.local.by.emphasis,
              "uniform" = uniform, "distance.local" = distance.local)
  
  inic <- lphom_inic[c(1L:6L, 11L)]
  names(inic) <- paste0(names(inic), "_init")
  
  # Caso de filas o columnas con cero votos
  filas0 <- which(rowSums(VTM_votos) == 0)
  colum0 <- which(colSums(VTM_votos) == 0)
  VTM[filas0, ] <- 0
  VTM.complete[filas0, ] <- 0
  OTM[colum0, ] <- 0
  det.bounds <- bounds_compound(origin = lphom_inic$origin,
                                destination = lphom_inic$destination, zeros = zeros)
  
  output <- list("VTM" = VTM, "VTM.votes" = VTM.votes, "OTM" = OTM, "HETe" = HETe,
                 "VTM.complete" = VTM.complete, "VTM.complete.votes" = VTM_votos,
                 "VTM.prop.units" = VTM_units, "VTM.votes.units" = votos_units,
                 "VTM.sequence" = VTM.sequence, "zeros" = zeros,
                 "VTM.prop.units.by.emphasis" = VTM.prop.units.by.emphasis, "errors" = errors,
                 "deterministic.bounds" = det.bounds, "inputs" = inputs, "origin" = lphom_inic$origin,
                 "destination" = lphom_inic$destination, "EHet" = EHet, "solution_init" = inic)
  class(output) <- c("rslphom", "ei_lp", "lphom")
  return(output)
}
