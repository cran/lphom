#' Implements lphom algorithm
#'
#' @description  Estimates RxC (JxK) vote transfer matrices (ecological contingency tables) with lphom
#'
#' @author Jose M. Pavia, \email{pavia@@uv.es}
#' @author Rafael Romero \email{rromero@@eio.upv.es}
#' @references Romero, R, Pavia, JM, Martin, J and Romero G (2020). Assessing uncertainty of voter transitions estimated from aggregated data. Application to the 2017 French presidential election. *Journal of Applied Statistics*, 47(13-15), 2711-2736. \doi{10.1080/02664763.2020.1804842}
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
#' @param lambda A number between 0 and 1, informing the relative weight the user assigns to the 
#'               `apriori` information. Setting `lambda = 0` is equivalent to not having a priori
#'               information (i.e., `apriori = NULL`). Default, `0.5`.
#' 
#' @param uniform A `TRUE/FALSE` value that informs whether census exits affect all the electoral options
#'                in a (relatively) similar fashion. If `uniform = TRUE` typically at least one of the equations  
#'                among equations (6) to (11) of Pavia (2022) is included in the underlying model.
#'                This parameter has never effect in `simultaneous` scenarios. It also has not impact 
#'                in `raw` and `regular` scenarios when no net exits are estimated by the function
#'                from the provided information. Default, `TRUE`.
#'
#' @param structural_zeros Default `NULL`. A list of vectors of length two, indicating the election options
#'                         for which no transfer of votes are allowed between election 1 and election 2.
#'                         For instance, when new_and_exit_voters is set to `"semifull"`,
#'                         lphom implicitly states `structural_zeros = list(c(J1, K2))`.
#'
#' @param integers A `TRUE/FALSE` value that indicates whether the LP solution of counts (votes) must be 
#'                 approximate to the closest integer solution using ILP to generate the final solution.
#'                 Default, `FALSE`.
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
#'  \item{VTM}{ A matrix of order J'xK' (where J'=J-1 or J and K'=K-1 or K) with the estimated percentages of row-standardized vote transitions from election 1 to election 2. 
#'             In `raw`, `regular`, `ordinary` and `enriched` scenarios when the percentage of net entries is small, less than 1% of the census in all units, 
#'             net entries are omitted (i.e., the number of rows of `VTM` is equal to J1) even when estimates for net entries different from zero are obtained. Likewise, in the same scenarios when the percentage of net exits is small, less than 1%
#'             of the census in all units, net exits are omitted (i.e., the number of rows of `VTM` is equal to K2) even when estimates for net exits different from zero are obtained.}
#'  \item{VTM.votes}{ A matrix of order J'xK' (where J'=J-1 or J and K'=K-1 or K) with the estimated vote transitions from election 1 to election 2.
#'             In `raw`, `regular`, `ordinary` and `enriched` scenarios when the percentage of net entries is small, less than 1% of the census, 
#'             net entries are omitted (i.e., J = J1) even when estimates for net entries different from zero are obtained. Likewise, in the same scenarios when the percentage of net exits is small, less than 1%
#'             of the census, net exits are omitted (i.e., K = K2) even when estimates for net exits different from zero are obtained.}
#'  \item{OTM}{ A matrix of order KxJ with the estimated percentages of the origin of the votes obtained for the different options of election 2.}
#'  \item{HETe}{ The estimated heterogeneity index defined in equation (11) of Romero et al. (2020).}
#'  \item{VTM.complete}{ A matrix of order JxK with the estimated proportions of row-standardized vote transitions from election 1 to election 2.
#'                      In `raw`, `regular`, `ordinary` and `enriched` scenarios, this matrix includes the row and the column corresponding to net entries 
#'                      and net exits (when they are present) even when they are really small, less than 1%.}
#'  \item{VTM.complete.votes}{ A matrix of order JxK with the estimated vote transitions from election 1 to election 2.
#'                      In `raw`, `regular`, `ordinary` and `enriched` scenarios, this matrix includes the row and the column corresponding to net entries 
#'                      and net exits (when they are present) even when they are really small, less than 1%.}
#'  \item{deterministic.bounds}{ A list of two matrices of order JxK containing for each vote transition the lower and upper proportions allowed given the observed aggregates.}
#'  \item{inputs}{ A list containing all the objects with the values used as arguments by the function.}
#'  \item{origin}{ A matrix with the final data used as votes of the origin election after taking into account the level of information available regarding to new entries and exits of the election censuses between the two elections.}
#'  \item{destination}{ A matrix with the final data used as votes of the origin election after taking into account the level of information available regarding to new entries and exits of the election censuses between the two elections.}
#'  \item{EHet}{ A matrix of order IxK measuring in each spatial unit a distance to the homogeneity hypothesis. That is, the differences under the homogeneity hypothesis between the actual recorded results and the expected results in each territorial unit for each option of election 2.}
#' @export
#'
#'
#' @family linear programing ecological inference functions
#' @seealso \code{\link{tslphom}} \code{\link{nslphom}} \code{\link{lclphom}} \code{\link{rslphom}}
#'
#' @examples
#' lphom(France2017P[, 1:8] , France2017P[, 9:12], new_and_exit_voters= "raw")
#
#' @importFrom lpSolve lp
#
lphom <- function(votes_election1,
                  votes_election2,
                  new_and_exit_voters = c("raw", "regular", "ordinary", "enriched", 
                                          "adjust1", "adjust2", "simultaneous", 
                                          "semifull", "full", "fullreverse", "gold"),
                  apriori = NULL,
                  lambda = 0.5,
                  uniform = TRUE,
                  structural_zeros = NULL,
                  integers = FALSE,
                  verbose = TRUE,
                  solver = "lp_solve",
                  integers.solver = "symphony",
                  ...){
  
  argg <- c(as.list(environment()), list(...))
  integers <- test_integers(argg)
  matrix.votes <- tests_inputs_lphom(argg)
  
  # inputs
  inputs <- list("votes_election1" = votes_election1, "votes_election2" = votes_election2,
                 "new_and_exit_voters" = new_and_exit_voters[1], "apriori" = apriori,
                 "lambda" = lambda, "uniform" = uniform, "structural_zeros" = structural_zeros,
                 "integers" = integers, "verbose" = verbose, "solver" = solver,
                 "integers.solver" = integers.solver)
  
  x0 <- matrix.votes$x
  y0 <- matrix.votes$y
  
  if (new_and_exit_voters[1] == "adjust1"){
    new_and_exit_voters[1] <- "simultaneous"
    x0 <- adjust_xy(x0, y0, integers)
  }

  if (new_and_exit_voters[1] == "adjust2"){
    new_and_exit_voters[1] <- "simultaneous"
    y0 <- adjust_xy(y0, x0, integers)
  }
    
  new_and_exit_voters <- scenario <- new_and_exit_voters[1]
  
  if (integers.solver == "lp_solve"){
    dec2counts <- dec2counts_lp
  } else {
    dec2counts <- dec2counts_symphony
  }
  
  # Data preparation
  net <- compute_net_voters(x0 = x0, y0 = y0)
  x <- net$x
  y <- net$y
  apriori <- completar_apriori(net = net, apriori = apriori) 
  
  # Parameters
  I <- nrow(x); J <- ncol(x); K <- ncol(y); J0 <- ncol(x0); K0 <- ncol(y0)
  JK <- J*K; IK <- I*K
  
  # Names of election options
  names1 <- colnames(x)
  names2 <- colnames(y)
  
  # Net entris and exits
  NET_ENTRIES <- NET_EXITS <- 0L
  if(J > J0) NET_ENTRIES <- x[, ncol(x)]  
  if(K > K0) NET_EXITS <- y[, ncol(y)]  
  d <- rowSums(y0) - rowSums(x0)
  net_entries <- sum(NET_ENTRIES)/sum(x)
  net_exits <- sum(NET_EXITS)/sum(x)
  
  if(verbose & (new_and_exit_voters %in% c("regular", "raw", "ordinary", "enriched") & 
                any(d != 0L))){
    message(paste0('*********************WARNING*********************\n',
                   'You are in a \"', new_and_exit_voters, '\" scenario.\n',
                   'The sums (by row) of origin and destination data differ. ',
                   'It is, for at least a unit, the total number ',
                   'of electors in both elections is not the same. \n',
                   'To guarantee the matching: A new category of census entries ',
                   '(NET_ENTRIES) has been included in the origin election and ',
                   'a new category of census exits (NET_EXITS) ',
                   'has been also included in the destination election.',
                   '\n Their aggregate importances, measured in percentage of the total census, are:\n ',
                   '    %NET_ENTRIES = ',round(100*net_entries, 4),'%\n ', 
                   '    %NET_EXITS = ',round(100*net_exits, 4),'%\n',
                   'If NET_ENTRIES and/or NET_EXITS are really small, less than 1% ',
                   'in all units, their corresponding results will not be displayed in the main ',
                   'output, VTM. They are anyway always included in VTM.complete. \n',
                   '*************************************************'))
  }
  
  # System of equations
  if (scenario == "simultaneous" |
      (scenario == "raw" & J == J0 & K == K0) |
      (scenario == "raw" & J > J0 & K == K0) |
      (scenario == "raw" & J == J0 & K > K0 & uniform == FALSE) |
      (scenario == "regular" & J == J0 & K == K0) |
      (scenario == "regular" & J > J0 & K == K0) |
      (scenario == "ordinary" & J == J0 & K == K0 & uniform == FALSE) |
      (scenario == "ordinary" & J == J0 & K > K0 & uniform == FALSE)) {
    sistema <- model_lphom_apriori_1_2(x, y, P0 = apriori, lambda = lambda, uniform = FALSE)
  }

  if ((scenario == "raw" & J == J0 & K > K0 & uniform == TRUE) |
      (scenario == "ordinary" & J == J0 & K == K0 & uniform == TRUE)) {
    sistema <- model_lphom_apriori_1_2(x, y, P0 = apriori, lambda = lambda, uniform = TRUE)
  }
  
  if (scenario == "ordinary" & J == J0 & K > K0  & uniform == TRUE) {
    sistema <- model_lphom_apriori_7(x, y, P0 = apriori, lambda = lambda, uniform = uniform)
  }
  
  if ((scenario == "raw" & J > J0 & K > K0) |
      (scenario == "regular" & J == J0 & K > K0) |
      (scenario == "ordinary" & J > J0 & K == K0) |
      (scenario == "enriched" & J == J0 & K == K0) |
      scenario == "semifull") {
    sistema <- model_lphom_apriori_3_4(x, y, P0 = apriori, lambda = lambda, uniform = uniform)
  }
  
  if ((scenario == "regular" & J > J0 & K > K0) |
      (scenario == "enriched" & J > J0 & K == K0) |
      scenario == "full") {
    sistema <- model_lphom_apriori_5_6(x, y, P0 = apriori, lambda = lambda, uniform = uniform)
  }
  
  if ((scenario == "ordinary" & J > J0 & K > K0) |
      (scenario == "enriched" & J == J0 & K > K0) |
      scenario == "fullreverse") {
    sistema <- model_lphom_apriori_8_9(x, y, P0 = apriori, lambda = lambda, uniform = uniform)
  }
  
  if ((scenario == "enriched" & J > J0 & K > K0) |
      scenario == "gold") {
    sistema <- model_lphom_apriori_10_11(x, y, P0 = apriori, lambda = lambda, uniform = uniform)
  }
  
  # Structural zero restrictions introduced by the user
  if (!is.null(structural_zeros)){
    sistema <- add_structural_zeros(sistema = sistema, 
                                    structural_zeros = structural_zeros, 
                                    K = K)
  }
  
  # Solution
  if (solver == "symphony"){
    sol <- Rsymphony::Rsymphony_solve_LP(obj = sistema$f,
                                         mat = sistema$A,
                                         dir = rep('==', length(sistema$b)),
                                         rhs = sistema$b)
  } else {
    sol <- suppressWarnings(lpSolve::lp(direction = 'min', 
                                        objective.in = sistema$f,
                                        const.mat = sistema$A,
                                        const.dir = rep('=', length(sistema$b)),
                                        const.rhs = sistema$b) )
  }
  
  z <- sol$solution
  pjk <- matrix(z[1:JK], J, K, TRUE, dimnames = list(names1, names2))
  eik <- y - x %*% pjk
  
  # Integer solution
  if (integers){
    vjk <- pjk*colSums(x)
    vjk <- dec2counts(vjk, colSums(x), colSums(y))
    pjk <- vjk/rowSums(vjk)
    dimnames(pjk) <- list(names1, names2)
    eik <- y - x %*% pjk
  }
  
  # Outputs
  colnames(eik) <- names2
  rownames(eik) <- rownames(x)
  EHet <- eik
  vjk <- pjk*colSums(x)
  vjk.complete <- vjk
  pkj <- t(vjk)/colSums(vjk)
  
  # Caso de filas o columnas con cero votos
  filas0 <- which(rowSums(vjk) == 0L)
  colum0 <- which(colSums(vjk) == 0L)
  pjk[filas0, ] <- 0L
  pkj[colum0, ] <- 0L
  
  pjk.complete <- pjk
  HIe <- 100*sum(abs(eik))/sum(vjk.complete)
  
  # Reduced solution
  if (net_entries > 0 & max(x[,ncol(x)]/rowSums(x)) < 0.01){
    pjk <- pjk[-J, ]; vjk <- vjk[-J, ]
  }
  if (net_exits > 0 & max(y[,ncol(y)]/rowSums(y)) < 0.01){
    pjk = pjk[, -K]; vjk <- vjk[, -K]
  }
  
  pjk <- round(100*pjk, 2)
  pkj <- round(100*pkj, 2)
  
  det.bounds <- bounds_compound(origin = x, destination = y, zeros = structural_zeros)[c(1,2)]
  
  output <- list("VTM" = pjk, "VTM.votes" = vjk, "OTM" = pkj, "HETe" = HIe,
                 "VTM.complete" = pjk.complete, "VTM.complete.votes" = vjk.complete,
                 "deterministic.bounds" = det.bounds, "inputs" = inputs, "origin" = x, 
                 "destination" = y, "EHet" = EHet)
  class(output) <- c("lphom", "ei_lp")
  return(output)
}





## Auxiliary functions

# Function that tests if all inputs of lphom (version apriori) are correct and converts
# votes_election1 and votes_election2 into matrix objects 
tests_inputs_lphom <- function(argg){
  
  argg.1 <- as.matrix(argg$votes_election1)
  argg.2 <- as.matrix(argg$votes_election2)
  
  # Test integers (incluido en la función test_integers()) 
  # condicion <- max(abs(argg.1 - round(argg.1))) + max(abs(argg.2 - round(argg.2)))
  # if(condicion > 0L & argg$integers)
  #  stop("Integer solutions cannot be computed. At least a value included in 'votes_election1' or 'votes_election2' is decimal.")
  
  # Test names 1
  if (!(argg$new_and_exit_voters[1L] %in% c("raw", "regular", "ordinary", "semifull", "enriched", "adjust1",
                                            "adjust2", "simultaneous", "full", "fullreverse", "gold")))
    stop('The value set for argument "new_and_exit_voters" is not allowed. The only allowed strings for "new_and_exit_voters" are "raw", "simultaneous", "regular", "ordinary", "enriched", "semifull", "full" and "gold".')
  
  # Test names 2
  if (!(argg$solver[1L] %in% c("lp_solve", "symphony")))
    stop('The value set for argument "solver" is not allowed. The only allowed strings for "solver" are "lp_solve" and "symphony".')
  
  # Test names 3
  if (!(argg$integers.solver[1L] %in% c("lp_solve", "symphony")))
    stop('The value set for argument "integers.solver" is not allowed. The only allowed strings for "integers.solver" are "lp_solve" and "symphony".')
  
  # Tests numeric data
  x <- as.matrix(argg$votes_election1)
  y <- as.matrix(argg$votes_election2)
  if (ncol(x) == 1L) x <- t(x)
  if (ncol(y) == 1L) y <- t(y)
  if (nrow(x) != nrow(y))
    stop('The number of spatial units is different in origin and destination.')
  if (argg$new_and_exit_voters[1L] %in% c("semifull", "simultaneous", "full", "fullreverse", "gold")){
    if (!identical(round(rowSums(x)), round(rowSums(y)))){
      texto <- paste0('The number of voters (electors) in Election 1 and ',
                      'Election 2 differ (in at least a territorial unit). This is not ',
                      'allowed in a \"', argg$new_and_exit_voters[1L], '\" scenario.')
      stop(texto)
    } 
  }
  if (min(x,y) < 0L) stop('Negative values are not allowed in arguments "votes_election1" and "votes_election2".')
  
  # Test lambda
  if(!is.null(argg$apriori)){
    if (min(argg$apriori, na.rm = T) < 0L) 
      stop('Negative values for a priori proportions are not allowed')
    if(argg$lambda > 1L | argg$lambda < 0L)
      stop('Only values in the interval [0, 1] are allowed for the "lambda" argument.')
  }
  
  return(list("x" = x, "y" = y))
}


# Function to compute net entries and exits and to generate the origin and destination matrices,
# including if it is the case net entries and exits
compute_net_voters <- function(x0, y0){
  NET_ENTRIES <- NET_EXITS <- rep(0L, nrow(x0))
  x <- x0
  y <- y0
  
  # Estimation of net entries/exits in the election census
  d <- rowSums(y) - rowSums(x)
  if (any(d != 0L)) {
    NET_ENTRIES[d > 0L] <- d[d > 0L]
    NET_EXITS[d < 0L] <- -d[d < 0L]
  }
  
  # Net entries and exits
  if (sum(NET_ENTRIES) > 0L){
    x <- cbind(x, NET_ENTRIES)
  }
  if (sum(NET_EXITS) > 0L){
    y <- cbind(y0, NET_EXITS)
  }
  
  return(list("x" = x, "y" = y))
}


# Funcion para completar la matriz apriori con NAs, en caso de que su tamanyo sea inferior al finalmente
# necesario despues de calcular el numero total de opciones de origen y de destino.
# Caso de que sea mayor al requerido tambien se soluciona.
completar_apriori <- function(net, apriori){
  # net is an object output of the function compute_net_voters()
  J0 <- nrow(apriori)
  K0 <- ncol(apriori)
  J <- ncol(net$x)
  K <- ncol(net$y)
  if (is.null(apriori)){
    apriori <- matrix(NA, J, K)
  } else{  
    if (J > J0) apriori <- rbind(apriori, matrix(NA, J - J0, ncol(apriori)))
    if (K > K0) apriori <- cbind(apriori, matrix(NA, nrow(apriori), K - K0))
  }
  apriori <- apriori[1L:J, 1L:K]
  return(apriori)
}


# Function para calcular todos los componentes necesarios del programa lineal para
# calcular matrices de transferencia de type I and type II structures, ver Pavia (2022).
model_lphom_apriori_1_2 <- function(X, Y, P0, lambda, uniform, ...){
  
  J <- ncol(X)
  K <- ncol(Y)
  I <- nrow(X)
  JK <- J*K
  IK <- I*K
  
  xt <- colSums(X)
  yt <- colSums(Y)
  
  # Constraints suma por filas unitaria de la matriz pjk
  A2 <- cbind(kronecker(diag(J), t(rep(1L, K))), matrix(0L, J, 2L*JK), matrix(0L, J, 2L*IK))
  # Coherencia de matriz global al aplicar la matriz pjk por filas obtener las columnas
  A3 <- cbind(kronecker(t(xt), diag(K)),  matrix(0L, K, 2L*JK), matrix(0L, K, 2L*IK))
  # Constraints of (aprox) matching with a priori results.
  #  A5 <- cbind(kronecker(diag(xt), diag(K)), matrix(0L, JK, 2L*IK), 
  #              kronecker(t(c(-1L, 1L)), diag(JK)))
  A4 <- cbind(diag(JK), kronecker(t(c(-1L, 1L)), diag(JK)), matrix(0L, JK, 2L*IK))
  # Constraints for uniformity in last column, complete
  A5 <- cbind(kronecker(cbind(rep(1L, J - 1L), -diag(J - 1L)), t(c(rep(0L, K - 1L), 1L))), 
              matrix(0L, J - 1L, 2*JK), matrix(0L, J - 1L, 2L*IK))
  # Constraints of (aprox) matching in the I spatial units 
  A9 <- cbind(kronecker(X, diag(K)), matrix(0L, IK, 2L*JK), t(kronecker(diag(IK), c(1L,-1L))))
  
  
  # Terminos independientes
  b2 <- c(rep(1L, J))
  b3 <- yt
  #  b5 <- as.vector(t(xt * P0))
  b4 <- as.vector(t(P0))
  b5 <- rep(0L, J - 1L)
  b9 <- as.vector(t(Y))
  
  fp <- rep(0L, JK) # pjk
  fe <- rep(1L - lambda, 2L*IK) # eik
  #  fs <- rep(lambda, 2L*JK) # epsilon_jk
  fs <- lambda*as.vector(t(matrix(rep(xt, length(yt)), length(xt), length(yt)))) # epsilon_jk  
  
  # Dealing with missing values in P0
  falta <- which(is.na(b4))
  A4[falta, ] <- 0L
  b4[falta] <- 0L
  fs[falta] <- 0L
  
  if (uniform){
    A <- rbind(A2, A3, A4, A5, A9)
    b <- c(b2, b3, b4, b5, b9)
    f <- c(fp, fs, fs, fe)
  } else {
    A <- rbind(A2, A3, A4, A9)
    b <- c(b2, b3, b4, b9)
    f <- c(fp, fs, fs, fe)
  }
  
  if (lambda == 1){
    b <- b[1L:(nrow(A) - IK)]
    if(min(xt[xt != 0]) > 1L) fs <- fs/min(xt[xt != 0])
    f <- c(fp, fs, fs)
    A <- A[1L:(nrow(A) - IK), 1L:(ncol(A) - 2L*IK)]
  }
  
  output <- list("A" = A, "b" = b, "f" = f)
  return(output)
}


# Function para calcular todos los componentes necesarios del programa lineal para
# calcular matrices de transferencia de type I and type II structures, ver Pavia (2022).
model_lphom_apriori_3_4 <- function(X, Y, P0, lambda, uniform, ...){
  
  J <- ncol(X)
  K <- ncol(Y)
  I <- nrow(X)
  JK <- J*K
  IK <- I*K
  
  xt <- colSums(X)
  yt <- colSums(Y)
  
  # Constraints suma por filas unitaria de la matriz pjk
  A2 <- cbind(kronecker(diag(J), t(rep(1L, K))), matrix(0L, J, 2L*JK), matrix(0L, J, 2L*IK))
  # Coherencia de matriz global al aplicar la matriz pjk por filas obtener las columnas
  A3 <- cbind(kronecker(t(xt), diag(K)),  matrix(0L, K, 2L*JK), matrix(0L, K, 2L*IK))
  
  # Constraints of (aprox) matching with a priori results.
  A4 <- cbind(diag(JK), kronecker(t(c(-1L, 1L)), diag(JK)), matrix(0L, JK, 2L*IK))
  # Constraints for uniformity in last column, except JK value
  A5 <- cbind(kronecker(cbind(rep(1L, J - 2L), -diag(J - 2L), rep(0L, J - 2L)), 
                        t(c(rep(0L, K - 1L), 1L))), 
              matrix(0L, J - 2L, 2L*JK), matrix(0L, J - 2L, 2L*IK))
  # Constraint of zero for JK value
  A6 <- t(c(rep(0L, JK - 1L), 1L, rep(0L, 2L*JK), rep(0L, 2L*IK)))
  # Constraints of (aprox) matching in the I spatial units 
  A9 <- cbind(kronecker(X, diag(K)), matrix(0L, IK, 2L*JK), t(kronecker(diag(IK), c(1L, -1L))))
  
  # Terminos independientes
  b2 <- c(rep(1L, J))
  b3 <- yt
  #  b5 <- as.vector(t(xt * P0))
  b4 <- as.vector(t(P0))
  b5 <- rep(0L, J - 2L)
  b6 <- 0
  b9 <- as.vector(t(Y))
  
  fp <- rep(0L, JK) # pjk
  fe <- rep(1L -lambda, 2L*IK) # eik
  #  fs <- rep(lambda, 2L*JK) # epsilon_jk
  fs <- lambda*as.vector(t(matrix(rep(xt, length(yt)), length(xt), length(yt)))) # epsilon_jk  
  
  # Dealing with missing values in P0
  falta <- which(is.na(b4))
  A4[falta, ] <- 0L
  b4[falta] <- 0L
  fs[falta] <- 0L
  
  if (uniform){
    A <- rbind(A2, A3, A4, A5, A6, A9)
    b <- c(b2, b3, b4, b5, b6, b9)
    f <- c(fp, fs, fs, fe)
  } else {
    A <- rbind(A2, A3, A4, A6, A9)
    b <- c(b2, b3, b4, b6, b9)
    f <- c(fp, fs, fs, fe)
  }
  
  if (lambda == 1){
    b <- b[1L:(nrow(A) - IK)]
    if(min(xt[xt != 0]) > 1L) fs <- fs/min(xt[xt != 0])
    f <- c(fp, fs, fs)
    A <- A[1L:(nrow(A) - IK), 1L:(ncol(A) - 2L*IK)]
  }
  
  output <- list("A" = A, "b" = b, "f" = f)
  return(output)
}


# Function para calcular todos los componentes necesarios del programa lineal para
# calcular matrices de transferencia de estructuras type V and type VI, ver Pavia (2022).
model_lphom_apriori_5_6 <- function(X, Y, P0, lambda, uniform, ...){
  
  J <- ncol(X)
  K <- ncol(Y)
  I <- nrow(X)
  JK <- J*K
  IK <- I*K
  
  xt <- colSums(X)
  yt <- colSums(Y)
  
  # Constraints suma por filas unitaria de la matriz pjk
  A2 <- cbind(kronecker(diag(J), t(rep(1L, K))), matrix(0L, J, 2L*JK), matrix(0L, J, 2L*IK))
  # Coherencia de matriz global al aplicar la matriz pjk por filas para obtener las columnas
  A3 <- cbind(kronecker(t(xt), diag(K)),  matrix(0L, K, 2L*JK), matrix(0L, K, 2L*IK))
  # Constraints of (aprox) matching with a priori results.
  A4 <- cbind(diag(JK), kronecker(t(c(-1L, 1L)), diag(JK)), matrix(0L, JK, 2L*IK))
  # Constraints for uniformity in last column, except JK and (J-1)K values
  A5 <- cbind(kronecker(cbind(rep(1L, J - 3L), -diag(J - 3L), matrix(0L, J - 3L, 2L)), 
                        t(c(rep(0L, K - 1L), 1L))), 
              matrix(0L, J - 3L, 2L*JK), matrix(0L, J - 3L, 2L*IK))
  # Constraint of zero for (J-1)K value
  A6 <- t(c(rep(0L, (J - 1L)*K - 1L), 1L, rep(0L, 2*JK + K), rep(0L, 2L*IK)))
  # Constraint of zero for JK value
  A7 <- t(c(rep(0L, JK - 1L), 1L, rep(0L, 2L*JK), rep(0L, 2L*IK)))
  
  # Constraints of (aprox) matching in the I spatial units 
  A9 <- cbind(kronecker(X, diag(K)), matrix(0L, IK, 2L*JK), t(kronecker(diag(IK), c(1L, -1L))))
  
  # Terminos independientes
  b2 <- c(rep(1L, J))
  b3 <- yt
  #  b5 <- as.vector(t(xt * P0))
  b4 <- as.vector(t(P0))
  b5 <- rep(0L, J - 3L)
  b6 <- b7 <- 0
  b9 <- as.vector(t(Y))
  
  fp <- rep(0L, JK) # pjk
  fe <- rep(1L - lambda, 2L*IK) # eik
  #  fs <- rep(lambda, 2L*JK) # epsilon_jk
  fs <- lambda*as.vector(t(matrix(rep(xt, length(yt)), length(xt), length(yt)))) # epsilon_jk  
  
  # Dealing with missing values in P0
  falta <- which(is.na(b4))
  A4[falta, ] <- 0L
  b4[falta] <- 0L
  fs[falta] <- 0L
  
  if (uniform){
    A <- rbind(A2, A3, A4, A5, A6, A7, A9)
    b <- c(b2, b3, b4, b5, b6, b7, b9)
    f <- c(fp, fs, fs, fe)
  } else {
    A <- rbind(A2, A3, A4, A6, A7, A9)
    b <- c(b2, b3, b4, b6, b7, b9)
    f <- c(fp, fs, fs, fe)
  }
  
  if (lambda == 1){
    b <- b[1L:(nrow(A) - IK)]
    if(min(xt[xt != 0]) > 1L) fs <- fs/min(xt[xt != 0])
    f <- c(fp, fs, fs)
    A <- A[1L:(nrow(A) - IK), 1L:(ncol(A) - 2L*IK)]
  }
  
  output <- list("A" = A, "b" = b, "f" = f)
  return(output)
}


# Function para calcular todos los componentes necesarios del programa lineal para
# calcular matrices de transferencia de estructuras type VII, ver Pavia (2022).
model_lphom_apriori_7 <- function(X, Y, P0, lambda, ...){
  
  J <- ncol(X)
  K <- ncol(Y)
  I <- nrow(X)
  JK <- J*K
  IK <- I*K
  
  xt <- colSums(X)
  yt <- colSums(Y)
  
  # Constraints suma por filas unitaria de la matriz pjk
  A2 <- cbind(kronecker(diag(J), t(rep(1L, K))), matrix(0L, J, 2L*JK), matrix(0L, J, 2L*IK))
  # Coherencia de matriz global al aplicar la matriz pjk por filas para obtener las columnas
  A3 <- cbind(kronecker(t(xt), diag(K)),  matrix(0L, K, 2L*JK), matrix(0L, K, 2L*IK))
  # Constraints of (aprox) matching with a priori results.
  A4 <- cbind(diag(JK), kronecker(t(c(-1L, 1L)), diag(JK)), matrix(0L, JK, 2L*IK))
  # Constraints for uniformity in penultimate column
  A5 <- cbind(kronecker(cbind(rep(1L, J - 1L), -diag(J - 1L)), t(c(rep(0L, K - 2L), 1L, 0L))), 
              matrix(0L, J - 1L, 2L*JK), matrix(0L, J - 1L, 2L*IK))
  # Constraints for uniformity in last column
  A6 <- cbind(kronecker(cbind(rep(1L, J - 1L), -diag(J - 1L)), t(c(rep(0L, K - 1L), 1L))), 
              matrix(0L, J - 1L, 2L*JK), matrix(0L, J - 1L, 2L*IK))  
  # Constraints of (aprox) matching in the I spatial units 
  A9 <- cbind(kronecker(X, diag(K)), matrix(0L, IK, 2L*JK), t(kronecker(diag(IK), c(1L, -1L))))
  
  # Terminos independientes
  b2 <- c(rep(1L, J))
  b3 <- yt
  #  b5 <- as.vector(t(xt * P0))
  b4 <- as.vector(t(P0))
  b5 <- rep(0L, J - 1L)
  b6 <- rep(0L, J - 1L)
  b9 <- as.vector(t(Y))
  
  fp <- rep(0L, JK) # pjk
  fe <- rep(1L - lambda, 2L*IK) # eik
  #  fs <- rep(lambda, 2L*JK) # epsilon_jk
  fs <- lambda*as.vector(t(matrix(rep(xt, length(yt)), length(xt), length(yt)))) # epsilon_jk  
  
  # Dealing with missing values in P0
  falta <- which(is.na(b4))
  A4[falta, ] <- 0L
  b4[falta] <- 0L
  fs[falta] <- 0L
  
  A <- rbind(A2, A3, A4, A5, A6, A9)
  b <- c(b2, b3, b4, b5, b6, b9)
  f <- c(fp, fs, fs, fe)
  
  if (lambda == 1){
    b <- b[1L:(nrow(A) - IK)]
    if(min(xt[xt != 0]) > 1L) fs <- fs/min(xt[xt != 0])
    f <- c(fp, fs, fs)
    A <- A[1L:(nrow(A) - IK), 1L:(ncol(A) - 2L*IK)]
  }
  
  output <- list("A" = A, "b" = b, "f" = f)
  return(output)
}


# Function para calcular todos los componentes necesarios del programa lineal para
# calcular matrices de transferencia de estructuras type VIII and type IX, ver Pavia (2022).
model_lphom_apriori_8_9 <- function(X, Y, P0, lambda, uniform, ...){
  
  J <- ncol(X)
  K <- ncol(Y)
  I <- nrow(X)
  JK <- J*K
  IK <- I*K
  
  xt <- colSums(X)
  yt <- colSums(Y)
  
  # Constraints suma por filas unitaria de la matriz pjk
  A2 <- cbind(kronecker(diag(J), t(rep(1L, K))), matrix(0L, J, 2L*JK), matrix(0L, J, 2L*IK))
  # Coherencia de matriz global al aplicar la matriz pjk por filas para obtener las columnas
  A3 <- cbind(kronecker(t(xt), diag(K)),  matrix(0L, K, 2L*JK), matrix(0L, K, 2L*IK))
  # Constraints of (aprox) matching with a priori results.
  A4 <- cbind(diag(JK), kronecker(t(c(-1L, 1L)), diag(JK)), matrix(0L, JK, 2L*IK))
  # Constraints for uniformity in penultimate column
  A5 <- cbind(kronecker(cbind(rep(1L, J - 2L), -diag(J - 2L), rep(0L, J-2)), 
                        t(c(rep(0L, K - 2L), 1L, 0L))), 
              matrix(0L, J - 2L, 2L*JK), matrix(0L, J - 2L, 2L*IK))
  # Constraints for uniformity in last column
  A6 <- cbind(kronecker(cbind(rep(1L, J - 2L), -diag(J - 2L), rep(0L, J-2)), 
                        t(c(rep(0L, K - 1L), 1L))), 
              matrix(0L, J - 2L, 2L*JK), matrix(0L, J - 2L, 2L*IK))
  # Constraint of zero for JK-1 value
  A7 <- t(c(rep(0L, JK - 2L), 1L, rep(0L, 2L*JK + 1L), rep(0L, 2L*IK)))
  # Constraint of zero for JK value
  A8 <- t(c(rep(0L, JK - 1L), 1L, rep(0L, 2L*JK), rep(0L, 2L*IK)))
  # Constraints of (aprox) matching in the I spatial units 
  A9 <- cbind(kronecker(X, diag(K)), matrix(0L, IK, 2L*JK), t(kronecker(diag(IK), c(1L, -1L))))
  
  # Terminos independientes
  b2 <- c(rep(1L, J))
  b3 <- yt
  #  b5 <- as.vector(t(xt * P0))
  b4 <- as.vector(t(P0))
  b5 <- rep(0L, J - 2L)
  b6 <- rep(0L, J - 2L)
  b7 <- b8 <- 0L
  b9 <- as.vector(t(Y))
  
  fp <- rep(0L, JK) # pjk
  fe <- rep(1L - lambda, 2L*IK) # eik
  #  fs <- rep(lambda, 2L*JK) # epsilon_jk
  fs <- lambda*as.vector(t(matrix(rep(xt, length(yt)), length(xt), length(yt)))) # epsilon_jk  
  
  # Dealing with missing values in P0
  falta <- which(is.na(b4))
  A4[falta, ] <- 0L
  b4[falta] <- 0L
  fs[falta] <- 0L
  
  if (uniform){
    A <- rbind(A2, A3, A4, A5, A6, A7, A8, A9)
    b <- c(b2, b3, b4, b5, b6, b7, b8, b9)
    f <- c(fp, fs, fs, fe)
  } else {
    A <- rbind(A2, A3, A4, A7, A8, A9)
    b <- c(b2, b3, b4, b7, b8, b9)
    f <- c(fp, fs, fs, fe)
  }
  
  if (lambda == 1){
    b <- b[1L:(nrow(A) - IK)]
    if(min(xt[xt != 0]) > 1L) fs <- fs/min(xt[xt != 0])
    f <- c(fp, fs, fs)
    A <- A[1L:(nrow(A) - IK), 1L:(ncol(A) - 2L*IK)]
  }
  
  output <- list("A" = A, "b" = b, "f" = f)
  return(output)
}


# Function para calcular todos los componentes necesarios del programa lineal para
# calcular matrices de transferencia de estructuras type X and type XI, ver Pavia (2022).
model_lphom_apriori_10_11 <- function(X, Y, P0, lambda, uniform, ...){
  
  J <- ncol(X)
  K <- ncol(Y)
  I <- nrow(X)
  JK <- J*K
  IK <- I*K
  
  xt <- colSums(X)
  yt <- colSums(Y)
  
  # Constraints suma por filas unitaria de la matriz pjk
  A2 <- cbind(kronecker(diag(J), t(rep(1L, K))), matrix(0L, J, 2L*JK), matrix(0L, J, 2L*IK))
  # Coherencia de matriz global al aplicar la matriz pjk por filas para obtener las columnas
  A3 <- cbind(kronecker(t(xt), diag(K)),  matrix(0L, K, 2L*JK), matrix(0L, K, 2L*IK))
  # Constraints of (aprox) matching with a priori results.
  A4 <- cbind(diag(JK), kronecker(t(c(-1L, 1L)), diag(JK)), matrix(0L, JK, 2L*IK))
  # Constraints for uniformity in penultimate column
  A5 <- cbind(kronecker(cbind(rep(1L, J - 3L), -diag(J - 3L),  matrix(0L, J - 3L, 2L)), 
                        t(c(rep(0L, K - 2L), 1L, 0L))), 
              matrix(0L, J - 3L, 2L*JK), matrix(0L, J - 3L, 2L*IK))
  # Constraints for uniformity in last column
  A6 <- cbind(kronecker(cbind(rep(1L, J - 3L), -diag(J - 3L),  matrix(0L, J - 3L, 2L)), 
                        t(c(rep(0L, K - 1L), 1L))), 
              matrix(0L, J - 3L, 2L*JK), matrix(0L, J - 3L, 2L*IK))
  # Constraint of zero for J(K-1)-1 value
  A7.1 <- t(c(rep(0L, (J - 1L)*K - 2L), 1L, rep(0L, 2L*JK + K + 1L), rep(0L, 2L*IK)))
  # Constraint of zero for J(K-1) value
  A7.2 <- t(c(rep(0L, (J - 1L)*K - 1L), 1L, rep(0L, 2L*JK + K), rep(0L, 2L*IK)))
  # Constraint of zero for JK-1 value
  A8.1 <- t(c(rep(0L, JK - 2L), 1L, rep(0L, 2L*JK + 1L), rep(0L, 2L*IK)))
  # Constraint of zero for JK value
  A8.2 <- t(c(rep(0L, JK - 1L), 1L, rep(0L, 2L*JK), rep(0L, 2L*IK)))
  # Constraints of (aprox) matching in the I spatial units 
  A9 <- cbind(kronecker(X, diag(K)), matrix(0L, IK, 2L*JK), t(kronecker(diag(IK), c(1L, -1L))))
  
  # Terminos independientes
  b2 <- c(rep(1L, J))
  b3 <- yt
  #  b5 <- as.vector(t(xt * P0))
  b4 <- as.vector(t(P0))
  b5 <- rep(0L, J - 3L)
  b6 <- rep(0L, J - 3L)
  b7 <- b8 <- rep(0L, 2)
  b9 <- as.vector(t(Y))
  
  fp <- rep(0L, JK) # pjk
  fe <- rep(1L - lambda, 2L*IK) # eik
  #  fs <- rep(lambda, 2L*JK) # epsilon_jk
  fs <- lambda*as.vector(t(matrix(rep(xt, length(yt)), length(xt), length(yt)))) # epsilon_jk  
  
  # Dealing with missing values in P0
  falta <- which(is.na(b4))
  A4[falta, ] <- 0L
  b4[falta] <- 0L
  fs[falta] <- 0L
  
  
  if (uniform){
    A <- rbind(A2, A3, A4, A5, A6, A7.1, A7.2, A8.1, A8.2, A9)
    b <- c(b2, b3, b4, b5, b6, b7, b8, b9)
    f <- c(fp, fs, fs, fe)
  } else {
    A <- rbind(A2, A3, A4, A7.1, A7.2, A8.1, A8.2, A9)
    b <- c(b2, b3, b4, b7, b8, b9)
    f <- c(fp, fs, fs, fe)
  }
  
  if (lambda == 1){
    b <- b[1L:(nrow(A) - IK)]
    if(min(xt[xt != 0]) > 1L) fs <- fs/min(xt[xt != 0])
    f <- c(fp, fs, fs)
    A <- A[1L:(nrow(A) - IK), 1L:(ncol(A) - 2L*IK)]
  }
  
  output <- list("A" = A, "b" = b, "f" = f)
  return(output)
}


# Function to add structural zeros to the system
add_structural_zeros <- function(sistema, structural_zeros, K){
  # sistema: una salida de la funcion model_lphom_apriori_
  
  # Structural zero restrictions introduced by the user
  Ast <- matrix(0L, length(structural_zeros), ncol(sistema$A))
  bst <- rep(0L, length(structural_zeros))
  for (i in 1L:length(structural_zeros)){
    Ast[i, K*(structural_zeros[[i]][1L] - 1L) + structural_zeros[[i]][2L]] <- 1L
  }
  sistema$A <- rbind(sistema$A, Ast)
  sistema$b <- c(sistema$b, bst)
  
  return(sistema)
}

