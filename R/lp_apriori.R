#' Implements lp_apriori models
#'
#' @description  Adjusts an initial J0xK0 vote transfer matrix (ecological contingency table) to guarantee (i) congruency with aggregate results and (ii) completeness.
#'
#' @author Jose M. Pavia, \email{pavia@@uv.es}
#' @references Pavia, JM (2022). Adjusting initial estimates of voter transition probabilities to guarantee consistency and completeness. The function lp_apriori of the R-package lphom.
#'
#' @param votes_election1 data.frame (or matrix) of order IxJ1 (or vector of length J1)
#'                        with the votes gained by (or the numbers corresponding to) the J1
#'                        political options competing on election 1 (or origin) in the I
#'                        territorial units considered.
#'
#' @param votes_election2 data.frame (or matrix) of order IxK2 (or vector of length K2)
#'                        with the votes gained by (or the numbers corresponding to) the K2
#'                        political options competing on election 2 (or destination) in the I
#'                        territorial units considered.
#'
#' @param apriori data.frame (or matrix) of order J0xK0 with an initial estimate of the 
#'                (row-standarized) voter transition proportions/factions, pjk0, between
#'                the first J0 election options of election 1 and the first K0 election options
#'                of election 2. It could be also a data.frame (matrix) of counts. 
#'                This matrix can contain some missing values. 
#'                        
#' @param weights Either a numeric matrix (or data.frame) of order J0xK0 of weights, wjk, or
#'                a character string indicating the structure of weights to be used. As character  
#'                string this argument admits seven different values: `constant`, `x`, `xy`, `expected`, `counts`,
#'                `sqrt`, or `sd`. Default, `constant` (i.e., wjk = 1). The wjk coefficients measure the
#'                (relative) degree of confidence we have in the a priori values pjk0. 
#'                Everything else constant, the greater a weight wjk the closer the estimated pjk and 
#'                the pjk0 proportions will be.
#'                As numeric matrix, this matrix can contain some missing values, usually located in the same cells than 
#'                the missing values of `apriori`.
#'
#' @param new_and_exit_voters A character string indicating the level of information available
#'                            regarding new entries and exits of the election censuses between the
#'                            two elections. This argument captures the different options discussed
#'                            in Pavia (2022). This argument admits eight values:
#'                            `raw`, `regular`, `ordinary`, `simultaneous`, `enriched`, 
#'                            `semifull`, `full` and `gold`. Default, `raw`.
#' 
#' @param uniform A TRUE/FALSE value that indicates if census exits affect all the electoral options in
#'                a (relatively) similar fashion; depending on the scenario any equation(s) among equations (6) to (11) of
#'                Pavia (2022) could be used in the underlying model. Default, `TRUE`.
#'
#' @param solver A character string indicating the linear programming solver to be used, only
#'               `lp_solve` and `symphony` are allowed. By default, `lp_solve`. The package `Rsymphony`
#'               needs to be installed for the option `symphony` to be used.
#'
#' @param integers A TRUE/FALSE value that indicates whether the LP solution of counts (votes) must be 
#'                 approximate to the closest integer solution using ILP to generate the final solution.
#'                 Default, TRUE.
#'  
#' @param integers.solver A character string indicating the linear programming solver to be used to 
#'                        to the closest integer solution, only `symphony` and `lp_solve` are allowed.
#'                        By default, `symphony`. The package `Rsymphony` needs to be installed for the
#'                        option `symphony` to be used. Only used when `integers = TRUE`. 
#'                 
#' @param ... Other arguments to be passed to the function. Not currently used.
#'  

#'
#' @details Description of the `new_and_exit_voters` argument in more detail.
#' \itemize{
#'  \item{`raw`: }{The default value. This argument accounts for the most plausible scenario when
#'                 adjusting vote transfer matrices. A scenario with two elections elapsed at least 
#'                 some months where only the raw election data recorded in I territorial units 
#'                 (where I can be equal to one), in which the area under study is divided, are  
#'                 available. In this scenario, net exits and net entries are estimated according to  
#'                 equation (7) of Romero et al. (2020). When both net entries and exits are no 
#'                 null, constraint (15) of Pavia (2022) applies. If `uniform = TRUE` constraints    
#'                 (7) of Pavia (2022) are also imposed. In this scenario, J could be equal to J1 or  
#'                 J1 + 1 and K equal to K2 or K2 + 1.}
#'  \item{`simultaneous`: }{This is the value to be used in simultaneous elections and when the user 
#'                is interested in adjusting other type of transfer matrices such as the one that typically arise in 
#'                classical ecological inference problems, In this scenario, the sum by rows of  
#'                `votes_election1` and `votes_election2` must coincide. In this case, the function  
#'                just implements the basic model defined by equations (1) to (5) of Pavia (2022).}
#'  \item{`regular`: }{This value accounts for a scenario with 
#'                 two elections elapsed at least some months where (i) the column J1  
#'                 of `votes_election1` corresponds to new young electors who have the right 
#'                 to vote for the first time, (ii) net exits and maybe other additional 
#'                 net entries are computed according to equation (7) of Romero et al. (2020), and 
#'                 (iii) we can (or not) assume that net exits affect equally all the first J1 - 1  
#'                 options of election 1. When both net entries and exits are no null, constraints 
#'                 (12) and (15) of Pavia (2022) apply and if `uniform = TRUE` constraints    
#'                 (11) of Pavia (2022) are also imposed. In this scenario, J could be equal   
#'                 to J1 or J1 + 1 and K equal to K2 or K2 + 1.}
#'  \item{`ordinary`: }{This value accounts for a scenario  
#'                 with two elections elapsed at least some months where (i) the column K1  
#'                 of `votes_election2` corresponds to electors who died in the interperiod  
#'                 election, (ii) net entries and maybe other additional net exits are 
#'                 computed according to equation (7) of Romero et al. (2020), and (iii) we can 
#'                 assume that exits affect equally all the J1 options of election 1.
#'                  When both net entries and exits are no null, constraints (13) and 
#'                 (15) of Pavia (2022) apply and if `uniform = TRUE` constraints    
#'                 (8) and (9) of Pavia (2022) are also imposed. In this scenario, J could be equal   
#'                 to J1 or J1 + 1 and K equal to K2 or K2 + 1.}
#'  \item{`enriched`: }{This value accounts for a scenario that somewhat combine `regular` and 
#'                 `ordinary` ecenarios. We consider two elections elapsed at least some months where   
#'                 (i) the column J1 of `votes_election1` corresponds to new young electors
#'                  who have the right to vote for the first time, (ii) the column K1 of  
#'                 `votes_election2` corresponds to electors who died in the interperiod  
#'                 election, (iii) other (net) entries and (net) exits are computed according 
#'                 to equation (7) of Romero et al. (2020), and (iv) we can assume 
#'                 (or not) that exits affect equally all the J1 - 1 options of election 1.
#'                 When both net entries and exits are no null, constraints (12) to 
#'                 (15) of Pavia (2022) apply and if `uniform = TRUE` constraints    
#'                 (10) and (11) of Pavia (2022) are also imposed. In this scenario, J could be equal   
#'                 to J1 or J1 + 1 and K equal to K2 or K2 + 1.}
#'  \item{`semifull`: }{This value accounts for a scenario with two elections elapsed at least some
#'                months, where: (i) the column J of `votes_election1` totals new
#'                electors (young and immigrants) that have the right to vote for the first time and
#'                (ii) the column K of `votes_election2` corresponds to total exits of the census
#'                lists (due to death or emigration). In this scenario, the sum by rows of
#'                `votes_election1` and `votes_election2` must agree and constraint (15)
#'                of Pavia (2022) apply. Additionally, if `uniform = TRUE` constraints   
#'                (8) of Pavia (2022) are also imposed.}
#'  \item{`full`: }{This value accounts for a scenario with two elections elapsed at least some
#'                months, where (i) the column J - 1 of `votes_election1` totals new young
#'                electors that have the right to vote for the first time, (ii) the column J
#'                of `votes_election1` measures new immigrants that have the right to vote and
#'                (iii) the column K of `votes_election2` corresponds to total exits of the census
#'                lists (due to death or emigration). In this scenario, the sum by rows of
#'                `votes_election1` and `votes_election2` must agree and constraints (13)
#'                and (15) of Pavia (2022) apply.  Additionally, if `uniform = TRUE` constraints   
#'                (11) of Pavia (2022) are also imposed.}
#'  \item{`gold`: }{This value accounts for a scenario similar to full, where total exits are
#'               separated out between exits due to emigration (column K - 1 of `votes_election2`)
#'               and death (column K of `votes_election2`). In this scenario, the sum by rows
#'               of `votes_election1` and `votes_election2` must agree. Constraints (12) to 
#'               (15) of Pavia (2022) apply and if `uniform = TRUE` constraints (10) and (11)   
#'               of Pavia (2022) are also imposed.}
#' }
#' 
#' @return
#' A list with the following components
#'  \item{VTM}{ A matrix of order JxK with the estimated percentages of row-standardized vote transitions from election 1 to election 2. 
#'             In `raw`, `regular`, `ordinary` and `enriched` scenarios when the percentage of net entries is small, less than 1% of the census, 
#'             net entries are omitted (i.e., J = J0) even if estimates for net entries different from zero are obtained. Likewise, in the same scenarios When the percentage of net exits is small, less than 1%
#'             of the census, net exits are omitted (i.e., K = K0) even if estimates for net exits different from zero are obtained.}
#'  \item{VTM.votes}{ A matrix of order JxK with the estimated vote transitions from election 1 to election 2.
#'             In `raw`, `regular`, `ordinary` and `enriched` scenarios when the percentage of net entries is small, less than 1% of the census, 
#'             net entries are omitted (i.e., J = J0) even if estimates for net entries different from zero are obtained. Likewise, in the same scenarios When the percentage of net exits is small, less than 1%
#'             of the census, net exits are omitted (i.e., K = K0) even if estimates for net exits different from zero are obtained.}
#'  \item{weights}{ A matrix of order JxK with the weights used to adjust the a priori vote transitions from  election 1 to election 2.}
#'  \item{OTM}{ A matrix of order KxJ with the estimated percentages of the origin of the votes obtained for the different options of election 2.}
#'  \item{VTM.complete}{ A matrix of order JxK with the estimated proportions of row-standardized vote transitions from election 1 to election 2.
#'                      In `raw`, `regular`, `ordinary` and `enriched` scenarios, this matrix includes the row and the column corresponding to net entries 
#'                      and net exits (when they are present) even when they are really small.}
#'  \item{VTM.complete.votes}{ A matrix of order JxK with the estimated vote transitions from election 1 to election 2.
#'                      In `raw`, `regular`, `ordinary` and `enriched` scenarios, this matrix includes the row and the column corresponding to net entries 
#'                      and net exits (when they are present) even when they are really small.}
#'  \item{inputs}{ A list containing all the objects with the values used as arguments by the function.}
#'  \item{origin}{ A vector with the final data used as votes of the origin election after taking into account the level of information available regarding to new entries and exits of the election censuses between the two elections.}
#'  \item{destination}{ A vector with the final data used as votes of the origin election after taking into account the level of information available regarding to new entries and exits of the election censuses between the two elections.}
#'
#' @export
#'
#' @family linear programing ecological inference functions
#' @seealso \code{\link{lphom}} \code{\link{tslphom}} \code{\link{nslphom}} \code{\link{lclphom}}
#'
#' @examples
#' P0 <- matrix(c(.75, .02, .15, .08, .01, .01, .97, .01,
#'                .01, .01, .01, .97, .01, .10, .80, .09,
#'                .20, .30, .30, .20, .10, .10, .50, .30,
#'                .10, .30, NA, NA, .25, .20, NA, NA), 
#'              byrow = TRUE, 8, 4)
#' mt <- lp_apriori(France2017P[, 1:8], France2017P[, 9:12], P0, integers = FALSE)
#
#' @importFrom lpSolve lp
#
lp_apriori <- function(votes_election1,
                       votes_election2,
                       apriori,
                       weights = "constant",
                       new_and_exit_voters ="raw",
                       uniform = TRUE, 
                       solver = "lp_solve",
                       integers = TRUE,
                       integers.solver = "symphony",
                       ...){
  
  argg <- c(as.list(environment()), list(...))
  datos <- test_inputs_lp_apriori(argg)
  # inputs
  inputs <- list("votes_election1" = votes_election1, "votes_election2" = votes_election2,
                 "apriori" = apriori, "weights" = weights , "new_and_exit_voters" = new_and_exit_voters[1], 
                 "uniform" = uniform, "solver" = solver, "integers" = integers,
                 "integers.solver" = integers.solver)
  
  
  
  # Aggregate results    
  x0 <- datos$x
  y0 <- datos$y
  
  # Tests apriori
  aprioris <- test_apriori_weights(x0 = x0, y0 = y0, P0 = apriori, 
                                   pesos0 = weights, scenario = new_and_exit_voters[1])
  apriori <- aprioris$apriori
  weights <- aprioris$weights
  
  # Aggregate results with net entries and exits
  datos <- lp_apriori_net(x0 = x0, y0 = y0)
  
  # Adjusting of data
  adjusting_function <- get(paste0("lp_apriori_", new_and_exit_voters[1]))
  solucion <- adjusting_function(x0 = datos$x0, y0 = datos$y0,
                            x = datos$x, y = datos$y,
                            P0 = apriori, weights = weights,
                            uniform = uniform, solver = solver)
  pjk <- solucion$VTM.complete
  # Pesos
  if (is.character(weights)){
    pesos <- calcular_weights(x = datos$x, y = datos$y, P0 = apriori, weights = weights)
  } else {
    pesos <- solucion$weights
  }
  
  # Votos   
  vjk <- pjk*datos$x
  if (integers){
    if (integers.solver == "lp_solve"){
      dec2counts <- dec2counts_lp
    } else {
      dec2counts <- dec2counts_symphony
    }
    vjk <- dec2counts(vjk, datos$x, datos$y)
  }
  
  rownames(vjk) <- rownames(pesos) <- names(datos$x)
  colnames(vjk) <- colnames(pesos) <- names(datos$y) 
  pjk <- vjk/rowSums(vjk)
  pkj <- t(vjk/colSums(vjk))
  
  # Caso de filas o columnas con cero votos
  filas0 <- which(rowSums(vjk) == 0L)
  colum0 <- which(colSums(vjk) == 0L)
  pjk[filas0, ] <- 0L
  pkj[colum0, ] <- 0L

  VTM <- VTM.reduced <- round(100*pjk, 2L)
  vjk.reduced <- vjk
  OTM <- round(100*pkj, 2L)
  
  # Presentation of the matrix when net_entries or net_exits are small
  net_entries <- net_exits <- Inf
  if(length(datos$x0) < length(datos$x)) 
    net_entries <- datos$x[length(datos$x)]/sum(datos$x)
  if(length(datos$y0) < length(datos$y)) 
    net_exits <- datos$y[length(datos$y)]/sum(datos$y)
  
  if(net_entries < 0.01){
    VTM.reduced <- VTM.reduced[-nrow(VTM), ]
    vjk.reduced <- vjk.reduced[-nrow(vjk), ]
  }
  
  if(net_exits < 0.01){
    VTM.reduced <- VTM.reduced[, -ncol(VTM)]
    vjk.reduced <- vjk.reduced[, -ncol(vjk)]
  }
  
  output <- list("VTM" = VTM.reduced, "VTM.votes" = vjk.reduced, "weights" = pesos, "OTM" = OTM,
                 "VTM.complete" = pjk, "VTM.complete.votes" = vjk,
                 "inputs" = inputs, "origin" = datos$x, "destination" = datos$y)
  class(output) <- c("lphom", "lp_apriori", "ei_lp")
  return(output)
}
