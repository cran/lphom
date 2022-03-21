#' Implements lphom algorithm
#'
#' @description  Estimates RxC vote transfer matrices (ecological contingency tables) with lphom
#'
#' @author Jose M. Pavia, \email{pavia@@uv.es}
#' @author Rafael Romero \email{rromero@@eio.upv.es}
#' @references Romero, R, Pavia, JM, Martin, J and Romero G (2020). Assessing uncertainty of voter transitions estimated from aggregated data. Application to the 2017 French presidential election. *Journal of Applied Statistics*, 47(13-15), 2711-2736. \doi{10.1080/02664763.2020.1804842}
#'
#' @param votes_election1 data.frame (or matrix) of order IxJ (likely of final order Ix(J-1)
#'                        in `regular` and `raw` scenarios when net entries are
#'                        estimated by the function) with the votes gained by the *J*
#'                        political options competing on election 1 (or origin) in the *I*
#'                        territorial units considered. In general, the row marginals 
#'                        of the *I* tables.
#'
#' @param votes_election2 data.frame (or matrix) of order IxK (likely of final order Ix(K-1)
#'                        in `regular` and `raw` scenarios  when net exits are
#'                        estimated by the function) with the votes gained by
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
#'                         lphom implicitly states `structural_zeros = list(c(J, K))` in case exits and/or
#'                         entries are computed because the sum by rows of `votes_election1` and
#'                         `votes_election2` does not coincide.
#'
#' @param integers A TRUE/FALSE value that indicates whether the LP solution of counts (votes) must be approximate
#'                 to the closest integer solution using ILP to generate the final solution. Default, FALSE.
#'
#' @param verbose A TRUE/FALSE value that indicates if the main outputs of the function should be
#'                printed on the screen. Default, FALSE.
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

#'
#' @details Description of the `new_and_exit_voters` argument in more detail.
#' \itemize{
#'  \item{`raw`: }{The default value. This argument accounts for the most plausible scenario when
#'                 estimating vote transfer matrices: A scenario with two elections elapsed at least some
#'                 months where only the raw election data recorded in the *I* territorial units, 
#'                 in which the area under study is divided, are available. 
#'                 In this scenario, net exits (basically deaths) and net entries (basically 
#'                 new young voters) are estimated according to equation (7) of Romero et al. (2020). 
#'                 Constraints defined by equations (8) and (9) of Romero et al. (2020) are imposed. 
#'                 In this scenario, when net exits and/or net entries are negligible (such as between 
#'                 the first- and second-round of French Presidential elections), they are omitted in 
#'                 the outputs.}
#'  \item{`regular`: }{For estimating vote transfer matrices, this value accounts for a scenario with 
#'                 two elections elapsed at least some months where (i) the column *J* of `votes_election1` 
#'                 corresponds to new young electors who have the right to vote for the first time and (ii)
#'                 net exits (basically a consequence of mortality), and maybe other additional net entries,
#'                 are computed according equation (7) of Romero et al. (2020), and (iii) we
#'                 assume net exits affect equally all the first *J-1* options of election 1,
#'                 hence (8) and (9) constraints of Romero et al. (2020) are imposed.}
#'  \item{`simultaneous`: }{This is the value to be used in a classical ecological inference problems, 
#'                such as for racial voting, and in a scenario with two simultaneous elections. 
#'                In this scenario, the sum by rows of `votes_election1` and `votes_election2` must coincide. 
#'                Constraints defined by equations (8) and (9) of Romero et al. (2020) are not included in 
#'                the model.}
#'  \item{`full`: }{This value accounts for a scenario with two elections elapsed at least some
#'                months, where: (i) the column *J-1* of `votes_election1` totals new young
#'                electors that have the right to vote for the first time; (ii) the column *J*
#'                of `votes_election1` measures new immigrants that have the right to vote; and
#'                (iii) the column *K* of `votes_election2` corresponds to total exits of the census
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
#'  \item{HETe}{ The estimated heterogeneity index defined in equation (11) of Romero et al. (2020).}
#'  \item{VTM.complete}{ A matrix of order J'xK' with the estimated proportions of row-standardized vote transitions from election 1 to election 2, including in `regular` and `raw` scenarios the row and the column corresponding to net_entries and net_exits even when they are really small, less than 1% in all units.}
#'  \item{VTM.complete.votes}{ A matrix of order J'xK' with the estimated vote transitions from election 1 to election 2, including in `regular` and `raw` scenarios the row and the column corresponding to net_entries and net_exits even when they are really small, less than 1% in all units.}
#'  \item{inputs}{ A list containing all the objects with the values used as arguments by the function.}
#'  \item{origin}{ A matrix with the final data used as votes of the origin election after taking into account the level of information available regarding to new entries and exits of the election censuses between the two elections.}
#'  \item{destination}{ A matrix with the final data used as votes of the origin election after taking into account the level of information available regarding to new entries and exits of the election censuses between the two elections.}
#'  \item{EHet}{ A matrix of order IxK measuring in each spatial unit a distance to the homogeneity hypothesis. That is, the differences under the homogeneity hypothesis between the actual recorded results and the expected results in each territorial unit for each option of election 2.}
#' @export
#'
#'
#' @family linear programing ecological inference functions
#' @seealso \code{\link{tslphom}} \code{\link{nslphom}} \code{\link{lclphom}}
#'
#' @examples
#' lphom(France2017P[, 1:8] , France2017P[, 9:12], new_and_exit_voters= "raw")
#
#' @importFrom lpSolve lp
#
lphom <- function(votes_election1,
                  votes_election2,
                  new_and_exit_voters = c("raw", "regular", "simultaneous", "full", "gold"),
                  structural_zeros = NULL,
                  integers = FALSE,
                  verbose = FALSE,
                  solver = "lp_solve",
                  integers.solver = "symphony",
                  ...){

# Loading package lpSolve
#  if (!require(lpSolve)) install.packages("lpSolve", repos = "http://cran.rstudio.com")
#  require(lpSolve)

  argg <- c(as.list(environment()), list(...))
  integers <- test_integers(argg)
  
  # inputs
  inputs <- list("votes_election1" = votes_election1, "votes_election2" = votes_election2,
                 "new_and_exit_voters" = new_and_exit_voters[1], "structural_zeros" = structural_zeros,
                 "integers" = integers, "verbose" = verbose, "solver" = solver,
                 "integers.solver" = integers.solver)

  # Data conditions
  x <- as.matrix(votes_election1)
  y <- as.matrix(votes_election2)
  if (nrow(x) != nrow(y))
    stop('The number of spatial units is different in origin and destination.')
  # new_and_exit_voters = match.arg(new_and_exit_voters)
  new_and_exit_voters = new_and_exit_voters[1]
  if (!(new_and_exit_voters %in% c("regular", "raw", "simultaneous", "full", "gold")))
    stop('Not allowed string for argument "new_and_exit_voters".
          The only allowed strings for "new_and_exit_voters" are "regular", "raw", "simultaneous", "full" and "gold".')
  if (new_and_exit_voters %in% c("simultaneous", "full", "gold")){
    if (!identical(round(rowSums(x)), round(rowSums(y)))){
      texto <- paste0('The number of voters (electors) in Election 1 and ',
                      'Election 2 differ in at least a territorial unit. This is not ',
                      'allowed in a \"', new_and_exit_voters, ' scenario \".')
      stop(texto)
    }
  }
  if (min(x,y) < 0) stop('Negative values for voters (electors) are not allowed')
  if (!(solver %in% c("symphony", "lp_solve")))
    stop('Only "symponhy", "lp_solve" are allowed as solvers')
  if (!(integers.solver %in% c("symphony", "lp_solve")))
    stop('Only "symponhy", "lp_solve" are allowed as solvers')
  
  if (integers.solver == "lp_solve"){
    dec2counts <- dec2counts_lp
  } else {
    dec2counts <- dec2counts_symphony
  }

  # Data preparation
  net_entries = net_exits = TRUE
  tt = sum(x)
  if (new_and_exit_voters %in% c("regular", "raw")){
    NET_ENTRIES = NET_EXITS = rep(0L, nrow(x))
    x = cbind(x, NET_ENTRIES); y = cbind(y, NET_EXITS)

    # Estimation of net entries/exits in the election census
    d = rowSums(y) - rowSums(x)
    if (any(d != 0L)) {
      th = sum(d[d > 0])
      te = -sum(d[d < 0])
      message(paste0('*********************WARNING*********************\n',
                     'You are in a \"', new_and_exit_voters, '\" scenario.\n',
                     'The sums (by row) of origin and destination data differ. ',
                     'It is, for at least a spatial unit, the total number ',
                     'of electors in both elections is not the same. \n',
                     '\n To guarantee the matching: A new category of census entries ',
                     '(NET_ENTRIES) is included in the origin election and ',
                     'a new category of census exits (NET_EXITS) ',
                     'is also included in the destination election.\n',
                     '\n %NET_ENTRIES=',100*th/tt,'%\n',
                     ' %NET_EXITS=',100*te/tt,'%\n',
                     '\n If NET_ENTRIES and/or NET_EXITS are really small, less than 1% ',
                     'in all units, their results will not be displayed in VTM, ',
                     'they are included in VTM.complete\n',
                     '\n See paper "A model to estimate voter transitions from aggregated',
                     ' data. Application to 2017 French presidential elections" for details. \n',
                     '*************************************************'))
      x[,ncol(x)] = d*(d > 0)
      y[,ncol(y)] = -d*(d < 0)
    }

    # Net entries and exits
    if (sum(x[,ncol(x)]) == 0L){
      net_entries = FALSE
      x = x[,-ncol(x)]
    }
    if (sum(y[,ncol(y)]) == 0L){
      net_exits = FALSE
      y = y[,-ncol(y)]
    }
  }
  # Parameters
  I = nrow(x); J = ncol(x); K = ncol(y); JK = J*K; IK = I*K
  # Names of election options
  names1 = colnames(x); names2 = colnames(y)
  # Constraints sum(pjk)=1
  a1 = cbind(kronecker(diag(J), t(rep(1L, K))), matrix(0L, J, 2L*IK))
#  a1 = matrix(0, J, JK)
#  for (j in 1:J) {
#    a1[j,((j-1)*K+1):(j*K)]=1
#  }
#  a1 = cbind(a1, matrix(0, J, 2*IK))
  b1 = rep(1L, J)
  # Constraints total match for K parties
  xt = colSums(x)
  yt = colSums(y)
  at = cbind(t(kronecker(xt, diag(K))), matrix(0L, K, 2L*IK))
#  at = cbind(kronecker(diag(1, J, J), t(rep(1,K))), matrix(0, J, 2*IK))
#  at = matrix(0, K, JK+2*IK)
#  for (j in 1:J){
#    at[,((j-1)*K+1):(j*K)] = diag(xt[j], K)
#  }
  bt = yt
  # Constraints to match votes in the I spatial units and K parties
#  ap = matrix(0, 0, JK+2*IK)
# bp = NA
#  for (i in 1:I) {
#    ai = matrix(0, K, 0)
#    for (j in 1:J) {ai = cbind(ai,diag(x[i,j],K))} # ai_0 = t(kronecker(x[i,], diag(1, K, K)))
#    ai = cbind(ai,matrix(0,K,2*IK))
#   bi= y[i,]
#    ap = rbind(ap,ai)
#    bp = c(bp,bi)
#  }
#  ap = cbind(kronecker(x, diag(K)), matrix(0, IK, 2*IK))
#  bp = bp[-1]
  bp = as.vector(t(y))
#  for (f in 1:IK) {ap[f,JK+c(2*f-1,2*f)] = c(1,-1)}
  ap = cbind(kronecker(x, diag(K)), t(kronecker(diag(IK), c(1L,-1L))))

  # Joining the three sets of constraints
  ae = rbind(a1,at,ap)
  be = c(b1,bt,bp)
  # Raw scenario. Constraints pjk(J, K)=0 & pjk(j,K) constant, related,
  # respectively, to net entries and net exits.
  if (new_and_exit_voters == "raw" & net_exits){
    if (net_entries){
      pb = yt[K]/sum(xt[1L:(J-1L)])
      ab = matrix(0L, J, JK+2L*IK)
      bb = rep(0L, J)
      for (j in 1L:J){
        ab[j,j*K] = 1L
        bb[j]=pb*(j<J)
      }
      ae = rbind(ae,ab)
      be = c(be,bb)
    } else {
      pb = yt[K]/sum(xt[1L:J])
      ab = matrix(0L, J, JK+2L*IK)
      bb = rep(0L, J)
      for (j in 1L:J){
        ab[j,j*K] = 1L
        bb[j]=pb
      }
      ae = rbind(ae,ab)
      be = c(be,bb)
    }
  }
  # Regular scenario. Constraint pjk(J, K)= pik(J-1,K) = 0 related to new young
  # voters and net entries and pjk(j,K) constant related to net exits.
  if (new_and_exit_voters == "regular" & net_exits){
    if (net_entries){
      pb = yt[K]/sum(xt[1L:(J-2L)])
      ab = matrix(0L, J , JK+2L*IK)
      bb = rep(0L, J)
      for (j in 1L:J){
        ab[j,j*K] = 1L
        bb[j]=pb*(j<(J-1L))
      }
      ae = rbind(ae,ab)
      be = c(be,bb)
    } else {
      pb = yt[K]/sum(xt[1L:(J-1L)])
      ab = matrix(0L, J, JK+2L*IK)
      bb = rep(0L, J)
      for (j in 1L:J){
        ab[j,j*K] = 1L
        bb[j]=pb*(j<J)
      }
      ae = rbind(ae,ab)
      be = c(be,bb)
    }
  }
  # Full scenario. Constraint pjk(J, K)= pik(J-1,K) = 0 related to new young
  # voters and immigrants and pjk(j,K) constant related to exits.
  if (new_and_exit_voters == "full"){
    pb = yt[K]/sum(xt[1L:(J-2L)])
    ab = matrix(0,J,JK + 2L*IK)
    bb = rep(0,J)
    for (j in 1L:J){
      ab[j,j*K] = 1L
      bb[j]=pb*(j<J-1L)
    }
    ae = rbind(ae,ab)
    be = c(be,bb)
  }
  # Gold scenario. Constraints pjk(J, K-1) = pik(J,K) = pjk(J-1, K-1) =
  # = pik(J-1,K) =0 related to new young voters and immigrants and
  # pjk(j,K-1) and pjk(j,K) constant related to exits.
  if (new_and_exit_voters == "gold"){
    # Column K-1
    pb = yt[K-1L]/sum(xt[1L:(J-2L)])
    ab = matrix(0L, J, JK+2L*IK)
    bb = rep(0,J)
    for (j in 1L:J){
      ab[j,(K-1L) + (j-1L)*K] = 1L
      bb[j]=pb*(j<(J-1L))
    }
    ae = rbind(ae,ab)
    be = c(be,bb)
    # Column K
    pb = yt[K]/sum(xt[1L:(J - 2L)])
    ab = matrix(0L, J, JK+2L*IK)
    bb = rep(0L, J)
    for (j in 1L:J){
      ab[j,j*K] = 1L
      bb[j]=pb*(j<(J-1L))
    }
    ae = rbind(ae,ab)
    be = c(be,bb)
  }
  # Structural zero restrictions introduced by the user
  if (length(structural_zeros) > 0){
    ast = matrix(0L, length(structural_zeros), JK+2L*IK)
    bst = rep(0L, length(structural_zeros))
    for (i in 1L:length(structural_zeros)){
      ast[i, K*(structural_zeros[[i]][1L]-1L)+structural_zeros[[i]][2L]] = 1L
    }
    ae = rbind(ae,ast)
    be = c(be,bst)
  }
  # Objective function, to minimize
  fun.obj = c(rep(0L, JK), rep(1L, 2L*IK))
  # Solution
  if (solver == "symphony"){
      sol = Rsymphony::Rsymphony_solve_LP(obj = fun.obj,
                                          mat = ae,
                                          dir = rep('==', length(be)),
                                          rhs = be)
  } else {
     sol = suppressWarnings(lpSolve::lp('min', fun.obj, ae, rep('=', length(be)), be) )

  }
  z = sol$solution
  pjk = matrix(z[1:JK], J, K, TRUE, dimnames = list(names1, names2))
  e = z[(JK+1L):(JK+2L*IK)]
  e = matrix(e, 2L, IK)
  e = e[1L,] - e[2L,]
  eik = t(matrix(e,K,I))
  if (integers){
    vjk <- pjk*colSums(x)
    vjk <- dec2counts(vjk, colSums(x), colSums(y))
    pjk <- vjk/rowSums(vjk)
    dimnames(pjk) <- list(names1, names2)
    eik = y - x %*% pjk
  }
  colnames(eik) = names2
  rownames(eik) = rownames(x)
  EHet = eik
  pkj = matrix(0L, K, J, dimnames=list(names2,names1))
  for (k in 1L:K) {
    for (j in 1L:J) {
      pkj[k,j]=xt[j]*pjk[j,k]/yt[k]
    }
  }
  pjk.complete <- pjk
  vjk <- pjk*colSums(x)
  vjk.complete <- vjk
  HIe = 100*sum(abs(eik))/sum(vjk.complete)
  if (new_and_exit_voters %in% c("regular", "raw")){
    if (net_entries & max(x[,ncol(x)]/rowSums(x)) < 0.01){
      pjk = pjk[-J,]; eik = eik[-J,]; pkj = pkj[,-J]; vjk = vjk[-J,]
    }
    if (net_exits & max(y[,ncol(y)]/rowSums(y)) < 0.01){
      pjk = pjk[,-K]; eik = eik[,-K]; pkj = pkj[-K,]; vjk = vjk[,-K]
    }
  }
  pjk = round(100*pjk, 2)
  pkj = round(100*pkj, 2)
  if (verbose){
    cat("\n\nEstimated Heterogeneity Index HETe:",round(HIe, 2),"%\n")
    cat("\n\n Matrix of vote transitions (in %) from Election 1 to Election 2\n\n")
    print(pjk)
    cat("\n\n Origin (%) of the votes obtained in Election 2\n\n")
    print(pkj)
  }
  
  # Caso de filas o columnas con cero votos
  filas0 <- which(rowSums(vjk.complete) == 0)
  colum0 <- which(colSums(vjk.complete) == 0)
  pjk[filas0, ] <- 0
  pjk.complete[filas0, ] <- 0
  pkj[colum0, ] <- 0
  
  output <- list("VTM" = pjk, "VTM.votes" = vjk, "OTM" = pkj, "HETe" = HIe,
                 "VTM.complete" = pjk.complete, "VTM.complete.votes" = vjk.complete,
                 "inputs" = inputs, "origin" = x, "destination" = y, "EHet" = EHet)
  class(output) <- c("lphom", "ei_lp")
  return(output)
}
