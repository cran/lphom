#' Implements lphom algorithm
#'
#' @description  Estimates RxC vote transfer matrices (ecological contingency tables) with lphom
#'
#' @author Jose M. Pavia, \email{pavia@@uv.es}
#' @author Rafael Romero \email{rromero@@eio.upv.es}
#' @references Romero, R, Pavia, JM, Martin, J and Romero G (2020). Assessing uncertainty of voter transitions estimated from aggregated data. Application to the 2017 French presidential election. *Journal of Applied Statistics*, 47(13-15), 2711-2736. \doi{10.1080/02664763.2020.1804842}
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
#' @param structural_zeros Default NULL. A list of vectors of length two, indicating the election options
#'                         for which no transfer of votes are allowed between election 1 and election 2.
#'                         For instance, when new_and_exit_voters is set to `"regular"`,
#'                         lphom implicitly `states structural_zeros = list(c(J, K))` in case exits and/or
#'                         entries are computed because the sum by rows of `votes_election1` and
#'                         `votes_election2` does not coincide.
#'
#' @param verbose A TRUE/FALSE value that indicates if the main outputs of the function should be
#'                printed on the screen. Default, FALSE.
#'
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
#'  \item{HETe}{ The estimated heterogeneity index defined in equation (11) of Romero et al. (2020).}
#'  \item{EHet}{ A matrix of order IxK measuring in each spatial unit a distance to the homogeneity hypothesis, that is, the differences under the homogeneity hypothesis between the actual recorded results and the expected results in each territorial unit for each option of election two.}
#'  \item{inputs}{ A list containing all the objects with the values used as arguments by the function.}
#'  \item{origin}{ A matrix with the final data used as votes of the origin election after taking into account the level of information available regarding to new entries and exits of the election censuses between the two elections.}
#'  \item{destination}{ A matrix with the final data used as votes of the origin election after taking into account the level of information available regarding to new entries and exits of the election censuses between the two elections.}
#'  \item{VTM.complete}{ A matrix of order J'xK' with the estimated proportions of vote transitions from election 1 to election 2, including in `regular` and `raw` scenarios the row and the column corresponding to net_entries and net_exits even when they are really small, less than 1% in all units.}
#' @export
#'
#'
#' @family linear programing ecological inference functions
#' @seealso \code{\link{tslphom}} \code{\link{nslphom}}
#'
#' @examples
#' lphom(France2017P[, 1:8] , France2017P[, 9:12], new_and_exit_voters= "raw", 
#'       structural_zeros = NULL, verbose = FALSE)
#
#' @importFrom lpSolve lp
#
lphom <- function(votes_election1, votes_election2,
                  new_and_exit_voters = c("regular", "raw", "simultaneous", "full", "gold"),
                  structural_zeros = NULL, verbose = FALSE){

  # Loading package lpSolve
#  if (!require(lpSolve)) install.packages("lpSolve", repos = "http://cran.rstudio.com")
#  require(lpSolve)

  # inputs
  inputs <- list("votes_election1" = votes_election1, "votes_election2" = votes_election2,
                 "new_and_exit_voters" = new_and_exit_voters[1], "structural_zeros" = structural_zeros,
                 "verbose" = verbose)

  # Data conditions
  x = as.matrix(votes_election1)
  y = as.matrix(votes_election2)
  if (nrow(x) != nrow(y))
    stop('The number of spatial units is different in origin and destination.')
  new_and_exit_voters = new_and_exit_voters[1]
  if (!(new_and_exit_voters %in% c("regular", "raw", "simultaneous", "full", "gold")))
    stop('Not allowed string for argument "new_and_exit_voters".
          The only allowed strings for "new_and_exit_voters" are "regular", "raw", "simultaneous", "full" and "gold".')
  if (new_and_exit_voters %in% c("simultaneous", "full", "gold")){
    if (!identical(round(rowSums(x)), round(rowSums(y)))){
      texto <- paste0('The number of voters (electors) in Election 1 and',
                      'Election 2 differ in at least a territorial unit. This is not ',
                      'allowed in a \"', new_and_exit_voters, ' scenario \".')
      stop(texto)
    }
  }
  if (min(x,y) < 0) stop('Negative values for voters (electors) are not allowed')

  # Data preparation
  net_entries = net_exits = TRUE
  tt = sum(x)
  if (new_and_exit_voters %in% c("regular", "raw")){
    NET_ENTRIES = NET_EXITS = rep(0,nrow(x))
    x = cbind(x, NET_ENTRIES); y = cbind(y, NET_EXITS)

    # Estimation of net entries/exits in the election census
    d = rowSums(y) - rowSums(x)
    if (any(d != 0)) {
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
    if (sum(x[,ncol(x)]) == 0){
      net_entries = FALSE
      x = x[,-ncol(x)]
    }
    if (sum(y[,ncol(y)]) == 0){
      net_exits = FALSE
      y = y[,-ncol(y)]
    }
  }
  # Parameters
  I = nrow(x); J = ncol(x); K = ncol(y); JK = J*K; IK = I*K
  # Names of election options
  names1 = colnames(x); names2 = colnames(y)
  # Constraints sum(pjk)=1
  a1 = matrix(0, J, JK)
  for (j in 1:J) {
    a1[j,((j-1)*K+1):(j*K)]=1
  }
  a1 = cbind(a1, matrix(0, J, 2*IK))
  b1 = rep(1,J)
  # Constraints total match for K parties
  xt = colSums(x)
  yt = colSums(y)
  at = matrix(0, K, JK+2*IK)
  for (j in 1:J){
    at[,((j-1)*K+1):(j*K)] = diag(xt[j], K)
  }
  bt = yt
  # Constraints to match votes in the I spatial units and K parties
  ap = matrix(0, 0, JK+2*IK)
  bp = NA
  for (i in 1:I) {
    ai = matrix(0, K, 0)
    for (j in 1:J) {ai = cbind(ai,diag(x[i,j],K))}
    ai = cbind(ai,matrix(0,K,2*IK))
    bi= y[i,]
    ap = rbind(ap,ai)
    bp = c(bp,bi)
  }
  bp = bp[-1]
  for (f in 1:IK) {ap[f,JK+c(2*f-1,2*f)] = c(1,-1)}
  # Joining the three sets of constraints
  ae = rbind(a1,at,ap)
  be = c(b1,bt,bp)
  # Raw scenario. Constraints pjk(J, K)=0 & pjk(j,K) constant, related,
  # respectively, to net entries and net exits.
  if (new_and_exit_voters == "raw" & net_exits){
    if (net_entries){
      pb = yt[K]/sum(xt[1:(J-1)])
      ab = matrix(0,J,JK+2*IK)
      bb = rep(0,J)
      for (j in 1:J){
        ab[j,j*K] = 1
        bb[j]=pb*(j<J)
      }
      ae = rbind(ae,ab)
      be = c(be,bb)
    } else {
      pb = yt[K]/sum(xt[1:J])
      ab = matrix(0,J,JK+2*IK)
      bb = rep(0,J)
      for (j in 1:J){
        ab[j,j*K] = 1
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
      pb = yt[K]/sum(xt[1:(J-2)])
      ab = matrix(0,J,JK+2*IK)
      bb = rep(0,J)
      for (j in 1:J){
        ab[j,j*K] = 1
        bb[j]=pb*(j<(J-1))
      }
      ae = rbind(ae,ab)
      be = c(be,bb)
    } else {
      pb = yt[K]/sum(xt[1:(J-1)])
      ab = matrix(0,J,JK+2*IK)
      bb = rep(0,J)
      for (j in 1:J){
        ab[j,j*K] = 1
        bb[j]=pb*(j<J)
      }
      ae = rbind(ae,ab)
      be = c(be,bb)
    }
  }
  # Full scenario. Constraint pjk(J, K)= pik(J-1,K) = 0 related to new young
  # voters and immigrants and pjk(j,K) constant related to exits.
  if (new_and_exit_voters == "full"){
    pb = yt[K]/sum(xt[1:(J-2)])
    ab = matrix(0,J,JK+2*IK)
    bb = rep(0,J)
    for (j in 1:J){
      ab[j,j*K] = 1
      bb[j]=pb*(j<J-1)
    }
    ae = rbind(ae,ab)
    be = c(be,bb)
  }
  # Gold scenario. Constraints pjk(J, K-1) = pik(J,K) = pjk(J-1, K-1) =
  # = pik(J-1,K) =0 related to new young voters and immigrants and
  # pjk(j,K-1) and pjk(j,K) constant related to exits.
  if (new_and_exit_voters == "gold"){
    # Column K-1
    pb = yt[K-1]/sum(xt[1:(J-2)])
    ab = matrix(0,J,JK+2*IK)
    bb = rep(0,J)
    for (j in 1:J){
      ab[j,(K-1)+(j-1)*K] = 1
      bb[j]=pb*(j<(J-1))
    }
    ae = rbind(ae,ab)
    be = c(be,bb)
    # Column K
    pb = yt[K]/sum(xt[1:(J-2)])
    ab = matrix(0,J,JK+2*IK)
    bb = rep(0,J)
    for (j in 1:J){
      ab[j,j*K] = 1
      bb[j]=pb*(j<(J-1))
    }
    ae = rbind(ae,ab)
    be = c(be,bb)
  }
  # Structural zero restrictions introduced by the user
  if (length(structural_zeros) > 0){
    ast = matrix(0, length(structural_zeros), JK+2*IK)
    bst = rep(0, length(structural_zeros))
    for (i in 1:length(structural_zeros)){
      ast[i, K*(structural_zeros[[i]][1]-1)+structural_zeros[[i]][2]] = 1
    }
    ae = rbind(ae,ast)
    be = c(be,bst)
  }
  # Objective function, to minimize
  fun.obj = c(rep(0,JK), rep(1,2*IK))
  # Solution
  sol=suppressWarnings(lpSolve::lp('min', fun.obj, ae, rep('=', length(be)), be))
  z = sol$solution
  pjk = matrix(z[1:JK], J, K, TRUE, dimnames = list(names1,names2))
  e = z[(JK+1):(JK+2*IK)]
  e = matrix(e,2,IK)
  e = e[1,] - e[2,]
  eik = t(matrix(e,K,I))
  colnames(eik) = names2
  rownames(eik) = rownames(x)
  EHet = eik
  pkj = matrix(0, K, J, dimnames=list(names2,names1))
  for (k in 1:K) {
    for (j in 1:J) {
      pkj[k,j]=xt[j]*pjk[j,k]/yt[k]
    }
  }
  pjk.complete <- pjk
  if (new_and_exit_voters %in% c("regular", "raw")){
    if (net_entries & max(x[,ncol(x)]/rowSums(x)) < 0.01){
      pjk = pjk[-J,]; eik = eik[-J,]; pkj = pkj[,-J]
    }
    if (net_exits & max(y[,ncol(y)]/rowSums(y)) < 0.01){
      pjk = pjk[,-K]; eik = eik[,-K]; pkj = pkj[-K,]
    }
  }
  HIe = 100*sum(abs(eik))/tt
  pjk = round(100*pjk,2)
  pkj = round(100*pkj,2)
  if (verbose){
    cat("\n\nEstimated Heterogeneity Index HIe:",round(HIe, 2),"%\n")
    cat("\n\n Matrix of vote transitions (in %) from Election 1 to Election 2\n\n")
    print(pjk)
    cat("\n\n Origin (%) of the votes obtained in Election 2\n\n")
    print(pkj)
  }
  output <- list("VTM" = pjk, "OTM" = pkj, "HETe" = HIe, "EHet" = EHet,
                 "inputs" = inputs, "origin" = x, "destination" = y, "VTM.complete" = pjk.complete)
  class(output) <- "lphom"
  return(output)
}
