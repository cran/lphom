#' Confidence Intervals for lphom estimates
#'
#' @description  Estimates confidence intervals for the (vote) transfer probabilities obtained with **lphom()**
#'
#' @author Jose M. Pavia, \email{pavia@@uv.es}
#' @author Rafael Romero \email{rromero@@eio.upv.es}
#' @references Romero, R, Pavia, JM, Martin, J and Romero G (2020). Assessing uncertainty of voter transitions estimated from aggregated data. Application to the 2017 French presidential election. *Journal of Applied Statistics*, 47(13-15), 2711-2736. \doi{10.1080/02664763.2020.1804842}
#' @references Martin, J (2020). Analisis de la incertidumbre en la estimacion de la movilidad electoral mediante el procedimiento plhom. PhD Dissertation.
#'
#' @param lphom.object An object output of the **lphom()** function.
#' @param level A number between 0 and 1 to be used as level of confidence for the intervals. By default 0.90
#' @param num.d Number maximum of different disturbances, `d`, to be initially considered. Positive integer greater than or equal to 5. By default, 11.
#' @param B Integer that determines the number of simulations to be performed for each disturbance value. By default, 30.
#'
#' @return
#' A list with the following components
#'    \item{TM.estimation}{ Transfer matrix of probability point estimates.}
#'    \item{TM.lower}{ Transfer matrix of lower values for the probability estimates.}
#'    \item{TM.upper}{ Transfer matrix of upper values for the probability estimates.}
#'    \item{level}{ Confidence level used when computing the confidence intervals.}

#' @export
#'
#' @seealso \code{\link{lphom}} \code{\link{error_lphom}}
#'
#' @examples
#' mt.lphom <- lphom(France2017P[, 1:8], France2017P[, 9:12], "raw", NULL, FALSE)
#' set.seed(533423)
#' confidence_intervals_pjk(mt.lphom, level = 0.90, num.d = 5, B = 8)
#
# @importFrom stats qnorm
#
confidence_intervals_pjk <- function(lphom.object,
                                     level = 0.90,
                                     num.d = 11,
                                     B = 30) {
  if(class(lphom.object) != "lphom"){
    stop("'lphom.object' must be output from 'lphom'")
  }
  if (level > 1 | level < 0)
    stop('level must be between 0 and 1')
  num.d <- as.integer(num.d); B <- as.integer(B)
  if (num.d < 5) stop('num.d must be a integer greater than 5')
  if (B < 1) stop('B must be a integer greater than 0')
  num.q <- (num.d + 1)/2
  d0.ref <- hallar_d0(lphom.object)
  id = d0.ref[1] + d0.ref[2]*((1:num.d)-num.q)
  while (id[1] < 0) {
    d0.ref[2] <- d0.ref[2]/2
    id = d0.ref[1] + d0.ref[2]*((1:num.d)-num.q)
  }
  id <- seq(id[1], id[1] + 0.25, d0.ref[2]) # Expansión de id
  escenarios <- simular_escenarios(lphom.object = lphom.object, M.d = id, B = B)
  J <- nrow(lphom.object$VTM.complete)
  K <- ncol(lphom.object$VTM.complete)
  ceros <- determinar_zeros_estructurales(lphom.object)
  MT.upper <- MT.lower <- array(NA, dim(lphom.object$VTM.complete))
  zeta <- stats::qnorm((1-level)/2)
  if (!is.null(ceros)){
    for (i in 1:length(ceros)){
      MT.upper[ceros[[i]][1], ceros[[i]][2]] <- MT.lower[ceros[[i]][1], ceros[[i]][2]] <- 0
    }
  }
  for (j in 1:J){
    for (k in 1:K){
      if (is.na(MT.upper[j,k])){
        sesgos <- simular_errores_estimacion_pjk(escenarios, j, k)
        ses.var <- modelos_ajuste_estimacion_pjk(sesgos, lphom.object$HETe)
        li <- ses.var[1] + zeta*ses.var[2]
        MT.upper[j,k] <- max(min(1, lphom.object$VTM.complete[j,k] - li), lphom.object$VTM.complete[j,k])
        ls <- ses.var[1] - zeta*ses.var[2]
        MT.lower[j,k] <- min(max(0, lphom.object$VTM.complete[j,k] - ls), lphom.object$VTM.complete[j,k])
      }
    }
  }
  dimnames(MT.upper) <- dimnames(MT.lower) <- dimnames(lphom.object$VTM.complete)
  output <- list("TM.estimation"  = lphom.object$VTM.complete,
                 "TM.lower" = MT.lower, "TM.upper" = MT.upper, "level" = level)
  return(output)
}
