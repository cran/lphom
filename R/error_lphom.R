#' 	Global error of a lphom estimated table
#'
#' @description  Estimation of the error index (EI) of a RxC vote transfer matrix obtained with **lphom()**
#'
#' @author Jose M. Pavia, \email{pavia@@uv.es}
#' @author Rafael Romero \email{rromero@@eio.upv.es}
#' @references Romero, R, Pavia, JM, Martin, J and Romero G (2020). Assessing uncertainty of voter transitions estimated from aggregated data. Application to the 2017 French presidential election. *Journal of Applied Statistics*, 47(13-15), 2711-2736. \doi{10.1080/02664763.2020.1804842}
#'
#'
#' @param lphom.object An object output of the **lphom()** function.
#' @param upper.alfa Upper bound that will not exceed by the EI estimate with a confidence 1 - alpha. By default, 0.10.
#' @param show.plot TRUE/FALSE. Indicates whether the graphical representation describing the relationship between EI and HETe estimated by simulation for the election under study should be displayed as a side-effect. By default, TRUE.
#' @param num.d Number maximum of different disturbances, `d`, to be initially considered. Positive integer greater than or equal to 5. By default, 11.
#' @param B Integer that determines the number of simulations to be performed for each disturbance value. By default, 30.
#'
#' @return
#' A list with the following components
#'    \item{EI.estimate}{ Point estimate for EI.}
#'    \item{EI.upper}{ Upper bound with confidence 1 - alpha of the EI estimate }
#'    \item{figure}{ ggplot2 object describing the graphical representation of the relationship between EI and HETe.}
#'    \item{equation}{ lm object of the adjustment between EI and HETe.}
#'    \item{statistics}{ A four column matrix with the values of HET, HETe, EI and d associated with each simulated scenario.}
#'    \item{TMs.real}{ Array with the simulated real transfer matrices associated with each scenario.}
#'    \item{TMs.estimate}{ Array with the estimated transfer matrices associated with each scenario.}
#'
#' @note ggplot2 is needed to be installed for this function to work.
#' @export
# @import ggplot2
#'
#' @seealso \code{\link{lphom}} \code{\link{confidence_intervals_pjk}}
#'
#' @examples
#' mt.lphom <- lphom(France2017P[, 1:8] , France2017P[, 9:12], new_and_exit_voters= "raw",
#'                   verbose = FALSE)
#' set.seed(253443)
#' example <- error_lphom(mt.lphom, show.plot = FALSE, num.d = 5, B = 8)
#' example$EI.estimate
#
#' @importFrom stats lm
#' @importFrom stats predict
#
error_lphom <- function(lphom.object, upper.alfa = 0.10,
                        show.plot = TRUE, num.d = 11, B = 30) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package ggplot2 needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if(class(lphom.object) != "lphom"){
    stop("'lphom.object' must be output from 'lphom'")
  }
  HETe <- EI <- medias <- maximos <- NULL
  if (upper.alfa > 0.5 | upper.alfa < 0)
    stop('upper.alfa must be between 0 and 0.5')
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
  log.HETe <- log(escenarios$estadisticos[,2])
  log.HETe2 <- log.HETe^2
  log.EI <- log(escenarios$estadisticos[,3])
  bbdd <- data.frame(log.EI, log.HETe, log.HETe2)
  equation <- stats::lm(log.EI ~ log.HETe + log.HETe2, data = bbdd)
  EI.pred <- exp(stats::predict(equation, newdata = bbdd, interval = "prediction",
                         level = 1 - upper.alfa))
  EI.estimation <- data.frame(log.HETe = log(lphom.object$HETe),
                              log.HETe2 = log(lphom.object$HETe)^2)
  EI.estimation <- exp(stats::predict(equation, newdata = EI.estimation,
                               interval = "prediction", level = 1 - upper.alfa))
  bbdd1 <- cbind(bbdd, medias = log(EI.pred[, 1]), maximos = log(EI.pred[, 3]))
  theme_own <- ggplot2::theme_bw() + ggplot2::theme(aspect.ratio = 1)
  #  p1 <- ggplot2::ggplot(bbdd1, aes(log.HETe, log.EI)) +
  #    geom_point(color = "blue", shape = 4) +
  #    geom_line(aes(y = medias), size = 1.5) +
  #    geom_line(aes(y = maximos), size = 1.5, color = "red", linetype="dashed") +
  #    geom_vline(xintercept=log(lphom.object$HETe), linetype="dotted", color = "green") +
  #    theme_own
  bbdd2 <- exp(bbdd1)
  names(bbdd2)[1:2] <- c("EI", "HETe")
  p2 <- ggplot2::ggplot(bbdd2, ggplot2::aes(HETe, EI)) +
    ggplot2::geom_point(color = "blue", shape = 4) +
    ggplot2::geom_line(ggplot2::aes(y = medias), size = 1.5) +
    ggplot2::geom_line(ggplot2::aes(y = maximos), size = 1.5, color = "red", linetype="dashed") +
    ggplot2::scale_y_log10("EI (axis distances in log-scale)") +
    ggplot2::scale_x_log10("HETe (axis distances in log-scale)") +
    #    geom_vline(xintercept = lphom.object$HETe, linetype="dotted", color = "green", size = 1.2) +
    #    geom_hline(yintercept = EI.estimation[c(1,3)], linetype="dotted", color = "green", size = 1.2) +
    ggplot2::geom_segment(ggplot2::aes(x = lphom.object$HETe, xend = lphom.object$HETe, y = 0,
                              yend = EI.estimation[1]), linetype="dotted", color = "green", size = 1.1) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = lphom.object$HETe, y = EI.estimation[1],
                              yend = EI.estimation[1]), linetype="dotted", color = "green", size = 1.1) +
    ggplot2::geom_segment(ggplot2::aes(x = lphom.object$HETe, xend = lphom.object$HETe, y = EI.estimation[1],
                              yend = EI.estimation[3]), linetype="dotted", color = "yellow", size = 1) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = lphom.object$HETe, y = EI.estimation[3],
                              yend = EI.estimation[3]), linetype="dotted", color = "yellow", size = 1) +
    theme_own
  if (show.plot) suppressWarnings(print(p2))
  output <-list("EI.estimate" = round(EI.estimation[1], 2),
                "EI.upper" = round(EI.estimation[3], 2),
                "figure" = p2, "equation" = equation,
                "statistics" = escenarios$estadisticos,
                "TMs.real" = escenarios$MT.reales, "TMs.estimates" = escenarios$MT.estimadas)
  return(output)
}
