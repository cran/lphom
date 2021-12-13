#' 	Graphical representation of a RxC ecological inference (vote transfer) matrix  
#'
#' @description Plot method for objects obtained with an algorithm of the lphom-family (lphom, tslphom, nslphom, tslphom_dual, nslphom_joint, ....).
#'
#' @author Jose M. Pavia, \email{pavia@@uv.es}
#'
#' @param x An object output of a **lphom** family algorithm.
#' @param complete A TRUE/FALSE argument informing if the complete matrix should be displayed. In `regular` and `raw` scenarios this matrix includes the row and the column corresponding to net_entries and net_exits even when they are really small, less than 1% in all units. Default, FALSE 
#' @param margins A TRUE/FALSE argument informing if the margins of the matrix should be displayed. Default TRUE.
#' @param digits Integer indicating the number of decimal places to be shown. Default, 2.
#' @param row.names Names to be used for the rows of the matrix.
#' @param col.names Names to be used for the columns of the matrix.
#' @param size.numbers A reference number indicating the average font size to be used for the transfer numbers. Default, 6.   
#' @param size.labels A number indicating the font size to be used for labels. Default, 4.
#' @param size.margins A number indicating the font size to be used for margin numbers. Default, 4. 
#' @param colour.cells Background base colour for cells.
#' @param colour.grid Colour to be used for grid lines. 
#' @param alpha A \[0,1\] number of colour transparency.
#' @param which A vector of integers informing the units for which the aggregate transfer matrix should be plotted. Default, NULL, the global matrix is shown.    
#' @param ... Other arguments passed on to methods. Not currently used.
#' @param type A character string indicating the solution (transfer matrix) to be plotted. Only valid for **_dual** algorithms. `type = "w"` stands for the weighted solution and `type = "a"` for the simple average solution. Default `w`.  
#' @param show.plot A TRUE/FALSE indicating if the plot should be displayed as a side-effect. By default, TRUE.
#' 
#' @return
#' Invisibly returns the (ggplot) description of the plot, which is a list with components that contain the plot itself, the data, information about the scales, panels etc.
#' 
#' @note ggplot2 is needed to be installed for this function to work.
#' 
# @import ggplot2
#'
#' @export
#' @method plot lphom
#' @examples
#' mt.ns <- nslphom(France2017P[, 1:8] , France2017P[, 9:12], new_and_exit_voters= "raw")
#' p <- plot(mt.ns, show.plot = FALSE)
#' p
#' 
plot.lphom <- function(x,
                       complete = FALSE,
                       margins = TRUE,
                       digits = 2,
                       row.names = NULL,
                       col.names = NULL,
                       size.numbers = 6,
                       size.labels = 4,
                       size.margins = 4,
                       colour.cells = "deeppink3",
                       colour.grid = "blanchedalmond",
                       alpha = 0.5,
                       which = NULL,
                       ...,
                       type = "w",
                       show.plot = TRUE){
  
  if (inherits(x, "ei_lp")){
    pintar_con <- plot_ei_lp
  } else if (inherits(x, "ei_dual")){
    pintar_con <- plot_ei_dual 
  } else if (inherits(x, "ei_joint")){
    pintar_con <- plot_ei_joint
  }
  
  p <- pintar_con(x = x, 
                  complete = complete,
                  margins = margins,
                  digits = digits,
                  row.names = row.names,
                  col.names = col.names,
                  size.numbers = size.numbers,
                  size.labels = size.labels,
                  size.margins = size.margins,
                  colour.cells = colour.cells,
                  colour.grid = colour.grid,
                  alpha = alpha,
                  which = which,
                  ...,
                  type = type)
  if (show.plot) print(p)
  return(p)
}

plot_ei_lp  <- function(x,
                        complete = FALSE,
                        margins = TRUE,
                        digits = 2,
                        row.names = NULL,
                        col.names = NULL,
                        size.numbers = 6,
                        size.labels = 4,
                        size.margins = 4,
                        colour.cells = "deeppink3",
                        colour.grid = "blanchedalmond",
                        alpha = 0.5,
                        which = NULL,
                        ...) 
{
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package ggplot2 needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (!complete){
    prop <- x$VTM
    votos <- x$VTM.votes
  } else {
    prop <- x$VTM.complete*100
    votos <- x$VTM.complete.votes
  }
  
  n.fil <- nrow(prop)
  n.col <- ncol(prop)

  if (!is.null(which)){
    if(substr(class(x)[1], 1, 5) == "lphom"){
      stop("Unit matrices are not available in the input object, please set 'which = NULL'")
    }
    if (max(which) > dim(x$VTM.votes.units)[3L] | min(which) < 1 | max(abs(which-round(which))) > 10^-5){
      stop("The 'which' argument that you are using is not valid. Please check it.")
    }
    x$VTM.complete.votes <- apply(x$VTM.votes.units[, , which], c(1, 2), sum)
    votos <- x$VTM.complete.votes[1L:n.fil, 1L:n.col]
    prop <- (x$VTM.complete.votes/rowSums(x$VTM.complete.votes)*100)[1L:n.fil, 1L:n.col]
  }  
  
  
  votos.filas <- round(rowSums(x$VTM.complete.votes)[1L:n.fil])
  votos.columnas <- round(colSums(x$VTM.complete.votes)[1L:n.col])
  
  prop2 <- as.vector(prop)

  ## base de trabajo  
  bbdd <- cbind(expand.grid(n.fil:1L, 1L:n.col), prop2,
                format(round(prop2, digits), n.small = digits))
  bbdd <- as.data.frame(bbdd)
  
  names(bbdd) <- c("y", "x", "coefficient", "label")
  bbdd$color <- paste0("gray", round((100 - round(bbdd$coefficient))/1.5))
  
  ## Tamanyos numeros
  factor.size <- log(votos/sum(votos)*100 + 1L)
  factor.size <- factor.size/max(max(factor.size)) + 0.5
  bbdd$size <- as.vector(factor.size*size.numbers)
  
  ## Se añaden marginales
  if (margins){
    suma.fila <- data.frame(y = n.fil:1L, x = n.col + 1L, coefficient = 0,
                            label = votos.filas, color = "gray27", 
                            size = size.margins)
    suma.columna <- data.frame(y = 0, x = 1L:n.col, coefficient = 0,
                               label = votos.columnas, color = "gray27",
                               size = size.margins)
    bbdd <- rbind(bbdd, suma.fila, suma.columna)
  }
  
  ## Se añaden nombres
  if (is.null(row.names)){
    row.names <- rownames(prop)
  }
  nombres.fila <- data.frame(y = n.fil:1L, x = 0, coefficient = 0,
                             label = row.names, color = "gray27", size = size.labels)
  
  if (is.null(col.names)){
    col.names <- colnames(prop)
  }
  nombres.columna <- data.frame(y = n.fil + 1L, x = 1L:n.col, coefficient = 0,
                                label = col.names, color = "gray27", size = size.labels)
  
  bbdd <- rbind(bbdd, nombres.fila, nombres.columna)
  
  p <- ggplot2::ggplot(bbdd, ggplot2::aes(x = !!quote(x), y = !!quote(y))) +
    ggplot2::geom_tile(ggplot2::aes(fill = !!quote(coefficient)), 
                                    color = colour.grid) + 
    ggplot2::scale_fill_continuous(high = scales::alpha(colour = colour.cells, alpha = alpha), 
                                   low = "white", trans = "sqrt") +
    ggplot2::geom_text(ggplot2::aes(label = !!quote(label)), 
                                     size = bbdd$size, colour = bbdd$color) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.position = "none",
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )
}

plot_ei_dual <- function(x,
                         complete = TRUE,
                         margins = TRUE,
                         digits = 2,
                         row.names = NULL,
                         col.names = NULL,
                         size.numbers = 6,
                         size.labels = 4,
                         size.margins = 4,
                         colour.cells = "deeppink3",
                         colour.grid = "blanchedalmond",
                         alpha = 0.5,
                         which = NULL,
                         ...,
                         type = "w"){
  
  if (type == "w"){
    x$VTM.complete<- x$VTM12.w
    x$VTM.complete.votes <- x$VTM.votes.w
  } else if (type =="a"){
    x$VTM.complete<- x$VTM12.a
    x$VTM.complete.votes <- x$VTM.votes.a
  } else {
    stop("Improper type argument, only 'w' and 'a' are allowed")
  }
  class(x)[2] <- "ei_lp"
  p <- plot_ei_lp(x = x, 
                  margins = margins,
                  digits = digits,
                  row.names = row.names,
                  col.names = col.names,
                  size.numbers = size.numbers,
                  size.labels = size.labels,
                  size.margins = size.margins,
                  colour.cells = colour.cells,
                  colour.grid = colour.grid,
                  alpha = alpha,
                  which = which,
                  ...)
}

plot_ei_joint <- function(x,
                          complete = TRUE,
                          margins = TRUE,
                          digits = 2,
                          row.names = NULL,
                          col.names = NULL,
                          size.numbers = 6,
                          size.labels = 4,
                          size.margins = 4,
                          colour.cells = "deeppink3",
                          colour.grid = "blanchedalmond",
                          alpha = 0.5,
                          which = NULL,
                          ...){
  
    x$VTM.complete<- x$VTM12
    x$VTM.complete.votes <- x$VTM.votes
    class(x)[2] <- "ei_lp"
    p <- plot_ei_lp(x = x, 
                    complete = TRUE,
                    margins = margins,
                    digits = digits,
                    row.names = row.names,
                    col.names = col.names,
                    size.numbers = size.numbers,
                    size.labels = size.labels,
                    size.margins = size.margins,
                    colour.cells = colour.cells,
                    colour.grid = colour.grid,
                    alpha = alpha,
                    which = which,
                    ...)
}

