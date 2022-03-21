#' 	Print a summary of a lphom-family object
#'
#' @description Print method for objects obtained with an algorithm of the lphom-family (lphom, tslphom, nslphom, tslphom_dual, nslphom_joint, ....).
#'
#' @author Jose M. Pavia, \email{pavia@@uv.es}
#'
#' @param x An object output of a **lphom** family algorithm.
#' @param ... Other arguments passed on to methods. Not currently used.
#' @param margins A TRUE/FALSE argument informing if the margins of the matrix should be displayed. Default TRUE.
#' @param digits Integer indicating the number of decimal places to be shown. Default, 2.
#' 
#'
#' @export
#' @method print lphom
#' @examples
#' mt.ns <- nslphom(France2017P[, 1:8] , France2017P[, 9:12], new_and_exit_voters= "raw")
#' print(mt.ns, digits = 2, margins = TRUE)
#'

#' @export
print.lphom  <- function(x, 
                         ...,
                         margins = TRUE,
                         digits = 2) 
{

  print.summary.lphom(x = summary.lphom(object = x),
                      margins = margins,
                      digits = digits)

}

#' 	Print a summary of a lphom-family object
#'
#' @description Print method for `summary.lphom` objects
#' @inheritParams print.lphom
#' @param x An `summary.lphom` class object.
#' @method print summary.lphom
#' @export
print.summary.lphom  <- function(x, 
                         ...,
                         margins = TRUE,
                         digits = 2) 
{
 
  if (length(x) == 5){ 
    tabla <- format(round(x$prop.matrix, digits), nsmall = digits)
    tabla <- apply(tabla, 2, as.character)
    rownames(tabla) <- rownames(x$prop.matrix)
  
    if (margins){
      nr <- nrow(tabla)
      tabla <- rbind(tabla, x$col.margins[1L:ncol(tabla)])
      tabla <- cbind(tabla,  c(format(x$row.margins[1L:nr], justify = "right"), ""))
#      colnames(tabla)[ncol(tabla)] <- " "
    }
    
    cat("Estimated global row-standardized transfer matrix \n")
    print(as.table(tabla))
    cat("Estimated Heterogeneity index (HETe):", x$HETe, "\n")
  
  } else if (length(x) == 4){ 
    tabla <- format(round(x$prop.matrix, digits), nsmall = digits)
    tabla <- apply(tabla, 2, as.character)
    rownames(tabla) <- rownames(x$prop.matrix)
    
    if (margins){
      nr <- nrow(tabla)
      tabla <- rbind(tabla, x$col.margins[1L:ncol(tabla)])
      tabla <- cbind(tabla,  c(format(x$row.margins[1L:nr], justify = "right"), ""))
    }
    
    cat("Adjusted row-standardized transfer matrix \n")
    print(as.table(tabla))
    
  } else {

    tabla.w <- format(round(x$prop.matrix.w, digits), nsmall = digits)
    tabla.w <- apply(tabla.w, 2, as.character)
    tabla.a <- format(round(x$prop.matrix.a, digits), nsmall = digits)
    tabla.a <- apply(tabla.a, 2, as.character)
    rownames(tabla.w) <- rownames(tabla.a) <- rownames(x$prop.matrix.w)
    
    if (margins){
      tabla.w <- rbind(tabla.w, x$col.margins)
      tabla.a <- rbind(tabla.a, x$col.margins)
      tabla.w <- cbind(tabla.w, c(format(x$row.margins, justify = "right"), ""))
      tabla.a <- cbind(tabla.a, c(format(x$row.margins, justify = "right"), ""))
    }
    
    cat("Estimated global row-standardized transfer matrix\n      (weigthed average of dual solutions)\n")
    print(as.table(tabla.w))
    cat("Estimated Heterogeneity index (HETe):", x$HETe.w, "\n")
    
    cat("\nEstimated global row-standardized transfer matrix\n        (simple average of dual solutions) \n")
    print(as.table(tabla.a))
    cat("Estimated Heterogeneity index (HETe):", x$HETe.a, "\n")
    
  }
  
}
