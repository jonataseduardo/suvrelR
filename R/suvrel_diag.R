
#' Estimate the metric that maximize the intergroup distance while minimizes
#' the intragroup distance
#'
#' \code{suvrel_diag()} Filter SNPs with MAF bellow a threshold
#
#' @import data.table
#' @export 
#'
#' @param datum data.table. Data with one allele per line. 
#'
#' @return data.table metric 

suvrel_diag <- 
  function(datum, value_var, group, feature){

     mean_dt <- 
       datum[, .(m = mean(get(value_var))), 
             by = c(group, feature)]

     K <- mean_dt[, .GRP, by = c(group)][, .N]

     e1 <- datum[, .(e_1 = var(get(value_var))), by = c(feature)]
     e2 <- mean_dt[, .(e_2 = (K - 1) * var(m)), by = c(feature)]
     erg <- e1[e2, on = c(feature)] 
     erg[, e := - e_1 - e_2]
     cte <- 1 / sqrt(erg[, sum(e ^ 2, na.rm = TRUE)])
     erg[, g := - e * cte]
     return(erg[!is.na(g), c(feature, 'g'), with = FALSE])
  }


