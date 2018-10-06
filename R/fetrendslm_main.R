## c:/Dropbox/Rpackages/fetrendslm/R/fetrendslm_main.R

##    Chandler Lutz
##    Questions/comments: cl.eco@cbs.dk
##    $Revisions:      1.0.0     $Date:  2018-10-04

#' Fixed Effects linear models with trends
#'
#' @import data.table Matrix
"_PACKAGE"

##For magrittr dot. From https://github.com/tidyverse/magrittr/issues/29
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))

#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

#' @importFrom R6 R6Class

