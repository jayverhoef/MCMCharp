#-------------------------------------------------------------------------------
#
#          logit
#
#-------------------------------------------------------------------------------

#' logit function
#'
#' logit function
#'
#' @param p a value between 0 and 1
#'
#' @return the logit of p
#'
#' @author Jay Ver Hoef
#' @export
  logit = function(p){log(p/(1 - p))}
  
#-------------------------------------------------------------------------------
#
#          expit
#
#-------------------------------------------------------------------------------

#' inverse logit function
#'
#' inverse logit function
#'
#' @param logp a real-valued number
#'
#' @return the inverse logit of logp, which will be between 0 and 1
#'
#' @author Jay Ver Hoef
#' @export
  expit = function(logp){exp(logp)/(1 + exp(logp))}
