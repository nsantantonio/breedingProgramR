#' getTrueQTLeff function
#'
#' function to (do something)
#'
#' @param simParam [value]
#' @param trait [value]. Default is 1
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
getTrueQTLeff <- function(simParam, trait = 1) { simParam$traits[[trait]]@addEff }
