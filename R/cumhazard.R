#' @import tibble
#' @import dplyr
#' @import tidyr
NULL

#' Get Cumulative Hazard at a Landmark Timepoint
#' @param eventRates A tibble containing period duration (\code{duration}) and event rate (\code{rate})
#' for specified periods.
#' @param landmark The landmark of interest to evaluate cumulative hazard.
#' @return A numeric which is the cumulative hazard at a landmark timepoint.
#' @examples
#' # Piecewise exponential event rates of 0.5 for time 0-3, 0.4 for time 3-6, and 0.5 after
#' cumhaz <- cumhazard(eventRates=tibble::tibble(duration = c(3,3,100),rate = c(0.5, 0.4, 0.3)),
#'                     landmark=12)
#'
#' @export

cumhazard <- function(eventRates = tibble::tibble(duration=c(3,100),rate=c(log(2)/5,log(2)/5*0.5)),landmark)
{
  finish <- cumsum(eventRates$duration)
  if(landmark>dplyr::last(finish)) {cat("Please enter a time within the enrollment duration \n"); return();}

  n_rates <- nrow(eventRates)
  if (n_rates == 1){
    if(eventRates$rate == 0) cumhz <- 0
    # generate exponential failure time if non-0 failure rate
    else cumhz <- eventRates$rate * landmark
  } else {
    time.int <- c(0,finish)
    time.sorted <- sort(c(time.int,landmark))
    time.int.cumhz <- diff(time.sorted)[1:(min(which(time.sorted==landmark))-1)]
    cumhz <- sum(eventRates$rate[1:length(time.int.cumhz)]*time.int.cumhz)
  }

  return(cumhz)
}
