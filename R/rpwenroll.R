#' @import tibble
#' @import dplyr
#' @import tidyr
NULL

#' Generate Piecewise Exponential Enrollment
#'
#'
#' @param n Number of observations.
#' Default of \code{NULL} yields random enrollment size.
#' @param enrollRates A tibble containing period duration (\code{duration}) and enrollment rate (\code{rate})
#' for specified enrollment periods.
#' If necessary, last period will be extended to ensure enrollment of specified \code{n}.
#' @return A vector of random enrollment times following piecewise exponential distribution.
#' @examples
#' # piecewise uniform (piecewise exponential inter-arrival times) for 10k patients enrollment
#' # enrollment rates of 5 for time 0-100, 15 for 100-300, and 30 thereafter
#' x <- rpwenroll(n=10000, enrollRates=tibble::tibble(rate = c(5, 15, 30), duration = c(100,200,100)))
#' plot(x,1:10000,
#'      main="Piecewise uniform enrollment simulation",xlab="Time",
#'      ylab="Enrollment")
#' # exponential enrollment
#' x <- rpwenroll(10000, enrollRates=tibble::tibble(rate = .03, duration = 1))
#' plot(x,1:10000,main="Simulated exponential inter-arrival times",
#'      xlab="Time",ylab="Enrollment")
#'
#' @export

rpwenroll <- function(n = NULL,
                      enrollRates = tibble(duration=c(1,2),rate=c(2,5))
){# take care of the simple case first
  if(nrow(enrollRates)==1) {
    # stop with error message if only 1 enrollment period and the enrollment rate is less or equal with 0
    if (enrollRates$rate <=0) stop("Please specify > 0 enrollment rate, otherwise enrollment cannot finish.")
    # otherwise, return inter-arrival exponential times
    else return(cumsum(stats::rexp(n=n,rate=enrollRates$rate)))
  }
  y <- enrollRates %>%
    dplyr::mutate(period=dplyr::row_number(),
                  finish=cumsum(duration),
                  lambda=duration*rate,
                  origin=dplyr::lag(finish,default=0)) %>%
    group_by(period) %>%
    dplyr::mutate(N=stats::rpois(n=1,lambda=lambda))
  # deal with extreme cases where none randomized in fixed intervals
  if (sum(y$N)==0){
    if (is.null(n)) return(NULL)
    # stop with error message if enrollment has not finished but enrollment rate for last period is less or equal with 0
    if (dplyr::last(enrollRates$rate) <=0) stop("Please specify > 0 enrollment rate for the last period; otherwise enrollment cannot finish.")
    # otherwise, return inter-arrival exponential times
    else return(cumsum(stats::rexp(n,rate=dplyr::last(enrollRates$rate)))+dplyr::last(y$finish))
  }
  # generate sorted uniform observations for Poisson count for each interval
  z <- tidyr::expand(y,enrollTime=sort(stats::runif(n=N,min=origin,max=finish)))
  if (is.null(n)) return(z$enrollTime) # if n not specified, return generated times
  if (nrow(z) >= n) return(z$enrollTime[1:n]) # if n already achieved, return first n observations
  # after specified finite intervals, add required additional observations with
  # exponential inter-arrival times
  nadd <- n-nrow(z)
  # stop with error message if enrollment has not finished but enrollment rate for last period is less or equal with 0
  if (dplyr::last(enrollRates$rate) <=0) stop("Please specify > 0 enrollment rate for the last period; otherwise enrollment cannot finish.")
  # Otherwise, return inter-arrival exponential times
  else return(c(z$enrollTime, cumsum(stats::rexp(nadd,rate=dplyr::last(enrollRates$rate)))+dplyr::last(y$finish)))
}
