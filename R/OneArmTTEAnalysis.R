#' @import tibble
#' @import dplyr
#' @import tidyr
#' @import survival
NULL

#' Perform analysis on the data of one-arm clinical trial with time-to-event endpoint
#'
#' This function can get analysis results on the input trial data using several approaches for
#' one-arm trial design with time-to-event endpoint. Default approaches include one-sample log-rank
#' test, Landmark Kaplen-Meier method and binary method which regards the survival of each subject
#' at a landmark is a binary variable. In addition, if \code{RWdata} is not \code{NULL}, the \code{RWdata}
#' input will be used as an external control and cox model will be used to evaluate the treatment effect
#' of input trial data (experimental arm) compared with the external control.

#' @param data Trial data. A tibble/data.frame containing \code{time} and \code{censor}, where \code{censor=1} indicates event.
#' @param eventRates.ctrl Event rates of historical control.
#' @param landmark The landmark of interest to evaluate the survival rate for Landmark Kaplan-Meier method and binary method.
#' @param RWdata The real world data to be used as external control; A tibble/data.frame containing \code{time} and \code{censor}, where \code{censor=1} indicates event; default is NULL.
#' @param RWSurvCal Indicator of whether to calculate historical cumulative hazard and survival rate at landmark from real world data; default is FALSE.
#' @param conf.type Type of confidence interval in the survival model; One of "\code{none}", "\code{plain}" (the default), "\code{log}", "\code{log-log}", "\code{logit}" or "\code{arcsin}".
#' @param alpha Type I error rate level.
#'
#' @return No visible return values.
#'
#' @details This function outputs a list of analysis results of each design, including:
#'          1) p-value of one-sample log-rank test,
#'          2) historical survival rate at landmark,
#'          3) survival rate estimate with confidence interval of landmark kaplan-meier method,
#'          4) survival rate estimate with confidence interval of binary method,
#'          5) p-value of binary method,
#'          6) hazard ratio estimate with confidence interval compared with real world data (if available),
#'          7) p-value of log-rank test compared with real world data (if available).
#'
#' @examples
#' \donttest{
#' library(survival)
#' data(example_data)
#' # Piecewise exponential of historical control
#' median.ctrl <- c(14.3, 1.5, 4.9)
#' eventRates.ctrl <- tibble::tibble(duration=c(4,2,100),rate=log(2)/median.ctrl)
#' OneArmTTEAnalysis(example_data, eventRates.ctrl, landmark=6)
#' }
#' @export


OneArmTTEAnalysis <- function(data, eventRates.ctrl, landmark, RWdata = NULL, RWSurvCal = FALSE, conf.type = 'plain', alpha = 0.05)
{
  if(sum(c("time","censor") %in% names(data))!=2) stop("Input trial data should be a tibble/data.frame containing the column names 'time' and 'censor'")
  if(!is.null(RWdata)){
    if(sum(c("time","censor") %in% names(RWdata))!=2) stop("Real world data should be a tibble/data.frame containing the column names 'time' and 'censor'")
  }

  rej.lr1 = rej.lm = rej.naive = rej.rwd = err = 0
  O=sum(data$censor)
  E = 0
  if(!is.null(RWdata) & RWSurvCal==TRUE) {
    fit.hist <- survfit(Surv(time,censor)~1, data=RWdata)
    s.null <- summary(fit.hist,times = landmark)$surv
    for(i in 1:length(data$time)) E = E + ifelse(data$time[i]<max(fit.hist$time), -log(summary(fit.hist,times = data$time[i])$surv), cumhazard(eventRates.ctrl,data$time[i]))
  } else {
    s.null <- exp(-cumhazard(eventRates.ctrl, landmark))
    for(i in 1:length(data$time)) E = E + cumhazard(eventRates.ctrl,data$time[i])
  }
  L1=(O-E)/sqrt(E)
  p.lr = stats::pnorm(L1)
  #L2=(O-E)/sqrt((O+E)/2)
  if(L1<(-stats::qnorm(1-alpha))) rej.lr1 = 1
  #if(L2<(-stats::qnorm(1-alpha))) rej.lr2 = rej.lr2+1

  fit <- survfit(Surv(time,censor)~1, data=data, conf.type = conf.type, conf.int=1-2*alpha)
  summ <- summary(fit)
  if(max(summ$time)<landmark) {
    lower = min(summ$lower[!is.na(summ$lower)])
    est = min(summ$surv[!is.na(summ$lower)])
    upper = min(summ$upper[!is.na(summ$lower)])
    err=err+1
    } else {
      lower = summary(fit,times = landmark)$lower
      est = summary(fit,times = landmark)$surv
      upper = summary(fit,times = landmark)$upper
    }
  if(lower>s.null) rej.lm = 1

  data <- data %>% mutate(naive_flag = ifelse(censor==0 & time<landmark,0,1)) %>% mutate(censor2 = ifelse(naive_flag==1 & time>landmark, 0, censor))
  bound <- stats::qbinom(1-alpha,sum(data$naive_flag),s.null)
  ne <- data %>% filter(naive_flag==1) %>% .$censor2 %>% sum()
  if(ne<sum(data$naive_flag)-bound) rej.naive = 1
  bin.test <- stats::binom.test(sum(data$naive_flag)-ne,sum(data$naive_flag),p = s.null, alternative = "greater", conf.level = 1-2*alpha)
  bin.est <- bin.test$estimate
  bin.lower <- bin.test$conf.int[1]
  bin.upper <- bin.test$conf.int[2]

  if(!is.null(RWdata)) {
    if(sum(c("time","censor","Trt") %in% names(RWdata))!=3) stop("Real world data should be a tibble/data.frame containing the column names 'time', 'censor' and 'Trt'")

    exp.data <- data %>% select(time, censor) %>% mutate(Trt="Experimental")
    ana.data <- dplyr::bind_rows(RWdata, exp.data)

    survfit.coxph <- coxph(Surv(time,censor)~Trt, data=ana.data)
    survfit.coxph_result <- summary(survfit.coxph)

    HR <- as.numeric(survfit.coxph_result$conf.int[1, 'exp(coef)'])

    HR_display <- paste0(
      sprintf(paste0("%.",2,"f"),
              round(HR, 2)
      ),
      " (",
      paste0(sprintf(paste0("%.",2,"f"),
                     round(survfit.coxph_result$conf.int[1,
                                                         c("lower .95", "upper .95")],
                           2)
      ), collapse = ", "
      ),")"
    )

    logrank_fit <- survdiff(Surv(time,censor)~Trt, data=ana.data)

    pval <- 1 - stats::pchisq(logrank_fit$chisq, df = 1)

    if(HR>1){
      Pval_test <- 1- pval/2
    }else{
      Pval_test <- pval/2
    }

    if(Pval_test<alpha) rej.rwd = 1

  }

  cat("One-sample log-rank test p-value:", formatC(p.lr, digits=4, format="f"), "\n");

  cat("Historical survival rate at landmark:", formatC(s.null, digits=3, format="f"), "\n");

  cat("Survival rate estimate (CI) of Landmark Kaplan-Meier Method:", paste0(round(est,3), " (", round(lower,3), ", ", round(upper,3), ")"), "\n")

  cat("Survival rate estimate (CI) of Binary Method:", paste0(round(bin.est,3), " (", round(bin.lower,3), ", ", round(bin.upper,3), ")"), "\n")

  cat("Binary method test p-value:", formatC(bin.test$p.value, digits=4, format="f"), "\n");

  if(!is.null(RWdata)) {
    cat("Hazard ratio estimate (CI) compared with real world data:", HR_display, "\n");
    cat("Log-rank test p-value compared with real world data:", formatC(Pval_test, digits=4, format="f"), "\n");
  }

}


