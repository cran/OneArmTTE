#' @import tibble
#' @import dplyr
#' @import tidyr
#' @import survival
NULL

#' Get operating characteristics of one-arm clinical trial design with time-to-event endpoint
#'
#' Using simulation, this function can get operating characterisitics of several approaches for
#' one-arm trial design with time-to-event endpoint. Default approaches include one-sample log-rank
#' test, Landmark Kaplen-Meier method and binary method which regards the survival of each subject
#' at a landmark is a binary variable. In addition, if \code{RWdata} is not \code{NULL}, the \code{RWdata}
#' input will be used as an external control and cox model will be used to evaluate the treatment effect
#' of simulated data (experimental arm) compared with the external control. The output includes
#' probability of rejecting null hypothesis of each design, average number of events at analysis,
#' and average analysis time after last patient in. When \code{eventRates} is same as \code{eventRates.ctrl},
#' the probability of rejecting null hypothesis is type I error; When \code{eventRates} is the alternative
#' hypothesis from desirable treatment effect, the probability of rejecting null hypothesis is power.
#'
#' @param n Number of subjects.
#' @param eventRates.ctrl Event rates of historical control.
#' @param eventRates Event rates of subjects in the trial.
#' @param enrollRates Enrollment rates of subjects in the trial.
#' @param dropoutRates Dropout rates of the subjects in the trial.
#' @param cutTime Analysis time after last patient in; not used if Event=TRUE.
#' @param landmark The landmark of interest to evaluate the survival rate for Landmark Kaplan-Meier method.
#' @param Event Indicator of whether the analysis is driven by number of events; default is FALSE.
#' @param n.event Number of events at analysis; not used if Event=FALSE.
#' @param RWdata The real world data to be used as external control; A tibble/data.frame containing \code{time} and \code{censor}, where \code{censor=1} indicates event; default is NULL.
#' @param RWSurvCal Indicator of whether to calculate historical cumulative hazard and survival rate at landmark from real world data as the null case; default is FALSE.
#' @param conf.type Type of confidence interval in the survival model; One of "\code{none}", "\code{plain}" (the default), "\code{log}", "\code{log-log}", "\code{logit}" or "\code{arcsin}".
#' @param alpha Type I error rate level.
#' @param nsim Number of simulations; default is 10000.
#' @param seed Seed for simulation.
#'
#' @return No visible return values.
#'
#' @details The function output a list of the operating characteristics including:
#'          1) probability of rejecting null hypothesis of each design,
#'           2) average number of events at analysis,
#'         3) average analysis time after last patient in.
#' @examples
#' \donttest{
#' library(survival)
#' # Piecewise exponential of historical control
#' median.ctrl <- c(14.3, 1.5, 4.9)
#' eventRates.ctrl <- tibble::tibble(duration=c(4,2,100),rate=log(2)/median.ctrl)
#' # Piecewise exponential assumption of treatment:
#' # Hazard ratio = 1 for time 0-3 and Hazard ratio = 0.47 after
#' eventRates.trt = tibble::tibble(duration=c(3,1,2,100),rate=log(2)/c(14.3, median.ctrl/0.47))
#' # Constant enrollment rates and dropout rates
#' enrollRates = tibble::tibble(duration=106, rate=14/3)
#' dropoutRates = tibble::tibble(duration=106, rate=0.2/12)
#' OneArmTTEDesign(n=40, eventRates.ctrl, eventRates.trt, enrollRates, dropoutRates, cutTime=3,
#'                 landmark=6, Event=FALSE, conf.type = 'plain', alpha=0.05, nsim=100, seed=43)
#' }
#' @export


OneArmTTEDesign <- function(n, eventRates.ctrl, eventRates, enrollRates, dropoutRates, cutTime, landmark, Event = FALSE, n.event, RWdata = NULL, RWSurvCal = FALSE, conf.type = 'plain', alpha = 0.05, nsim = 10000, seed = 43)
{
  if(!is.null(RWdata)){
    if(sum(c("time","censor") %in% names(RWdata))!=2) stop("Real world data should be a tibble/data.frame containing the column names 'time' and 'censor'")
  }

  if(!is.null(RWdata) & RWSurvCal==TRUE) {
    fit.hist <- survfit(Surv(time,censor)~1, data=RWdata)
    s.null <- summary(fit.hist,times = landmark)$surv
  } else s.null <- exp(-cumhazard(eventRates.ctrl, landmark))

  rej.lr1 = rej.lr2 = rej.lm = rej.naive = rej.rwd = 0
  nevent = tcut = tcut.lpi = rep(0,nsim)
  err=0

  set.seed(seed)

  for(sim in 1:nsim){
    simdata = tibble::tibble(enrollTime=rpwenroll(n,enrollRates),failTime=rpwexp(n,eventRates),dropoutTime=rpwexp(n,dropoutRates)) %>%
      mutate(cte=pmin(dropoutTime,failTime)+enrollTime, fail=(failTime <= dropoutTime)*1)

    if(!Event) {
      simdata <- simdata %>% mutate(censor=ifelse(cte<=max(enrollTime)+cutTime,fail,0)) %>%
        mutate(time=ifelse(censor==1,failTime,pmin(cte-enrollTime,max(enrollTime)+cutTime-enrollTime)))
      t.cut <- max(simdata$enrollTime)+cutTime
    } else {
      t.cut <- sort(simdata$cte[simdata$fail==1])[min(sum(simdata$fail),n.event)]
      simdata <- simdata %>% mutate(censor=ifelse(cte<=t.cut,fail,0)) %>%
        mutate(time=ifelse(censor==1,failTime,pmin(cte-enrollTime,t.cut-enrollTime)))

    }

    if(sum(simdata$censor)==0) {sim=sim-1; next} else
    {
      O=sum(simdata$censor)
      E = 0
      if(!is.null(RWdata) & RWSurvCal==TRUE) {
        for(i in 1:length(simdata$time)) E = E + ifelse(simdata$time[i]<max(fit.hist$time), -log(summary(fit.hist,times = simdata$time[i])$surv), cumhazard(eventRates.ctrl,simdata$time[i]))
      } else {
        for(i in 1:length(simdata$time)) E = E + cumhazard(eventRates.ctrl,simdata$time[i])
      }
      L1=(O-E)/sqrt(E)
      L2=(O-E)/sqrt((O+E)/2)
      if(L1<(-stats::qnorm(1-alpha))) rej.lr1 = rej.lr1+1
      if(L2<(-stats::qnorm(1-alpha))) rej.lr2 = rej.lr2+1

      fit <- survfit(Surv(time,censor)~1, data=simdata, conf.type = conf.type, conf.int=1-2*alpha)
      summ <- summary(fit)
      if(max(summ$time)<landmark) {lower = min(summ$lower[!is.na(summ$lower)]);err=err+1}
      else lower = summary(fit,times = landmark)$lower
      if(lower>s.null) {
        rej.lm = rej.lm + 1
        if(max(fit$time)<landmark) pest = min(summ$surv[!is.na(summ$surv)])
        else pest = summary(fit,times = landmark)$surv
      }
      nevent[sim] = O
      tcut[sim] = t.cut
      tcut.lpi[sim] = t.cut-max(simdata$enrollTime)

      simdata <- simdata %>% mutate(naive_flag = ifelse(t.cut-enrollTime<landmark-1e-7 | (censor==0 & time<landmark),0,1)) %>% mutate(censor2 = ifelse(naive_flag==1 & time>landmark, 0, censor))
      bound <- stats::qbinom(1-alpha,sum(simdata$naive_flag),s.null)
      ne <- simdata %>% filter(naive_flag==1) %>% .$censor2 %>% sum()
      if(ne<sum(simdata$naive_flag)-bound) rej.naive=rej.naive+1

      if(!is.null(RWdata)) {

        exp.data <- simdata %>% select(time, censor) %>% mutate(Trt="Experimental")
        ana.data <- dplyr::bind_rows(RWdata %>% mutate(Trt="Control"), exp.data)

        survfit.coxph <- coxph(Surv(time,censor)~Trt, data=ana.data)
        survfit.coxph_result <- summary(survfit.coxph)

        HR <- as.numeric(survfit.coxph_result$conf.int[1, 'exp(coef)'])

        logrank_fit <- survdiff(Surv(time,censor)~Trt, data=ana.data)

        pval <- 1 - stats::pchisq(logrank_fit$chisq, df = 1)

        if(HR>1){
          Pval_test <- 1- pval/2
        }else{
          Pval_test <- pval/2
        }

        if(Pval_test<alpha) rej.rwd = rej.rwd + 1

      }
    }


  }

  cat("Probability of rejecting the null hypothesis (One-Sample Log-rank Test, %):", formatC(rej.lr1/nsim*100, digits=1, format="f"), "\n");

  cat("Probability of rejecting the null hypothesis (Landmark Kaplen-Meier Method, %):", formatC(rej.lm/nsim*100, digits=1, format="f"), "\n");

  cat("Probability of rejecting the null hypothesis (Binary Method, %):", formatC(rej.naive/nsim*100, digits=1, format="f"), "\n");

  if(!is.null(RWdata)) cat("Probability of rejecting the null hypothesis (Compare with Real World Data, %):", formatC(rej.rwd/nsim*100, digits=1, format="f"), "\n");

  cat("Average number of events at analysis:", formatC(mean(nevent), digits=1, format="f"), "\n");

  cat("Average analysis time after last patient in:", formatC(mean(tcut.lpi), digits=1, format="f"), "\n");

  invisible(list(rej.lr1=rej.lr1/nsim, rej.lr2=rej.lr2/nsim, rej.lm=rej.lm/nsim, rej.naive=rej.naive/nsim, rej.rwd/nsim, n.event=nevent/nsim, t.ana=mean(tcut), t.ana.lpi=mean(tcut.lpi), err=err/nsim))
}
