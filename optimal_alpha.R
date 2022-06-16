#Optimal alpha t-test R code version 1.1, updated for compatibility with R version 3.x.
#Authored by Joe Mudge (questions or comments? contact: joe.mudge83@gmail.com).

library(ggplot2)

#' Calculate a beta t-test
#' @param n1 Numeric. The sample size for group 1
#' @param n2 Numeric. The sample size for group 2. For a one sample test, enter any value >= 3 for \code{n2}. \code{n2} will be ignored.
#' @param d Numeric. The 'Cohen's d' standardized critical effect size. Cohen's d is the difference between group means divided by the pooled within-group standard deviation.
#' @param sig.level Numeric. The significance level. Must be between \code{0} and \code{1}. Defaults to \code{0.05}.
#' @param type Character string. The type of t-test being undertaken. Must be \code{"two.sample"}, \code{"one.sample"}, or \code{"paired"}. If ignored, \code{"two.sample"} is the default.
#' @param tails Character string. The number of tails being examined. Must be either \code{"two.tailed"} or \code{"one.tailed"}. If ignored, \code{"two.tailed"} is the default.
beta.t.test <- function(n1 = NULL,
                        n2 = NULL,
                        d = NULL,
                        sig.level = 0.05,
                        type = c("two.sample", "one.sample", "paired"),
                        tails = c("two.tailed","one.tailed")) {
  if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > sig.level | sig.level > 1)) 
    stop("sig.level must be a numeric value between 0 and 1")
  if (!is.null(n1) && n1 < 2) 
    stop("The value of n1 (number of observations in the first group) must be at least 2")
  
  type <- match.arg(type)
  
  tails <- match.arg(tails)
  
  d <- abs(d)
  
  tsample <- switch(type,
                    one.sample = 1,
                    two.sample = 2,
                    paired = 1)
  
  tside <- switch(tails,
                  one.tailed = 1,
                  two.tailed = 2)
  
  if (tside == 1) {
    p.body <- quote({
      nu <- switch(type,
                   one.sample = n1 - 1,
                   two.sample = n1 + n2 - 2,
                   paired = n1 - 1)
      delta <- d * switch(type,
                          one.sample = sqrt(n1),
                          two.sample = (1 / sqrt(1 / n1 + 1 / n2)),
                          paired = sqrt(n1))
      
      pt(q = qt(p = sig.level / tside,
                df = nu,
                lower.tail = FALSE),
         df = nu,
         ncp = delta,
         lower.tail = FALSE)
    })
  }
  if (tside == 2) {
    p.body <- quote({
      nu <- switch(type,
                   one.sample = n1 - 1,
                   two.sample = n1 + n2 - 2,
                   paired = n1 - 1)
      qu <- qt(p = sig.level / tside,
               df = nu,
               lower.tail = FALSE)
      delta <- d * switch(type,
                          one.sample = sqrt(n1),
                          two.sample = (1 / sqrt(1 / n1 + 1 / n2)),
                          paired = sqrt(n1))
      pt(q = qu,
         df = nu,
         ncp = delta,
         lower.tail = FALSE) +
        pt(q = -qu,
           df = nu,
           ncp = delta,
           lower.tail = TRUE)
    })
  }
  1 - eval(p.body)
}

#' Calculate average error
#' @param alpha Numeric. The alpha value to use.
#' @param n1 Numeric. The sample size for group 1
#' @param n2 Numeric. The sample size for group 2. For a one sample test, enter any value >= 3 for \code{n2}. \code{n2} will be ignored.
#' @param d Numeric. The 'Cohen's d' standardized critical effect size. Cohen's d is the difference between group means divided by the pooled within-group standard deviation.
#' @param T1T2cratio Numeric. The cost ratio of Type I errors relative to Type II errors. T1T2cratio is set at \code{1} as a default, making Type I and Type II errors equally serious.
#' @param HaHopratio Numeric. The prior probability of the alternate hypothesis relative to the prior probability of the null hypothesis. \code{HaHopratio} is set at \code{1} as a default, to not weight alpha and beta by their prior probabilities (assuming they are unknown).
#' @param type Character string. The type of t-test being undertaken. Must be \code{"two.sample"}, \code{"one.sample"}, or \code{"paired"}. If ignored, \code{"two.sample"} is the default.
#' @param tails Character string. The number of tails being examined. Must be either \code{"two.tailed"} or \code{"one.tailed"}. If ignored, \code{"two.tailed"} is the default.
w.average.error <- function(alpha = NULL,
                            n1 = NULL,
                            n2 = NULL,
                            d = NULL,
                            T1T2cratio = 1,
                            HaHopratio = 1,
                            type = c("two.sample", "one.sample", "paired"),
                            tails = c("two.tailed","one.tailed")) {
  numerator <- alpha * T1T2cratio + HaHopratio * beta.t.test(n1 = n1,
                                                             n2 = n2,
                                                             d = d,
                                                             sig.level = alpha,
                                                             type = type,
                                                             tails  = tails)
  denominator <- HaHopratio + T1T2cratio
  
  numerator / denominator
} 

#' Calculate the minimum average error
#' @param n1 Numeric. The sample size for group 1
#' @param n2 Numeric. The sample size for group 2. For a one sample test, enter any value >= 3 for \code{n2}. \code{n2} will be ignored.
#' @param d Numeric. The 'Cohen's d' standardized critical effect size. Cohen's d is the difference between group means divided by the pooled within-group standard deviation.
#' @param T1T2cratio Numeric. The cost ratio of Type I errors relative to Type II errors. T1T2cratio is set at \code{1} as a default, making Type I and Type II errors equally serious.
#' @param HaHopratio Numeric. The prior probability of the alternate hypothesis relative to the prior probability of the null hypothesis. \code{HaHopratio} is set at \code{1} as a default, to not weight alpha and beta by their prior probabilities (assuming they are unknown).
#' @param type Character string. The type of t-test being undertaken. Must be \code{"two.sample"}, \code{"one.sample"}, or \code{"paired"}. If ignored, \code{"two.sample"} is the default.
#' @param tails Character string. The number of tails being examined. Must be either \code{"two.tailed"} or \code{"one.tailed"}. If ignored, \code{"two.tailed"} is the default.
min.average.error <- function(n1 = NULL,
                              n2 = NULL,
                              d = NULL,
                              T1T2cratio = 1,
                              HaHopratio = 1,
                              type = c("two.sample", "one.sample", "paired"),
                              tails = c("two.tailed","one.tailed")) {
  # Note that only f, interval, and tol are arguments for optimize()
  # The rest are passed to the function f (w.average.error())
  unlist(optimize(f = w.average.error,
                  interval = c(0, 1),
                  tol = 0.0000000000001,
                  n1 = n1,
                  n2 = n2,
                  d = d,
                  T1T2cratio = T1T2cratio,
                  HaHopratio = HaHopratio,
                  type = type,
                  tails = tails))[2]
}

#' Calculate the optimal alpha
#' @param n1 Numeric. The sample size for group 1
#' @param n2 Numeric. The sample size for group 2. For a one sample test, enter any value >= 3 for \code{n2}. \code{n2} will be ignored.
#' @param d Numeric. The 'Cohen's d' standardized critical effect size. Cohen's d is the difference between group means divided by the pooled within-group standard deviation.
#' @param T1T2cratio Numeric. The cost ratio of Type I errors relative to Type II errors. T1T2cratio is set at \code{1} as a default, making Type I and Type II errors equally serious.
#' @param HaHopratio Numeric. The prior probability of the alternate hypothesis relative to the prior probability of the null hypothesis. \code{HaHopratio} is set at \code{1} as a default, to not weight alpha and beta by their prior probabilities (assuming they are unknown).
#' @param type Character string. The type of t-test being undertaken. Must be \code{"two.sample"}, \code{"one.sample"}, or \code{"paired"}. If ignored, \code{"two.sample"} is the default.
#' @param tails Character string. The number of tails being examined. Must be either \code{"two.tailed"} or \code{"one.tailed"}. If ignored, \code{"two.tailed"} is the default.
alpha <- function(n1 = NULL,
                  n2 = NULL,
                  d = NULL,
                  T1T2cratio = 1,
                  HaHopratio = 1,
                  type = c("two.sample", "one.sample", "paired"),
                  tails = c("two.tailed","one.tailed")) {
  # Note that only f, interval, and tol are arguments for optimize()
  # The rest are passed to the function f (w.average.error())
  unlist(optimize(w.average.error,
                  c(0, 1),
                  tol = 0.000000000001,
                  n1 = n1,
                  n2 = n2,
                  d = d,
                  T1T2cratio = T1T2cratio,
                  HaHopratio = HaHopratio,
                  type = type,
                  tails = tails))[1]
}

beta<-function (n1=NULL,n2=NULL,d=NULL,T1T2cratio=1,HaHopratio=1,type = c("two.sample", "one.sample", "paired"),tails = c("two.tailed","one.tailed")) 
  ((T1T2cratio+HaHopratio)*min.average.error(n1=n1,n2=n2,d=d,T1T2cratio=T1T2cratio,HaHopratio=HaHopratio,type=type,tails=tails)-T1T2cratio*alpha(n1=n1,n2=n2,d=d,T1T2cratio=T1T2cratio,HaHopratio=HaHopratio,type=type,tails=tails))/HaHopratio
#' Calculate the optimal beta
#' @param n1 Numeric. The sample size for group 1
#' @param n2 Numeric. The sample size for group 2. For a one sample test, enter any value >= 3 for \code{n2}. \code{n2} will be ignored.
#' @param d Numeric. The 'Cohen's d' standardized critical effect size. Cohen's d is the difference between group means divided by the pooled within-group standard deviation.
#' @param T1T2cratio Numeric. The cost ratio of Type I errors relative to Type II errors. T1T2cratio is set at \code{1} as a default, making Type I and Type II errors equally serious.
#' @param HaHopratio Numeric. The prior probability of the alternate hypothesis relative to the prior probability of the null hypothesis. \code{HaHopratio} is set at \code{1} as a default, to not weight alpha and beta by their prior probabilities (assuming they are unknown).
#' @param type Character string. The type of t-test being undertaken. Must be \code{"two.sample"}, \code{"one.sample"}, or \code{"paired"}. If ignored, \code{"two.sample"} is the default.
#' @param tails Character string. The number of tails being examined. Must be either \code{"two.tailed"} or \code{"one.tailed"}. If ignored, \code{"two.tailed"} is the default.


optab<-function (n1=NULL,n2=NULL,d=NULL,T1T2cratio=1,HaHopratio=1,type = c("two.sample", "one.sample", "paired"),tails = c("two.tailed","one.tailed")) {  
#' Produce a table with statistics including optimal alpha, optimal beta, overall probability of error, and cost-weighted probability of error
#' @param n1 Numeric. The sample size for group 1
#' @param n2 Numeric. The sample size for group 2. For a one sample test, enter any value >= 3 for \code{n2}. \code{n2} will be ignored.
#' @param d Numeric. The 'Cohen's d' standardized critical effect size. Cohen's d is the difference between group means divided by the pooled within-group standard deviation.
#' @param T1T2cratio Numeric. The cost ratio of Type I errors relative to Type II errors. T1T2cratio is set at \code{1} as a default, making Type I and Type II errors equally serious.
#' @param HaHopratio Numeric. The prior probability of the alternate hypothesis relative to the prior probability of the null hypothesis. \code{HaHopratio} is set at \code{1} as a default, to not weight alpha and beta by their prior probabilities (assuming they are unknown).
#' @param type Character string. The type of t-test being undertaken. Must be \code{"two.sample"}, \code{"one.sample"}, or \code{"paired"}. If ignored, \code{"two.sample"} is the default.
#' @param tails Character string. The number of tails being examined. Must be either \code{"two.tailed"} or \code{"one.tailed"}. If ignored, \code{"two.tailed"} is the default.
  
  list(
    "test type"=match.arg(type),
    "tails"=match.arg(tails),
    "output"=t(data.frame(
      "sample size 1"=n1,
      "sample size 2"=n2,
      "Cohen's d effect size"=d,
      "Type I/II error cost ratio"=T1T2cratio,
      "Ha/Ho prior probability ratio"=HaHopratio,
      "overall probability of error"=(alpha(n1=n1,n2=n2,d=d,T1T2cratio=T1T2cratio,HaHopratio=HaHopratio,type=type,tails=tails)+HaHopratio*beta(n1=n1,n2=n2,d=d,T1T2cratio=T1T2cratio,HaHopratio=HaHopratio,type=type,tails=tails))/(1+HaHopratio),
      "cost-weighted probability of error"=min.average.error(n1=n1,n2=n2,d=d,T1T2cratio=T1T2cratio,HaHopratio=HaHopratio,type=type,tails=tails), 
      "optimal alpha"=alpha(n1=n1,n2=n2,d=d,T1T2cratio=T1T2cratio,HaHopratio=HaHopratio,type=type,tails=tails),
      "optimal beta"=beta(n1=n1,n2=n2,d=d,T1T2cratio=T1T2cratio,HaHopratio=HaHopratio,type=type,tails=tails),row.names="values"))
  )
  
}

optab.plot <- function(n1=NULL,n2=NULL,d=NULL,a_range=c(0,1),step=0.001,T1T2cratio=1,HaHopratio=1,type = c("two.sample", "one.sample", "paired"),tails = c("two.tailed","one.tailed"),xlim=c(0,0.5),ylim=c(0,0.5)) {  
  rng <- seq(a_range[1],a_range[2],by=step)
  out <- data.frame("alpha"=rng,"w"=numeric(length(rng)))
#' Create a plot for a range of alpha values, marking the optimal alpha and the minimum average error,
#' @param n1 Numeric. The sample size for group 1
#' @param n2 Numeric. The sample size for group 2. For a one sample test, enter any value >= 3 for \code{n2}. \code{n2} will be ignored.
#' @param d Numeric. The 'Cohen's d' standardized critical effect size. Cohen's d is the difference between group means divided by the pooled within-group standard deviation.
#' @param a_range Numeric vector. The minimum and maximum values for alpha. Defaults to \code{c(0, 1)}.
#' @param T1T2cratio Numeric. The cost ratio of Type I errors relative to Type II errors. T1T2cratio is set at \code{1} as a default, making Type I and Type II errors equally serious.
#' @param HaHopratio Numeric. The prior probability of the alternate hypothesis relative to the prior probability of the null hypothesis. \code{HaHopratio} is set at \code{1} as a default, to not weight alpha and beta by their prior probabilities (assuming they are unknown).
#' @param type Character string. The type of t-test being undertaken. Must be \code{"two.sample"}, \code{"one.sample"}, or \code{"paired"}. If ignored, \code{"two.sample"} is the default.
#' @param tails Character string. The number of tails being examined. Must be either \code{"two.tailed"} or \code{"one.tailed"}. If ignored, \code{"two.tailed"} is the default.
#' @param xlim Numeric vector. The minimum and maximum values for the x axis of the plot. Defaults to \code{c(0, 0.5)}.
#' @param ylim Numeric vector. The minimum and maximum values for the y axis of the plot. Defaults to \code{c(0, 0.5)}.
  for (a in rng) {
    out[out$alpha==a,2]=w.average.error(a,n1,n2,d,T1T2cratio,HaHopratio,type,tails)
  }
  p<-ggplot(data=out,aes(x=alpha,y=w))+geom_line(size=1.5)+coord_cartesian(xlim=xlim,ylim=ylim)+
    xlab("Alpha")+ylab("Average Error Rate (alpha + beta)")+
    #ggtitle("Average error by alpha level")+
    geom_hline(yintercept=min.average.error(n1=n1,n2=n2,d=d,T1T2cratio=T1T2cratio,HaHopratio=HaHopratio,type=type,tails=tails),linetype=2,color="steelblue")+
    geom_vline(xintercept=alpha(n1=n1,n2=n2,d=d,T1T2cratio=T1T2cratio,HaHopratio=HaHopratio,type=type,tails=tails),linetype=2,color="steelblue")
    p+theme(text=element_text(size=18),axis.text.y=element_text(angle=90))
}

effect.size <- function(x1,n1,s1,type,pctChg,x2,n2,s2,paired,threshold=0){
  rho = 0.6
  if(type=="pct"){
    xd <- x1*pctChg/100
#' Calculate effect size
#' @param x1 Numeric. The initial value for group 1.
#' @param n1 Numeric. The sample size for group 1.
#' @param s1 Numeric.
#' @param type Character string. The type of effect size calculation. Must be \code{"pct"}, \code{"two.sample"}, or \code{"threshold"}. If ignored, \code{"pct"} is the default.
#' @param x2 Optional, numeric. The initial value for group 2. Only used if \code{type} is \code{"two.sample"}.
#' @param n2 Optional, numeric. The sample size for group 2. Only used if \code{type} is \code{"two.sample"}.
#' @param s2 Optional, numeric. Only used if \code{type} is \code{"two.sample"}.
#' @param paired Optional, logical. Should the test be paired? Only used if \code{type} is \code{"pct"} or \code{"threshold"}. Defaults to \code{TRUE}
#' @param threshold Optional, numeric. The difference threshold. Only used if \code{type} is \code{"threshold"}.
    if(paired) {
      sp = sqrt( (n1-1)/n1*s1*2*(1-rho) )
    } else {
      sp <- sqrt( (n1-1)*s1/n1  )
    }
  } else if (type=="two.sample") {
    xd <- abs(x1-x2)
    sp <- sqrt(((n1-1)*s1+(n2-1)*s2)/(n1+n2-2))
  } else { # Threshold
    xd <- abs(x1-threshold)
    if(paired) {
      sp = sqrt( (n1-1)/n1*s1*2*(1-rho) )
    } else {
      sp <- sqrt( (n1-1)*s1/n1  )
    }
  }
