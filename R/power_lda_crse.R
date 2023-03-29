#' Generate and analyze longitudinal data
#' with a correctly specified CRSE model
#'
#' Simulates data for \strong{time}, a time (in)variant 
#' predictor \strong{x}, and a \strong{response} dependent 
#' on both \strong{x}, \strong{time}, and their interaction. 
#' Then applies various analytical methods, and 
#' applies CRSE to a model of correct fixed effects. 
#'
#' \strong{Size of data:} Determined by number of samples 
#' (clusters) and the number of waves (timepoints). 
#'
#' \strong{Time intervals:} Time is spaced in 1-unit
#' increments
#' 
#' \strong{Predictor variable:} The predictor \strong{x} 
#' comes from a standard normal distribution.
#' 
#' \strong{Random effects:} Random effects define group 
#' structure/dependency. Group-level/cluster-level 
#' effects are Level 2 effects. Observation-level effects
#' are Level 1 effects. In longitudinal data, the "group"
#' is a single individual. The "observations" are the unique
#' measurements from a given individual at a given time. 
#' This simulation only generates a random intercept (RI) per cluster
#' drawn from a normal distribution centered at 0 with 
#' user-specified standard deviation. The greater the variance of 
#' the Level 2 effects, relative to the Level 1 variation,
#' the more similar the values within a group.
#' 
#' \strong{Data generation:} Data is generated as a linear 
#' combination of \code{x}, \code{time}, a random intercept (RI)
#' and individual level error (E). The relationship
#' between \code{y} and predictors \code{x} and \code{time}
#' can be linear or exponential.
#'
#' \deqn{y = b0 + b1x + b2time + b3x*time + RI + E}
#' \deqn{y = b0 + e^{b1x + b2time + b3x*time} + RI + E}
#'
#' The following methods are then applied to the data:
#'
#' \strong{Naive linear regression (NLR):} Ignores clustering and time
#' when regressing outcome on predictor(s). The tested model 
#' is \eqn{y = x}
#'
#' \strong{Cluster robust standard error (CRSE):} Applies a CRSE correction
#' to the standard errors of GLM output. Clustering is the group
#' membership. The tested model is \eqn{y = x + time + x*time}
#'
#' \strong{Linear mixed model (LMM):} Regresses outcome on the 
#' predictor(s) and time, while accounting for group dependency
#' with random effects. The tested model is 
#' \eqn{y = x + time + x*time + (1|group)}
#'
#' To fit models, CRSE uses \code{lmtest} and the clustered variance
#' estimator from \code{sandwich} using type "HC1".  
#' LMM uses \code{lmer} and \code{lmerTest}.
#'
#' In randomly generated \code{x} has no variation instead of
#' reporting results for each method, the function returns a 
#' data.frame with an error code \code{99999} as all entries
#' in the return value.
#'
#' @param clus Numeric. Number of clusters. Default is 50.
#' @param waves Numeric. Number of timepoints per sample. Default is 3.
#' @param gvar Numeric. Group-level variance of random 
#' effect(s). Default is 1.
#' @param ivar Numeric. Individual-level variance. Default
#' is 1. 
#' @param xVary Logical. Does predictor \code{x} vary over time?
#' Default is \code{F}.
#' @param yType Numeric. Relationship between predictors
#' and outcome. Linear=1, exponential=2
#' Default is Linear=1. 
#' @param b0 Numeric. True intercept of the model. Default
#' is 1.
#' @param b1 Numeric. True effect of \code{x} on \code{y}. 
#' Default is 1.
#' @param b2 Numeric. True effect of \code{time} on \code{y}. 
#' Default is 1.
#' @param b3 Numeric. True effect of \code{x*time} on \code{y}. 
#' Default is 1.
#'
#' @returns Dataframe of analysis results
#' \itemize{
#' 	\item \strong{Estimate} - Coefficient estimate
#'	\item \strong{Std.Err} - Coefficient standard error
#'	\item \strong{Test.Stat} - Coefficient test statistic
#'	\item \strong{p.value} - Coefficient p-value
#'	\item \strong{Variable} - Predictor term in model
#'	\item \strong{Method} - Statistical analysis method
#'	\item \strong{Conv} - Convergence/error codes or messages. 
#' }
#' @examples
#' power_lda_crse()
#' @export

power_lda_crse = function(clus=50, waves=3, gvar=1, ivar=1, xVary=F, yType=1, b0=1, b1=1, b2=1, b3=1){
 
  # Generate clusters
  clusters = rep(1:clus, each=waves) # cluster number
  n = length(clusters) # total samples

  # Generate time - in 1 unit increments
  times = rep(1:waves, times=clus)
 
  # Random Effects: Group level
  ri = rnorm(clus, mean=0, sd=sqrt(gvar))

  # Individual Noise
  e = rnorm(n, mean=0, sd=sqrt(ivar))
  
  # Generate X: Standard Normal predictor
  if(xVary==F){ x = rnorm(n=clus, mean=0, sd=1) }
  if(xVary==T){ x = rnorm(n=n, mean=0, sd=1) }

  # Generate outcome: yType {1:linear, 2:exponential}	
  if(yType==1){
      # Linear outcomes
      if(xVary==F){ y = b0 + b1*x[clusters] + b2*times + b3*x[clusters]*times + ri[clusters] + e }
      if(xVary==T){ y = b0 + b1*x + b2*times + b3*x*times + ri[clusters] + e }	
    }
  if(yType==2){
      # Exponential outcome 
      if(xVary==F){ y = b0 + exp(1)^(b1*x[clusters] + b2*times + b3*x[clusters]*times) + ri[clusters] + e }
      if(xVary==T){ y = b0 + exp(1)^(b1*x + b2*times + b3*x*times) + ri[clusters] + e }
    }

  # All Simulated Data
  if(xVary==F){ mydata = data.frame(id=clusters, y=y, x=x[clusters], time=times, ri=ri[clusters], e=e) }
  if(xVary==T){ mydata = data.frame(id=clusters, y=y, x=x, time=times, ri=ri[clusters], e=e) }

  
  # ANALYSIS
  # NLR model: y~x
  # CRSE model: y~x+time+x:time
  # LMM model: y~x+time+x:time+(1|group)
  
  if(length(unique(mydata$x))==1){
    
    return(data.frame(Estimate=99999, Std.Err=99999, Test.Stat=99999, p.value=99999, Variable=99999, Method=99999, Conv=99999))
    
  }else{ 
    
    # Base Formulas
    f1 = 'y~x' # Formula NLR
    f2 = 'y~x*time' # Formula CRSE
    f3 = 'y~x*time + (1|id)' # Formula LMM
    
    # NLR
    m1 = glm(as.formula(f1), data=mydata)
    nlr = as.data.frame(coef(summary(m1)))
    nlr = nlr[-1,] 
    nlr$Variable = rownames(nlr)
    nlr$Method = 'NLR'
    nlr$Conv = m1$converged
    rownames(nlr) = 1:nrow(nlr)
    colnames(nlr)[1:4] = c('Estimate','Std.Err','Test.Stat','p.value')
    
    # CRSE
    m2f = glm(as.formula(f2), data=mydata) 
    m2 = lmtest::coeftest(m2f, vcov. = sandwich::vcovCL(m2f, cluster=mydata[,'id'], type='HC1'))
    c_col = colnames(m2); c_row = rownames(m2); c_dim = dim(m2)
    crse = data.frame(matrix(m2, nrow=c_dim[1], ncol=c_dim[2]))
    rownames(crse) = c_row
    colnames(crse) = c_col
    crse = crse[-1,]
    crse$Variable = rownames(crse)
    crse$Method = 'CRSE'
    crse$Conv = m2f$converged
    rownames(crse) = 1:nrow(crse)
    colnames(crse)[1:4] = c('Estimate','Std.Err','Test.Stat','p.value')
    
    # LMM
    m3 = lmerTest::lmer(as.formula(f3), data=mydata)
    lme = as.data.frame(coef(summary(m3))[,-3]) # remove df column
    lme = lme[-1,] 
    lme$Variable = rownames(lme)
    lme$Method = 'LMM'
    lme$Conv = ifelse(length(m3@optinfo$conv$lme4) > 0, m3@optinfo$conv$lme4$messages, 0)
    rownames(lme) = 1:nrow(lme)
    colnames(lme)[1:4] = c('Estimate','Std.Err','Test.Stat','p.value')
    
    out = rbind(nlr, crse, lme)
    return(out)
  }
  
}

