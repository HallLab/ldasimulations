#' Generate and analyze longitudinal data with 
#' a nonvarying predictor effect over time
#'
#' Simulates data for \strong{time}, a time (in)variant 
#' predictor \strong{x}, and a \strong{response} dependent 
#' on the marginal effects of \strong{x} and \strong{time}. 
#' Then applies various analytical methods. 
#'
#' \strong{Size of data:} Determined by number of samples 
#' (clusters) and the number of waves (timepoints). 
#'
#' \strong{Time intervals:} Time is spaced in 1-unit
#' increments
#'
#' \strong{Predictor variable:} The predictor \strong{x} 
#' is fixed to a standard normal distribution. 
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
#' \deqn{y = b0 + b1x + b2time + RI + E}
#' \deqn{y = b0 + b1x + e^{b2time} + RI + E}
#'
#' The following methods are then applied to the data:
#'
#' \strong{Naive linear regression (NLR):} Ignores clustering and time
#' when regressing outcome on predictor(s).
#'
#' \strong{Cluster robust standard error (CRSE):} Applies a CRSE correction
#' to the standard errors of the aforementioned NLR output. Clustering is the group
#' membership.
#'
#' \strong{Aggregate linear regression (AGG):} Averages all variables over 
#' time within each cluster. Regresses the average outcome against the 
#' average predictor(s).
#'
#' \strong{Linear mixed model (LMM):} Regresses outcome on the 
#' predictor(s) and time, while accounting for group dependency
#' with random effects. Two models - LMM1 and LMM2 - are implemented.
#' LMM1 is the correctly specified model, while LMM2 
#' evaluates the interaction between \strong{x} and
#' \strong{time}
#'
#' The models tested are 
#' \eqn{y = x} for NLR, CRSE, and AGG, 
#' \eqn{y = x + time + (1|group)} for LMM1, and
#' \eqn{y = x + time + x*time + (1|group)} for LMM2
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
#' and outcome. Linear=1, exponential=2.
#' Default is Linear=1. 
#' @param b0 Numeric. True intercept of the model. Default
#' is 1.
#' @param b1 Numeric. True effect of \code{x} on \code{y}. 
#' Default is 1.
#' @param b2 Numeric. True effect of \code{time} on \code{y}. 
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
#' power_lda_noInt()
#' @export

power_lda_noInt = function(clus=50, waves=3, gvar=1, ivar=1, xVary=F, yType=1, b0=1, b1=1, b2=1){
	
	# Generate clusters
	clusters = rep(1:clus, each=waves)
	
	# Generate time - in 1 unit increments
	times = rep(1:waves, times=clus)
		
	# Random Effects: Group level
	ri = rnorm(n=clus, mean=0, sd=sqrt(gvar))
	
	# Individual Noise
	e = rnorm(n=length(clusters), mean=0, sd=sqrt(ivar))
	
	# Generate X: Standard Normal predictor
	if(xVary){
		#x = rbinom(length(ids), size=1, prob=0.5)
		x = rnorm(n=length(clusters), mean=0, sd=1)
	}else{
		#x = rbinom(clus, size=1, prob=0.5)
		x = rnorm(n=clus, mean=0, sd=1)

	}
	
	# Generate outcome: yType {1=linear, 2=exponential}
	if(xVary){
			if(yType==1){ y = b0 + (b1*x) + (b2*times) + ri[clusters] + e } # linear
			if(yType==2){ y = b0 + (b1*x) + exp((b2*times)) + ri[clusters] + e } # exponential
		}else{	
			if(yType==1){ y = b0 + (b1*x[clusters]) + (b2*times) + ri[clusters] + e } # linear
			if(yType==2){ y = b0 + (b1*x[clusters]) + exp((b2*times)) + ri[clusters] + e } # exponential
		}
	
	# Format data
	if(xVary){
		dat = data.frame(id=clusters, y=y, x=x, time=times, ri=ri[clusters], e=e)
	}else{
		dat = data.frame(id=clusters, y=y, x=x[clusters], time=times, ri=ri[clusters], e=e)
	}
	
	if(length(unique(dat$x))==1){
	
	# Return an error value if x does not vary	
	return(data.frame(Estimate=99999, Std.Err=99999, Test.Stat=99999, p.value=99999, Variable=99999, Method=99999, Conv=99999))
		
	}else{
		
	# Base Formulas
    # Formula NLR, CRSE, AGG 'y~x' 
	# Formula LMM1 'y~x + time + (1|id)' 
	# Formula LMM2 'y~x + time + x*time + (1|id)' 
	
	# Analysis
	m1 = glm(y~x, data=dat)
	nlr = as.data.frame(coef(summary(m1)))
	nlr = nlr[-1,]
	nlr$Variable = rownames(nlr)
	nlr$Method='NLR'
	nlr$Conv = m1$converged # Convergence
	rownames(nlr) = 1:nrow(nlr)
	colnames(nlr)[1:4] = c('Estimate','Std.Err','Test.Stat','p.value')
	
	m2 = lmtest::coeftest(m1, vcov. = sandwich::vcovCL(m1, cluster=dat[,'id'], type='HC1'))
	crse = data.frame(matrix(m2, nrow=dim(m2)[1], ncol=dim(m2)[2]))
	rownames(crse) = rownames(m2)
	colnames(crse) = colnames(m2)
	crse = crse[-1,]
	crse$Variable = rownames(crse)
	crse$Method='CRSE'
	crse$Conv = m1$converged
	rownames(crse) = 1:nrow(crse)
	colnames(crse)[1:4] = c('Estimate','Std.Err','Test.Stat','p.value')
	
	col_name = colnames(dat); col_num = dim(dat)[2]; row_num = length(unique(dat[,'id']))
	d_g = by(dat, dat['id'], function(x) colMeans(x, na.rm=T)) # Average vars per group
	d_g_dat = as.data.frame( matrix(unlist(d_g), nrow=row_num, ncol=col_num, byrow=T) ) # Reformat into dataframe
	colnames(d_g_dat) = col_name
	m3 = glm(y~x, data=d_g_dat)
	agg = as.data.frame(coef(summary(m3)))
	agg = agg[-1,] 
	agg$Variable = rownames(agg)
	agg$Method = 'AGG'
	agg$Conv = m3$converged
	rownames(agg) = 1:nrow(agg)
	colnames(agg)[1:4] = c('Estimate','Std.Err','Test.Stat','p.value')
	
	m4 = lmerTest::lmer(y~x+time+(1|id), data=dat)
	lmm1 = as.data.frame(coef(summary(m4))[,-3]) # Drop DF column
	lmm1 = lmm1[-1,] 
	lmm1$Variable = rownames(lmm1)
	lmm1$Method = 'LMM1'
	lmm1$Conv = ifelse(length(m4@optinfo$conv$lme4) > 0, m4@optinfo$conv$lme4$messages, 0) # 0=no messages
	rownames(lmm1) = 1:nrow(lmm1)
	colnames(lmm1)[1:4] = c('Estimate','Std.Err','Test.Stat','p.value')
	
	m5 = lmerTest::lmer(y~x*time+(1|id), data=dat)
	lmm2 = as.data.frame(coef(summary(m5))[,-3]) # Drop DF column
	lmm2 = lmm2[-1,] 
	lmm2$Variable = rownames(lmm2)
	lmm2$Method = 'LMM2'
	lmm2$Conv = ifelse(length(m5@optinfo$conv$lme4) > 0, m5@optinfo$conv$lme4$messages, 0) # 0=no messages
	rownames(lmm2) = 1:nrow(lmm2)
	colnames(lmm2)[1:4] = c('Estimate','Std.Err','Test.Stat','p.value')
	
	out = rbind(nlr, crse, agg, lmm1, lmm2)
	return(out)
	
	}
}

