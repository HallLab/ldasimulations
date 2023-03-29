#' Generate longitudinal/repeated measures data
#'
#' Simulates data for \strong{time}, a time (in)variant 
#' predictor \strong{x}, and a \strong{response} dependent 
#' on both \strong{x}, \strong{time} and their interaction. 
#'
#' \strong{Size of data:} Determined by number of samples 
#' (clusters) and the number of waves (timepoints). 
#'
#' \strong{Time intervals:} There are three options for 
#' time interval spacing. 
#' \itemize{
#'	\item \emph{Equal} - Time intervals are equal and even.
#' All individuals have the same timepoints.
#'	\item \emph{Increasing} - Time intervals increase with
#' each measure. All individuals have the same timepoints.
#'	\item \emph{Random} - Time intervals are from a random
#' uniform distribution. Each individual has different 
#' timepoints. 
#' }
#' 
#' \strong{Predictor variable:} The predictor \strong{x} 
#' can come from a binomial, normal, chi-square, or 
#' uniform distribution.  If binary, the user determines 
#' P(succcess). If normal, \strong{x} is centered 
#' at 0 with user-specified standard deviation. If
#' chi-square, the user specifies degrees-of-freedom. If 
#' uniform, the user sets the minimum and maximum possible values.
#' \strong{x} can be fixed within or vary within a cluster. 
#' 
#' \strong{Random effects:} Random effects define group 
#' structure/dependency. Group-level/cluster-level 
#' effects are Level 2 effects. Observation-level effects
#' are Level 1 effects. In longitudinal data, the "group"
#' is a single individual. The "observations" are the unique
#' measurements from a given individual at a given time. 
#' A cluster can have a random intercept (RI) and, if desired,
#' a random slope (RS). Random effects are drawn from normal 
#' distributions centered at 0 with user-specified standard
#' deviations. If both a RI and RS are designated, 
#' the user can specify the covariance between
#' these two random effects. The greater the variance of 
#' the Level 2 effects, relative to the Level 1 variation,
#' the more similar the values within a group.
#' That is, the response values within the same
#' group will be more similar than those without. 
#' 
#' \strong{Data generation:} Data is generated as a linear 
#' or nonlinear combination of \code{x}, \code{time}, 
#' Level 2 effects (RI and RS), and Level 1 effects 
#' (error; E). RI, RS, and E are random error terms 
#' which add variability to the effects of \code{x} and
#' \code{time} on the outcome \code{y}. The relationship
#' between \code{y} and predictors \code{x} and \code{time}
#' can be linear, exponential, or parabolic.
#'
#' \strong{Linear relationship:}
#' \deqn{y = b0 + b1x + b2time + b3x*time + Error}
#'
#' \strong{Exponential relationship:}
#' \deqn{y = b0 + e^{b1x + b2time + b3x*time} + Error}
#'
#' \strong{Parabolic relationship:}
#' \deqn{y = b0 + b1x + b2time^2 + b3x*time^2 + Error}
#'
#' \strong{Error terms:}
#' \itemize{
#'	\item \emph{RI error} - RI + E
#'	\item \emph{RI + RS error} - RI + RS*time + E
#' }
#'
#' @param clus Numeric. Number of clusters. Default is 100.
#' @param waves Numeric. Number of timepoints per sample. Default is 3.
#' @param tType Numeric. Type of interval spacing. If tType=1,
#' equal intervals. If tType=2, increasing intervals.
#' If tType=3, random intervals. Default is 1.
#' @param incr Numeric. The increment of time spacing when 
#' tType=1 or tType=2. Default is 1. 
#' @param mint Numeric. Lower limit of time, if tType=3.
#' Must be greater than 0, but less than maxt. Default is 1.
#' @param maxt Numeric. Upper limit of time, if tType=3.
#' Must be greater than mint. Default is 10.
#' @param rEff String. Type of random effect(s). "RI" for
#' a random intercept-only model. "RIS" for a model with both
#' a random intercept and random slope. Default is "RI".
#' @param xType Numeric. Predictor type. Binary=1, normal=2, 
#' chi-square=3, or uniform=4. Default is 1.
#' @param xVary Logical. Does predictor \code{x} vary over time?
#' Default is \code{F}.
#' @param xsd Numeric. Standard deviation of \code{x}, if normal.
#' Default is 1.
#' @param xdf Numeric. Degrees-of-freedom of \code{x}, if 
#' chi-square. Default is 1.
#' @param xmin Numeric. Min possible value of \code{x}, 
#' if uniform. Default is 1.
#' @param xmax Numeric. Max possible value of \code{x}, 
#' if uniform. Default is 10.
#' @param xp Numeric. P(success) if \code{x} is binary. Range
#' (0,1). Default is 0.5.
#' @param yType Numeric. Relationship between predictors
#' and outcome. Linear=1, exponential=2, parabolic=3. 
#' Default is Linear=1. 
#' @param gvar Numeric. Group-level variance of random 
#' effect(s). Default is 1.
#' @param ivar Numeric. Observation-level variance. Default
#' is 1. 
#' @param cov Numeric. Covariance of random effects, if 
#' \code{xEff="RIS"}. Default is 0.5.
#' @param b0 Numeric. True intercept of the model. Default
#' is 1.
#' @param b1 Numeric. True effect of \code{x} on \code{y}. 
#' Default is 1.
#' @param b2 Numeric. True effect of \code{time} on \code{y}. 
#' Default is 1.
#' @param b3 Numeric. True effect of \code{x*time} on \code{y}. 
#' Default is 1.
#' @param mType Numeric. Type of missingness. 1 = random 
#' missingness; 2 = non-random missingness
#' @param mPer Numeric. Fraction of missing data. Range 
#' (0, 1). If mType = 1, then randomly deletes mPer*100%
#' of rows from who data frame. If mType = 2, deletes
#' mPer*100% rows from the groups selected to 
#' have missingness.
#' @param propDel Numeric. Range (0,1). Only applies
#' if mType=2. Denotes proportion of groups from which 
#' data is to be missing.
#'
#' @returns Dataframe of simulated data.
#' \itemize{
#'	\item \strong{id} - Sample (cluster) labels
#'	\item \strong{y} - Response variable
#'	\item \strong{x} - Predictor variable
#'	\item \strong{time} - Time variable
#'	\item \strong{ri} - Random intercept value
#'	\item \strong{rs} - Random slope value, if applicable
#'	\item \strong{e} - Individual error 
#' }
#' @examples
#' rpm_data(clus=5, tType=1)
#' rpm_data(clus=5, tType=2)
#' rpm_data(clus=5, tType=3)
#' rpm_data(clus=5, ype=2)
#' rpm_data(clus=5, yType=3)
#' rpm_data(clus=5, mType=1)
#' rpm_data(clus=5, mType=2)
#' @export

# X predictor: xType {1:binary, 2:normal, 3:Chisq, 4:Uniform}
# Outcome: yType {1:linear, 2:exponential, 3:parabolic}	
# Missingness: mType {0:complete data, 1:random missingness, 2:targeted missingness}

rpm_data = function(clus=100, waves=3, tType=1, incr=1, mint=1, maxt=10, rEff='RI', gvar=1, ivar=1, cov=0.5, xType=1, xVary=F, xsd=1, xdf=1, xmin=1, xmax=10, xp=0.5, yType=1, b0=1, b1=1, b2=1, b3=1, mType=0, mPer=0.5, propDel=0.5){
	
	# Group IDs
	clusters = rep(1:clus, each=waves) # cluster number
	n = length(clusters) # total samples
	
	# Time Waves: tType {1:balanced, 2:increasing, 3:random}
	times = make_time(clus=clus, waves=waves, tType=tType, incr=incr, mint=mint, maxt=maxt)
	
	# Random Effects: Group level
	if(rEff=='RI'){ ri = rnorm(clus, mean=0, sd=sqrt(gvar)) }
	if(rEff=='RIS'){sig = matrix(c(gvar, cov, cov, gvar), 2, 2); mat = MASS::mvrnorm(n = clus, mu = c(0,0), Sigma = sig) } # col1=r int, col2=r slope}
	
	# Individual Noise
	e = rnorm(n, mean=0, sd=sqrt(ivar))

	# X predictor: xType {1:binary, 2:normal, 3:chisq, 4:uniform}	
	if(xType==1){
		# Binary predictor
		if(xVary==F){ 
			uniq=1
			while(uniq == 1){
				# Generate binary predictor
				x = rbinom(n=clus, size=1, prob=xp) 
				# Stop loop if there is variation
				if(length(unique(x))==2){break}
				# If no variation, loop continues
			}		
		} 
		
		if(xVary==T){ 
			uniq=1
			while(uniq == 1){
				# Generate binary predictor
				x = rbinom(n=n, size=1, prob=xp)
				# Stop loop if there is variation
				if(length(unique(x))==2){break}
				# If no variation, loop continues
			}
		} 
			 
	}
	
	if(xType==2){
		# Normal predictor
		if(xVary==F){ x = rnorm(n=clus, mean=0, sd=xsd) }
		if(xVary==T){ x = rnorm(n=n, mean=0, sd=xsd) }
	}
	
	if(xType==3){
		# Chisq predictor
		if(xVary==F){ x = rchisq(n=clus, df=xdf) }
		if(xVary==T){ x = rchisq(n=n, df=xdf) }	
	}
	
	if(xType==4){
		# Uniform predictor
		if(xVary==F){ x = runif(n=clus, min=xmin, max=xmax)}	
		if(xVary==T){ x = runif(n=n, min=xmin, max=xmax)}	
	}
	
	
	# Outcome: yType {1:linear, 2:exponential, 3:parabolic}	
	if(rEff=='RI'){ # RI only
		#print('RI only')
		if(yType==1){
			#print('RI & yType==1')
			# Linear outcomes
			if(xVary==F){ y = b0 + b1*x[clusters] + b2*times + b3*x[clusters]*times + ri[clusters] + e }
			if(xVary==T){ y = b0 + b1*x + b2*times + b3*x*times + ri[clusters] + e }	
		}
		if(yType==2){
			#print('RI & yType==2')
			# Exponential outcome 
			# y = time*exp(1)^(x) + ri + e
			if(xVary==F){ y = b0 + exp(1)^(b1*x[clusters] + b2*times + b3*x[clusters]*times) + ri[clusters] + e }
			if(xVary==T){ y = b0 + exp(1)^(b1*x + b2*times + b3*x*times) + ri[clusters] + e }
		}
		if(yType==3){
			#print('RI & yType==3')
			# Parabolic outcome 
			# y = x + time^2 + x*time^2 + ri + e
			if(xVary==F){ y = b0 + b1*x[clusters] + b2*times^2 + b3*x[clusters]*times^2 + ri[clusters] + e }
			if(xVary==T){ y = b0 + b1*x + b2*times^2 + b3*x*times^2 + ri[clusters] + e }	
		}
	}
	if(rEff=='RIS'){ # RI + RS
		#print('RI & RS')
		if(yType==1){
			#print('RIS & yType==1')
			# Linear outcomes
			if(xVary==F){ y = b0 + b1*x[clusters] + b2*times + b3*x[clusters]*times + (mat[,1][clusters] + mat[,2][clusters]*times + e) }
			if(xVary==T){ y = b0 + b1*x + b2*times + b3*x*times + (mat[,1][clusters] + mat[,2][clusters]*times + e) }	
		}
		if(yType==2){
			#print('RIS & yType==2')
			# Exponential outcome 
			# y = time*exp(1)^(x) + ri + e
			if(xVary==F){ y = b0 + exp(1)^(b1*x[clusters] + b2*times + b3*x[clusters]*times) + (mat[,1][clusters] + mat[,2][clusters]*times + e) }
			if(xVary==T){ y = b0 + exp(1)^(b1*x + b2*times + b3*x*times) + (mat[,1][clusters] + mat[,2][clusters]*times + e) }
		}
		if(yType==3){
			#print('RIS & yType==3')
			# Parabolic outcome 
			# y = x + time^2 + x*time^2 + ri + e
			if(xVary==F){ y = b0 + b1*x[clusters] + b2*times^2 + b3*x[clusters]*times^2 + (mat[,1][clusters] + mat[,2][clusters]*times + e) }
			if(xVary==T){ y = b0 + b1*x + b2*times^2 + b3*x*times^2 + (mat[,1][clusters] + mat[,2][clusters]*times + e) }	
		}
	}
	
	#y = b0 + b1*x[clusters] + b2*times + b3*x[clusters]*times + ri[clusters] + e
	# All Simulated Data
	if(rEff=='RI'){ # RI only
		if(xVary==F){ mydata = data.frame(id=clusters, y=y, x=x[clusters], time=times, ri=ri[clusters], e=e) }
		if(xVary==T){ mydata = data.frame(id=clusters, y=y, x=x, time=times, ri=ri[clusters], e=e) }
	}
	if(rEff=='RIS'){ # RI + RS
		if(xVary==F){ mydata = data.frame(id=clusters, y=y, x=x[clusters], time=times, ri=mat[,1][clusters], rs=mat[,2][clusters], e=e) }
		if(xVary==T){ mydata = data.frame(id=clusters, y=y, x=x, time=times, ri=mat[,1][clusters], rs=mat[,2][clusters], e=e) }
	}
	
	
		
	# Missingness: mType {0:complete data, 1:random missingness, 2:targeted missingness}
	
	if(mType==0){ return(mydata) }
	if(mType==1|mType==2){
		newdata = del_data(df=mydata, mType=mType, mPer=mPer, propDel=propDel) 	
		newdata = newdata[complete.cases(newdata),]
		return(newdata)
	}
	
}


