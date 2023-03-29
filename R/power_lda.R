#' Apply multiple analytic methods to longitudinal data
#'
#' Outputs the results of six analytical methods applied to longitudinal
#' input data. If the data was simulated with known 
#' parameters, the results can be used to compare power.
#' The input must be a dataframe with, at minimum, variables 
#' corresponding to cluster, a numerical predictor, a time variable,
#' and a numerical outcome variable. Cluster IDs must be numbers.
#' As this is longitudinal data, a "cluster" refers to a given sample 
#' in the data, with repeated measurements. 
#' 
#' The following methods are applied to the data:
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
#' \strong{Fixed effects model (FE):} Regresses the outcome on the predictor(s), 
#' time, and group. Explicitly models time and group, but doesn't account for 
#' within group dependency.
#'
#' \strong{Generalized estimating equation (GEE):} Regresses outcome on
#' the predictor(s) and time, while accounting for group dependency in
#' the assumed residual correlation structure.
#'
#' \strong{Linear mixed model (LMM):} Regresses outcome on the 
#' predictor(s) and time, while accounting for group dependency
#' with random effects.
#'
#' If no genotype is specified, the models tested are 
#' \eqn{y = x} for NLR, CRSE, and AGG methods, 
#' \eqn{y = x + time + x*time + group} for FE, and 
#' \eqn{y = x + time + x*time} for GEE and LMM. For LMM, the 
#' random effects are specified with an additional
#' \eqn{ + (1|group)} for a random intercept or 
#' \eqn{ + (1 + time|group)} for a random intercept and slope.
#'
#' If genotype is specified, the models are  
#' \eqn{y = x + geno + x*geno} for NLR, CRSE, and AGG methods, 
#' \eqn{y = x + geno + time + x*geno + x*time + geno*time + x*geno*time + group}
#' for FE, and 
#' \eqn{y = x + geno + time + x*geno + x*time + geno*time + x*geno*time} 
#' for FE, GEE, and LMM. The LMM random effects are specified
#' as described previously.
#'
#' To fit models, CRSE uses \code{lmtest} and the clustered variance
#' estimator from \code{sandwich} using type "HC1". GEE uses \code{geepack}.  
#' LMM uses \code{lmer} and \code{lmerTest}.
#'
#' In randomly generated \code{x} has no 
#' variation (e.g. can happen if there are few clusters and a binary
#' predictor with a high probability of success), instead of
#' reporting results for each method, the function returns a 
#' data.frame with an error code \code{99999} as all entries
#' in the return value.
#'
#' @param df Dataframe. Dataframe object with data to be analyzed.
#' @param group String. Column name of grouping/cluster variable. Must be numerical.
#' @param y String. Column name of outcome variable aka "y".
#' @param x String. Column name of predictor variable aka "x".  
#' @param time String. Column name of time variable aka "time".
#' @param geno String. Column name of genotype data aka "geno", if applicable.
#' @param rEff String. Random effect(s) for  the LMM. "RI" for
#' a random intercept-only model. "RIS" for a model with both
#' a random intercept and random slope. Default is "RI".
#' @param GEEcor String. Residual correlation structure for the GEE. 
#' Must match one of the correlation structures implemented by 
#' \code{geepack}. Default is 'independence'.
#'
#' @returns Dataframe of analysis results
#' \itemize{
#' 	\item \strong{Estimate} - Coefficient estimate
#'	\item \strong{Std.Err} - Coefficient standard error
#'	\item \strong{Test.Stat} - Coefficient test statistic
#'	\item \strong{p.value} - Coefficient p-value
#'	\item \strong{Variable} - Predictor term in model
#'	\item \strong{Outcome} - Outcome term in model
#'	\item \strong{Method} - Statistical analysis method
#'	\item \strong{Conv} - Convergence/error codes or messages. 
#' }
#'
#'  In the results table, for NLR, AGG, CRSE,
#' 	and FE the convergence code is 1 if converged and 0 if not. For GEE,
#'  it prints 0 if there are no errors, and prints the error message if
#'  not For LMM, it likewise reports 0 for no errors, and prints the 
#'	warning message if not. Note that for GEE and LMM, the reported 
#'  warning may deal with something other than convergence 
#' 
#"
#' @examples
#' df=data.frame(ids=rep(1:5, each=10), time=rep(1:10, times=5), x=rnorm(50), y=rnorm(50))
#' power_lda(df=df, group='ids', y='y', x='x', time='time', rEff='RI', GEEcor='exchangeable')
#' @export


power_lda = function(df=NULL, group=NULL, y=NULL, x=NULL, time=NULL, geno=NULL, 
					rEff='RI', GEEcor='independence'){
						
	# True variable names
	cols = c(group, y, x, time, geno) # geno NULL if absent
	# Desired varnames
	col_n = c('group','y','x','time','geno')
	# Select vars
	d = df[,cols]
	# Rename vars
	colnames(d) = col_n[1:length(cols)] # length cols 5 if geno present, 4 if geno absent
	
	
	if(length(unique(d$x))==1){
		
		return(data.frame(Estimate=99999, Std.Err=99999, Test.Stat=99999, p.value=99999, Variable=99999, Outcome=99999, Method=99999, Conv=99999))
		
	}else{
			
	# Define model formulas
	if(is.null(geno)){
		
		# Base Formulas
		f1 = 'y~x' # Formula NLR, CRSE, Agg
		f2 = 'y~x*time' # Formula FE, GEE, LMM
	
		# LMM formulas
		if(rEff=='RI') {f3 = paste0(f2, '+ (1|group)')}
		if(rEff=='RIS') {f3 = paste0(f2, '+ (1 + time|group)')}
		
	} else{
		
		# Base Formulas
		f1 = 'y~x*geno' # Formula NLR, CRSE, Agg
		f2 = 'y~x*geno*time' # Formula FE, GEE, LMM
	
		# LMM formulas
		if(rEff=='RI') {f3 = paste0(f2, '+ (1|group)')}
		if(rEff=='RIS') {f3 = paste0(f2, '+ (1 + time|group)')}
		
	}
	
	
		# NLR
		m1 = glm(as.formula(f1), data=d)
		nlr = as.data.frame(coef(summary(m1)))
		nlr = nlr[-1,] 
		nlr$Variable = rownames(nlr)
		nlr$Outcome = y
		nlr$Method = 'NLR'
		nlr$Conv = m1$converged
		rownames(nlr) = 1:nrow(nlr)
		colnames(nlr)[1:4] = c('Estimate','Std.Err','Test.Stat','p.value')
		
		# CRSE
		m2 = lmtest::coeftest(m1, vcov. = sandwich::vcovCL(m1, cluster=d[,'group'], type='HC1'))
		c_col = colnames(m2); c_row = rownames(m2); c_dim = dim(m2)
		crse = data.frame(matrix(m2, nrow=c_dim[1], ncol=c_dim[2]))
		rownames(crse) = c_row
		colnames(crse) = c_col
		crse = crse[-1,]
		crse$Variable = rownames(crse)
		crse$Outcome = y
		crse$Method = 'CRSE'
		crse$Conv = m1$converged
		rownames(crse) = 1:nrow(crse)
		colnames(crse)[1:4] = c('Estimate','Std.Err','Test.Stat','p.value')
		
		# ALR
		col_name = colnames(d)
		col_num = dim(d)[2]
		row_num = length(unique(d[,'group']))
		d_g = by(d, d['group'], function(x) colMeans(x, na.rm=T)) # Average vars per group; "averages" group ID too, so make sure IDs are numerical
		d_g_dat = as.data.frame( matrix(unlist(d_g), nrow=row_num, ncol=col_num, byrow=T) ) # Reformat into dataframe
		colnames(d_g_dat) = col_name
		
		m3 = glm(as.formula(f1), data=d_g_dat)
		agg = as.data.frame(coef(summary(m3)))
		agg = agg[-1,] 
		agg$Variable = rownames(agg)
		agg$Outcome = y
		agg$Method = 'AGG'
		agg$Conv = m3$converged
		rownames(agg) = 1:nrow(agg)
		colnames(agg)[1:4] = c('Estimate','Std.Err','Test.Stat','p.value')
		
		# FE
		f_fe = paste0(f2, '+ as.factor(group)')
		
		if(is.null(geno)){
			cols_fe = c('x','time','x:time')
		}else{cols_fe = c('x','time','geno','x:geno','x:time','geno:time','x:geno:time')}

		m4 = glm(as.formula(f_fe), data=d)
		fe = as.data.frame(coef(summary(m4))[cols_fe,]) 
		fe$Variable = rownames(fe)
		fe$Outcome = y
		fe$Method = 'FE'
		fe$Conv = m4$converged
		rownames(fe) = 1:nrow(fe)
		colnames(fe)[1:4] = c('Estimate','Std.Err','Test.Stat','p.value')
		
		# GEE
		m5 = geepack::geeglm(as.formula(f2), data=d[complete.cases(d),], id=d[complete.cases(d),'group'], corstr=GEEcor) # Only automatically handles missing Y data
		gee = as.data.frame(coef(summary(m5)))
		gee = gee[-1,] 
		gee$Variable = rownames(gee)
		gee$Outcome = y
		gee$Method = 'GEE'
		gee$Conv = summary(m5)$error
		rownames(gee) = 1:nrow(gee)
		colnames(gee)[1:4] = c('Estimate','Std.Err','Test.Stat','p.value')
		
		# LMM
		m6 = lmerTest::lmer(as.formula(f3), data=d)
		lme = as.data.frame(coef(summary(m6))[,-3]) # remove df column
		lme = lme[-1,] 
		lme$Variable = rownames(lme)
		lme$Outcome = y
		lme$Method = 'LMM'
		lme$Conv = ifelse(length(m6@optinfo$conv$lme4) > 0, m6@optinfo$conv$lme4$messages, 0)
		rownames(lme) = 1:nrow(lme)
		colnames(lme)[1:4] = c('Estimate','Std.Err','Test.Stat','p.value')
		
		out = rbind(nlr, crse, agg, fe, gee, lme)
		return(out)
		

	}
	
}
















