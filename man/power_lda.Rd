% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/power_lda.R
\name{power_lda}
\alias{power_lda}
\title{Apply multiple analytic methods to longitudinal data}
\usage{
power_lda(
  df = NULL,
  group = NULL,
  y = NULL,
  x = NULL,
  time = NULL,
  geno = NULL,
  rEff = "RI",
  GEEcor = "independence"
)
}
\arguments{
\item{df}{Dataframe. Dataframe object with data to be analyzed.}

\item{group}{String. Column name of grouping/cluster variable. Must be numerical.}

\item{y}{String. Column name of outcome variable aka "y".}

\item{x}{String. Column name of predictor variable aka "x".}

\item{time}{String. Column name of time variable aka "time".}

\item{geno}{String. Column name of genotype data aka "geno", if applicable.}

\item{rEff}{String. Random effect(s) for  the LMM. "RI" for
a random intercept-only model. "RIS" for a model with both
a random intercept and random slope. Default is "RI".}

\item{GEEcor}{String. Residual correlation structure for the GEE.
Must match one of the correlation structures implemented by
\code{geepack}. Default is 'independence'.}
}
\value{
Dataframe of analysis results
\itemize{
\item \strong{Estimate} - Coefficient estimate
\item \strong{Std.Err} - Coefficient standard error
\item \strong{Test.Stat} - Coefficient test statistic
\item \strong{p.value} - Coefficient p-value
\item \strong{Variable} - Predictor term in model
\item \strong{Outcome} - Outcome term in model
\item \strong{Method} - Statistical analysis method
\item \strong{Conv} - Convergence/error codes or messages.
}

In the results table, for NLR, AGG, CRSE,
and FE the convergence code is 1 if converged and 0 if not. For GEE,
it prints 0 if there are no errors, and prints the error message if
not For LMM, it likewise reports 0 for no errors, and prints the
warning message if not. Note that for GEE and LMM, the reported
warning may deal with something other than convergence
}
\description{
Outputs the results of six analytical methods applied to longitudinal
input data. If the data was simulated with known
parameters, the results can be used to compare power.
The input must be a dataframe with, at minimum, variables
corresponding to cluster, a numerical predictor, a time variable,
and a numerical outcome variable. Cluster IDs must be numbers.
As this is longitudinal data, a "cluster" refers to a given sample
in the data, with repeated measurements.
}
\details{
The following methods are applied to the data:

\strong{Naive linear regression (NLR):} Ignores clustering and time
when regressing outcome on predictor(s).

\strong{Cluster robust standard error (CRSE):} Applies a CRSE correction
to the standard errors of the aforementioned NLR output. Clustering is the group
membership.

\strong{Aggregate linear regression (AGG):} Averages all variables over
time within each cluster. Regresses the average outcome against the
average predictor(s).

\strong{Fixed effects model (FE):} Regresses the outcome on the predictor(s),
time, and group. Explicitly models time and group, but doesn't account for
within group dependency.

\strong{Generalized estimating equation (GEE):} Regresses outcome on
the predictor(s) and time, while accounting for group dependency in
the assumed residual correlation structure.

\strong{Linear mixed model (LMM):} Regresses outcome on the
predictor(s) and time, while accounting for group dependency
with random effects.

If no genotype is specified, the models tested are
\eqn{y = x} for NLR, CRSE, and AGG methods,
\eqn{y = x + time + x*time + group} for FE, and
\eqn{y = x + time + x*time} for GEE and LMM. For LMM, the
random effects are specified with an additional
\eqn{ + (1|group)} for a random intercept or
\eqn{ + (1 + time|group)} for a random intercept and slope.

If genotype is specified, the models are
\eqn{y = x + geno + x*geno} for NLR, CRSE, and AGG methods,
\eqn{y = x + geno + time + x*geno + x*time + geno*time + x*geno*time + group}
for FE, and
\eqn{y = x + geno + time + x*geno + x*time + geno*time + x*geno*time}
for FE, GEE, and LMM. The LMM random effects are specified
as described previously.

To fit models, CRSE uses \code{lmtest} and the clustered variance
estimator from \code{sandwich} using type "HC1". GEE uses \code{geepack}.
LMM uses \code{lmer} and \code{lmerTest}.

In randomly generated \code{x} has no
variation (e.g. can happen if there are few clusters and a binary
predictor with a high probability of success), instead of
reporting results for each method, the function returns a
data.frame with an error code \code{99999} as all entries
in the return value.
}
\examples{
df=data.frame(ids=rep(1:5, each=10), time=rep(1:10, times=5), x=rnorm(50), y=rnorm(50))
power_lda(df=df, group='ids', y='y', x='x', time='time', rEff='RI', GEEcor='exchangeable')
}
