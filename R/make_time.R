#' Define timepoints
#'
#' Creates the timepoints (e.g. waves, measurement occasions)
#' for simulating longitudinal data.
#' 
#' There are three options for time interval spacing.
#' \itemize{
#'	\item \emph{Equal} - Time intervals are equal and even.
#' All individuals have the same timepoints.
#'	\item \emph{Increasing} - Time intervals increase with
#' each wave. All individuals have the same timepoints.
#'	\item \emph{Random} - Time intervals are random, selected
#' from a uniform distribution. Each individual has 
#' different timepoints 
#' }
#'
#' @param clus Numeric. Number of samples/clusters.
#' @param waves Numeric. Number of timepoints
#' @param tType Numeric. Type of interval spacing. If tType=1,
#' equal intervals. If tType=2, increasing intervals.
#' If tType=3, random intervals. Default is 1.
#' @param incr Numeric. The increment of time spacing when 
#' tType=1 or tType=3. Default is 1.
#' @param mint Numeric. Lower limit of time, if tType=3.
#' Must be greater than 0, but less than maxt. Default is 1. 
#' @param maxt Numeric. Upper limit of time, if tType=3.
#' Must be greater than mint. Default is 10.
#' 
#' @returns Vector of timepoints, grouped by cluster.
#' @examples
#' make_time(clus=10, waves=3, tType=1)
#' make_time(clus=10, waves=3, tType=2, incr=5)
#' make_time(clus=10, waves=3, tType=3)
#' @export

# type of time: {1:balanced, 2:increasing, 3:random}
make_time = function(clus, waves, tType=1, incr=1, mint=1, maxt=10){
	
	if(tType==1){ # balanced intervals, same for all clusters
		t = seq(by=incr, length.out=waves)
		all_tpts = rep(t, times = clus)
		return(all_tpts)
	}
	
	if(tType==2){ # increasing intervals, same for all clusters
		w = waves-2
		its = 1 + (0:w)*incr
		t = sapply(1:(waves-1), function(x) 1+sum(its[1:x]))
		t = c(1, t)
		all_tpts = rep(t, times = clus)
		return(all_tpts)
	}
		
	if(tType==3){ # random intervals
		all_tpts  = runif(n=waves*clus, min=mint, max=maxt)
		return(all_tpts)
	}
	
		
}


