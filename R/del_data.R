#' Create missing data
#'
#' Deletes records from longitudinal data frame 
#' in a random or non-random manner. For random deletion,
#' function deletes a designated percent of rows at random.
#' If deletion is random, the function specifies a subset
#' of groups to have missing data, and then deletes a 
#' designated percent of rows from only those groups. 
#' (Note: In longitudinal data, a "group" is an individual
#' with repeated measurements). Function designed so that
#' at least one record (row) will be deleted.
#'
#' @param df Data frame. Data frame of longitudinal
#' data.
#' @param mType Numeric. Type of missingness. 1 = random 
#' missingness; 2 = non-random missingness. Default is 1.
#' @param mPer Numeric. Proportion of missing data. Range 
#' (0, 1). If mType = 1, then randomly deletes mPer*100%
#' of rows from who data frame. If mType = 2, deletes
#' mPer*100% rows from the groups selected to 
#' have missingness. Default is 0.5.
#' @param propDel Numeric. Range (0,1). Only applies
#' if mType=2. Denotes proportion of groups from which 
#' data is to be missing. Default is 0.1.
#' 
#' @returns Data frame with rows deleted
#' @examples
#' del_data(df, mType=1, mPer=0.05)
#' del_data(df, mType=2, mPer=0.1, propDel=0.5)

del_data = function(df, mType=1, mPer=0.5, propDel=0.1)	{
	
	if(mType==1){
		num_miss = round(nrow(df)*mPer) # detemine total number of missing records 
		if(num_miss==0){num_miss=1} # ensure at least 1 missing record
		index = sample(rownames(df), size = num_miss, replace=F) # select rows to be deleted
 		index = index[order(index)] # delete rows
		df[index,] = NA
		return(df)
		
	}
	
	if(mType==2){
		ids = unique(df$id) # get list of cluster ids
		ids_miss = round(length(ids)*propDel) # get number of ids to have missingness
		if(ids_miss==0){ids_miss=1} # ensure at least 1 id with missing data
		ids_list = sample(ids, size=ids_miss, replace=F) # select ids w/missingness
		index = which(df$id%in%ids_list) # row indices of ids w/missingness
		num_miss = round(length(index)*mPer) # determine total missing records(rows) w/in subset
		if(num_miss==0){num_miss=1} # ensure at least 1 missing record(row)
		index2 = sample(index, size = num_miss, replace=F) # select rows to be deleted
		df[index2,] = NA # delete rows
		return(df)	
	}
	
}











