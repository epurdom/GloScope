#' @rdname getMetrics
#' @details The function `bootGloscope` is a wrapper to the
#'   \code{\link[boot]{boot}} function for creating bootstraps of one of the
#'   metrics calculated by `getMetrics`. Most users will probably prefer `bootCI`
#' @param R number of bootstrap replicates. See \code{\link[boot]{boot}}.
#' @param ... arguments passed to \code{\link[boot]{boot}}
#' @return `bootGloscope` returns an object of class `boot` created by \code{\link[boot]{boot}}.
#' @seealso \code{\link[boot]{boot}}
#' @examples
#' # single bootstrap of anosim
#' bootout<-bootGloscope(dist_result,sample_metadata,"sample_id",
#'   metric="anosim",group_var="phenotype")
#' #work with the boot object using functions in boot package:
#' library(boot)
#' print(bootout)
#' boot.ci(bootout)
#' 
#' @importFrom boot boot
#' @export
#' 
bootGloscope<-function(dist_mat, metadata_df, metrics="anosim",sample_id, group_vars,R=1000,...){
  metadata_df<-.testDistMeta(dist_mat,metadata_df,sample_id)
  if(!c(group_vars) %in% names(metadata_df)) stop("group_vars does not define a variable in metadata_df")
  combinedData<-cbind(dist_mat,metadata_df)
  nsamples<-ncol(dist_mat)
  if(length(group_vars)>1) stop("group_vars must be single character value")
  if(length(metrics)>1) stop("metrics must be single character value")
  bootfun<-function(df,i){#x are ids
    d<-df[,1:nsamples]
    m<-df[,(nsamples+1):ncol(df)]
    getMetrics(d[i,i], m[i,], metrics=metrics,sample_id=sample_id, 
                group_vars=group_vars,checkData=FALSE, permuteTest=FALSE)$statistic
  }
  return(boot::boot(data=combinedData,statistic=bootfun,R=R,...))
    
}

#' @rdname getMetrics
#' @order 2
#' @param ci_type Single character value. The type of confidence interval to
#'   compute. Passed to argument `type` in \code{\link[boot]{boot.ci}}
#' @param ci_conf Scalar value between 0 and 1. The confidence level requested.
#'   Passed to argument `conf` in \code{\link[boot]{boot.ci}}.
#' @details `bootCI` is a wrapper function to
#'   \code{\link[boot]{boot.ci}}. `boot.ci` can be called directly on the output
#'   of `bootGloscope`. The main advantage of `bootCI` is to calculate
#'   bootstrap CI over multiple choices of metrics, variables, and/or distance
#'   matrices. Unlike `boot.ci`, `bootCI` does not allow different
#'   choices of confidence interval types or levels, so `ci_type` and `ci_level`
#'   must be of length 1. For this kind of multiplicity, call `boot.ci` directly
#'   on the output of `bootGloscope`.
#' @return `bootCI` creates a data frame containing the statistic for
#'   each combination of metric and grouping variable with columns with the
#'   upper and lower bounds of the requested confidence intervals
#' \itemize{
#'   \item metric 
#'   \item grouping
#'   \item statistic
#'   \item lower
#'   \item upper
#' }
#' @seealso \code{\link[boot]{boot.ci}}
#' @examples
#' # calculate many bootstraps -- for speed up we set R ridiculously low
#' manyboot<-bootCI(list("Distance 1"=dist_result,"Another distance"=dist_result),
#'   sample_metadata,"sample_id",
#'   metrics=c("anosim","silhouette"),group_vars=c("phenotype","grouping"),R=20)
#' 
#' @importFrom boot boot.ci
#' @export
bootCI<-function(dist_mat, metadata_df, metrics="anosim",sample_id, group_vars,R=1000,ci_type=c("perc","norm","basic", "stud",  "bca"),ci_conf=0.95,...){
  ci_type=match.arg(ci_type)
  if(length(ci_conf)!=1) stop("only single value of ci_conf is allowed")
  if(!inherits(dist_mat,"list")){
    dist_mat<-list(dist_mat)
  }
  if(inherits(dist_mat,"list") & is.null(names(dist_mat))){
    names(dist_mat)<-paste("Distance",1:length(dist_mat),sep="")
  }
  if(length(ci_type)!=1) stop("only a single type of confidence interval can be done with this wrapper function")
  comb<-expand.grid(distName=names(dist_mat),g=group_vars,m=metrics)
  perComboFunction<-function(distName,g,m){
    distName<-as.character(distName)
    g<-as.character(g)
    m<-as.character(m)
    out<-bootGloscope(dist_mat[[distName]],metadata_df,sample_id, metrics=m,group_vars=g,...)
    out.ci<-boot.ci(out,type=ci_type,conf=ci_conf)
    ci_type<-match.arg(ci_type,names(out.ci)) #make sure have the full name
    nCols<-ncol(out.ci[[ci_type]])
    return(data.frame(distance=distName,grouping=g,metric=m,statistic=out.ci$t0,lower=out.ci[[ci_type]][,nCols-1],upper=out.ci[[ci_type]][,nCols]))
    
  }
  do.call("rbind",do.call("mapply",c(list(FUN=perComboFunction,SIMPLIFY=FALSE),comb)))
}
