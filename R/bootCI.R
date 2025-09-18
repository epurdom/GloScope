#' @rdname get_metrics
#' @details The function `boot_gloscope` is a wrapper to the
#'   \code{\link[boot]{boot}} function for creating bootstraps of one of the
#'   metrics calculated by `get_metrics`.
#' @param R number of bootstrap replicates. See \code{\link[boot]{boot}}.
#' @param ... arguments passed to \code{\link[boot]{boot}}
#' @return `boot_gloscope` returns an object created by \code{\link[boot]{boot}}.
#' @seealso \code{\link[boot]{boot}}
#' @examples
#' # bootstrap anosim (nonsensical for only 5 samples!) 
#' library(boot)
#' bootout<-boot_gloscope(dist_result,sample_metadata,"sample_id",
#'   metric="anosim",group_var="phenotype")
#' print(bootout)
#' boot.ci(bootout)
#' manyboot<-bootCI_gloscope(dist_result,sample_metadata,"sample_id",
#'   metric=c("anosim","silhouette"),group_var="phenotype")
#' @importFrom boot boot
#' @export
#' 
boot_gloscope<-function(dist_mat, metadata_df, metrics="anosim",sample_id, group_vars,R=1000,...){
  metadata_df<-.testDistMeta(dist_mat,metadata_df,sample_id)
  if(!c(group_vars) %in% names(metadata_df)) stop("group_vars does not define a variable in metadata_df")
  combinedData<-cbind(dist_mat,metadata_df)
  nsamples<-ncol(dist_mat)
  if(length(group_vars)>1) stop("group_vars must be single character value")
  if(length(metrics)>1) stop("metrics must be single character value")
  bootfun<-function(df,i){#x are ids
    d<-df[,1:nsamples]
    m<-df[,(nsamples+1):ncol(df)]
    get_metrics(d[i,i], m[i,], metrics=metrics,sample_id=sample_id, 
                group_vars=group_vars,checkData=FALSE, permuteTest=FALSE)$statistic
  }
  return(boot::boot(data=combinedData,statistic=bootfun,R=R,...))
    
}

bootCI_gloscope<-function(dist_mat, metadata_df, metrics="anosim",sample_id, group_vars,R=1000,ci_type="percent",ci_conf=0.95,...){
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
    out<-boot_gloscope(dist_mat[[distName]],metadata_df,sample_id, metric=m,group_var=g,...)
    out.ci<-boot.ci(out,type=ci_type,conf=ci_conf)
    ci_type<-match.arg(ci_type,names(out.ci)) #make sure have the full name
    nCols<-ncol(out.ci[[ci_type]])
    return(data.frame(distance=distName,grouping=g,metric=m,statistic=out.ci$t0,lower=out.ci[[ci_type]][,nCols-1],upper=out.ci[[ci_type]][,nCols]))
    
  }
  do.call("rbind",do.call("mapply",c(list(FUN=perComboFunction,SIMPLIFY=FALSE),comb)))
}
