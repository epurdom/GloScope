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
#' 
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


