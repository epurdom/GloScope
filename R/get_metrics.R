

#' @title Calculate metrics on separation in distance matrix due to grouping variable
#'
#' @description These functions are wrappers for calculating common metrics for
#'   the amount of separation in a distance matrix due to a grouping (factor)
#'   variable and creating bootstrap confidence intervals and permutation tests.
#' 
#' @order 1
#' 
#' @param dist_mat The divergence matrix output of `gloscope()`. Should be a
#'   symmetric, square matrix. For `bootCI` the argument can be a list
#'   of distance matrices.
#' @param metadata_df A data frame contains each sample's metadata. Note this is
#'   NOT at the cell-level, and should have the same number of rows as dist_mat.
#' @param sample_id The column name or index in metadata_df that contains the
#'   sample ID. This is for ensuring alignment between the dist_mat and the
#'   metadata_df. The rownames of dist_mat are expected to match the sample_id
#'   values.
#' @param metrics vector of statistics to calculate. For `bootstrap_gloscope`
#'   must be single value.
#' @param group_vars vector of names of grouping variables in metadata_df for
#'   which to calculate metrics. For `bootstrap_gloscope` must be single value.
#' @param checkData Whether to check whether dist_mat, metadata_df, and
#'   sample_id match, for example in terms of dimensions and rownames. Mainly
#'   used internally.
#' @param permuteTest whether to run permutation tests on each of the metrics
#' @param permutations if `permuteTest=TRUE`, an integer value defines the
#'   number of permutations. Can also except output of
#'   \code{\link[permute]{how}} for more fine control of the permutation
#'   mechanisms.
#' @return `getMetrics` creates a data frame containing the statistic for each
#'   combination of metric and grouping variable with columns
#' \itemize{
#'   \item metric 
#'   \item grouping
#'   \item statistic
#'   \item pval (if `permuteTest=TRUE`)
#' }
#'
#' @details The function `getMetrics` is a simple wrapper for calculating statistics that
#'   summarize the difference between distances within and between groupings. If
#'   the variable defined by group_var does not have at least two groupings, the
#'   function will return a NA. 
#' @details The options "anosim" and "adonis2" are wrappers to the functions of
#'   that name in the package `vegan`; we have turned off the permutation
#'   testing option of those functions. The functions in `vegan` have greater
#'   capability, and in particular \code{\link[vegan]{adonis2}} has capability
#'   to handle more complicated testing paradigms than a simple grouping factor.
#'   Permutation tests for these statistics are handled by the functions of
#'   `vegan`.
#' @details The option "silhouette" calls the function of that name from the
#'   package `cluster` and calculates the silhouette width of each group and
#'   then averages them across groups. The permutation test is coded making use
#'   of the package `permute`, similar to `vegan`, so that control of the
#'   permutation mechanism is possible in the same way.
#' @seealso \code{\link[vegan]{anosim}}, \code{\link[vegan]{adonis2}}, 
#'  \code{\link[cluster]{silhouette}}, \code{\link[permute]{how}}
#' @examples
#' data(example_SCE_small)
#' sample_ids <- SingleCellExperiment::colData(example_SCE_small)$sample_id
#' # Run gloscope on first 10 PCA embeddings
#' # We use 'KNN' option for speed ('GMM' is slightly slower)
#' pca_embeddings <- SingleCellExperiment::reducedDim(example_SCE_small,"PCA")
#' pca_embeddings_subset <- pca_embeddings[,seq_len(10)] # select the first 10 PCs
#' dist_result <- gloscope(pca_embeddings_subset, sample_ids,
#'    dens="KNN",
#'    BPPARAM = BiocParallel::SerialParam(RNGseed=2))
#' # make a per-sample metadata
#' sample_metadata <- as.data.frame(unique(SingleCellExperiment::colData(example_SCE_small)[,c(1,2)]))
#' # make another variable
#' sample_metadata$grouping<-c(rep(c("A","B"),each=2),"A")
#' getMetrics(dist_result,metadata_df=sample_metadata, sample_id="sample_id", 
#'   group_vars="phenotype")
#' # run permutation tests:
#' getMetrics(dist_result,metadata_df=sample_metadata, sample_id="sample_id", 
#'   group_vars=c("phenotype","grouping"), permuteTest=TRUE)
#'   
#' @importFrom cluster silhouette
#' @importFrom vegan anosim
#' @importFrom permute how 
#' @importFrom permute permute
#' @importFrom vegan adonis2
#' @export
#'
getMetrics <- function(dist_mat, metadata_df, metrics=c("anosim","adonis2","silhouette"),
                        sample_id, group_vars, checkData=TRUE, permuteTest=FALSE,
                        permutations=100) {
  metrics<-match.arg(metrics,several.ok = TRUE)
  if(any(!c(sample_id,group_vars) %in% names(metadata_df))) stop("sample_id or group_vars does not define a variable in metadata_df")
  # test that at least 2 levels of each factor
  testnlevels<-sapply(group_vars,function(x){nlevels(factor(x))>1})
  # turn into dist for methods; align order using labels
  m <- .as_full_matrix(dist_mat)
  labs <- rownames(m)
  stopifnot(!is.null(labs))
  if(checkData) metadata_df<-.testDistMeta(dist_mat,metadata_df,sample_id)
  metadata_df[,sample_id] <- as.character(metadata_df[,sample_id])
  d <- as.dist(m)
  
  # to turn off permutation calculations:
  if(!permuteTest){
    permControl<-permute::how(
      minperm=0, maxperm=2,nperm=0
    )
  }else{
    permControl=permutations
  }
  combinedData<-cbind(dist_mat,metadata_df)
  nsamples<-ncol(dist_mat)
  
  getanosim<-function(varname){
    f<-factor(metadata_df[,varname])
    if(nlevels(f)<2) return(NA)
    else{
      out<-vegan::anosim(d,f ,permutations=permControl)
      if(!permuteTest) return(out$statistic)
      else return(c(statistic=out$statistic, pval=out$signif))
    }
  }
  getadonis2<-function(varname){
    f<-factor(metadata_df[,varname])
    if(nlevels(f)<2) return(NA)
    else{
      out<-vegan::adonis2(d~f,perm=permControl)
      if(!permuteTest) return(out$F[1])
      else{ return(c(statistic=out$F[1],pval=out$`Pr(>F)`[1]))}
    }      
  }
  
  getsil<-function(varname){
    if(nlevels(factor(metadata_df[,varname]))<2) return(NA)
    else{
      meansil<-function(data,ind){
        d<-data[,1:nsamples]
        m<-data[,(nsamples+1):ncol(data)]
        fac<-factor(m[,varname])
        sil<-cluster::silhouette(as.integer(factor(fac)), d)
        return(mean(sil[,"sil_width"]))
      }
      if(!permuteTest){
        return(meansil(combinedData,ind=1:nsamples))
      }
      else{
        if(!inherits(permControl,"how")) permControl<-how(nperm=permControl)
        return(.permGeneric(data=combinedData,nobs=nsamples,statistic=meansil,control=permControl))
      }
    }
  }
  
  combs<-expand.grid(metric=metrics,grouping=group_vars)
  out<-apply(as.matrix(combs),1,function(x){
    metric_x<-as.character(x[[1]])
    var_x<-as.character(x[[2]])
    if(metric_x=="anosim") return(getanosim(var_x))
    if(metric_x=="adonis2") return(getadonis2(var_x))
    if(metric_x=="silhouette") return(getsil(var_x))
  })
  if(!permuteTest) combs$statistic<-out
  else combs<-cbind(combs,t(out))
  return(combs)
}


.permGeneric<-function(data,nobs,statistic,control){
  ## generate randomisation distribution 
  stat.obs<-statistic(data,ind=1:nobs)
  stat.permu<-numeric(length = control$nperm)
  for(i in seq_along(control$nperm)){
    ## return a permutation
    want <- permute::permute(i, nobs, control=control)
    ## calculate permute stat
    stat.permu[i] <-  statistic(data,ind=want)
  }
  ## pval from permutation test
  pval <- sum(abs(c(stat.permu,stat.obs)) >= abs(stat.obs))/(control$nperm+1)
  ## return value
  return(c(statistic = stat.obs, pval = pval))
}