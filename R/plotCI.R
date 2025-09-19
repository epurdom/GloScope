#' @title Plot confidence intervals
#'
#' @description This function creates a `ggplot` object that plots the
#'   confidence intervals created by `bootCI_gloscope`
#'
#' @param ci_df A data.frame that is the output of \code{\link{bootCI_gloscope}}
#' @param ci_df A data frame contains each sample's metadata. Note this is
#'   NOT at the cell-level.
#' @param color_by The column name or index in ci_df that should be used
#'   to color the confidence intervals by. 
#' @param group_by The column name or index in ci_df that should be used to
#'   determine how to group the confidence intervals. If missing all confidence
#'   intervals will be plotted in an order determined internally.
#' @param dodge_width value passed to `width` argument of
#'   \code{\link[ggplot2]{position_dodge}} if `group_by` variable is given.
#'   Controls separation between confidence intervals with the same grouping
#'   value.
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
#' sample_metadata$fakeGroup<-c(rep(c("A","B"),each=2),"A")
#' manyboot<-bootCI_gloscope(dist_result,
#'   sample_metadata,"sample_id",
#'   metrics=c("anosim","silhouette"),group_vars=c("phenotype","fakeGroup"),R=20)
#' plotCI(manyboot,group_by="metric",color_by="grouping")
#' 
#' @export
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_pointrange
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 position_dodge

plotCI<-function(ci_df,color_by,group_by,dodge_width=.5){

  if(!all(c('distance','grouping','metric','statistic','lower','upper') %in% names(ci_df)))
    stop("names of ci_df must match the expected output of bootCI_gloscope")
  whCols<-c("metric","distance","grouping")
  listCols<-paste(whCols,collapse=",")
  whMulti<-sapply(ci_df[,whCols,drop=FALSE],function(x){nlevels(factor(x))>1})
  if(!missing(color_by)){
    if(!c(color_by) %in% whCols) stop("color_by must be one of: ",listCols)
    if(!color_by %in% whCols[whMulti]) stop("color_by must have more than one distinct value")
  }
  if(!missing(group_by)){
    if(!c(group_by) %in% whCols) stop("group_by must be one of: ",listCols)
    if(!group_by %in% whCols[whMulti]) stop("group_by must have more than one distinct value")
  }
  if(missing(group_by)){
    pos="identity"
    if(sum(whMulti)>0){
      ci_df$ConfInt <- apply(ci_df[,whCols][,whMulti,drop=FALSE],1,paste,collapse=":")
      ci_df$ConfInt <- factor(ci_df$ConfInt)
    }
    else{ci_df$ConfInt<-""}
    group_by="ConfInt"
    if(missing(color_by)){
      ggp <- ggplot(ci_df,
                    aes(x = .data[[group_by]], y = .data$statistic)
      )
    }
    else{
      ggp <- ggplot(ci_df,
                      aes(x = .data[[group_by]], y = .data$statistic,
                          color = .data[[color_by]])
        )
    }
  }
  else{
    if(sum(whMulti)>2) stop("Cannot accurate group confidence intervals if they vary by more than two of:",listCols)
    pos <- position_dodge(width = dodge_width)
    if(missing(color_by)){
      ggp <- ggplot(ci_df,
                    aes(x = .data[[group_by]], y = .data$statistic)
      )
    }
    else{
      ggp <- ggplot(ci_df,
                    aes(x = .data[[group_by]], y = .data$statistic,
                        color = .data[[color_by]])
      )
    }
  }  
  ggp <- ggp + geom_pointrange(aes(ymin = .data$lower, ymax = .data$upper),
                  position = pos, size = 0.9)
  return(ggp)
}

