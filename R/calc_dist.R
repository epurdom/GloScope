#' @title Compute statistical divergences between GloScope representations of
#'   samples
#'
#' @description The `.calc_dist` function calculates a statistical divergences
#'   between two specified samples of an input list of sample GloScope
#'   representations. Different subroutines are called depending on the
#'   divergence specified (symmetric KL or Jensen-Shannon) and density
#'   estimation used for the GloScope representations (GMM or KNN).
#'
#' @param mod_list A named list with each sample's estimated density
#' @param s1 The name or index of the first sample in a pair (must be a key in
#'   mod_list)
#' @param s2 The name or index of the second sample in a pair
#' @param df_list A named list with each sample's reduced dimension embedding
#' @param dist_mat The distance metric to use (KL or JS)
#' @param dens The density estimation method (GMM or KNN)
#' @param r Number of Monte Carlo simulations to generate
#' @param k Number of k nearest neighbours for KNN density estimation, default k
#'   = 50.
#' @param KNN_params a list of arguments for either `FNN::KL.dist` (KL) or `RANN::nn2` (JS)
#' @param varapp Boolean for using variation approximation of KL divergence;
#'   NOTE: Currently disabled
#' @param epapp Boolean for applying an epsilon perturbation to MC calculated
#'   KL, default = FALSE
#' @param ep Epsilon perturbation size to add to MC KL divergence calculation,
#'   default = NA
#' @return The estimated statistical divergence between two GloScope
#'   represenations
#' @noRd
.calc_dist <- function(mod_list, s1, s2, df_list,
                dist_mat = c("KL","JS"), dens = c("GMM","KNN"), r, k,
                KNN_params, varapp = FALSE, epapp = FALSE, ep = NA){
    dens<-match.arg(dens)
    dist_mat<-match.arg(dist_mat)
    if(dist_mat == "KL"){
        if(dens == "KNN"){
        if (typeof(KNN_params) == "list" & "k" %in% names(KNN_params)){
            stop("k cannot be specified in `KNN_params`. This is specified by `k` instead.")
        }
        KNN_params <- c(list(k=k),KNN_params)
            mydist <- .calc_kl(mod_list = mod_list, df_list = df_list, sample1 = s1, sample2 = s2,
                dens = dens, KNN_params = KNN_params)
        } else{
            # symmeterize by hand
            mydist <- .calc_kl(mod_list = mod_list, df_list = df_list, sample1 = s1, sample2 = s2,
                dens = dens, r = r, varapp = varapp, epapp = epapp, ep = ep) +
                .calc_kl(mod_list = mod_list, df_list = df_list, sample1 = s2, sample2 = s1,
                    dens = dens, r = r, varapp = varapp, epapp = epapp, ep = ep)
        }
    }  else if(dist_mat == "JS"){
        mydist <- .calc_JS(mod_list = mod_list, df_list = df_list, sample1 = s1, sample2 = s2,
            dens = dens, r = r, KNN_params = KNN_params)
    }

    return(mydist)
}

#' @title Calculate the KL divergence between a single pair of samples
#'
#' @description The `.calc_kl` function calculates the symmetric KL divergence
#'   between two GloScope representations. This is implemented with Monte Carlo
#'   approximation (with an optional epsilon perturbation term) if GMM is used
#'   for the density estimate of each cell or a plug-in formula if KNN is used
#'   for the density.
#'
#' @param df_list A named list with each sample's reduced dimension embedding
#' @param mod_list A named list with each sample's estimated density
#' @param sample1 The name or index of the first sample in a pair (must be a key
#'   in mod_list)
#' @param sample2 The name or index of the second sample in a pair
#' @param dens The density estimation method (GMM or KNN)
#' @param r Number of Monte Carlo simulations to generate
#' @param KNN_params a list of arguments to the FNN:KL.dist function
#' @param varapp Boolean for using variation approximation of KL divergence;
#'   NOTE: Currently disabled'
#' @param epapp Boolean for applying an epsilon perturbation to MC calculated
#'   KL, default = FALSE
#' @param ep Epsilon perturbation size to add to MC KL divergence calculation,
#'   default = NA
#' @return a numeric value of distance between sample1 and sample2's
#'   distribution.
#' @importFrom FNN KL.dist
#' @importFrom stats predict
#' @noRd
.calc_kl <- function(mod_list, df_list, sample1, sample2, dens, r = 10000 ,
            KNN_params = NULL, varapp=FALSE, epapp = FALSE, ep = NA){
    if(dens == "GMM"){
        mclust_mod1 <- mod_list[[sample1]]
        mclust_mod2 <- mod_list[[sample2]]
        s <- .sample_mclust(mclust_mod1, r=r)
        pi_1 <- mclust_mod1$parameters$pro
        pi_2 <- mclust_mod2$parameters$pro
        cov_1 <- mclust_mod1$parameters$variance$sigma
        cov_2 <- mclust_mod2$parameters$variance$sigma

        mu_1 <- mclust_mod1$parameters$mean
        mu_2 <- mclust_mod2$parameters$mean

        # Deprecated option for approximating KL
        if(varapp) {
            stop("A variational approximation to the KL divergence is not available at this time")
        } else {
            dens1 <- stats::predict(mclust_mod1, s, what = "dens", logarithm = TRUE)
            dens2 <- stats::predict(mclust_mod2, s, what = "dens", logarithm = TRUE)
            if(epapp) {
                kl <- sum(dens1 - (dens2+ep))/r
            } else {
                kl <- sum(dens1 - dens2) / r
            }
        }
    } else if(dens == "KNN"){
        KNN_params$X <- as.matrix(mod_list[[sample1]])
        KNN_params$Y <- as.matrix(mod_list[[sample2]])
        kl <- do.call(FNN::KL.dist,KNN_params)[KNN_params$k]
    }
    return(kl)
}

#' @title Calculate the Jensen-Shannon distance between a single pair of samples
#'
#' @description The `.calc_JS` function calculates the Jensen-Shannon distance
#'   between two GloScope representations. This is implemented with Monte Carlo
#'   approximation (with an optional epsilon perturbation term) if GMM is used
#'   for the density estimate of each cell or a plug-in formula if KNN is used
#'   for the density.
#'
#' @param df_list A named list with each sample's reduced dimension embedding
#' @param mod_list A named list with each sample's estimated density
#' @param sample1 The name or index of the first sample in a pair (must be a key
#'   in mod_list)
#' @param sample2 The name or index of the second sample in a pair
#' @param dens The density estimation method (GMM or KNN)
#' @param r Number of Monte Carlo simulations to generate
#' @param KNN_params a list of arguments for the `RANN::nn2` function
#' @return a numeric value of distance between sample1 and sample2's
#'   distribution.
#' @importFrom stats predict
#' @noRd
.calc_JS <- function(mod_list, df_list, sample1, sample2, dens,
            r , KNN_params){
    if(dens == "GMM"){
        mclust_mod1 <- mod_list[[sample1]]
        mclust_mod2 <- mod_list[[sample2]]
        s1 <- .sample_mclust(mclust_mod1, r=r)
        s2 <- .sample_mclust(mclust_mod2, r=r)
        dens1_1 <- stats::predict(mclust_mod1, s1, what = "dens", logarithm = TRUE)
        dens2_1 <- stats::predict(mclust_mod2, s1, what = "dens", logarithm = TRUE)
        mixture_1 <- log(1/2*exp(dens1_1) + 1/2*exp(dens2_1))
        dens1_2 <- stats::predict(mclust_mod1, s2, what = "dens", logarithm = TRUE)
        dens2_2 <- stats::predict(mclust_mod2, s2, what = "dens", logarithm = TRUE)
        mixture_2 <- log(1/2*exp(dens1_2) + 1/2*exp(dens2_2))
        js <- sum(dens1_1 - mixture_1)/(2*r) + sum(dens2_2 - mixture_2)/(2*r)
    }else if(dens == "KNN"){
        k <- KNN_params$k
        knn1_1 <- .knn_query(df_list, input = sample1, query = sample1, KNN_params = KNN_params)
        knn2_2 <- .knn_query(df_list, input = sample2, query = sample2, KNN_params = KNN_params)
        knn2_1 <- .knn_query(df_list, input = sample2, query = sample1, KNN_params = KNN_params)
        knn1_2 <- .knn_query(df_list, input = sample1, query = sample2, KNN_params = KNN_params)
        knn1 <- mod_list[[sample1]]
        knn2 <- mod_list[[sample2]]
        if(dim(knn1)[2] != dim(knn2)[2]){stop("This method assumes both densities are of equal dimension.")}
        ndim <- dim(knn1)[2]
        js <- 1/(2*dim(knn1)[1])*sum(log(2*dim(knn2)[1] * knn2_1^ndim/(dim(knn2)[1] * knn2_1^ndim +
                (dim(knn1)[1]-1) * knn1_1^ndim))) +
            1/(2*dim(knn2)[1])*sum(log(2*dim(knn1)[1] * knn1_2^ndim/(dim(knn1)[1] * knn1_2^ndim +
                (dim(knn2)[1]-1) * knn2_2^ndim)))
    }
    return(js)
}
