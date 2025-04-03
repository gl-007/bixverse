# preprocessing ----------------------------------------------------------------

#' @title Prepare class for ICA
#'
#' @description
#' This is the generic function for doing the necessary preprocessing for
#' running independent component analysis.
#'
#' @param object The class, see [bixverse::bulk_coexp()]. Ideally, you
#' should run [bixverse::preprocess_bulk_coexp()] before applying this function.
#' @param fast_svd Boolean. Shall randomised SVD be used for the whitening.
#' This is faster and usually causes little precision loss.
#' @param random_seed Integer. Seed for the randomised SVD. Only relevant, if
#' fast_svd = `TRUE`.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return `bulk_coexp` with the needed data for ICA in the
#' properties of the class.
#'
#' @export
ica_processing <- S7::new_generic(
  name = "ica_processing",
  dispatch_args = "object",
  fun = function(object,
                 fast_svd = TRUE,
                 random_seed = 123L,
                 .verbose = TRUE) {
    S7::S7_dispatch()
  }
)


#' @export
#'
#' @importFrom magrittr `%>%`
#' @importFrom magrittr `%$%`
#' @importFrom zeallot `%<-%`
#' @import data.table
#'
#' @method ica_processing bulk_coexp
S7::method(ica_processing, bulk_coexp) <- function(object,
                                                   fast_svd = TRUE,
                                                   random_seed = 123L,
                                                   .verbose = TRUE) {
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  checkmate::qassert(fast_svd, "B1")
  checkmate::qassert(random_seed, "I1")
  checkmate::qassert(.verbose, "B1")

  # Function body
  if (purrr::is_empty(S7::prop(object, "processed_data")[['processed_data']])) {
    warning("No pre-processed data found. Defaulting to the raw data")
    target_mat <- S7::prop(object, "raw_data")
  } else {
    target_mat <- S7::prop(object, "processed_data")[['processed_data']]
  }

  # Whiten the data
  if (.verbose)
    message("Preparing the whitening of the data for ICA.")
  if (fast_svd && .verbose)
    message("Using randomised SVD for whitening (faster, but less precise)")
  if (!fast_svd && .verbose)
    message("Using full SVD for whitening (slowed, but more precise)")
  c(X1, K) %<-% rs_prepare_whitening(
    x = target_mat,
    fast_svd = fast_svd,
    seed = random_seed,
    rank = nrow(target_mat),
    oversampling = NULL,
    n_power_iter = NULL
  )

  S7::prop(object, "processed_data")[["X1"]] <- X1
  S7::prop(object, "processed_data")[["K"]] <- K
  S7::prop(object, "params")["detection_method"] <- "ICA-based"

  return(object)
}

# component identification -----------------------------------------------------

#' @title Iterate over different ncomp parameters for ICA
#'
#' @description
#' This function allows to iterate over a vector of ncomp to identify which
#' ncomp parameter to choose for your data set. The idea is to generate stability
#' profiles over the different ncomps and identify a 'sweet spot' of good
#' stability of the identified independent components
#'
#' @param object The class, see [bixverse::bulk_coexp()]. You need to apply
#' [bixverse::ica_processing()] before running this function.
#' @param ica_type String, element of `c("logcosh", "exp")`.
#' @param iter_params List. This list controls the randomisation parameters for
#' the ICA runs, see [bixverse::params_ica_randomisation()] for estimating
#' stability. Has the following elements:
#' \itemize{
#'  \item cross_validate - Boolean. Shall the data be split into different
#'  chunks on which ICA is run. This will slow down the function substantially,
#'  as every chunk needs to whitened again.
#'  \item random_init - Integer. How many random initialisations shall be used
#'  for the ICA runs.
#'  \item folds - If `cross_validate` is set to `TRUE` how many chunks shall be
#'  used. To note, you will run per ncomp random_init * fold ICA runs which
#'  can quickly increase.
#' }
#' @param ncomp_params List. Parameters for the ncomp to iterate through, see
#' [bixverse::params_ica_ncomp()]. In the standard setting, `c(2, 3, 4, 5)` will
#' be tested and then in steps until max_no_comp will be tested, i.e.,
#' `c(2, 3, 4, 5, 10, 15, ..., max_no_comp - 5, max_no_comp)`.
#' \itemize{
#'  \item max_no_comp - Maximum number of ncomp to test.
#'  \item steps - Integer. In which steps to move from 5 onwards.
#'  \item custom_seq - An integer vector. If you wish to provide a custom version
#'  of no_comp to iterate through.
#' }
#' @param ica_params List. The ICA parameters, see [bixverse::params_ica_general()]
#' wrapper function. This function generates a list containing:
#' \itemize{
#'  \item maxit - Integer. Maximum number of iterations for ICA.
#'  \item alpha - Float. The alpha parameter for the logcosh version of ICA.
#'  Should be between 1 to 2.
#'  \item max_tol - Maximum tolerance of the algorithm.
#'  \item verbose - Controls verbosity of the function.
#' }
#' @param random_seed Integer. For reproducibility.
#' @param .verbose Boolean. Controls verbosity.
#'
#' @return `bulk_coexp` with the added information of stability of the components
#' and other data to plot to choose the right `ncomp`.
#'
#' @export
ica_evaluate_comp <- S7::new_generic(
  name = "ica_evaluate_comp",
  dispatch_args = "object",
  fun = function(object,
                 ica_type = c("logcosh", "exp"),
                 iter_params = params_ica_randomisation(),
                 ncomp_params = params_ica_ncomp(),
                 ica_params = params_ica_general(),
                 random_seed = 42L,
                 .verbose = TRUE) {
    S7::S7_dispatch()
  }
)

#' @import data.table
#' @importFrom magrittr `%>%`
#' @importFrom zeallot `%<-%`
#'
#' @export
#'
#' @method ica_evaluate_comp bulk_coexp
S7::method(ica_evaluate_comp, bulk_coexp) <- function(object,
                                                      ica_type = c("logcosh", "exp"),
                                                      iter_params = params_ica_randomisation(),
                                                      ncomp_params = params_ica_ncomp(),
                                                      ica_params = params_ica_general(),
                                                      random_seed = 42L,
                                                      .verbose = TRUE) {
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  checkmate::assertChoice(ica_type, c("logcosh", "exp"))
  assertIcaParams(ica_params)
  assertIcaNcomps(ncomp_params)
  assertIcaIterParams(iter_params)
  checkmate::qassert(.verbose, "B1")

  detection_method <- S7::prop(object, "params")[["detection_method"]]

  # Early return
  if (is.null(detection_method) &&
      detection_method != "ICA-based") {
    warning(
      paste(
        "This class does not seem to be set for ICA-based module detection",
        "Returning class as is."
      )
    )
    return(object)
  }

  # Get data
  X <- S7::prop(object, "processed_data")[['processed_data']]
  X1 <- S7::prop(object, "processed_data")[["X1"]]
  K <- S7::prop(object, "processed_data")[["K"]]

  # Prepare the n_comp vector
  n_comp_vector <- if (is.null(ncomp_params$custom_seq)) {
    with(ncomp_params, c(2, 3, 4, seq(
      from = 5, to = max_no_comp, by = steps
    )))
  } else {
    ncomp_params$custom_seq
  }

  # TODO Control the max_no_comp rate pending the number of samples.

  if (.verbose)
    message(sprintf(
      "Using a total of %i different n_comp parameters",
      length(n_comp_vector)
    ))

  # Set up the loop
  if (iter_params$cross_validate) {
    total_randomisations <- iter_params$random_init * iter_params$folds
    no_ica_runs <- total_randomisations * length(n_comp_vector)
    if (.verbose)
      message(
        sprintf(
          "Using a CV-like approach with %i folds and %i random initialisations for a total of %i ICA runs",
          iter_params$folds,
          iter_params$random_init,
          no_ica_runs
        )
      )
  } else {
    total_randomisations <- iter_params$random_init
    no_ica_runs <- iter_params$random_init * length(n_comp_vector)
    if (.verbose)
      message(
        sprintf(
          "Using %i random initialisations for a total of %i ICA runs",
          iter_params$random_init,
          no_ica_runs
        )
      )
  }

  all_scores <- c()
  all_convergence <- c()

  pb = txtProgressBar(initial = 0,
                      max = length(n_comp_vector),
                      style = 3)

  for (i in seq_along(n_comp_vector)) {
    no_comp <- n_comp_vector[[i]]
    # Get the combined S matrix and convergence information
    c(s_combined, converged) %<-% with(iter_params, switch(
      as.integer(iter_params$cross_validate) + 1,
      rs_ica_iters(
        x1 = X1,
        k = K,
        no_comp = no_comp,
        no_random_init = random_init,
        ica_type = ica_type,
        random_seed = random_seed,
        ica_params = ica_params
      ),
      rs_ica_iters_cv(
        x = X_raw,
        no_comp = no_comp,
        no_folds = folds,
        no_random_init = random_init,
        ica_type = ica_type,
        random_seed = random_seed,
        ica_params = ica_params
      )
    ))

    c(stability_scores, centrotype) %<-% .community_stability(
      no_comp = as.integer(no_comp),
      s = s_combined,
      return_centrotype = FALSE
    )

    all_scores <- append(all_scores, sort(stability_scores, decreasing = TRUE))
    all_convergence <- append(all_convergence, converged)

    setTxtProgressBar(pb, i)
  }

  close(pb)

  ica_comps_rep <- purrr::map(n_comp_vector, \(x) {
    rep(x, x)
  })
  ica_comps_rep <- do.call(c, ica_comps_rep)
  ica_comps_no <- purrr::map(n_comp_vector, \(x) {
    seq_len(x)
  })
  ica_comps_no <- do.call(c, ica_comps_no)

  convergence_split <- split(all_convergence, ceiling(seq_along(all_convergence) /
                                                        total_randomisations))
  names(convergence_split) <- n_comp_vector

  ica_stability_res <- list(
    component_rank = ica_comps_no,
    no_components = ica_comps_rep,
    stability = all_scores
  ) %>% data.table::setDT()

  prop_converged <- purrr::imap_dfr(convergence_split, \(bool, x) {
    total_converged <- sum(bool) / length(bool)
    data.table(no_components = as.integer(x), converged = total_converged)
  })

  ica_stability_res_sum = ica_stability_res[, .(
    median_stability = median(stability),
    max_stability = max(stability),
    percentile_75 = quantile(stability, 0.75)
  ), .(no_components)] %>%
    merge(., prop_converged, by = 'no_components')

  stability_params <- list(
    ica_type = ica_type,
    tested_ncomps = n_comp_vector,
    iter_params = iter_params,
    ica_params = ica_params,
    random_seed = random_seed
  )

  # Assign stuff
  S7::prop(object, "outputs")[['ica_stability_res']] <- ica_stability_res
  S7::prop(object, "outputs")[['ica_stability_res_sum']] <- ica_stability_res_sum
  S7::prop(object, "params")[["ica_stability_assessment"]] <- stability_params

  return(object)
}


#' @title Run stabilised ICA with a given number of components
#'
#' @description
#' This function runs stabilised ICA with the defined number of components.
#'
#' @param object The class, see [bixverse::bulk_coexp()]. You need to apply
#' [bixverse::ica_processing()] before running this function.
#' @param no_comp Integer. Number of components you wish to use for the ICA
#' run. [bixverse::ica_evaluate_comp()] can give you an idea in terms of
#' stability of the number of components at different levels.
#' @param ica_type String, element of `c("logcosh", "exp")`.
#' @param iter_params List. This list controls the randomisation parameters for
#' the ICA runs, see [bixverse::params_ica_randomisation()] for estimating
#' stability. Has the following elements:
#' \itemize{
#'  \item cross_validate - Boolean. Shall the data be split into different
#'  chunks on which ICA is run. This will slow down the function substantially,
#'  as every chunk needs to whitened again.
#'  \item random_init - Integer. How many random initialisations shall be used
#'  for the ICA runs.
#'  \item folds - If `cross_validate` is set to `TRUE` how many chunks shall be
#'  used. To note, you will run per ncomp random_init * fold ICA runs which
#'  can quickly increase.
#' }
#' @param ica_params List. The ICA parameters, see [bixverse::params_ica_general()]
#' wrapper function. This function generates a list containing:
#' \itemize{
#'  \item maxit - Integer. Maximum number of iterations for ICA.
#'  \item alpha - Float. The alpha parameter for the logcosh version of ICA.
#'  Should be between 1 to 2.
#'  \item max_tol - Maximum tolerance of the algorithm.
#'  \item verbose - Controls verbosity of the function.
#' }
#' @param random_seed Integer. For reproducibility.
#' @param consistent_sign Boolean. If set to `TRUE`, for each source the absolute
#' maximum value will be positive, i.e., the sign will be inverted so that the
#' absolute bigger tail is set to positive floats.
#' @param .verbose Boolean. Controls verbosity.
#'
#' @return `bulk_coexp` with the the source matrix S, mixing matrix A and other
#' parameters added to the slots.
#'
#' @export
ica_stabilised_results <- S7::new_generic(
  name = "ica_stabilised_results",
  dispatch_args = "object",
  fun = function(object,
                 no_comp,
                 ica_type = c("logcosh", "exp"),
                 iter_params = params_ica_randomisation(),
                 ica_params = params_ica_general(),
                 random_seed = 42L,
                 consistent_sign = TRUE,
                 .verbose = TRUE) {
    S7::S7_dispatch()
  }
)


#' @export
#'
#' @import data.table
#' @importFrom magrittr `%>%`
#' @importFrom zeallot `%<-%`
#'
#' @method ica_stabilised_results bulk_coexp
S7::method(ica_stabilised_results, bulk_coexp) <- function(object,
                                                           no_comp,
                                                           ica_type = c("logcosh", "exp"),
                                                           iter_params = params_ica_randomisation(),
                                                           ica_params = params_ica_general(),
                                                           random_seed = 42L,
                                                           consistent_sign = TRUE,
                                                           .verbose = TRUE) {
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  checkmate::qassert(no_comp, "I1")
  checkmate::assertChoice(ica_type, c("logcosh", "exp"))
  assertIcaIterParams(iter_params)
  assertIcaParams(ica_params)
  checkmate::qassert(random_seed, "I1")
  checkmate::qassert(consistent_sign, "B1")
  checkmate::qassert(.verbose, "B1")

  detection_method <- S7::prop(object, "params")[["detection_method"]]

  # Early return
  if (is.null(detection_method) &&
      detection_method != "ICA-based") {
    warning(
      paste(
        "This class does not seem to be set for ICA-based module detection",
        "Returning class as is."
      )
    )
    return(object)
  }

  # Get the attributes
  X <- S7::prop(object, "processed_data")[['processed_data']]
  X1 <- S7::prop(object, "processed_data")[["X1"]]
  K <- S7::prop(object, "processed_data")[["K"]]

  # Get the combined S matrix and convergence information
  c(s_combined, converged) %<-% with(iter_params, switch(
    as.integer(iter_params$cross_validate) + 1,
    rs_ica_iters(
      x1 = X1,
      k = K,
      no_comp = no_comp,
      no_random_init = random_init,
      ica_type = ica_type,
      random_seed = random_seed,
      ica_params = ica_params
    ),
    rs_ica_iters_cv(
      x = X_raw,
      no_comp = no_comp,
      no_folds = folds,
      no_random_init = random_init,
      ica_type = ica_type,
      random_seed = random_seed,
      ica_params = ica_params
    )
  ))

  # Get the component stability and centrotypes
  c(stability_scores, centrotype) %<-% .community_stability(
    no_comp = as.integer(no_comp),
    s = s_combined,
    return_centrotype = TRUE
  )

  colnames(centrotype) <- sprintf("IC_%i", 1:no_comp)
  rownames(centrotype) <- colnames(X)

  if (consistent_sign) {
    centrotype <- apply(centrotype, 2, .flip_ica_loading_signs)
  }

  S <- t(centrotype)
  A <- t(X1) %*% MASS::ginv(S)
  rownames(A) <- rownames(X)
  colnames(A) <- rownames(S)

  ica_meta <- list(component = sprintf("IC_%i", 1:no_comp),
                   stability = stability_scores) %>% data.table::setDT()

  result <- list(S = S, A = A, ica_meta = ica_meta)

  result_params <- list(
    no_comp = no_comp,
    ica_type = ica_type,
    iter_params = iter_params,
    ica_params = ica_params,
    random_seed = random_seed,
    consistent_sign = consistent_sign,
    converged = converged
  )

  S7::prop(object, "final_results") <- result
  S7::prop(object, "params")[["ica_final_gen"]] <- result_params

  return(object)
}


# general ICA function --------------------------------------------------------

#' Fast ICA via Rust
#'
#' @description
#' This functions is a wrapper over the Rust implementation of fastICA. It has
#' the same two options `c("logcosh", "exp")` to run ICA in parallel modus. This
#' function expected
#'
#' @param X Numeric matrix. Processed data. Output of
#' [bixverse::rs_prepare_whitening()].
#' @param K The K matrix. Whitening matrix. Output of
#' [bixverse::rs_prepare_whitening()].
#' @param n_icas Integer. Number of independent components to recover.
#' @param ica_fun String, element of `c("logcosh", "exp")`.
#' @param seed Integer. Seed to ensure reproducible results.
#' @param ica_params List. The ICA parameters, see [bixverse::params_ica_general()]
#' wrapper function. This function generates a list containing:
#' \itemize{
#'  \item maxit - Integer. Maximum number of iterations for ICA.
#'  \item alpha - Float. The alpha parameter for the logcosh version of ICA.
#'  Should be between 1 to 2.
#'  \item max_tol - Maximum tolerance of the algorithm.
#'  \item verbose - Controls verbosity of the function.
#' }
#'
#' @returns A list containing:
#' \itemize{
#'  \item w The mixing matrix w.
#'  \item A ICA results matrix A.
#'  \item S ICA results matrix S.
#'  \item converged Boolean indicating if algorithm converged.
#' }
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
fast_ica_rust <- function(X,
                          K,
                          n_icas,
                          ica_fun = c("logcosh", "exp"),
                          seed = NULL,
                          ica_params = params_ica_general()) {
  # Checks
  checkmate::assertMatrix(X, mode = "numeric")
  checkmate::assertMatrix(K, mode = "numeric")
  checkmate::qassert(n_icas, sprintf("I1[2,%i]", ncol(X)))
  checkmate::assertChoice(ica_fun, c("logcosh", "exp"))
  checkmate::qassert(seed, c("I1", "0"))
  assertIcaParams(ica_params)
  # Function
  K <- matrix(K[1:n_icas, ], n_icas, dim(X)[1])
  X1 <- K %*% X
  set.seed(seed)
  w_init <- matrix(rnorm(n_icas * n_icas), nrow = n_icas, ncol = n_icas)

  c(a, converged) %<-% rs_fast_ica(X1, w_init, ica_type = ica_fun, ica_params = ica_params)

  w <- a %*% K
  S <- w %*% X
  A <- t(w) %*% solve(w %*% t(w))

  res <- list(
    w = w,
    A = A,
    S = S,
    converged = converged
  )

  return(res)
}

# helpers ----------------------------------------------------------------------

## plotting --------------------------------------------------------------------

#' @title Plot the stability of the ICA components
#'
#' @description
#' Helper function to plot the individual stability profiles over the tested
#' ncomps.
#'
#' @param object The class, see [bixverse::bulk_coexp()]. You need to apply
#' [bixverse::ica_evaluate_comp()] before running this function.
#'
#' @export
plot_ica_stability_individual <- S7::new_generic(
  name = "plot_ica_stability_individual",
  dispatch_args = "object",
  fun = function(object) {
    S7::S7_dispatch()
  }
)

#' @method plot_ica_stability_individual bulk_coexp
S7::method(plot_ica_stability_individual, bulk_coexp) <- function(object) {
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  # Function body
  plot_df <- S7::prop(object, "outputs")[['ica_stability_res']]
  if (is.null(plot_df)) {
    warning("No ICA stability results found. Did you run ica_evaluate_comp()? Returning NULL.")
    return(NULL)
  }

  p1 <- ggplot(data = plot_df,
               mapping = aes(x = component_rank, y = stability)) +
    geom_line(mapping = aes(color = factor(no_components)), linewidth = 1) +
    scale_color_viridis_d(option = "C") +
    theme_minimal() +
    labs(color = "No ICAs") +
    ylim(0, 1) +
    xlab("Component rank") +
    ylab("Stability index")


  p2 <- ggplot(data = plot_df, aes(x = component_rank, y = stability)) +
    geom_point(mapping = aes(colour = factor(no_components))) +
    scale_color_viridis_d(option = "C") +
    theme_minimal() +
    labs(color = "No ICAs") +
    ylim(0, 1) +
    xlab("Component rank") +
    ylab("Stability index")

  p1 + p2 + patchwork::plot_annotation(title = "Stability of independent components", subtitle = "Over different ncomps and randomisations")
}

#' @title Plot the maximum and 75 %ile of the stability.
#'
#' @description
#' Helper function to plot the individual maximum and 75th percentile in terms
#' of stability over the tested ncomp to identify inflection points in the
#' stability profiles.
#'
#' @param object The class, see [bixverse::bulk_coexp()]. You need to apply
#' [bixverse::ica_evaluate_comp()] before running this function.
#'
#' @export
plot_ica_stability_summarised <- S7::new_generic(
  name = "plot_ica_stability_summarised",
  dispatch_args = "object",
  fun = function(object) {
    S7::S7_dispatch()
  }
)

#' @method plot_ica_stability_summarised bulk_coexp
S7::method(plot_ica_stability_summarised, bulk_coexp) <- function(object) {
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  # Function body
  plot_df <- S7::prop(object, "outputs")[['ica_stability_res_sum']]
  if (is.null(plot_df)) {
    warning("No ICA stability results found. Did you run ica_evaluate_comp()? Returning NULL.")
    return(NULL)
  }

  p <- ggplot(data = plot_df,
              mapping = aes(x = no_components, y = max_stability)) +
    geom_point(
      mapping = aes(fill = percentile_75),
      shape = 21,
      size = 3,
      alpha = .75
    ) +
    ylim(0, 1) +
    xlab("Component rank") +
    ylab("Max stability index") +
    scale_fill_viridis_c(option = "G") +
    labs(fill = "Stability 75%-ile:") +
    theme_minimal() +
    ggtitle(label = "Maximum and 75th %-ile component stability") +
    stat_smooth(
      method = "loess",
      se = FALSE,
      span = 0.5,
      color = "darkred",
      linetype = "dashed",
      linewidth = 0.5
    )

  p
}

## component stability ---------------------------------------------------------

#' Assess stability of ICA components
#'
#' @description
#' Assesses the stability of a component given some random iterations (in terms
#' of random initialisations and potentially additional cross-validation). It
#' will calculate the (absolute) Pearson's correlation between the different
#' generated independent components and run hierarchical clustering on top of
#' these and assess the cluster stability. If desired, it can also return the
#' centrotype of each cluster.
#'
#' @param no_comp Integer. Number of components that were set for this run.
#' @param s Numeric matrix. Combined S matrices from the ICA runs. Assumes
#' rows = features and column individual independent components.
#' @param return_centrotype Boolean. Shall the centrotype be returned for this
#' run.
#'
#' @returns A list containing:
#' \itemize{
#'  \item stability_scores - The stability scores for the given runs.
#'  \item centrotype - The centrotype component of each cluster.
#' }
.community_stability <- function(no_comp, s, return_centrotype) {
  # Visible global function stuff...
  as.dist <- hclust <- cutree <- NULL
  # Checks
  checkmate::qassert(no_comp, "I1")
  checkmate::assertMatrix(s, mode = "numeric")
  checkmate::qassert(return_centrotype, "B1")
  # Body
  abs_cor <- abs(rs_cor(s, spearman = FALSE))
  dist <- as.dist(1 - abs_cor)

  clusters <- hclust(dist)

  clusterCut <- cutree(tree = clusters, k = no_comp)

  scores <- vector(mode = 'double', length = no_comp)
  if (return_centrotype) {
    centrotypes <- vector(mode = "list", length = no_comp)
  }

  for (cluster in seq_len(no_comp)) {
    cluster_indx <- which(clusterCut == cluster)
    not_cluster_indx <- which(clusterCut != cluster)
    within_cluster <- sum(abs_cor[cluster_indx, cluster_indx]) / length(cluster_indx)^2
    outside_cluster <- sum(abs_cor[cluster_indx, not_cluster_indx]) / (length(cluster_indx) * length(not_cluster_indx))
    scores[cluster] <- within_cluster - outside_cluster
    if (return_centrotype) {
      max_similarity <- max(rowSums(abs_cor[cluster_indx, cluster_indx]))
      temp <- which(rowSums(abs_cor[cluster_indx, cluster_indx]) == max_similarity)
      centrotypes[[cluster]] <- s[, cluster_indx[temp], drop = FALSE]
    }
  }

  res <- list(stability_scores = scores, centrotype = if (return_centrotype) {
    do.call(cbind, centrotypes)
  } else {
    NULL
  })

  res
}

## sign flipping ---------------------------------------------------------------


#' Flips the ICA source sign
#'
#' @description
#' Makes sure that the largest absolute value in the ICA source matrix S will
#' be positive.
#'
#' @param x Numeric vector. The source signal for that independent component.
#'
#' @returns Returns a consistent ICA feature loading
.flip_ica_loading_signs = function(x) {
  feature_sign <- sign(x)
  max_val_sign <- feature_sign[which(abs(x) == max(abs(x)))]
  y <- if (max_val_sign == 1) {
    x
  } else {
    -x
  }
  y
}
