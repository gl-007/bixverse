# methods - simple correlations ------------------------------------------------

#' @title Prepare correlation-based module detection
#'
#' @description
#' This function will calculate the correlation coefficients between the genes,
#' using the highly variable genes (if available, otherwise the function will
#' use the raw data). The data will be stored in a memory-efficient format
#' in the properties of the class.
#'
#' @param object The class, see [bixverse::bulk_coexp()]. Ideally, you
#' should run [bixverse::preprocess_bulk_coexp()] before applying this function.
#' @param correlation_method String. Option of `c("pearson", "spearman")`.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return The class with added data to the properties for subsequent usage.
#'
#' @export
cor_module_processing <- S7::new_generic(
  name = "cor_module_processing",
  dispatch_args = "object",
  fun = function(object,
                 correlation_method = c("pearson", "spearman"),
                 .verbose = TRUE) {
    S7::S7_dispatch()
  }
)


#' @export
#'
#' @method cor_module_processing bulk_coexp
S7::method(cor_module_processing, bulk_coexp) <- function(object,
                                                          correlation_method = c("pearson", "spearman"),
                                                          .verbose = TRUE) {
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  checkmate::assertChoice(correlation_method, c("pearson", "spearman"))
  checkmate::qassert(.verbose, "B1")

  # Function body
  if (purrr::is_empty(S7::prop(object, "processed_data")[['processed_data']])) {
    warning("No pre-processed data found. Defaulting to the raw data")
    target_mat <- S7::prop(object, "raw_data")
  } else {
    target_mat <- S7::prop(object, "processed_data")[['processed_data']]
  }

  spearman <- if (correlation_method == 'pearson') {
    if (.verbose)
      message("Using Pearson correlations.")
    FALSE
  } else {
    if (.verbose)
      message("Using Spearman correlations.")
    TRUE
  }

  # Calculate the upper triangle of correlation matrix
  cor_diagonal <- rs_cor_upper_triangle(target_mat, spearman = spearman, shift = 1L)

  # Save data to memory friendly R6 class
  cor_data <- upper_triangular_cor_mat$new(
    cor_coef = cor_diagonal,
    features = colnames(target_mat),
    shift = 1L
  )

  correlation_params <- list(spearman = spearman, type = 'simple')

  S7::prop(object, "processed_data")[["correlation_res"]] <- cor_data
  S7::prop(object, "params")[["correlation_params"]] <- correlation_params
  S7::prop(object, "params")["detection_method"] <- "correlation-based"

  return(object)
}

# methods - differential correlations ------------------------------------------

#' @title Prepare differential correlation-based module detection
#'
#' @description
#' This function will calculate the differential correlation between the stored
#' data set in the class and another background data set. To do so, it uses a
#' Fisher transformation of the correlation coefficients and calculates a Z
#' score based on the delta. The function will automatically subset into shared
#' features between the two data sets.
#'
#' @param object The class, see [bixverse::bulk_coexp()]. Ideally, you
#' should run [bixverse::preprocess_bulk_coexp()] before applying this function.
#' @param background_mat Numerical matrix. The background data set.
#' @param correlation_method String. Option of `c("pearson", "spearman")`.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return The class with added data to the properties for subsequent usage.
#'
#' @export
diffcor_module_processing <- S7::new_generic(
  name = "diffcor_module_processing",
  dispatch_args = "object",
  fun = function(object,
                 background_mat,
                 correlation_method = c("pearson", "spearman"),
                 .verbose = TRUE) {
    S7::S7_dispatch()
  }
)

#' @export
#'
#' @method diffcor_module_processing bulk_coexp
S7::method(diffcor_module_processing, bulk_coexp) <- function(object,
                                                              background_mat,
                                                              correlation_method = c("pearson", "spearman"),
                                                              .verbose = TRUE) {
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  checkmate::assertMatrix(background_mat, mode = 'numeric')
  checkmate::assertChoice(correlation_method, c("pearson", "spearman"))
  checkmate::qassert(.verbose, "B1")

  # Function
  if (purrr::is_empty(S7::prop(object, "processed_data")[['processed_data']])) {
    warning("No pre-processed data found. Defaulting to the raw data")
    target_mat <- S7::prop(object, "raw_data")
  } else {
    target_mat <- S7::prop(object, "processed_data")[['processed_data']]
  }

  spearman <- if (correlation_method == 'pearson') {
    if (.verbose)
      message("Using Pearson correlations.")
    FALSE
  } else {
    if (.verbose)
      message("Using Spearman correlations.")
    TRUE
  }

  features <- colnames(target_mat)
  shared_features <- intersect(colnames(target_mat), colnames(background_mat))

  # Early return if there are no shared features
  if (length(shared_features) == 0) {
    warning("No shared features identified. Returning class as is.")
    return(object)
  }
  if (.verbose)
    message(
      sprintf(
        "A total of %i shared features were identified and used for differential correlation.",
        length(shared_features)
      )

    )

  target_mat <- target_mat[, shared_features]
  background_mat <- background_mat[, shared_features]

  combined_mad_df <- list(
    feature = shared_features,
    MAD_target = matrixStats::colMads(target_mat),
    MAD_background = matrixStats::colMads(background_mat)
  ) %>% data.table::setDT()

  diff_cor_res <- rs_differential_cor(target_mat, background_mat, spearman = spearman)

  cor_data <- upper_triangle_diffcor_mat$new(diff_cor_res = diff_cor_res, features = shared_features)

  correlation_params <- list(
    spearman = spearman,
    type = 'differential correlation',
    no_intersecting_features = length(shared_features)
  )

  S7::prop(object, "processed_data")[["differential_cor_res"]] <- cor_data
  S7::prop(object, "processed_data")[["differential_cor_feature_meta"]] <- combined_mad_df
  S7::prop(object, "params")[["correlation_params"]] <- correlation_params
  S7::prop(object, "params")["detection_method"] <- "differential correlation-based"

  object
}


# methods - graph-based gene module detection ----------------------------------

#' @title Iterate through different epsilon parameters
#'
#' @description
#' This functions iterates through a set of provided epsilons and checks for
#' each one to what extend the underlying affinity matrix will follow a power
#' law distribution.
#'
#' @param object The class, see [bixverse::bulk_coexp()]. You need to run
#' [bixverse::diffcor_module_processing()] before running this function.
#' @param epsilons Vector of floats. The different epsilon parameters you
#' would like to run.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return The class with added data to the properties for subsequent usage.
#'
#' @export
cor_module_check_epsilon <- S7::new_generic(
  name = "cor_module_check_epsilon",
  dispatch_args = "object",
  fun = function(object,
                 epsilons = c(0.25, seq(from = 0.5, to = 5, by = 0.5)),
                 .verbose = TRUE) {
    S7::S7_dispatch()
  }
)

#' @export
#'
#' @importFrom magrittr `%>%`
#' @importFrom zeallot `%->%`
#' @import data.table
#'
#' @method cor_module_check_epsilon bulk_coexp
S7::method(cor_module_check_epsilon, bulk_coexp) <- function(object,
                                                             epsilons = c(0.25, seq(from = 0.5, to = 5, by = 0.5)),
                                                             .verbose = TRUE) {
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  checkmate::qassert(epsilons, "R+")
  checkmate::qassert(.verbose, "B1")

  detection_method <- S7::prop(object, "params")[["detection_method"]]

  # Early return
  if (is.null(detection_method) ||
      detection_method != "correlation-based") {
    warning(
      paste(
        "This class does not seem to be set for correlation-based module detection",
        "Returning class as is."
      )
    )
    return(object)
  }

  # Pull out the correlation results
  cor_res <- S7::prop(object, "processed_data")$correlation_res
  c(cor_vector, n_features, shift) %<-% cor_res$get_cor_vector()

  # Prepare everything for iterating through the epsilons
  epsilons <- sort(epsilons, decreasing = TRUE)
  dist_vec <- 1 - abs(cor_vector)
  dist_vec <- data.table::fifelse(dist_vec < 0, 0, dist_vec)

  if (.verbose)
    message(sprintf("Testing %i epsilons.", length(epsilons)))

  epsilon_data <- rs_rbf_iterate_epsilons(
    dist = dist_vec,
    epsilon_vec = epsilons,
    original_dim = n_features,
    shift = shift,
    rbf_type = "bump"
  )

  r_square_vals <- apply(epsilon_data, 2, .scale_free_fit)

  r_square_data <- list(epsilon = epsilons, r2_vals = r_square_vals) %>%
    data.table::setDT()

  S7::prop(object, "outputs")[['epsilon_data']] <- r_square_data

  return(object)

}


#' @title Iterate through Leiden resolutions for graph-based community detection.
#'
#' @description
#' This function will identify gene modules based on affinity graphs from the
#' single correlation or differential correlation methods. Briefly, in the case
#' of single correlation, the graph is generated based on the absolute
#' correlation coefficients that are subjected to a Gaussian affinity kernel.
#' This reduces spurious correlations and leaves a sparsely connected graph.
#' In the case of differential correlations, the graph is generated based on
#' significant differential correlations if one of the two correlations reached
#' the defined minimum thresholds.\cr
#' Subsequently, Leiden community detection is applied on the respective graph
#' through a range of resolutions that the user can define. The function then
#' returns meta information about the resolutions (which can also be plotted) to
#' identify the best suitable resolution parameter to identify co-expression modules.
#'
#' @param object The class, see [bixverse::bulk_coexp()].
#' @param resolution_params List. Parameters for the resolution search, see
#' [bixverse::params_graph_resolution()]. Contains:
#' \itemize{
#'  \item min_res - Float. Minimum resolution to test.
#'  \item max_res - Float. Maximum resolution to test.
#'  \item number_res - Integer. Number of resolutions to test between the
#'  `max_res` and `min_res.`
#' }
#' @param graph_params List. Parameters for the generation of the (differential)
#' correlation graph, see [bixverse::params_cor_graph()]. Contains:
#' \itemize{
#'  \item Epsilon - Defines the epsilon parameter for the radial basis
#'  function. Defaults to 1, but should be ideally optimised.
#'  \item min_cor - Float. Minimum absolute correlation that needs to be
#'  observed in either data set. Only relevant for differential correlation-based
#'  graphs.
#'  \item fdr_threshold - Float. Maximum FDR for the differential correlation
#'  p-value.
#'  \item verbose - Boolean. Controls verbosity of the graph generation.
#' }
#' @param random_seed Integer. Random seed.
#' @param min_genes Integer. Minimum number of genes that should be in a
#' community.
#' @param parallel Boolean. Parallelise the Leiden clustering.
#' @param max_workers Integer. Maximum number of workers to use if parallel is
#' set to `TRUE`.
#' @param .verbose Controls the verbosity of the function.
#'
#' @return The class with added data to the properties.
#'
#' @export
cor_module_check_res <- S7::new_generic(
  name = "cor_module_check_res",
  dispatch_args = "object",
  fun = function(object,
                 resolution_params = params_graph_resolution(),
                 graph_params = params_cor_graph(),
                 random_seed = 123L,
                 min_genes = 10L,
                 parallel = TRUE,
                 max_workers = as.integer(parallel::detectCores() / 2),
                 .verbose = TRUE) {
    S7::S7_dispatch()
  }
)


#' @export
#'
#' @importFrom magrittr `%>%`
#' @importFrom zeallot `%->%`
#' @importFrom future plan multisession sequential
#' @import data.table
#'
#' @method cor_module_check_res bulk_coexp
S7::method(cor_module_check_res, bulk_coexp) <- function(object,
                                                         resolution_params = params_graph_resolution(),
                                                         graph_params = params_cor_graph(),
                                                         random_seed = 123L,
                                                         min_genes = 10L,
                                                         parallel = TRUE,
                                                         max_workers = as.integer(parallel::detectCores() / 2),
                                                         .verbose = TRUE) {

  # Checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  assertCorGraphParams(graph_params)
  assertGraphResParams(resolution_params)
  checkmate::qassert(min_genes, "I1")
  checkmate::qassert(parallel, "B1")
  checkmate::qassert(max_workers, "I1")
  checkmate::qassert(.verbose, "B1")

  detection_method <- S7::prop(object, "params")[["detection_method"]]

  # Early return
  if (is.null(detection_method) ||
      !detection_method %in% c("correlation-based", "differential correlation-based")) {
    warning(
      paste(
        "This class does not seem to be set for correlation-based module detection",
        "Returning class as is."
      )
    )
    return(object)
  }

  c(graph, graph_params) %<-% with(graph_params, switch(
    detection_method,
    "correlation-based" = get_cor_graph(
      object = object,
      epsilon = epsilon,
      .verbose = verbose
    ),
    "differential correlation-based" = get_diffcor_graph(
      object = object,
      min_cor = min_cor,
      fdr_threshold = fdr_threshold,
      .verbose = verbose
    )
  ))

  resolutions <- with(resolution_params, exp(seq(log(min_res), log(max_res), length.out = number_res)))

  if (.verbose)
    message(sprintf("Iterating through %i resolutions", length(resolutions)))

  if (parallel) {
    if (.verbose)
      message(sprintf("Using parallel computation over %i cores.", max_workers))

    future::plan(future::multisession(workers = max_workers))
  } else {
    if (.verbose)
      message("Using sequential computation.")
    future::plan(future::sequential())
  }

  community_df_res <- furrr::future_map(
    resolutions,
    \(res) {
      set.seed(random_seed)
      community <- igraph::cluster_leiden(
        graph,
        objective_function = 'modularity',
        resolution = res,
        n_iterations = 5L
      )

      modularity <- igraph::modularity(x = graph, membership = community$membership)

      community_df <- data.table::data.table(
        resolution = res,
        node_name = community$names,
        membership = community$membership,
        modularity = modularity
      )
    },
    .progress = .verbose,
    .options = furrr::furrr_options(seed = TRUE)
  ) %>% data.table::rbindlist(.)

  # To make the message trace prettier, if set to verbose
  if (.verbose)
    message("")

  future::plan(future::sequential())
  gc()

  community_df_res[, combined_id := sprintf("id_%s_%s", resolution, membership)]

  cluster_summary <- community_df_res[, .N, combined_id] %>%
    .[, good_clusters := N >= min_genes] %>%
    data.table::merge.data.table(.,
                                 unique(community_df_res[, c('resolution', 'combined_id')]),
                                 by.x = 'combined_id',
                                 by.y = 'combined_id',
    ) %>%
    .[, .(
      good_clusters = sum(good_clusters),
      avg_size = mean(N),
      max_size = max(N)
    ), resolution]

  resolution_results <- community_df_res[, .(no_clusters = length(unique(membership)),
                                             modularity = unique(modularity)), resolution] %>%
    data.table::merge.data.table(., cluster_summary, by.x = 'resolution', by.y = 'resolution')



  # Assign stuff
  S7::prop(object, "params")[["correlation_graph"]] <- graph_params
  S7::prop(object, "outputs")[['resolution_results']] <- resolution_results
  S7::prop(object, "outputs")[['cor_graph']] <- graph

  return(object)
}

#' @title Identify correlation-based gene modules via graphs.
#'
#' @description
#' This function leverages graph-based clustering to identify gene co-expression
#' modules. The class has the option to sub-cluster large communities within
#' their respective sub graphs, akin to the approach taken by Barrio-Hernandez,
#' et al.
#'
#' @param object The class, see [bixverse::bulk_coexp()].
#' @param resolution The Leiden resolution parameter you wish to use. If NULL,
#' it will use the optimal one identified by [bixverse::cor_module_check_res()].
#' If nothing can be found, will default to 1.
#' @param min_size Integer. Minimum size of the communities.
#' @param max_size Integer. Maximum size of the communities.
#' @param subclustering Boolean. Shall after a first clustering communities that
#' are too large be further sub clustered. Defaults to `TRUE`.
#' @param random_seed Integer. Random seed.
#' @param .graph_params List. Parameters for the generation of the (differential)
#' correlation graph, see [bixverse::params_cor_graph()]. Contains:
#' \itemize{
#'  \item Epsilon - Defines the epsilon parameter for the radial basis
#'  function. Defaults to 2, but should be ideally optimised.
#'  \item min_cor - Float. Minimum absolute correlation that needs to be
#'  observed in either data set. Only relevant for differential correlation-based
#'  graphs.
#'  \item fdr_threshold - Float. Maximum FDR for the differential correlation
#'  p-value.
#'  \item verbose - Boolean. Controls verbosity of the graph generation.
#' }
#' This parameter is only relevant if you did *not* run
#' [bixverse::cor_module_check_res()].
#' @param .max_iters Integer. If sub clustering is set to `TRUE`, what shall be the
#' maximum number of iterations. Defaults to 100L.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @return The class with added data to the properties.
#'
#' @references Barrio-Hernandez, et al., Nat Genet, 2023.
#'
#' @export
cor_module_final_modules <- S7::new_generic(
  name = "cor_module_final_modules",
  dispatch_args = "object",
  fun = function(object,
                 resolution = NULL,
                 min_size = 10L,
                 max_size = 500L,
                 subclustering = TRUE,
                 random_seed = 123L,
                 .graph_params = params_cor_graph(),
                 .max_iters = 100L,
                 .verbose = TRUE) {
    S7::S7_dispatch()
  }
)

#' @export
#'
#' @importFrom magrittr `%>%`
#' @importFrom zeallot `%->%`
#' @import data.table
#'
#' @method cor_module_final_modules bulk_coexp
S7::method(cor_module_final_modules, bulk_coexp) <- function(object,
                                                             resolution = NULL,
                                                             min_size = 10L,
                                                             max_size = 500L,
                                                             subclustering = TRUE,
                                                             random_seed = 123L,
                                                             .graph_params = params_cor_graph(),
                                                             .max_iters = 100L,
                                                             .verbose = TRUE) {
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  checkmate::qassert(resolution, c("0", "N1"))
  checkmate::qassert(min_size, "I1")
  checkmate::qassert(max_size, "I1")
  checkmate::qassert(subclustering, "B1")
  checkmate::qassert(random_seed, "I1")
  assertCorGraphParams(.graph_params)
  checkmate::qassert(.max_iters, "I1")
  checkmate::qassert(.verbose, "B1")

  detection_method <- S7::prop(object, "params")[["detection_method"]]

  # Early return
  if (is.null(detection_method) &&
      detection_method %in% c("correlation-based", "differential correlation-based")) {
    warning(
      paste(
        "This class does not seem to be set for correlation-based module detection",
        "Returning class as is."
      )
    )
    return(object)
  }

  # Get the graph
  if (is.null(S7::prop(object, "outputs")[["cor_graph"]])) {
    # Deal with the case a graph was not yet generated...
    warning(
      paste(
        "No correlation graph found. Did you run cor_module_check_res()?",
        "Generating correlation graph based on standard parameters.",
        sep = "\n"
      )
    )

    c(graph, graph_params) %<-% with(.graph_params, switch(
      detection_method,
      "correlation-based" = get_cor_graph(
        object = object,
        kernel_bandwidth = kernel_bandwidth,
        min_affinity = min_affinity,
        .verbose = verbose
      ),
      "differential correlation-based" = get_diffcor_graph(
        object = object,
        min_cor = min_cor,
        fdr_threshold = fdr_threshold,
        .verbose = verbose
      )
    ))

    S7::prop(object, "params")[["correlation_graph"]] <- graph_params
    S7::prop(object, "outputs")[['cor_graph']] <- graph

  } else {
    graph <- S7::prop(object, "outputs")[["cor_graph"]]
  }

  # Final resolution
  if (is.null(resolution)) {
    resolution_results <- S7::prop(object, "outputs")[["resolution_results"]]
    final_resolution <- if (!is.null(resolution_results)) {
      if (.verbose)
        message("Using resolution with best modularity.")
      resolution_results[modularity == max(modularity), resolution]
    } else {
      warning("No resolution results found and none provided. Will default to a resolution of 1.")
      1
    }
  }

  # Do a first clustering
  set.seed(random_seed)

  final_gene_communities <- igraph::cluster_leiden(
    graph,
    objective_function = 'modularity',
    resolution = final_resolution,
    n_iterations = 5L
  )

  clusters_df <- data.table::data.table(node_id = final_gene_communities$names,
                                        cluster_id = final_gene_communities$membership)

  node_frequency <- clusters_df[, .N, .(cluster_id)]

  # If sub clustering is active, do that
  if (subclustering) {
    if (.verbose)
      message("Sub-clustering larger communities until they are below max_size.")
    clusters_with_too_many_nodes <- node_frequency[N > max_size, cluster_id]
    final_clusters <- clusters_df[!cluster_id %in% clusters_with_too_many_nodes]

    for (i in seq_along(clusters_with_too_many_nodes)) {
      cluster_i <- clusters_with_too_many_nodes[i]
      nodes_in_cluster <- clusters_df[cluster_id == cluster_i, node_id]
      finalised_clusters <- data.table()
      # Loop through, until all clusters are below the minimum genes or max
      # iterations is hit
      l <- 1
      while (length(nodes_in_cluster) != 0) {
        set.seed(random_seed + l)

        sub_graph_l <- igraph::subgraph(graph,
                                        data.table::chmatch(nodes_in_cluster, igraph::V(graph)$name))

        # Restarting at a very small resolution
        clusters_red <- igraph::cluster_leiden(sub_graph_l, resolution = 0.1 + l * 0.05)

        subclusters <- data.table(node_id = clusters_red$names,
                                  cluster_id = clusters_red$membership)

        subclusters_frequency <- subclusters[, .N, .(cluster_id)]
        clusters_small_enough <- subclusters_frequency[N <= max_size, cluster_id]

        good_clusters <- subclusters[cluster_id %in% clusters_small_enough] %>%
          dplyr::mutate(cluster_id = paste0(i, paste(rep("sub", l), collapse = ""), cluster_id))

        finalised_clusters <- rbind(finalised_clusters, good_clusters)

        l <- l + 1
        if (l == .max_iters) {
          break
        }

        nodes_in_cluster <- setdiff(nodes_in_cluster, good_clusters$node_id)
      }

      final_clusters <- rbind(final_clusters, finalised_clusters)
    }

    node_frequency_updated <- final_clusters[, .N, .(cluster_id)]
  } else {
    node_frequency_updated <- node_frequency
    final_clusters <- clusters_df
  }

  # Finalise the data
  final_communities <- node_frequency_updated[N <= max_size &
                                                N >= min_size, cluster_id]
  final_clusters_filtered <- final_clusters[cluster_id %in% final_communities]

  cluster_name_prettifier <- setNames(
    paste("cluster", seq_along(
      unique(final_clusters_filtered$cluster_id)
    ), sep = "_"),
    unique(final_clusters_filtered$cluster_id)
  )

  final_clusters_filtered[, cluster_id := cluster_name_prettifier[cluster_id]]

  results_param <- list(
    resolution = resolution,
    seed = random_seed,
    min_size = min_size,
    max_size = max_size,
    max_iters = .max_iters
  )

  S7::prop(object, "final_results") <- final_clusters_filtered
  S7::prop(object, "params")[["module_final_gen"]] <- results_param

  return(object)
}


# methods - helpers ------------------------------------------------------------

## power law calculations ------------------------------------------------------

#' Calculate the goodness of fit for a power law distribution.
#'
#' @param k Numeric vector. The vector of
#' @param breaks Integer. Number of breaks for fitting the data.
#' @param plot Boolean. Shall the log-log plot be generated.
#'
#' @returns The R2 value of of the goodness of fit.
.scale_free_fit <- function(k, breaks = 50L, plot = FALSE) {
  # Visible global function stuff...
  lm <- NULL
  # Checks
  checkmate::qassert(k, "R>=50")
  checkmate::qassert(breaks, "I1")
  checkmate::qassert(plot, "B1")
  # Deal with stupid case that someone supplies something small here
  if (breaks > length(k)) {
    breaks <- ceiling(k / 10L)
  }
  k_discrete <- cut(k, breaks)
  dk <- tapply(k, k_discrete, mean)
  dk_p <- tapply(k, k_discrete, length) / length(k)
  log_dk <- log10(dk)
  log_dk_p <- log10(dk_p)

  if (plot) {
    plot(x = log_dk, y = log_dk_p)
  }

  summary(lm(log_dk ~ log_dk_p))$r.squared
}

## graph generation ------------------------------------------------------------

#' @title Get correlation-based graph
#'
#' @description
#' Helper function to get a correlation-based igraph from the class
#'
#' @param object The class, see [bixverse::bulk_coexp()].
#' @param epsilon Float. The epsilon parameter for the RBF function, in this
#' case the bump function.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return A list with the following elements:
#' \itemize{
#'  \item graph - The igraph
#'  \item params - A list that contains the parameters of the graph generation
#'  and general graph information (node, edge numbers).
#' }
#'
#' @export
get_cor_graph <- S7::new_generic(
  name = 'get_cor_graph',
  dispatch_args = 'object',
  fun = function(object,
                 epsilon,
                 .verbose) {
    S7::S7_dispatch()
  }
)


#' @export
#' @method get_cor_graph bulk_coexp
S7::method(get_cor_graph, bulk_coexp) <- function(object, epsilon, .verbose) {
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  checkmate::qassert(epsilon, "R1")
  checkmate::qassert(.verbose, "B1")
  # Function body
  cor_res <- S7::prop(object, "processed_data")$correlation_res
  graph_df <- cor_res$get_data_table(.verbose = .verbose) %>%
    .[, cor_abs := abs(cor)] %>%
    .[, dist := 1 - cor_abs] %>%
    .[, dist := data.table::fifelse(dist < 0, 0, dist)] %>%
    .[, affinity := rs_rbf_function(x = dist,
                                    epsilon = epsilon,
                                    rbf_type = "bump")] %>%
    .[affinity > 0] %>%
    .[, c("feature_a", "feature_b", "affinity")] %>%
    data.table::setnames(
      .,
      old = c("feature_a", "feature_b", "affinity"),
      new = c("from", "to", "weight")
    )

  if (.verbose)
    message(sprintf(
      "Generating correlation-based graph with %i edges.",
      nrow(graph_df)
    ))

  graph <- igraph::graph_from_data_frame(graph_df, directed = FALSE)
  graph <- igraph::simplify(graph)

  graph_params <- list(
    epsilon = epsilon,
    no_nodes = length(igraph::V(graph)),
    no_edges = length(igraph::E(graph))
  )

  list(graph = graph, params = graph_params)
}

#' @title Get differential correlation-based graph
#'
#' @description
#' Helper function to get a differential correlation-based igraph from the class
#'
#' @param object The class, see [bixverse::bulk_coexp()].
#' @param min_cor Float. The minimum absolute correlation that needs to be
#' present in either data set.
#' @param fdr_threshold Float. The maximum FDR that is tolerated for the
#' generation of the graph.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @return A list with the following elements:
#' \itemize{
#'  \item graph - The igraph
#'  \item params - A list that contains the parameters of the graph generation
#'  and general graph information (node, edge numbers).
#' }
#'
#' @export
get_diffcor_graph <- S7::new_generic(
  name = 'get_diffcor_graph',
  dispatch_args = 'object',
  fun = function(object,
                 min_cor = 0.2,
                 fdr_threshold = 0.05,
                 .verbose = TRUE) {
    S7::S7_dispatch()
  }
)

#' @export
#' @method get_diffcor_graph bulk_coexp
S7::method(get_diffcor_graph, bulk_coexp) <- function(object,
                                                      min_cor = 0.2,
                                                      fdr_threshold = 0.05,
                                                      .verbose = TRUE) {
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  checkmate::qassert(min_cor, "R1[0,1]")
  checkmate::qassert(fdr_threshold, "R1[0,1]")
  checkmate::qassert(.verbose, "B1")
  # Function body
  cor_res <- S7::prop(object, "processed_data")[["differential_cor_res"]]
  graph_df <- cor_res$get_data_table(.verbose = .verbose) %>%
    .[, delta_cor := cor_a - cor_b] %>%
    .[, `:=`(
      cor_a = abs(cor_a),
      cor_b = abs(cor_b),
      fdr = rs_fdr_adjustment(p_val)
    )] %>%
    .[fdr <= fdr_threshold & (cor_a >= min_cor | cor_b >= min_cor)] %>%
    .[, weight := rs_range_norm(abs(delta_cor), max_val = 1, min_val = 0.05)] %>%
    .[, c("feature_a", "feature_b")] %>%
    data.table::setnames(
      .,
      old = c("feature_a", "feature_b"),
      new = c("from", "to")
    )

  if (.verbose)
    message(sprintf(
      "Generating differential correlation-based graph with %i edges.",
      nrow(graph_df)
    ))

  graph <- igraph::graph_from_data_frame(graph_df, directed = FALSE)
  graph <- igraph::simplify(graph)

  graph_params <- list(
    min_cor = min_cor,
    fdr_threshold = fdr_threshold,
    no_nodes = length(igraph::V(graph)),
    no_edges = length(igraph::E(graph))
  )

  list(
    graph = graph,
    params = graph_params
  )
}

## getters ---------------------------------------------------------------------

#' @title Return the resolution results
#'
#' @description
#' Getter function to get the resolution results (if available).
#'
#' @param object The class, see [bixverse::bulk_coexp()].
#'
#' @return If resolution results were found, returns the data.table. Otherwise,
#' throws a warning and returns NULL.
#'
#' @export
get_resolution_res <- S7::new_generic(
  name = 'get_resolution_res',
  dispatch_args = 'object',
  fun = function(object) {
    S7::S7_dispatch()
  }
)

#' @export
#' @method get_resolution_res bulk_coexp
S7::method(get_resolution_res, bulk_coexp) <- function(object) {
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  # Body
  resolution_results <- S7::prop(object, "outputs")[['resolution_results']]
  if (is.null(resolution_results))
    warning("No resolution results found. Did you run cor_module_check_res()? Returning NULL.")

  resolution_results
}

## plotting --------------------------------------------------------------------

# This one has a shared generic...

#' @export
#'
#' @import ggplot2
#'
#' @method plot_resolution_res bulk_coexp
S7::method(plot_resolution_res, bulk_coexp) <- function(object, print_head = TRUE, ...) {
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  checkmate::qassert(print_head, "B1")
  # Body
  plot_df <- S7::prop(object, "outputs")[['resolution_results']]
  if (is.null(plot_df)) {
    warning("No resolution results found. Did you run cor_module_check_res()? Returning NULL.")
    return(NULL)
  }
  plot_df <- data.table::setorder(plot_df, -modularity)
  if (print_head)
    print(head(plot_df))
  p <- ggplot(data = plot_df,
              mapping =  aes(x = resolution, y = modularity)) +
    geom_point(
      mapping = aes(size = log10(good_clusters), fill = log10(avg_size)),
      shape = 21,
      alpha = .7
    ) +
    xlab("Leiden cluster resolution") +
    ylab("Modularity") +
    theme_minimal() +
    scale_fill_viridis_c() +
    scale_size_continuous(range = c(2, 8)) +
    labs(size = "Number of good clusters (log10)", fill = "Average cluster size (log10)") +
    ggtitle("Resolution vs. modularity", subtitle = 'With cluster number and size')
  p
}


#' @title Plot the epsilon vs. power law goodness of fit result
#'
#' @description
#' Plots the epsilon results (if they can be found in the class). The x-axis
#' reflects the different epsilon parameters for the radial basis function,
#' and the y-axis the R2 value that the resulting networks follows a power law
#' distribution (i.e., scale free topology).
#'
#' @param object The class, see [bixverse::bulk_coexp()].
#'
#' @return If epsilon results were found, returns the ggplot. Otherwise, throws
#' a warning and returns NULL.
#'
#' @export
plot_epsilon_res <- S7::new_generic(
  name = 'plot_epsilon_res',
  dispatch_args = 'object',
  fun = function(object) {
    S7::S7_dispatch()
  }
)

#' @export
#'
#' @import ggplot2
#'
#' @method plot_epsilon_res bulk_coexp
S7::method(plot_epsilon_res, bulk_coexp) <- function(object) {
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  # Body
  plot_df <- S7::prop(object, "outputs")[['epsilon_data']]
  if (is.null(plot_df)) {
    warning(
      "No resolution results found. Did you run cor_module_check_epsilon()? Returning NULL."
    )
    return(NULL)
  }

  p <- ggplot(data = plot_df, aes(x = epsilon, y = r2_vals)) +
    geom_point(size = 3, shape = 21) +
    geom_line() +
    theme_minimal() +
    ylim(0, 1) +
    xlab("Epsilon") +
    ylab("Goodness of fit (R2)") +
    ggtitle("Epsilon vs. scale free topology",
            subtitle = "Goodness of fit for log(connectivity) ~ log(p(connectivity)")

  p
}

