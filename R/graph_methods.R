# network_diffusions -----------------------------------------------------------

## diffusion methods -----------------------------------------------------------

#' Diffuse seed genes over a network
#'
#' @description
#' This function takes a diffusion vector and leverages personalised page-rank
#' diffusion to identify influential nodes. These can be used subsequently for
#' community detection or check AUROC values given a set of genes.
#'
#' @param object The underlying class [bixverse::network_diffusions()].
#' @param diffusion_vector A named vector with values to use for the reset
#' parameter in the personalised page-rank diffusion. Names should represent
#' node names of the graph.
#' @param summarisation If there are duplicated names in the `diffusion_vector`
#' how to summarise these.
#'
#' @return The class with added diffusion score based on a single set of seed
#' genes. Additionally, the seed genes are stored in the class.
#'
#' @export
diffuse_seed_nodes <- S7::new_generic(
  name = "diffuse_seed_nodes",
  dispatch_args = "object",
  fun = function(object,
                 diffusion_vector,
                 summarisation = c("max", "mean", "harmonic_sum")) {
    S7::S7_dispatch()
  }
)


#' @export
#'
#' @importFrom magrittr `%>%`
#'
#' @method diffuse_seed_nodes network_diffusions
S7::method(diffuse_seed_nodes, network_diffusions) <-
  function(object,
           diffusion_vector,
           summarisation = c("max", "mean", "harmonic_sum")) {
    # Checks
    checkmate::assertClass(object, "bixverse::network_diffusions")
    checkmate::assertNumeric(diffusion_vector)
    checkmate::assertNamed(diffusion_vector, .var.name = "diffusion_vector")
    checkmate::assertChoice(summarisation, c("max", "mean", "harmonic_sum"))

    # Body
    ## Create the diffusion vector
    diffusion_vector <- .summarise_scores(diffusion_vector, summarisation = summarisation)
    nodes_names <- igraph::V(S7::prop(object, "graph"))$name
    seed_nodes <- intersect(names(diffusion_vector), nodes_names)
    diff_vec <- rep(0, length(nodes_names)) %>% `names<-`(nodes_names)
    for (node in seed_nodes) {
      diff_vec[node] <- diffusion_vector[node]
    }

    if (sum(diff_vec) == 0) {
      stop("No scores found to diffuse over the network. Please check the names
      and/or values.")
    }

    ## Create the page-rank diffusion
    page_rank_score <- igraph::page_rank(S7::prop(object, "graph"), personalized = diff_vec)

    ## Assign and return
    S7::prop(object, "diffusion_res") <- page_rank_score$vector
    S7::prop(object, "params")["diffusion_type"] <- "single"
    S7::prop(object, "params")[["seed_nodes"]] <- seed_nodes

    return(object)
  }

#' Diffuse seed genes in a tied manner over a network
#'
#' @description
#' This function takes two sets of diffusion vector and leverages tied diffusion
#' to identify an intersection of influential nodes.
#'
#' @param object The underlying class [bixverse::network_diffusions()].
#' @param diffusion_vector_1 The first named vector with values to use for the
#' reset parameter in the personalised page-rank diffusion. Names should
#' represent node names of the graph.
#' @param diffusion_vector_2 The second named vector with values to use for the
#' reset parameter in the personalised page-rank diffusion. Names should
#' represent node names of the graph.
#' @param summarisation If there are duplicated names in the `diffusion_vector`
#' how to summarise
#' these.
#' @param score_aggregation How to summarise the tied scores.
#' @param .verbose Controls verbosity of the function.
#'
#' @return The class with added diffusion score based on a two sets of seed
#' genes. Additionally, the seed genes are stored in the class.
#'
#' @export
tied_diffusion <- S7::new_generic(
  name = "tied_diffusion",
  dispatch_args = "object",
  fun = function(object,
                 diffusion_vector_1,
                 diffusion_vector_2,
                 summarisation = c("max", "mean", "harmonic_sum"),
                 score_aggregation = c("min", "max", "mean"),
                 .verbose = FALSE) {
    S7::S7_dispatch()
  }
)


#' @export
#'
#' @importFrom magrittr `%>%`
#'
#' @method tied_diffusion network_diffusions
S7::method(tied_diffusion, network_diffusions) <-
  function(object,
           diffusion_vector_1,
           diffusion_vector_2,
           summarisation = c("max", "mean", "harmonic_sum"),
           score_aggregation = c("min", "max", "mean"),
           .verbose = FALSE) {
    # Checks
    checkmate::assertClass(object, "bixverse::network_diffusions")
    checkmate::assertNumeric(diffusion_vector_1)
    checkmate::assertNamed(diffusion_vector_1, .var.name = "diffusion_vector_1")
    checkmate::assertNumeric(diffusion_vector_2)
    checkmate::assertNamed(diffusion_vector_2, .var.name = "diffusion_vector_2")
    checkmate::assertChoice(summarisation, c("max", "mean", "harmonic_sum"))
    checkmate::assertChoice(score_aggregation, c("min", "max", "mean"))
    checkmate::qassert(.verbose, "B1")

    # Body
    ## Create the diffusion vectors
    diffusion_vector_1 <- .summarise_scores(diffusion_vector_1, summarisation = summarisation)
    diffusion_vector_2 <- .summarise_scores(diffusion_vector_2, summarisation = summarisation)
    nodes_names <- igraph::V(S7::prop(object, "graph"))$name
    seed_nodes_1 <- intersect(names(diffusion_vector_1), nodes_names)
    seed_nodes_2 <- intersect(names(diffusion_vector_2), nodes_names)
    diff_vec_1 <- diff_vec_2 <- rep(0, length(nodes_names)) %>%
      `names<-`(nodes_names)
    for (node in seed_nodes_1) {
      diff_vec_1[node] <- diffusion_vector_1[node]
    }
    for (node in seed_nodes_2) {
      diff_vec_2[node] <- diffusion_vector_2[node]
    }
    if ((sum(diff_vec_1) == 0) || (sum(diff_vec_1) == 0)) {
      stop(
        "No scores found on first and/or second of the diffusion vectors.
        Please check the names and/or values."
      )
    }

    ## First diffusion
    score_1 <- igraph::page_rank(S7::prop(object, "graph"), personalized = diff_vec_1)$vector

    ## Second diffusion
    directed <- S7::prop(object, "params")[["directed_graph"]]
    score_2 <- if (directed) {
      if (.verbose) {
        message(
          "Directed graph found. Function will use transpose of adjacency for
          second diffusion."
        )
      }
      adj <- igraph::as_adjacency_matrix(S7::prop(object, "graph"))
      adj_t <- Matrix::t(adj)
      igraph_obj_t <- igraph::graph_from_adjacency_matrix(adj_t)
      igraph::page_rank(igraph_obj_t, personalized = diff_vec_2)$vector
    } else {
      if (.verbose) {
        message("Undirected graph found. Using graph as is for second
        diffusion.")
      }
      igraph::page_rank(S7::prop(object, "graph"), personalized = diff_vec_2)$vector
    }

    ## Summarise the scores
    final_tiedie_diffusion <- switch(
      score_aggregation,
      "min" = pmin(score_1, score_2),
      "max" = pmax(score_1, score_2),
      rowMeans(cbind(score_1, score_2))
    )

    ## Assign and return
    S7::prop(object, "diffusion_res") <- final_tiedie_diffusion
    S7::prop(object, "params")["diffusion_type"] <- "tied"
    S7::prop(object, "params")[["seed_nodes"]] <- list("set_1" = seed_nodes_1, "set_2" = seed_nodes_2)

    return(object)
  }

## community detection ---------------------------------------------------------

#' Identify privileged communities based on a given diffusion vector
#'
#' @description Detects privileged communities after a diffusion based on seed
#' nodes.
#'
#' @param object The underlying class [bixverse::network_diffusions()].
#' @param diffusion_threshold Float. How much of the network to keep based on
#' the diffusion values. 0.25 for example would keep the 25% nodes with the
#' highest scores. This was the default in the original paper.
#' @param community_params List. Parameters for the community detection within
#' the reduced network, see [bixverse::params_community_detection()]. A list
#' with the following items:
#' \itemize{
#'  \item max_nodes - Integer. Number of maximum nodes per community. Larger
#'  communities will be recursively subclustered.
#'  \item min_nodes - Integer. Minimum number of nodes per community.
#'  \item min_seed_nodes - Integer. Minimum number of seed genes that have to
#'  be found in a given community.
#'  \item initial_res - Float. Initial resolution parameter for the Leiden
#'  clustering.
#' }
#' @param seed Random seed.
#' @param .verbose Controls the verbosity of the function.
#' @param .max_iters Controls how many iterations shall be tried for the
#' sub-clustering. To note, in each iteration of the sub-clustering, the
#' resolution parameter is increased by 0.05, to identify more granular
#' communities within the sub communities.
#'
#' @return The class with added diffusion community detection results (if any
#' could be identified with the provided parameters).
#'
#' @export
community_detection <- S7::new_generic(
  name = "community_detection",
  dispatch_args = "object",
  fun = function(object,
                 diffusion_threshold = 0.25,
                 community_params = params_community_detection(),
                 seed = 42L,
                 .verbose = FALSE,
                 .max_iters = 100L) {
    S7::S7_dispatch()
  }
)


#' @export
#'
#' @import data.table
#' @importFrom magrittr `%>%`
#'
#' @method community_detection network_diffusions
S7::method(community_detection, network_diffusions) <- function(object,
                                                                diffusion_threshold = 0.25,
                                                                community_params = params_community_detection(),
                                                                seed = 42L,
                                                                .verbose = FALSE,
                                                                .max_iters = 100L) {
  # Bindings
  `.` <- N <- cluster_id <- node_id <- cluster_size <- seed_nodes_no <-
    seed_nodes_no <- seed_nodes_1 <- seed_nodes_2 <- NULL
  # Checks
  checkmate::assertClass(object, "bixverse::network_diffusions")
  checkmate::qassert(diffusion_threshold, "R1[0,1]")
  assertCommunityParams(community_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")
  checkmate::qassert(.max_iters, "I1")

  # Body
  ## Reduce the graph
  diffusion_score <- S7::prop(object, "diffusion_res")
  # Early return
  if (length(diffusion_score) == 0) {
    warning(
      "The diffusion score has length 0. Likely you did not run the diffusion
      methods. Returning class as is."
    )
    return(object)
  }
  nodes_to_include <- diffusion_score %>%
    sort(decreasing = TRUE) %>%
    .[1:ceiling(diffusion_threshold * length(diffusion_score))]

  red_graph <- igraph::subgraph(S7::prop(object, "graph"), names(nodes_to_include))

  ## First clustering
  set.seed(seed)

  if (S7::prop(object, "params")$directed_graph) {
    red_graph <- igraph::as_undirected(red_graph, mode = "each")
  }

  final_clusters <- with(community_params, {
    first_clusters <- igraph::cluster_leiden(
      red_graph,
      resolution = intial_res,
      n_iterations = 5,
      objective_function = "modularity"
    )

    clusters_df <- data.table::data.table(node_id = first_clusters$names,
                                          cluster_id = first_clusters$membership)

    node_frequency <- clusters_df[, .N, .(cluster_id)]

    ## Subclustering
    clusters_with_too_many_nodes <- node_frequency[N > max_nodes, cluster_id]
    final_clusters <- clusters_df[!cluster_id %in% clusters_with_too_many_nodes]

    for (i in seq_along(clusters_with_too_many_nodes)) {
      cluster_i <- clusters_with_too_many_nodes[i]
      nodes_in_cluster <- clusters_df[cluster_id == cluster_i, node_id]
      finalised_clusters <- data.table()
      # Loop through, until all clusters are below the minimum genes or max
      # iterations is hit
      l <- 1
      while (length(nodes_in_cluster) != 0) {
        set.seed(seed + l)
        if (.verbose) {
          message("Cluster ", i, " gets subclustered. Iter: ", l)
        }
        red_graph_l <- igraph::subgraph(red_graph,
                                        data.table::chmatch(nodes_in_cluster, igraph::V(red_graph)$name))

        clusters_red <- igraph::cluster_leiden(
          red_graph_l,
          resolution = intial_res + l * 0.05,
          n_iterations = 5,
          objective_function = "modularity"
        )

        subclusters <- data.table(node_id = clusters_red$names,
                                  cluster_id = clusters_red$membership)
        subclusters_frequency <- subclusters[, .N, .(cluster_id)]
        clusters_small_enough <- subclusters_frequency[N <= max_nodes, cluster_id]

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

    return(final_clusters)
  })

  ## Add the seed node information based on diffusion type
  diffusion_type <- S7::prop(object, "params")$diffusion_type

  final_node_frequency <- with(community_params, {
    if (diffusion_type == "single") {
      seed_nodes <- S7::prop(object, "params")$seed_nodes

      final_clusters[, .(
        cluster_size = length(node_id),
        seed_nodes_no = sum(node_id %in% seed_nodes)
      ), .(cluster_id)]
    } else {
      seed_nodes_set_1 <- S7::prop(object, "params")$seed_nodes$set_1
      seed_nodes_set_2 <- S7::prop(object, "params")$seed_nodes$set_2

      final_clusters[, .(
        cluster_size = length(node_id),
        seed_nodes_1 = sum(node_id %in% seed_nodes_set_1),
        seed_nodes_2 = sum(node_id %in% seed_nodes_set_2)
      ), .(cluster_id)]
    }
  })

  ## Finalise the clusters
  clusters_to_take <- with(community_params, {
    if (diffusion_type == "single") {
      final_node_frequency[cluster_size >= min_nodes &
                             seed_nodes_no >= min_seed_nodes, cluster_id]
    } else {
      final_node_frequency[cluster_size >= min_nodes &
                             seed_nodes_1 >= min_seed_nodes &
                             seed_nodes_2 >= min_seed_nodes, cluster_id]
    }
  })

  # Early return
  if (length(clusters_to_take) == 0) {
    warning("No communities found with the given parameters.
    Returning class as is.")
    return(object)
  }

  finalised_clusters_clean <- final_clusters[cluster_id %in% clusters_to_take]

  ks_vals <- vector(mode = "numeric", length = length(clusters_to_take))

  for (i in seq_along(clusters_to_take)) {
    cluster <- clusters_to_take[i]
    cluster_nodes <- finalised_clusters_clean[cluster_id == cluster, node_id]
    ks <- suppressWarnings(ks.test(diffusion_score[cluster_nodes], diffusion_score[which(!names(diffusion_score) %in% cluster_nodes)], alternative = "less"))
    ks_vals[i] <- ks$p.value
  }

  ks_val_df <- data.table(cluster_id = clusters_to_take, ks_pval = ks_vals)

  final_result <- purrr::reduce(list(finalised_clusters_clean, ks_val_df, final_node_frequency),
                                merge,
                                by = "cluster_id") %>%
    .[, `:=`(diffusion_score = diffusion_score[node_id],
             cluster_id = as.character(cluster_id))]


  cluster_name_prettifier <- setNames(paste("cluster", seq_along(unique(
    final_result$cluster_id
  )), sep = "_"),
  unique(final_result$cluster_id))

  final_result[, cluster_id := cluster_name_prettifier[cluster_id]]

  ## Assign and return
  S7::prop(object, "final_results") <- final_result
  S7::prop(object, "params")[["community_params"]] <- with(
    community_params,
    list(
      diffusion_threshold = diffusion_threshold,
      max_nodes = max_nodes,
      min_nodes = min_seed_nodes,
      min_seed_nodes = min_seed_nodes
    )
  )

  return(object)
}

## utils ----

#' Calculate the AUROC for a diffusion score
#'
#' @description
#' This functions can take a given `network_diffusions` class and calculates an
#' AUC and generates a Z-score based on random permutation of `random_aucs` for
#' test for statistical significance if desired.
#'
#' @param object The underlying class [bixverse::network_diffusions()].
#' @param hit_nodes String vector. Which nodes in the graph are considered a
#' 'hit'.
#' @param auc_iters Integer. How many iterations to run to approximate the
#' AUROC.
#' @param random_aucs Integer. How many random AUROCs to calculate to estimate
#' the Z-score. Only of relevance if permutation test is set to `TRUE`.
#' @param permutation_test Boolean. Shall a permutation based Z-score be
#' calculated.
#' @param seed Integer. Random seed.
#'
#' @return List with AUC and Z-score as the two named elements if permutations
#' test set to TRUE; otherwise just the AUC.
#'
#' @export
calculate_diffusion_auc <- S7::new_generic(
  name = "calculate_diffusion_auc",
  dispatch_args = "object",
  fun = function(object,
                 hit_nodes,
                 auc_iters = 10000L,
                 random_aucs = 1000L,
                 permutation_test = FALSE,
                 seed = 42L) {
    S7::S7_dispatch()
  }
)


#' @export
#'
#' @importFrom magrittr `%>%`
#'
#' @method tied_diffusion network_diffusions
S7::method(calculate_diffusion_auc, network_diffusions) <-
  function(object,
           hit_nodes,
           auc_iters = 10000L,
           random_aucs = 1000L,
           permutation_test = FALSE,
           seed = 42L) {
    # Checks
    checkmate::assertClass(object, "bixverse::network_diffusions")
    checkmate::qassert(hit_nodes, "S+")
    checkmate::qassert(auc_iters, "I1")
    checkmate::qassert(seed, "I1")
    checkmate::qassert(permutation_test, "B1")
    if (permutation_test) {
      checkmate::qassert(random_aucs, "I1")
    }
    # Body
    diffusion_score <- S7::prop(object, "diffusion_res")
    if (length(diffusion_score) == 0) {
      warning(
        "The diffusion score has length 0. Likely you did not run the diffusion
        methods. Returning NULL."
      )
      return(NULL)
    }
    pos_scores <- diffusion_score[hit_nodes]
    neg_scores <- diffusion_score[which(!names(diffusion_score) %in% hit_nodes)]
    auc <- rs_fast_auc(
      pos_scores = pos_scores,
      neg_scores = neg_scores,
      iters = auc_iters,
      seed = seed
    )
    if (permutation_test) {
      random_aucs <- rs_create_random_aucs(
        score_vec = diffusion_score,
        size_pos = length(pos_scores),
        random_iters = random_aucs,
        auc_iters = auc_iters,
        seed = seed
      )

      z <- (auc - mean(random_aucs)) / sd(random_aucs)

      to_ret <- list(auc = auc, z = z)
    } else {
      to_ret <- auc
    }

    return(to_ret)
  }

# rbh_graph ----

## graph generation ----

#' Generate an RBH graph.
#'
#' @description
#' This function will generate an RBH graph based on set similarity between
#' gene modules. You have the option to use an overlap coefficient instead of
#' Jaccard similarity and to specify a minimum similarity.
#'
#' @param object The underlying class, see [bixverse::rbh_graph()].
#' @param minimum_similarity The minimum similarity to create an edge.
#' @param overlap_coefficient Shall the overlap coefficient be used instead of
#' Jaccard similarity.
#' @param .debug Debug flat that will create print messages from Rust.
#'
#' @return The class with added properties.
#'
#' @export
generate_rbh_graph <- S7::new_generic(
  name = "generate_rbh_graph",
  dispatch_args = "object",
  fun = function(object,
                 minimum_similarity,
                 overlap_coefficient = FALSE,
                 .debug = FALSE) {
    S7::S7_dispatch()
  }
)


#' @export
#'
#' @importFrom magrittr `%>%`
#'
#' @method generate_rbh_graph rbh_graph
S7::method(generate_rbh_graph, rbh_graph) <-
  function(object,
           minimum_similarity,
           overlap_coefficient = FALSE,
           .debug = FALSE) {
    # Assigns
    origin_modules <- `.` <- similiarity <- origin <- target <-
      target_modules <- NULL
    # Checks
    checkmate::assertClass(object, "bixverse::rbh_graph")
    checkmate::qassert(minimum_similarity, "R[0, 1]")
    checkmate::qassert(overlap_coefficient, "B1")
    checkmate::qassert(.debug, "B1")

    # Body
    list_of_list <- S7::prop(object, "module_data")

    rbh_results <- rs_rbh_sets(
      module_list = list_of_list,
      overlap_coefficient = overlap_coefficient,
      min_similarity = minimum_similarity,
      debug = .debug
    )

    rbh_results$origin_modules[rbh_results$origin_modules == "NA"] <- NA
    rbh_results$target_modules[rbh_results$target_modules == "NA"] <- NA

    origin_vector <- unlist(purrr::map2(rbh_results$origin, rbh_results$comparisons, ~ {
      rep(.x, each = .y)
    }))

    target_vector <- unlist(purrr::map2(rbh_results$target, rbh_results$comparisons, ~ {
      rep(.x, each = .y)
    }))

    rbh_results_dt <- data.table::data.table(
      origin = origin_vector,
      target = target_vector,
      origin_modules = rbh_results$origin_modules,
      target_modules = rbh_results$target_modules,
      similiarity = rbh_results$similarity
    ) %>%
      .[!is.na(origin_modules) &
          similiarity >= minimum_similarity] %>%
      .[, `:=`(
        combined_origin = paste(origin, origin_modules, sep = "_"),
        combined_target = paste(target, target_modules, sep = "_")
      )]

    edge_dt <- rbh_results_dt[, c("combined_origin", "combined_target", "similiarity")] %>%
      data.table::setnames(
        .,
        old = c("combined_origin", "combined_target", "similiarity"),
        new = c("from", "to", "weight")
      )

    rbh_igraph <- igraph::graph_from_data_frame(edge_dt, directed = FALSE)

    ## Assign and return
    S7::prop(object, "rbh_edge_df") <- rbh_results_dt
    S7::prop(object, "rbh_graph") <- rbh_igraph
    S7::prop(object, "params")[["rbh_graph_gen"]] <- list(minimum_similarity = minimum_similarity,
                                                          overlap_coefficient = overlap_coefficient)

    return(object)
  }


#' Find RBH communities
#'
#' @description
#' This function will identify communities in the reciprocal best hit (RBH)
#' graph. It will iterate through resolutions and add the results to the
#' class. Additionally, a column will be added that signifies the resolution
#' with the best modularity.
#'
#' @param object The underlying class, see [bixverse::rbh_graph()].
#' @param resolution_params List. Parameters for the resolution search, see
#' [bixverse::params_graph_resolution()]. Contains:
#' \itemize{
#'  \item min_res - Float. Minimum resolution to test.
#'  \item max_res - Float. Maximum resolution to test.
#'  \item number_res - Integer. Number of resolutions to test between the
#'  `max_res` and `min_res.`
#' }
#' @param parallel Boolean. Shall the resolution search be in parallel.
#' @param max_workers Integer. Number of maximum cores to use. Defaults to half
#' of the identified cores.
#' @param random_seed Integer. Random seed for reproducibility.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return The class with added community detection results.
#'
#' @export
find_rbh_communities <- S7::new_generic(
  name = "find_rbh_communities",
  dispatch_args = "object",
  fun = function(object,
                 resolution_params = params_graph_resolution(),
                 max_workers = as.integer(parallel::detectCores() / 2),
                 parallel = TRUE,
                 random_seed = 42L,
                 .verbose = TRUE) {
    S7::S7_dispatch()
  }
)


#' @export
#'
#' @importFrom magrittr `%>%`
#' @importFrom future plan multisession sequential
#' @import data.table
#'
#' @method find_rbh_communities rbh_graph
S7::method(find_rbh_communities, rbh_graph) <- function(object,
                                                        resolution_params = params_graph_resolution(),
                                                        max_workers = as.integer(parallel::detectCores() / 2),
                                                        parallel = TRUE,
                                                        random_seed = 42L,
                                                        .verbose = TRUE) {
  # Checks
  checkmate::assertClass(object, "bixverse::rbh_graph")
  assertGraphResParams(resolution_params)
  checkmate::qassert(parallel, "B1")
  checkmate::qassert(max_workers, "I1")
  checkmate::qassert(.verbose, "B1")

  if (is.null(S7::prop(object, "rbh_graph"))) {
    warning("No RBH graph yet generated. Returning class as is.")
    return(object)
  }

  graph <- S7::prop(object, "rbh_graph")

  if (.verbose)
    message(sprintf("Iterating through %i resolutions", length(resolutions)))

  resolutions <- with(resolution_params, exp(seq(log(min_res), log(max_res), length.out = number_res)))

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

  community_df_res[, best_modularity := modularity == max(modularity)]

  S7::prop(object, "final_results") <- community_df_res

  return(object)
}


# helpers ----------------------------------------------------------------------

## utils -----------------------------------------------------------------------

#' Summarise gene scores if they are duplicates.
#'
#' @param x Named numeric.
#' @param summarisation String. Which summary function to use.
#'
#' @return Named numeric.
#'
#' @export
#'
#' @importFrom magrittr `%$%`
.summarise_scores <- function(x,
                              summarisation = c("max", "mean", "harmonic_sum")) {
  # devtools::check() stuff
  value <- . <- node_name <- setNames <- NULL
  # Checks
  checkmate::assertNumeric(x)
  checkmate::assertNamed(x, .var.name = "x")
  checkmate::assertChoice(summarisation, c("max", "mean", "harmonic_sum"))
  # Body
  dt <- data.table::data.table(node_name = names(x), value = x)
  summary_fun <- switch(
    summarisation,
    "mean" = rlang::expr(mean(value)),
    "max" = rlang::expr(max(value)),
    rlang::expr(bixverse::ot_harmonic_score(value)) # Default case
  )
  res <-
    rlang::eval_tidy(rlang::quo(dt[, .(value = !!summary_fun), .(node_name)])) %$%
    setNames(value, node_name)
  res
}

## plots -----------------------------------------------------------------------

#' @export
#'
#' @import ggplot2
#'
#' @method plot_resolution_res rbh_graph
S7::method(plot_resolution_res, rbh_graph) <- function(object, print_head = TRUE, ...) {
  checkmate::assertClass(object, "bixverse::rbh_graph")
  # Ignoring print_head for this class

  plot_df <- S7::prop(object, "final_results")
  if (is.null(plot_df)) {
    warning("No resolution results found. Did you run cor_module_check_res()? Returning NULL.")
    return(NULL)
  }
  plot_df <- plot_df[, c("resolution", "modularity")] %>%
    unique()
  p <- ggplot(data = plot_df,
              mapping =  aes(x = resolution, y = modularity)) +
    geom_point(size = 3,
               shape = 21,
               alpha = .7) +
    xlab("Leiden cluster resolution") +
    ylab("Modularity") +
    theme_minimal() +
    ggtitle("Resolution vs. modularity")
  p
}
