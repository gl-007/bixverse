# class ------------------------------------------------------------------------

#' @title Bulk RNAseq co-expression modules
#'
#' @description
#' Class for applying various co-expression module detection methods on top of
#' bulk RNAseq data.
#'
#' @param raw_data The raw count matrix. Rows = samples, columns = features.
#' @param meta_data Metadata information on the samples. It expects to have a
#' column sample_id and case_control column.
#'
#' @section Properties:
#' \describe{
#'   \item{raw_data}{A numerical matrix of the provided raw data.}
#'   \item{meta_data}{A data.table with the meta-information about the samples.}
#'   \item{processed_data}{A list in which various types of processed data will
#'   be stored.}
#'   \item{outputs}{A list in which key outputs will be stored.}
#'   \item{params}{A (nested) list that will store all the parameters of the
#'   applied function.}
#'   \item{final_results}{A data.table that will contain the final results.}
#' }
#'
#' @return Returns the `bulk_coexp` class for further operations.
#'
#' @export
bulk_coexp <- S7::new_class(
  # Names, parents
  name = "bulk_coexp",
  parent = bixverse_base_class,

  # Properties
  properties = list(
    raw_data = S7::class_double,
    meta_data = S7::class_data.frame,
    processed_data = S7::class_list,
    outputs = S7::class_list,
    params = S7::class_list,
    final_results = S7::class_any
  ),


  constructor = function(raw_data, meta_data) {
    # Checks
    checkmate::assertMatrix(raw_data, mode = "numeric")
    checkmate::assertDataFrame(meta_data)
    checkmate::assertNames(names(meta_data),
                           must.include = c("sample_id", "case_control"))
    checkmate::assertTRUE(all(rownames(raw_data) %in% meta_data$sample_id))

    params <- list(
      "original_dim" = dim(raw_data)
    )

    S7::new_object(
      S7::S7_object(),
      raw_data = raw_data,
      meta_data = data.table::as.data.table(meta_data),
      processed_data = list(),
      outputs = list(),
      params = params,
      final_results = data.table::data.table()
    )
  }
)

# utils ------------------------------------------------------------------------

## getters ---------------------------------------------------------------------

#' Return the outputs from bulk_coexp
#'
#' @description
#' Getter function to extract the outputs from the [bixverse::bulk_coexp()]
#' class.
#'
#' @param object The underlying `bulk_coexp` class, see
#' [bixverse::bulk_coexp()].
#'
#' @return Returns the outputs stored in the class.
#'
#' @export
get_outputs <- S7::new_generic(
  name = "get_outputs",
  dispatch_args = "object",
  fun = function(object) {
    S7::S7_dispatch()
  }
)


#' @method get_outputs bulk_coexp
#'
#' @export
S7::method(get_outputs, bulk_coexp) <-
  function(object) {
    # Checks
    checkmate::assertClass(
      object, "bixverse::bulk_coexp"
    )

    # Return
    return(S7::prop(object, "outputs"))
  }

## print -----------------------------------------------------------------------

#' @name print.bulk_coexp
#' @title print Method for bulk_coexp object
#'
#' @description
#' Print a bulk_coexp object.
#'
#' @param x An object of class `bulk_coexp`.
#' @param ... Additional arguments (currently not used).
#'
#' @returns Invisibly returns `x`.
#'
#' @method print bulk_coexp
S7::method(print, bulk_coexp) <- function(x, ...) {
  # Pre-processing
  preprocessed <- !is.null(S7::prop(x, "processed_data")[['processed_data']])
  features <- if(preprocessed) {
    no_features <- S7::prop(x, "processed_data")[['feature_meta']] %>%
      .[(hvg)] %>%
      nrow()
    sprintf('  Number of HVG: %i.\n', no_features)
  } else {
    ''
  }
  # Prints for different methods
  method_info <- if (is.null(S7::prop(x, "params")[["detection_method"]])) {
    # Case of nothing has been applied
    ''
  } else if (S7::prop(x, "params")[["detection_method"]] == "cPCA") {
    # Contrastive PCA
    no_intersecting_features <- length(S7::prop(x, "params")[["cPCA_params"]][['intersecting_features']])
    paste0(
      " Detection method: cPCA.\n",
      sprintf("  No of intersecting features: %i.\n", no_intersecting_features)
    )
  } else if (S7::prop(x, "params")[["detection_method"]] == "correlation-based") {
    # For simple correlations
    non_parametric <- S7::prop(x, "params")[['correlation_params']][['spearman']]
    graph_generated <- !is.null(S7::prop(x, "params")[['correlation_graph']][['no_nodes']])
    paste0(
      " Detection method: correlation based.\n",
      sprintf("  Non-parametric correlation applied: %s.\n", non_parametric),
      sprintf("  Graph generated: %s.\n", graph_generated)
    )
  }

  cat(
    "Bulk co-expression module class (bulk_coexp).\n",
    " Pre-processing done: ",
    preprocessed,
    ".\n",
    features,
    method_info,
    sep = ''
  )

  invisible(x)
}

# methods ----------------------------------------------------------------------

#' Process the raw data
#'
#' @description
#' Function to do general pre-processing on top of the [bixverse::bulk_coexp()].
#' Options to do scaling, HVG selection, etc.
#'
#' @param object The underlying class, see [bixverse::bulk_coexp()].
#' @param hvg Integer or float. If an integer, the top `hvg` genes will be
#' included; if float, the float has to be between 0 and 1, representing the
#' percentage of genes to include.
#' @param mad_threshold Float. Instead of of selecting number or proportion of
#' genes, you can also provide a mad_threshold.
#' @param scaling Boolean. Shall the data be scaled.
#' @param scaling_type String. You have the option to use normal scaling or
#' robust scaling.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @return Returns the class with the `processed_data` data slot populated and
#' applied parameters added to the `params` slot.
#'
#' @export
#'
#' @import data.table
#' @importFrom magrittr `%>%`
#' @importFrom magrittr `%$%`
preprocess_bulk_coexp <- S7::new_generic(
  "preprocess_bulk_coexp",
  "object",
  fun = function(object,
                 hvg = NULL,
                 mad_threshold = NULL,
                 scaling = FALSE,
                 scaling_type = c("normal", "robust"),
                 .verbose = TRUE) {
    S7::S7_dispatch()
  }
)

#' @method preprocess_bulk_coexp bulk_coexp
#' @export
S7::method(preprocess_bulk_coexp, bulk_coexp) <- function(object,
                                                          hvg = NULL,
                                                          mad_threshold = NULL,
                                                          scaling = FALSE,
                                                          scaling_type = c("normal", "robust"),
                                                          .verbose = TRUE) {
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  checkmate::qassert(mad_threshold, c("R1", "0"))
  nfeatures <- S7::prop(object, "params")[['original_dim']][2]
  checkmate::qassert(hvg, c("R1[0,1]", sprintf("I1[0,%i]", nfeatures), "0"))
  checkmate::qassert(scaling, "B1")
  if (scaling) {
    checkmate::assertChoice(scaling_type, c("normal", "robust"))
  }

  mat <- S7::prop(object, "raw_data")

  feature_meta <- data.table::data.table(
    feature_name = colnames(mat),
    mean_exp = colMeans(mat),
    MAD = matrixStats::colMads(mat),
    var_exp = matrixStats::colVars(mat)
  ) %>%
    data.table::setorder(-MAD)

  if (!is.null(hvg) & !is.null(mad_threshold)) {
    choice <- menu(
      choices = c("MAD threshold", "HVG threshold"),
      title = paste(
        "You have provided both a MAD and HVG",
        "threshold. Which one do you want to use?"
      )
    )
    if (choice == 1) {
      hvg <- NULL
    } else {
      mad_threshold <- NULL
    }
  }

  if (is.null(hvg) & is.null(mad_threshold)) {
    hvg = 1
  }

  if (!is.null(hvg)) {
    no_genes_to_take <-
      ifelse(is.integer(hvg), hvg, ceiling(hvg * ncol(mat)))
    hvg_genes <- feature_meta[1:no_genes_to_take, feature_name]
  } else {
    hvg_genes <- feature_meta[MAD >= mad_threshold, feature_name]
  }

  if (.verbose)
    message(sprintf("A total of %i genes will be included.", length(hvg_genes)))

  feature_meta[, hvg := feature_name %in% hvg_genes]

  # Process the matrix
  matrix_processed <- mat[, hvg_genes]

  if (scaling) {
    fun <-
      ifelse(scaling_type == "normal",
             "scale",
             "bixverse::robust_scale")
    matrix_processed <- rlang::eval_tidy(rlang::quo(apply(
      matrix_processed, 1, !!!rlang::parse_exprs(fun)
    )))
  }

  processing_params <- list(
    mad_threshold = if (is.null(mad_threshold)) {
      "not applicable"
    } else {
      mad_threshold
    },
    hvg = if (is.null(hvg)) {
      "not applicable"
    } else {
      hvg
    },
    scaling = scaling,
    scaling_type = if (!scaling) {
      "not applicable"
    } else {
      scaling_type
    }
  )

  S7::prop(object, "params")[['preprocessing']] <-
    processing_params
  S7::prop(object, "processed_data")[['processed_data']] <-
    matrix_processed
  S7::prop(object, "processed_data")[['feature_meta']] <-
    feature_meta

  # Return
  object
}

# plots ------------------------------------------------------------------------

#' @title Plot the highly variable genes
#'
#' @description
#' Plots the median-absolute deviation of the genes and applied thresholds.
#' Expects that [bixverse::preprocess_bulk_coexp()] was run and will throw an
#' error otherwise.
#'
#' @param object The underlying class, see [bixverse::bulk_coexp()].
#' @param bins Integer. Number of bins to plot.
#'
plot_hvgs <- S7::new_generic(
  "plot_hvgs",
  "object",
  fun = function(object, bins = 50L) {
    S7::S7_dispatch()
  }
)

#' @method plot_hvgs bulk_coexp
#'
#' @import ggplot2
#'
#' @export
S7::method(plot_hvgs, bulk_coexp) <- function(object, bins = 50L) {
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  checkmate::qassert(bins, "I1")
  # Early return
  if(is.null(S7::prop(object, "params")[['preprocessing']])) {
    warning("No pre-processing data found. Returning NULL.")
    return(NULL)
  }
  plot_df <- S7::prop(object, "processed_data")[['feature_meta']]

  p <- ggplot(data = plot_df, mapping = aes(x = MAD)) +
    geom_histogram(mapping = aes(fill = hvg),
                   bins = 50L,
                   color = 'black',
                   alpha = 0.7) +
    scale_fill_manual(values = setNames(c('orange', 'grey'), c(TRUE, FALSE))) +
    ggtitle("Distribution of MAD across the genes", subtitle = "And included genes") +
    theme_minimal()

  mad_threshold <- S7::prop(object, "params")[['preprocessing']][['mad_threshold']]

  if (mad_threshold != "not applicable") {
    p <- p + geom_vline(xintercept = mad_threshold,
                        linetype = 'dashed',
                        color = 'darkred')
  }

  return(p)
}
