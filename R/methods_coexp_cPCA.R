# processing -------------------------------------------------------------------

#' Prepare class for contrastive PCA
#'
#' @description
#' This function will prepare the `bulk_coexp` for subsequent usage of the
#' contrastive PCA functions. This is based on the work of Abid, et al.
#'
#' @param object The underlying class, see [bixverse::bulk_coexp()].
#' @param background_matrix Numeric matrix. The background matrix you wish to
#' remove. You should apply any data transformation to this matrix, too!
#' @param scale Boolean. Shall the data be scaled. Defaults to FALSE.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return `bulk_coexp` with the needed data for contrastive PCA in the
#' properties of the class.
#'
#' @references Abid, et al., Nature Communications, 2018
#'
#' @export
contrastive_pca_processing <- S7::new_generic(
  name = "contrastive_pca_processing",
  dispatch_args = "object",
  fun = function(object,
                 background_matrix,
                 scale = FALSE,
                 .verbose = TRUE) {
    S7::S7_dispatch()
  }
)

#' @importFrom magrittr `%>%`
#' @importFrom magrittr `%$%`
#' @import data.table
#'
#' @export
#'
#' @method contrastive_pca_processing bulk_coexp
S7::method(contrastive_pca_processing, bulk_coexp) <-
  function(object,
           background_matrix,
           scale = FALSE,
           .verbose = TRUE) {
    # Checks
    checkmate::assertClass(object, "bixverse::bulk_coexp")
    checkmate::assertMatrix(background_matrix, mode = "numeric")
    checkmate::qassert(scale, "B1")
    checkmate::qassert(.verbose, "B1")

    # Function body
    if(purrr::is_empty(S7::prop(object, "processed_data")[['processed_data']])) {
      warning("No pre-processed data found. Defaulting to the raw data")
      target_mat <- S7::prop(object, "raw_data")
    } else {
      target_mat <- S7::prop(object, "processed_data")[['processed_data']]
    }

    background_mat <- background_mat

    intersecting_features <- intersect(
      colnames(target_mat),
      colnames(background_mat)
    )
    if (.verbose) {
      message(sprintf(
        "A total of %i features/genes were identified",
        length(intersecting_features)
      ))
    }

    target_mat <- target_mat[, intersecting_features]
    background_mat <- background_mat[, intersecting_features]

    if (scale) {
      target_mat <- scale(target_mat, scale = scale)
      background_mat <- scale(background_mat, scale = scale)
    }

    target_covar <- rs_covariance(target_mat)
    background_covar <- rs_covariance(background_mat)

    rownames(target_covar) <- rownames(background_covar) <- intersecting_features
    colnames(target_covar) <- colnames(background_covar) <- intersecting_features

    internal_params <- list(
      "intersecting_features" = intersecting_features,
      "scaled_data" = scale
    )

    # Data
    S7::prop(object, "processed_data")[["target_mat"]] <-
      target_mat
    S7::prop(object, "processed_data")[["background_mat"]] <-
      background_mat
    # Covariance matrices
    S7::prop(object, "processed_data")[["target_covar"]] <-
      target_covar
    S7::prop(object, "processed_data")[["background_covar"]] <-
      background_covar

    # Set the object to a cPCA analysis
    S7::prop(object, "params")["detection_method"] <- "cPCA"
    S7::prop(object, "params")[["cPCA_params"]] <- internal_params

    # Return
    object
  }


# methods ----------------------------------------------------------------------

#' Apply contrastive PCA.
#'
#' @description
#' Applies the contrastive PCA algorithm given a specified alpha and a number of
#' contrastive principal components to extract.
#'
#' @param object The underlying class, see [bixverse::bulk_coexp()].
#' @param alpha Alpha parameter to use.
#' @param no_pcs Number of contrastive PCs to generate.
#'
#' @return `bulk_coexp` with additional data in the slots
#'
#' @references Abid, et al., Nature Communications, 2018
#'
#' @export
apply_contrastive_pca <- S7::new_generic(
  name = "apply_contrastive_pca",
  dispatch_args = "object",
  fun = function(object, alpha, no_pcs) {
    S7::S7_dispatch()
  }
)


#' @importFrom magrittr `%>%`
#' @importFrom magrittr `%$%`
#' @import data.table
#'
#' @export
#'
#' @method apply_contrastive_pca bulk_coexp
S7::method(apply_contrastive_pca, bulk_coexp) <-
  function(object, alpha, no_pcs) {
    # Checks
    checkmate::assertClass(object, "bixverse::bulk_coexp")
    checkmate::qassert(alpha, "N1")
    checkmate::qassert(no_pcs, "I1")
    detection_method <- S7::prop(object, "params")["detection_method"]
    checkmate::assertTRUE(detection_method == "cPCA")

    # Extract data
    target_covar <- S7::prop(object, "processed_data")[["target_covar"]]
    background_covar <- S7::prop(object, "processed_data")[["background_covar"]]
    target_mat <- S7::prop(object, "processed_data")[["target_mat"]]

    # Run cPCA
    c(factors, loadings) %<-% rs_contrastive_pca(
      target_covar = target_covar,
      background_covar = background_covar,
      target_mat = target_mat,
      alpha = alpha,
      n_pcs = no_pcs,
      return_loadings = TRUE
    )

    colnames(factors) <- colnames(loadings) <- sprintf("cPC_%i", seq(1:10))
    rownames(factors) <- rownames(target_mat)
    rownames(loadings) <- S7::prop(object, "params")[["cPCA_params"]][["intersecting_features"]]

    S7::prop(object, "params")[["cPCA_params"]]['final_alpha'] <- alpha
    S7::prop(object, "params")[["cPCA_params"]]['n_pcs'] <- no_pcs

    S7::prop(object, "outputs")[['cPCA_factors']] <- factors
    S7::prop(object, "outputs")[['cPCA_loadings']] <- loadings

    return(object)
  }

# plotting ---------------------------------------------------------------------

#' Plot various alphas for the contrastive PCA
#'
#' @description
#' This function will plot various alphas to highlight the most interesting
#' alpha parameters akin to the implementation of contrastive PCA in Python.
#'
#' @param object The underlying class, see [bixverse::bulk_coexp()]. You need
#' to apply [bixverse::contrastive_pca_processing()] to the function for this
#' method to work. Checkmate will raise errors otherwise.
#' @param label_column An optional sample label column. Needs to exist in the
#' meta_data of the `bulk_coexp` class.
#' @param min_alpha Minimum alpha to test.
#' @param max_alpha Maximum alpha to test.
#' @param n_alphas Number of alphas to test. The function will generate a series
#' of alphas from log(min_alpha) to log(max_alpha) to test out.
#' @param .verbose Controls verbosity of function.
#'
#' @return A ggplot showing the impact of various alpha parameters on the
#' samples in form of 2D plots.
#'
#' @references Abid, et al., Nature Communications, 2018
#'
#' @export
c_pca_plot_alphas <- S7::new_generic(
  name = "c_pca_plot_alphas",
  dispatch_args = "object",
  fun = function(object,
                 label_column = NULL,
                 min_alpha = .1,
                 max_alpha = 100,
                 n_alphas = 10L,
                 .verbose = TRUE) {
    S7::S7_dispatch()
  })



#' @export
#'
#' @importFrom magrittr `%>%`
#' @importFrom magrittr `%$%`
#' @importFrom zeallot `%<-%`
#' @import data.table
#' @import ggplot2
#'
#' @method c_pca_plot_alphas bulk_coexp
S7::method(c_pca_plot_alphas, bulk_coexp) <- function(object,
                                                      label_column = NULL,
                                                      min_alpha = .1,
                                                      max_alpha = 100,
                                                      n_alphas = 10L,
                                                      .verbose = TRUE) {
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  checkmate::qassert(label_column, c("S1", "0"))
  checkmate::qassert(min_alpha, "N1")
  checkmate::qassert(max_alpha, sprintf("N1(%f,]", min_alpha))
  checkmate::qassert(n_alphas, "I1")
  detection_method <- S7::prop(object, "params")["detection_method"]
  checkmate::assertTRUE(detection_method == "cPCA")

  # Get data
  target_covar <- S7::prop(object, "processed_data")[["target_covar"]]
  background_covar <- S7::prop(object, "processed_data")[["background_covar"]]
  target_mat <- S7::prop(object, "processed_data")[["target_mat"]]
  meta_data <- S7::prop(object, "meta_data")

  # Calculate the sequence of alphas
  alpha_seq <- c(0, exp(seq(
    log(min_alpha),
    log(max_alpha),
    length.out = n_alphas - 1
  )))

  # Calculate the contrastive PCs...
  cpcs_list <- purrr::map(alpha_seq, ~ {
    # Rust baby...
    c(factors, loadings) %<-% rs_contrastive_pca(
      target_covar = target_covar,
      background_covar = background_covar,
      target_mat = target_mat,
      alpha = .x,
      n_pcs = 2,
      return_loadings = FALSE
    )

    colnames(factors) <- c("cPC1", "cPC2")

    factors
  }) %>%
    `names<-`(paste0("Alpha: ", round(alpha_seq, digits = 2)))

  plot_df <- purrr::imap_dfr(
    cpcs_list,
    ~ {
      alpha <- .y
      cPCs <- .x
      res <- data.table::data.table(
        cPCs,
        alpha = as.character(alpha)
      )
    }
  )
  data.table::setDT(plot_df)

  if (!is.null(label_column) &&
    label_column %in% names(meta_data)) {
    if (.verbose) {
      message(
        "Found the ",
        label_column,
        " in the meta data. Adding labels to the graph."
      )
    }
    labels <- unlist(meta_data[, ..label_column])
    plot_df[, label := rep(as.character(labels), n_alphas)]
    add_labels <- T
  } else {
    if (.verbose) {
      message("Either no label column was provided or could not be found.
              No labels added to the graph")
    }
    add_labels <- F
  }

  plot_df[, alpha := factor(alpha, levels = paste0("Alpha: ", round(alpha_seq, digits = 2)))]

  plot <- ggplot(
    plot_df,
    aes(
      x = cPC1,
      y = cPC2
    )
  ) +
    facet_wrap(~alpha,
      scales = "free"
    ) +
    ggtitle("Alpha parameter iteration.",
      subtitle = "Impact on first 2 cPCs:"
    ) +
    theme_minimal() +
    theme(legend.position = "right")
  if (add_labels) {
    plot <- plot +
      geom_point(aes(col = label)) +
      scale_color_brewer(palette = "Spectral")
  } else {
    plot <- plot +
      geom_point()
  }

  plot
}
