# Gene Set Enrichment versions ----

#' Run gene ontology enrichment with elimination method.
#'
#' @description
#' This method takes the gene_ontology_data and a target gene set and performs
#' an GSE enrichment leveraging ontological information. It starts at the lowest
#' levels of the ontology and tests if there is significant enrichment for any
#' GO terms. If the threshold of the p-value is below the elimination threshold,
#' the genes from this term will be removed from all its ancestors. The function
#' then proceeds to the next level of the ontology and repeats the process.
#'
#' @param object The underlying class, see [bixverse::gene_ontology_data()].
#' @param target_genes String. The target genes you wish to apply the GSEA over.
#' @param minimum_overlap Integer. Threshold for the minimal overlap.
#' @param fdr_threshold Float. Threshold for maximum fdr to include in the
#' output.
#' @param elim_threshold Float. Threshold from which p-value onwards the
#' elimination on the ancestors shall be conducted.
#' @param min_genes Integer. Minimum number of genes that have to be included in
#' the gene ontology term. If NULL, it will default to the number of minimum
#' genes stored in `gene_ontology_data`.
#' @param .debug Boolean. Shall information from the Rust function be displayed.
#' For debugging purposes.
#'
#' @return data.table with enrichment results.
#'
#' @export
gse_go_elim_method <- S7::new_generic(
  name = "gse_go_elim_method",
  dispatch_args = "object",
  fun = function(object,
                 target_genes,
                 minimum_overlap = 3L,
                 fdr_threshold = 0.05,
                 elim_threshold = 0.05,
                 min_genes = NULL,
                 .debug = FALSE) {
    S7::S7_dispatch()
  }
)


#' @export
#'
#' @import data.table
#' @importFrom magrittr `%>%`
#'
#' @method gse_go_elim_method gene_ontology_data
S7::method(gse_go_elim_method, gene_ontology_data) <-
  function(object,
           target_genes,
           minimum_overlap = 3L,
           fdr_threshold = 0.05,
           elim_threshold = 0.05,
           min_genes = NULL,
           .debug = FALSE) {
    # Initial assignment
    `.` <- pvals <- fdr <- hits <- NULL
    # First check
    checkmate::assertClass(object, "bixverse::gene_ontology_data")
    checkmate::qassert(target_genes, "S+")
    checkmate::qassert(fdr_threshold, "R+[0,1]")
    checkmate::qassert(elim_threshold, "R+[0,1]")
    checkmate::qassert(minimum_overlap, "I1")
    checkmate::qassert(min_genes, c("0", "I1"))
    checkmate::qassert(.debug, "B1")
    # Extract relevant data from the S7 object
    if (is.null(min_genes)) {
      min_genes <- S7::prop(object, "min_genes")
    }
    levels <- names(S7::prop(object, "levels"))
    go_info <- S7::prop(object, "go_info")

    gene_universe_length <- length(unique(unlist(S7::prop(object, "go_to_genes"))))

    results_go <- rs_gse_geom_elim(
      target_genes = target_genes,
      levels = levels,
      go_obj = object,
      gene_universe_length = gene_universe_length,
      min_genes = min_genes,
      elim_threshold = elim_threshold,
      debug = .debug
    )

    results_go_dt <- data.table(do.call(cbind, results_go[-1])) %>%
      .[, `:=`(
        go_id = results_go$go_ids,
        fdr = p.adjust(pvals, method = "BH")
      )] %>%
      data.table::setcolorder(
        .,
        c(
          "go_id",
          "odds_ratios",
          "pvals",
          "fdr",
          "hits",
          "gene_set_lengths"
        )
      )

    results_go_dt <- merge(go_info, results_go_dt, by = "go_id") %>%
      data.table::setorder(., pvals) %>%
      .[(fdr <= fdr_threshold) & (hits >= minimum_overlap)]

    results_go_dt
  }


#' Run gene ontology enrichment with elimination method over a list.
#'
#' @description
#' This method takes the gene_ontology_data and a list of target gene sets and
#' performs an GSE enrichment leveraging ontological information. It starts at
#' the lowest levels of the ontology and tests if there is significant
#' enrichment for any GO terms. If the threshold of the p-value is below the
#' elimination threshold, the genes from this term will be removed from all its
#' ancestors. The function then proceeds to the next level of the ontology and
#' repeats the process. The class will leverage Rust threading to parallelise
#' the process.
#'
#' @param object The underlying class, see [bixverse::gene_ontology_data()].
#' @param target_gene_list List. The target genes list you wish to apply the
#' gene set enrichment analysis over.
#' @param minimum_overlap Integer. Threshold for the minimal overlap.
#' @param fdr_threshold Float. Threshold for maximum fdr to include in the
#' output.
#' @param elim_threshold Float. Threshold from which p-value onwards the
#' elimination on the ancestors shall be conducted.
#' @param min_genes Integer. Minimum number of genes that have to be included in
#' the gene ontology term. If NULL, it will default to the number of minimum
#' genes stored in `gene_ontology_data`.
#' @param .debug Boolean. Shall information from the Rust function be displayed.
#' For debugging purposes. Warning: should you run this command over a large
#' list, you will have a large print output!
#'
#' @return data.table with enrichment results.
#'
#' @export
gse_go_elim_method_list <- S7::new_generic(
  name = "gse_go_elim_method_list",
  dispatch_args = "object",
  fun = function(object,
                 target_gene_list,
                 minimum_overlap = 3L,
                 fdr_threshold = 0.05,
                 elim_threshold = 0.05,
                 min_genes = NULL,
                 .debug = FALSE) {
    S7::S7_dispatch()
  }
)


#' @export
#'
#' @import data.table
#' @importFrom magrittr `%>%`
#'
#' @method gse_go_elim_method_list gene_ontology_data
S7::method(gse_go_elim_method_list, gene_ontology_data) <-
  function(object,
           target_gene_list,
           minimum_overlap = 3L,
           fdr_threshold = 0.05,
           elim_threshold = 0.05,
           min_genes = NULL,
           .debug = FALSE) {
    # Binding checks
    `.` <- pvals <- fdr <- hits <- target_set_name <- NULL
    # First check
    checkmate::assertClass(object, "bixverse::gene_ontology_data")
    checkmate::assertList(target_gene_list, types = "character")
    checkmate::qassert(fdr_threshold, "R+[0,1]")
    checkmate::qassert(elim_threshold, "R+[0,1]")
    checkmate::qassert(minimum_overlap, "I1")
    checkmate::qassert(min_genes, c("0", "I1"))
    checkmate::qassert(.debug, "B1")
    # Extract relevant data from the S7 object
    if (is.null(min_genes)) {
      min_genes <- S7::prop(object, "min_genes")
    }
    levels <- names(S7::prop(object, "levels"))
    go_info <- S7::prop(object, "go_info")

    gene_universe_length <- length(unique(unlist(S7::prop(object, "go_to_genes"))))

    results_go <- rs_gse_geom_elim_list(
      target_genes_list = target_gene_list,
      levels = levels,
      go_obj = object,
      gene_universe_length = gene_universe_length,
      min_genes = min_genes,
      elim_threshold = elim_threshold,
      debug = .debug
    )

    target_set_names <- purrr::map2(
      names(target_gene_list), results_go$no_test,
      ~ {
        rep(.x, each = .y)
      }
    )

    target_set_names <- do.call(c, target_set_names)

    cols_to_select <- c("pvals", "odds_ratios", "hits", "gene_set_lengths")

    results_go_dt <- data.table(do.call(cbind, results_go[cols_to_select])) %>%
      .[, `:=`(
        go_id = results_go$go_ids,
        target_set_name = target_set_names
      )] %>%
      .[, fdr := p.adjust(pvals, method = "BH"), by = target_set_name]

    results_go_dt <- merge(results_go_dt, go_info, by = "go_id") %>%
      data.table::setorder(., pvals) %>%
      data.table::setcolorder(
        .,
        c(
          "target_set_name",
          "go_name",
          "go_id",
          "odds_ratios",
          "pvals",
          "fdr",
          "hits",
          "gene_set_lengths"
        )
      ) %>%
      .[(fdr <= fdr_threshold) & (hits >= minimum_overlap)]

    results_go_dt
  }
