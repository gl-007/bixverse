# methods ----

#' Calculate the Resnik and Lin semantic similarity for an ontology.
#'
#' @description This function calculates the semantic similarities based on
#' Resnik and Lin similarity for a given ontology. In this case, the ontology
#' is stored within an [bixverse::ontology()] class.
#'
#' @param object `ontology class`. See [bixverse::ontology()].
#' @param terms Character vector or NULL. The function will calculate the
#' the semantic similarities for these terms. If NULL, the semantic similarity
#' for everything will be calculated.
#' @param .verbose
#'
#' @return The class with added semantic similarities to the properties.
#'
#' @export
calculate_semantic_sim_onto <- S7::new_generic(
  name = "calculate_semantic_sim_onto",
  dispatch_args = "object",
  fun = function(object,
                 terms = NULL,
                 .verbose = TRUE) {
    S7::S7_dispatch()
  }
)

#' @export
#'
#' @import data.table
#' @importFrom magrittr `%>%`
#'
#' @method calculate_semantic_sim_onto ontology
S7::method(calculate_semantic_sim_onto, ontology) <-
  function(object,
           terms = NULL,
           .verbose = TRUE) {
    # Checks
    checkmate::assertClass(object, "bixverse::ontology")
    checkmate::qassert(terms, c("0", "S+"))

    if (is.null(terms)) {
      if (.verbose)
        message(
          paste0(
            "No terms provided. Function will calculate semantic ",
            "similarities for all terms."
          )
        )
      terms <- unique(unlist(S7::prop(object, "parent_child_dt")))
    }

    ancestor_list <- S7::prop(object, "ancestor_list")
    information_content_list <- S7::prop(object, "information_content_list")
    max_ic <- max(unlist(information_content_list))

    onto_similarities <- rs_onto_similarity(terms = terms,
                                            ancestor_list = ancestor_list,
                                            ic_list = information_content_list) %>%
      data.table::setDT() %>%
      .[, resnik_sim_norm := resnik_sim / max_ic]

    S7::prop(object, "semantic_similarities") <- onto_similarities

    return(object)
  }

# helpers ----

#' Calculate the Resnik and Lin semantic similarity
#'
#' @description This function calculates the semantic similarities based on
#' Resnik and Lin similarity for a given ontology.
#'
#' @param terms Vector of strings. The terms in the ontology you wish to screen.
#' @param ancestor_list List. Names being the term and the elements in the
#' list the names of the ancestors, see [bixverse::get_ontology_ancestors()].
#' @param ic_list List. The names being the term and the elements the
#' information content of this given term. Needs to be a single float! See
#' [bixverse::calculate_information_content()].
#' @param max_ic Double. The maximum information content observed in the data.
#' This value will be used to calculate the normalised Resnik similarity.
#'
#' @return data.table with the following columns:
#' \itemize{
#'  \item term1 - String, the first term.
#'  \item term2 - String, the second term.
#'  \item resnik_sim - Float, the unnormalised Resnik similarity.
#'  \item lin_sim - Float, the Lin similarity.
#'  \item resnik_sim_norm - Float, the normalised Resnik similarity (i.e.,
#'  divided by max information content observed.)
#' }
#'
#' @export
#'
#' @import data.table
calculate_semantic_sim <- function(terms, ancestor_list, ic_list, max_ic) {
  # Checks
  checkmate::qassert(terms, "S+")
  checkmate::assertList(ancestor_list, types = "character")
  checkmate::assert_named(ancestor_list)
  checkmate::assertList(ic_list, types =  "double")
  checkmate::assert_named(ic_list)
  checkmate::qassert(max_ic, "N1")

  onto_similarities <- rs_onto_similarity(terms = terms,
                                          ancestor_list = ancestor_list,
                                          ic_list = ic_list) %>%
    data.table::setDT() %>%
    .[, resnik_sim_norm := resnik_sim / max_ic]

  return(onto_similarities)
}


#' Return ancestor terms from an ontology
#'
#' @description this function will return all ancestor terms based on a provided
#' data.table with parent-child terms
#'
#' @param parent_child_dt data.table. The data.table with column parent and
#' child.
#'
#' @return A named list with ancestor terms as values
#'
#' @export
get_ontology_ancestors <- function(parent_child_dt) {
  checkmate::assertDataTable(parent_child_dt)
  checkmate::assert(all(c("parent", "child") %in% colnames(parent_child_dt)))

  # Deep copy to avoid side effects
  edge_df <- data.table::copy(parent_child_dt)
  data.table::setnames(edge_df, c("parent", "child"), c("to", "from"))

  graph <- igraph::graph_from_data_frame(edge_df[, c("from", "to")])
  ancestor_DT <- graph %>%
    igraph::ego(order = igraph::vcount(graph),
                mode = "out") %>%
    setNames(igraph::V(graph)$name) %>%
    Map(f = names) %>%
    stack() %>%
    rev() %>%
    setNames(names(edge_df)) %>%
    data.table::as.data.table() %>%
    .[from != to] %>%
    .[, lapply(.SD, as.character)]

  result <- split(ancestor_DT$from, ancestor_DT$to)

  return(result)
}


#' Calculate the information content for each ontology term
#'
#' @description this function will calculate the information content of each
#' provided term based on a list of ancestors, which is a named list of terms
#' with ancestor identifiers as their values. Can be calculated using
#' [bixverse::get_ancestors()]. The information content is calculated as
#' `-log2(number descendant/total terms in the ontology)`. More information
#' can be found [here](https://yulab-smu.top/biomedical-knowledge-mining-book/semantic-similarity-overview.html).
#'
#' @param ancestor_list List. Named list of terms with ancestor identifiers as
#' their values
#'
#' @return A named list of each term and their information content as values.
#'
#' @export
#'
#' @import data.table
#' @importFrom magrittr `%>%`
calculate_information_content <- function(ancestor_list) {
  checkmate::assertList(ancestor_list)

  total_terms = length(unique(c(
    unlist(ancestor_list), names(ancestor_list)
  )))
  ic = as.data.table(stack(ancestor_list)) %>%
    .[, `:=`(values = as.character(values), ind = as.character(ind))] %>%
    setnames(., c("values", "ind"), c("ancestor", "descendant")) %>%
    unique() %>%
    .[, nb_desc := .N, by = ancestor] %>%
    .[, ic := -log2(nb_desc / total_terms)] %>%
    setnames(., "ancestor", "id")

  result <- as.list(setNames(ic$ic, ic$id))

  return(result)
}
