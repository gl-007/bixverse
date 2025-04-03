# checks -----------------------------------------------------------------------

## correlation params ----------------------------------------------------------

#' Check correlation graph parameters
#'
#' @description Checkmate extension for checking the graph parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
checkCorGraphParams <- function(x) {
  # Checkmate extension
  res <- checkmate::checkList(x)
  if (!isTRUE(res))
    return(res)
  res <- checkmate::checkNames(names(x),
                               must.include = c("epsilon", "min_cor", "fdr_threshold", "verbose"))
  if (!isTRUE(res))
    return(res)
  rules <- list(
    "epsilon" = "R1",
    "min_cor" = "R1[0, 1]",
    "fdr_threshold" = "R1[0, 1]",
    "verbose" = "B1"
  )
  res <- purrr::imap_lgl(x, \(x, name) {
    checkmate::qtest(x, rules[[name]])
  })
  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        "The following element `%s` in graph params does not conform the expected format. \
        min_cor and fdr_threshold need to be between 0 and 1, epsilon a double. and .verbose \
        a boolean.",
        broken_elem
      )
    )
  }
  return(TRUE)
}

#' Check resolution graph parameters
#'
#' @description Checkmate extension for checking the resolution parameters for
#' community detection with Leiden.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
checkGraphResParams <- function(x) {
  # Checkmate extension
  res <- checkmate::checkList(x)
  if (!isTRUE(res))
    return(res)
  res <- checkmate::checkNames(names(x), must.include = c("min_res", "max_res", "number_res"))
  if (!isTRUE(res))
    return(res)
  rules <- list("min_res" = "R1",
                "max_res" = "R1",
                "number_res" = "I1")
  res <- purrr::imap_lgl(x, \(x, name) {
    checkmate::qtest(x, rules[[name]])
  })
  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        "The following element `%s` in resolution params does not conform the expected format. \
        min_res and max_res need to be doubles and number res needs to be an integer",
        broken_elem
      )
    )
  }
  return(TRUE)
}

## ica params ------------------------------------------------------------------

#' Check ICA parameters
#'
#' @description Checkmate extension for checking the ICA parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
checkIcaParams <- function(x) {
  # Checkmate extension
  res <- checkmate::checkList(x)
  if (!isTRUE(res))
    return(res)
  res <- checkmate::checkNames(names(x),
                               must.include = c("maxit", "alpha", "max_tol", "verbose"))
  if (!isTRUE(res))
    return(res)
  rules <- list(
    "maxit" = "I1",
    "alpha" = "R1[1, 2]",
    "max_tol" = "R1(0, 1)",
    "verbose" = "B1"
  )
  res <- purrr::imap_lgl(x, \(x, name) {
    checkmate::qtest(x, rules[[name]])
  })
  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        "The following element `%s` in resolution params does not conform the expected format. \
        maxit needs to be an integer, alpha between 1 and 2, 0 < max_tol < 1 and verbose a \
        boolean",
        broken_elem
      )
    )
  }
  return(TRUE)
}


#' Check ICA no of component parameters
#'
#' @description Checkmate extension for checking the ICA number of component
#' parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
checkIcaNcomps <- function(x) {
  # Checkmate extension
  res <- checkmate::checkList(x)
  if (!isTRUE(res))
    return(res)
  res <- checkmate::checkNames(names(x), must.include = c("max_no_comp", "steps", "custom_seq"))
  if (!isTRUE(res))
    return(res)
  rules <- list(
    "max_no_comp" = "I1",
    "steps" = "I1",
    "custom_seq" = c("0", "I+")
  )
  res <- purrr::imap_lgl(x, \(x, name) {
    checkmate::qtest(x, rules[[name]])
  })
  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        "The following element `%s` in resolution params does not conform the expected format.
        max_no_comp and steps need to be an integer, and custom sequence either NULL or a vector \
        of integers",
        broken_elem
      )
    )
  }
  return(TRUE)
}

#' Check ICA randomisation parameters
#'
#' @description Checkmate extension for checking the ICA randomisation parameters
#' for a version of stabilised ICA.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
checkIcaIterParams <- function(x) {
  # Checkmate extension
  res <- checkmate::checkList(x)
  if (!isTRUE(res))
    return(res)
  res <- checkmate::checkNames(names(x),
                               must.include = c("cross_validate", "random_init", "folds"))
  if (!isTRUE(res))
    return(res)
  rules <- list(
    "cross_validate" = "B1",
    "random_init" = "I1",
    "folds" = "I1"
  )
  res <- purrr::imap_lgl(x, \(x, name) {
    checkmate::qtest(x, rules[[name]])
  })
  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        "The following element `%s` in resolution params does not conform the expected format. \
        random_init and steps folds to be an integer, and cross_validate a boolean",
        broken_elem
      )
    )
  }
  return(TRUE)
}

## community detections --------------------------------------------------------

#' Check community detection parameters
#'
#' @description Checkmate extension for checking the community detection
#' parameters for identifying genetically privileged communities.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
checkCommunityParams <- function(x) {
  # Checkmate extension
  res <- checkmate::checkList(x)
  if (!isTRUE(res))
    return(res)
  res <- checkmate::checkNames(
    names(x),
    must.include = c("max_nodes", "min_nodes", "min_seed_nodes", "initial_res")
  )
  if (!isTRUE(res))
    return(res)
  rules <- list(
    "max_nodes" = sprintf("I1[%i,)", x$min_nodes),
    "min_nodes" = "I1",
    "min_seed_nodes" = "I1",
    "initial_res" = "N1"
  )
  res <- purrr::imap_lgl(x, \(x, name) {
    checkmate::qtest(x, rules[[name]])
  })
  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        "The following element `%s` in community params does not conform the expected format. \
        min_nodes, max_nodes and min_seed_genes need to be integers (with max_nodes > min_nodes) \
        and initial resolution a double.",
        broken_elem
      )
    )
  }
  return(TRUE)
}

# asserts ----------------------------------------------------------------------

## correlation params ----------------------------------------------------------

#' Assert correlation graph parameters
#'
#' @description Checkmate extension for asserting the graph parameters.
#'
#' @inheritParams checkCorGraphParams
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
assertCorGraphParams <- checkmate::makeAssertionFunction(checkCorGraphParams)


#' Assert resolution graph parameters
#'
#' @description Checkmate extension for asserting the resolution parameters for
#' community detection with Leiden.
#'
#' @inheritParams checkGraphResParams
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
assertGraphResParams <- checkmate::makeAssertionFunction(checkGraphResParams)

## ica params ------------------------------------------------------------------

#' Assert ICA parameters
#'
#' @description Checkmate extension for asserting the ICA parameters
#'
#' @inheritParams checkIcaParams
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
assertIcaParams <- checkmate::makeAssertionFunction(checkIcaParams)

#' Assert ICA no of component parameters
#'
#' @description Checkmate extension for asserting the ICA number of component
#' parameters.
#'
#' @inheritParams checkIcaNcomps
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
assertIcaNcomps <- checkmate::makeAssertionFunction(checkIcaNcomps)

#' Assert ICA randomisation parameters
#'
#' @description Checkmate extension for asserting the ICA randomisation
#' parameters for a version of stabilised ICA.
#'
#' @inheritParams checkIcaIterParams
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
assertIcaIterParams <- checkmate::makeAssertionFunction(checkIcaIterParams)

## community detections --------------------------------------------------------

#' Assert community detection parameter
#'
#' @description Checkmate extension for asserting the community detection
#' parameters for identifying genetically privileged communities.
#'
#' @inheritParams checkCommunityParams
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
assertCommunityParams <- checkmate::makeAssertionFunction(checkCommunityParams)

