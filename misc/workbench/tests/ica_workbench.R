# Check ICA stuff ---

## Test real data ----

library(devtools)
library(ggplot2)
library(magrittr)
library(zeallot)

devtools::document()

gtex_brain <- recount3::create_rse_manual(
  project = "BRAIN",
  project_home = "data_sources/gtex",
  organism = "human",
  annotation = "gencode_v29",
  type = "gene"
)

coldata <- SummarizedExperiment::colData(gtex_brain) |> as.data.frame()
rowdata <- SummarizedExperiment::rowData(gtex_brain) |> as.data.frame()

d <- edgeR::DGEList(SummarizedExperiment::assay(gtex_brain))
d <- edgeR::calcNormFactors(d, method = 'upperquartile')
to_keep <- suppressWarnings(edgeR::filterByExpr(d))
d <- d[to_keep, ]
d <- edgeR::cpm(d, log = TRUE)

d <- as.matrix(d)

rextendr::document()
devtools::document()
devtools::load_all()
# devtools::check()

new_meta_data <- data.table::data.table(
  sample_id = rownames(coldata),
  case_control = 'case',
  gtex_subgrp = coldata$gtex.smtsd
)

samples_to_keep <- new_meta_data[gtex_subgrp == "Brain - Putamen (basal ganglia)", sample_id]
data_1 = t(d)[samples_to_keep, ]
meta_data = new_meta_data[gtex_subgrp == "Brain - Putamen (basal ganglia)"]

ica_test = bulk_coexp(raw_data = data_1, meta_data = meta_data)
ica_test = preprocess_bulk_coexp(ica_test, mad_threshold = 1)
ica_test = ica_processing(ica_test)

ica_test = ica_evaluate_comp(
  ica_test,
  ica_type = 'logcosh',
  ncomp_params = params_ica_ncomp(max_no_comp = 75L)
)

plot_ica_stability_individual(ica_test)

plot_ica_stability_summarised(ica_test)

devtools::load_all()

ica_test <- ica_stabilised_results(ica_test, no_comp = 40L, ica_type = "logcosh")

outputs <- get_results(ica_test)

outputs$S[1:5, 1:5]

outputs$A[1:5, 1:5]

outputs$ica_meta

get_results(ica_test)

# Write a final component function

?rs_prepare_whitening

?ica_evaluate_comp

object = ica_test
no_comp = 50L
ica_type = "logcosh"
iter_params = params_ica_randomisation()
ica_params = params_ica_general()
random_seed = 42L
consistent_sign = TRUE
.verbose = TRUE

X <- S7::prop(object, "processed_data")[['processed_data']]
X1 <- S7::prop(object, "processed_data")[["X1"]]
K <- S7::prop(object, "processed_data")[["K"]]

do.call(c, list(1, 2,3))

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
  return_centrotype = TRUE
)

colnames(centrotype) <- sprintf("IC_%i", 1:no_comp)
rownames(centrotype) <- colnames(X_raw)

centrotype[1:5, 1:5]



centrotype <- apply(centrotype, 2, .flip_ica_loading_signs)

S <- t(centrotype)

A <- t(X1) %*% MASS::ginv(S)
rownames(A) <- rownames(X)
colnames(A) <- rownames(S)

ica_meta <- list(
  component = sprintf("IC_%i", 1:no_comp),
  stability = stability_scores
) %>% data.table::setDT()


# Is my ICA implementation correct ... ? ---------------------------------------

## FastICA ---------------------------------------------------------------------

S <- cbind(sin((1:1000) / 20), rep((((
  1:200
) - 100) / 100), 5))
A <- matrix(c(0.291, 0.6557, -0.5439, 0.5572), 2, 2)
X <- S %*% A

a <- fastICA::fastICA(
  X,
  2,
  alg.typ = "parallel",
  fun = "logcosh",
  alpha = 1,
  method = "R",
  row.norm = FALSE,
  maxit = 200,
  tol = 0.0001,
  verbose = TRUE
)

par(mfcol = c(2, 3))
plot(
  1:1000,
  S[, 1],
  type = "l",
  main = "Original Signals",
  xlab = "",
  ylab = ""
)
plot(1:1000,
     S[, 2],
     type = "l",
     xlab = "",
     ylab = "")
plot(
  1:1000,
  X[, 1],
  type = "l",
  main = "Mixed Signals",
  xlab = "",
  ylab = ""
)
plot(1:1000,
     X[, 2],
     type = "l",
     xlab = "",
     ylab = "")
plot(
  1:1000,
  a$S[, 1],
  type = "l",
  main = "ICA source estimates",
  xlab = "",
  ylab = ""
)
plot(1:1000,
     a$S[, 2],
     type = "l",
     xlab = "",
     ylab = "")

## Rust ------------------------------------------------------------------------

S <- cbind(sin((1:1000) / 20), rep((((
  1:200
) - 100) / 100), 5))
A <- matrix(c(0.291, 0.6557, -0.5439, 0.5572), 2, 2)
X <- S %*% A

c(X_norm, K) %<-% rs_prepare_whitening(X, TRUE, 123L, NULL, NULL, NULL)

rextendr::document()
rextendr::clean()
devtools::load_all()

?fast_ica_rust

ica_res_rs <- fast_ica_rust(
  X_norm,
  K,
  n_icas = 2L,
  ica_fun = "logcosh",
  seed = 42L
)

par(mfcol = c(2, 3))
plot(
  1:1000,
  S[, 1],
  type = "l",
  main = "Original Signals",
  xlab = "",
  ylab = ""
)
plot(1:1000,
     S[, 2],
     type = "l",
     xlab = "",
     ylab = "")
plot(
  1:1000,
  X[, 1],
  type = "l",
  main = "Mixed Signals",
  xlab = "",
  ylab = ""
)
plot(1:1000,
     X[, 2],
     type = "l",
     xlab = "",
     ylab = "")
plot(
  1:1000,
  ica_res_rs$S[1, ],
  type = "l",
  main = "ICA source estimates",
  xlab = "",
  ylab = ""
)
plot(1:1000,
     ica_res_rs$S[2, ],
     type = "l",
     xlab = "",
     ylab = "")

