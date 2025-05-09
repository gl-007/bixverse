# Load in and extend the Rust documentation...
rextendr::document()
devtools::document()
devtools::load_all()
devtools::install()
# devtools::check()

# devtools::install()

system.file("extdata", package = "bixverse") |> list.files()

library(magrittr)

# General hypergeom tests ------------------------------------------------------

go_data_dt <- get_go_human_data()

go_data_s7 <- gene_ontology_data(go_data_dt, min_genes = 3L)

protein_coding_genes <- unique(unlist(go_data_s7@go_to_genes))

seed <- 123
set.seed(seed)

universe <- protein_coding_genes
gene_sets_no <- 5000
target_gene_sets_no <- 150

gene_sets <- purrr::map(
  1:gene_sets_no,
  ~ {
    set.seed(seed + .x + 1)
    size <- sample(20:100, 1)
    sample(universe, size, replace = FALSE)
  }
)

names(gene_sets) <- purrr::map_chr(
  1:gene_sets_no,
  ~ {
    set.seed(seed + .x + 1)
    paste(sample(LETTERS, 3), collapse = "")
  }
)

target_gene_sets <- purrr::map(
  1:target_gene_sets_no,
  ~ {
    set.seed(.x * seed)
    size <- sample(50:100, 1)
    sample(universe, size, replace = FALSE)
  }
)

names(target_gene_sets) <- purrr::map_chr(
  1:target_gene_sets_no,
  ~ {
    set.seed(seed + .x + 1)
    paste(sample(letters, 3), collapse = "")
  }
)

target_genes <- target_gene_sets[[1]]

tictoc::tic()
t1 <- gse_hypergeometric(
  target_genes = target_genes,
  gene_set_list = gene_sets,
  gene_universe = universe,
  minimum_overlap = 0L,
  threshold = 1
)
tictoc::toc()

tictoc::tic()
t2 <- gse_hypergeometric_list(
  target_gene_sets,
  gene_set_list = gene_sets,
  minimum_overlap = 0L,
  threshold = 1
)
tictoc::toc()

tictoc::tic()
t3 <- gse_go_elim_method(
  go_data_s7,
  target_genes,
  minimum_overlap = 0L,
  fdr_threshold = 1
)
tictoc::toc()

length(target_gene_sets)

tictoc::tic()
t4 <- gse_go_elim_method_list(
  go_data_s7,
  target_gene_sets,
  minimum_overlap = 0L,
  fdr_threshold = 1
)
tictoc::toc()
