library(magrittr)
library(data.table)

devtools::document()
rextendr::document()
devtools::load_all()
devtools::check()

# Community detections ----

edge_data = arrow::read_parquet("~/Desktop/string_clean.parquet") %>%
  as.data.table()
#
# edge_data_clean = edge_data %>%
#   .[interactionScore >= 0.85] %>%
#   dplyr::select(from = `:START_ID`, to = `:END_ID`)

devtools::load_all()

test_class = network_diffusions(edge_data, weighted = FALSE, directed = FALSE)



set.seed(123)
genes = sample(igraph::V(test_class@graph)$name, 10)
genes.2 = sample(igraph::V(test_class@graph)$name, 10)
genes.3 = sample(igraph::V(test_class@graph)$name, 25)
diffusion_vector = rep(1, 10) %>% `names<-`(genes)
diffusion_vector.2 = rep(1, 10) %>% `names<-`(genes.2)

test_class <- diffuse_seed_nodes(test_class, diffusion_vector, 'max')

get_params(test_class, TRUE, TRUE)

devtools::load_all()

test_class <- tied_diffusion(
  object = test_class,
  diffusion_vector_1 = diffusion_vector,
  diffusion_vector_2 = diffusion_vector.2,
  summarisation = 'max',
  score_aggregation = 'min'
)

get_results(test_class)

get_params(test_class, TRUE, TRUE)

test_class <- community_detection(
  test_class,
  community_params = params_community_detection(
    min_seed_nodes = 1L
  ),
  diffusion_threshold = .5
)

x <- get_results(test_class)

x

# RBH graph ----

protein_coding_genes <- data.table::fread("~/Desktop/protein_coding_genes.csv")

universe <- protein_coding_genes$id[1:500]

sets_per_origin <- 100
gene_sets_no <- 100

gene_sets_no <- sets_per_origin * length(LETTERS)

seed <- 123
set.seed(seed)

gene_sets <- purrr::map(1:gene_sets_no, ~ {
  set.seed(seed + .x + 1)
  size <- sample(20:100, 1)
  sample(universe, size, replace = FALSE)
})

names(gene_sets) <- purrr::map_chr(1:gene_sets_no, ~ {
  set.seed(seed + .x + 1)
  paste(sample(LETTERS, 5), collapse = "")
})

gene_sets_df <- purrr::imap(gene_sets, ~{
  data.table::data.table(
    genes = .x,
    name = .y
  )
})

origins <- rep(LETTERS, each = sets_per_origin)

gene_sets_df <- purrr::map2(gene_sets_df, origins, ~{
  .x[, origin := .y]
}) %>% rbindlist

module_df = gene_sets_df

rbh_class = rbh_graph(
  gene_sets_df,
  dataset_col = 'origin',
  module_col = 'name',
  value_col = 'genes'
)


rbh_class = generate_rbh_graph(
  rbh_class,
  minimum_similarity = 0.2,
  overlap_coefficient = T,
  .debug = FALSE
)

rbh_class = find_rbh_communities(rbh_class)

plot_resolution_res(rbh_class)
