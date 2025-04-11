
library(Matrix)
library(tictoc)
library(igraph)
library(reshape2)
library(polars)
library(data.table)
library(dplyr)
library(igraph)
library(purrr)
library(furrr)
library(bixverse)
library(rlang)

library(bixverse)
library(polars)
library(data.table)
data_path = "/Users/liesbeth/Datascience/data/processed/OpenTargets_platform"

ontology_df = pl$read_parquet(file.path(data_path,"OT_edges_disease_hierarchy.parquet"))$
  filter(pl$col("source:string") == "efo")$
  rename(":END_ID" = "to",
         ":START_ID" = "from",
         "relation_type:string" = "relation")$
  select("from", "to", "relation")

parent_child = as.data.table(ontology_df$
                               filter(pl$col("relation") == "parent_of")$
                               rename("from" = "parent", "to" = "child")$
                               select("parent", "child"))

tic()
ontology_obj = bixverse::ontology(parent_child)
ontology_obj = calculate_semantic_sim_onto(ontology_obj)
similarities = get_semantic_similarities(ontology_obj)
toc()




ancestor_list = get_ancestors(ontology)
ic_list = calculate_information_content(ancestor_list)
terms = names(ancestor_list)
max_ic = -log2(1/length(ancestor_list))

sim_matrix <- get_similarity_matrix(terms, ancestor_list, ic_list, similarity_type = c("resnik", "lin"), return_matrix = T)

## what happens if mac_ic != 1 to the lin similarity? should we not give a third matrix back, resnik, lin and norm_resnik?
## the message needs to be changed: "Generating the full matrix format of the correlation matrix."



## testing ontologysimilarity package


library(ontologyIndex)
library(ontologySimilarity)
data(hpo)
set.seed(1)


tic()
information_content <- descendants_IC(hpo)
sim_mat <-
  get_sim_grid(
    ontology = hpo,
    information_content = information_content,
    term_sim_method = "resnik",
    term_sets = as.list(hpo$id)
  )
toc()

library(bixverse)
library(polars)
terms = hpo$id
ancestor_list = hpo$ancestors
ontology_ancestors = pl$DataFrame(stack(ancestor_list) %>%
  mutate_if(is.factor, as.character) %>%
  rename(from = values,
         to = ind))
information_content = calculate_information_content(ontology_ancestors)
ic_list = split(as.data.table(information_content)$ic, as.data.table(information_content)$id)
tic()
bxv_sim <- rs_onto_similarity_both(terms = terms,
                                   ancestor_list = ancestor_list,
                                   ic_list = ic_list,
                                   max_ic = max(unlist(ic_list)))

bxv_sim_mat <- bixverse:::upper_triangular_cor_mat$new(
  cor_coef = bxv_sim$resnik_sim,
  features = bxv_sim$terms,
  shift = 1L
)$get_cor_matrix()
toc()



#
# ## some functions for the R function
# get_index <- function(name, graph){
#   match(name, V(graph)$name)
# }
# LCA = function(distance, n1, n2){
#   d = rowSums(distance[, c(n1, n2)])
#   d = d[!is.infinite(d)]
#   if(length(d) == 0){
#     NA
#   }else{
#     names(d[which.min(d)])
#   }
# }
# get_distances <- function(ontology){
#   graph = igraph::graph_from_data_frame(ontology )
#   distance = igraph::distances(graph, V(graph), mode="out")
# }
# get_similarity <- function(information_content, ontology, term1, term2){
#   checkmate::assert(all(c(term1, term2) %in% information_content$select("id")$to_series()$to_list()))
#   ## calculate LCA
#   distance <- get_distances(ontology)
#   lca = LCA(distance, n1, n2)
#
#   resnik = information_content$filter(pl$col("id") == names(lca))$select("ic", "norm_ic")$to_data_frame()
#   lin = (2*resnik$ic)/(sum(information_content$filter(pl$col("id")$is_in(c(term1, term2)))$select("ic")$to_data_frame()))
#
#   return(data.table(term1 = term1, term2 = term2, resnik = resnik, lin.ic = lin))
# }
# # run it in R only
# ontology = ontology$filter(pl$col("relation") == "parent_of") %>% as.data.table()
# nodes = nodes %>% as.data.table()
#
# tic()
# d = get_distances(ontology)
# md = reshape2::melt(d)
#
# seq.int <- seq(1, nrow(md), 2000000)
# results = list()
# n = 1
# for(i in 1:length(seq.int)){
#   start = seq.int[i]
#   end = seq.int[i+1]-1
#   if(is.na(end)){
#     end = nrow(md)
#   }
#   tmp = md[start:end,] %>%
#     as.data.table() %>%
#     .[!is.infinite(value) & Var1 != Var2]
#   results[[n]] <- tmp
#   n = n + 1
# }
# results <- do.call(rbind, results)
#
# ## so now we have all disease pairs that connect somewhere in the tree
# plan(multisession, workers = 5)
# options(future.globals.maxSize = 6 * 1024^3)
#
# lca_all = future_map2(results$Var1, results$Var2,
#                       ~{
#                         LCA(d, .x, .y)
#                       })
# results <- results[, lca := unlist(lca_all)]
#
# results_pl <- pl$DataFrame(results)$
#   with_columns(pl$col("Var1")$cast(pl$String)$alias("node1"), pl$col("Var2")$cast(pl$String)$alias("node2"))$
#   join(information_content$select("id", "ic", "norm_ic")$rename("ic" = "resnik", "norm_ic" = "norm_resnik"),
#        left_on = "lca", right_on = "id")$
#   join(information_content$select("id", "ic")$rename("ic" = "ic_var1"), left_on = "node1", right_on = "id")$
#   join(information_content$select("id", "ic")$rename("ic" = "ic_var2"), left_on = "node2", right_on = "id")$
#   with_columns(((2*pl$col("resnik"))$div(pl$col("ic_var1") + pl$col("ic_var2")))$alias("lin"))$
#   select("node1", "node2", "lca", "resnik", "norm_resnik", "lin")
# toc()
#
#
# n1 = "MONDO_0008648"
# n2 = "Orphanet_7"
# lca = LCA(d, n1, n2)
# lca
# information_content$filter(pl$col("id") == lca)
