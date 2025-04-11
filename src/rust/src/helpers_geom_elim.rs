use extendr_api::prelude::*;
use std::collections::{HashMap, HashSet};

use crate::helpers_hypergeom::*;
use crate::utils_r_rust::{r_list_to_hashmap, r_list_to_hashmap_set};

///////////////////////
// Types & Structure //
///////////////////////

/// Type alias for the go identifier to gene Hashmap
type GeneMap = HashMap<String, HashSet<String>>;

/// Type alias for the ancestor to go identifier HashMap
type AncestorMap = HashMap<String, Vec<String>>;

/// Type alias for the ontology level to go identifier HashMap
type LevelMap = HashMap<String, Vec<String>>;

/// Return structure of the `process_ontology_level()` ontology function.
pub struct GoElimLevelResults {
    pub go_ids: Vec<String>,
    pub pvals: Vec<f64>,
    pub odds_ratios: Vec<f64>,
    pub hits: Vec<u64>,
    pub gene_set_lengths: Vec<u64>,
}

/// Structure that contains the gene ontology and key functions to do apply the
/// elimination method.
pub struct GeneOntology {
    pub go_to_gene: GeneMap,
    pub ancestors: AncestorMap,
    pub levels: LevelMap,
}

impl GeneOntology {
    /// Returns the ancestors of a given gene ontology term identifier
    pub fn get_ancestors(&self, id: &String) -> Option<&Vec<String>> {
        self.ancestors.get(id)
    }

    /// Returns the gene ontology term identifiers for a given level of the
    /// ontology.
    pub fn get_level_ids(&self, id: &String) -> Option<&Vec<String>> {
        self.levels.get(id)
    }

    /// Remove genes from defined sets of genes
    pub fn remove_genes(&mut self, ids: &[String], genes_to_remove: &HashSet<String>) {
        for id in ids.iter() {
            if let Some(gene_set) = self.go_to_gene.get_mut(id) {
                gene_set.retain(|gene| !genes_to_remove.contains(gene));
            }
        }
    }

    /// Get the genes based on an array of Strings.
    pub fn get_genes_list(&self, ids: &[String]) -> (Vec<String>, Vec<&HashSet<String>>) {
        let keys: HashSet<String> = self.go_to_gene.clone().into_keys().collect();

        let id_keys: HashSet<_> = ids.iter().cloned().collect();

        let ids_final: Vec<String> = id_keys.intersection(&keys).cloned().collect();

        let gene_sets: Vec<_> = ids_final
            .iter()
            .filter_map(|s| self.go_to_gene.get(s))
            .collect();

        (ids_final, gene_sets)
    }

    /// Get the genes for one specific ID
    pub fn get_genes(&self, id: &String) -> Option<&HashSet<String>> {
        self.go_to_gene.get(id)
    }
}

///////////////
// Functions //
///////////////

/// Take the S7 go_data_class and return the necessary Rust types for further
/// processing.
pub fn prepare_go_data(go_obj: Robj) -> extendr_api::Result<(GeneMap, AncestorMap, LevelMap)> {
    let go_to_genes = go_obj.get_attrib("go_to_genes").unwrap().as_list().unwrap();
    let ancestors = go_obj.get_attrib("ancestry").unwrap().as_list().unwrap();
    let levels = go_obj.get_attrib("levels").unwrap().as_list().unwrap();

    let go_to_genes = r_list_to_hashmap_set(go_to_genes)?;
    let ancestors = r_list_to_hashmap(ancestors)?;
    let levels = r_list_to_hashmap(levels)?;

    Ok((go_to_genes, ancestors, levels))
}

/// Process a given ontology level
pub fn process_ontology_level(
    target_genes: &[String],
    level: &String,
    go_obj: &mut GeneOntology,
    min_genes: i64,
    gene_universe_length: u64,
    elim_threshold: f64,
    debug: bool,
) -> GoElimLevelResults {
    // Get the identfiers of that level and clean everything up
    let go_ids = go_obj.get_level_ids(level);
    let go_ids_final: &Vec<String> = go_ids.unwrap();
    let level_data = go_obj.get_genes_list(go_ids_final);
    let mut go_identifiers = level_data.0;
    let mut go_gene_sets = level_data.1;

    // Remove gene ontology terms that do not have a minimum of genes
    let filtered_sets: Vec<(String, &HashSet<String>)> = go_identifiers
        .iter()
        .zip(go_gene_sets.iter())
        .filter(|(_, set)| set.len() as i64 >= min_genes)
        .map(|(string, set)| (string.clone(), *set))
        .collect();

    go_identifiers.clear();
    go_gene_sets.clear();

    for (string, set) in filtered_sets {
        go_identifiers.push(string);
        go_gene_sets.push(set);
    }

    let trials = target_genes.iter().collect::<Vec<_>>().len() as u64;
    let gene_set_lengths = go_gene_sets
        .clone()
        .into_iter()
        .map(|s| s.len() as u64)
        .collect::<Vec<u64>>();
    let hits = count_hits_hash(go_gene_sets, target_genes);

    // Calculate p-values and odds ratios
    let pvals: Vec<f64> = hits
        .iter()
        .zip(gene_set_lengths.iter())
        .map(|(hit, gene_set_length)| {
            hypergeom_pval(
                *hit,
                *gene_set_length,
                gene_universe_length - *gene_set_length,
                trials,
            )
        })
        .collect();
    let odds_ratios: Vec<f64> = hits
        .iter()
        .zip(gene_set_lengths.iter())
        .map(|(hit, gene_set_length)| {
            hypergeom_odds_ratio(
                *hit,
                *gene_set_length - *hit,
                trials - *hit,
                gene_universe_length - *gene_set_length - trials + *hit,
            )
        })
        .collect();

    // Identify the GO terms were to apply the elimination on (if any)
    let go_to_remove: Vec<String> = go_identifiers
        .iter()
        .zip(pvals.iter())
        .filter(|(_, pval)| pval <= &&elim_threshold)
        .map(|(string, _)| string.clone())
        .collect();

    if debug {
        let no_terms = go_to_remove.len();
        println!(
            "At level {} a total of {} gene ontology terms will be affected by elimination.",
            level, no_terms
        );
    }

    for term in go_to_remove.iter() {
        let ancestors = go_obj.get_ancestors(term);
        let ancestors_final: Vec<String> = ancestors.cloned().unwrap_or_else(Vec::new);

        if let Some(genes_to_remove) = go_obj.get_genes(term) {
            let genes_to_remove = genes_to_remove.clone();
            go_obj.remove_genes(&ancestors_final, &genes_to_remove);
        }
    }

    GoElimLevelResults {
        go_ids: go_identifiers,
        pvals,
        odds_ratios,
        hits,
        gene_set_lengths,
    }
}
