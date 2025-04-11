use extendr_api::prelude::*;
use rayon::prelude::*;

use crate::helpers_ontology::*;
use crate::utils_r_rust::r_list_to_hashmap_set;
use crate::utils_rust::flatten_vector;

/// Calculate the semantic similarity in an ontology
///
/// @description This function calculates the Resnik and Lin similarity for a given ontology.
///
/// @param terms Vector of strings. The terms in the ontology you wish to screen.
/// @param ancestor_list R list with names being the term and the elements in the list the names
/// of the ancestors.
/// @param ic_list R list with the names being the term and the elements the information content
/// of this given term. Needs to be a single float!
///
/// @return A list with:
/// \itemize{
///   \item term1 - String, the first term.
///   \item term2 - String, the second term.
///   \item resnik_sim - Float, the unnormalised Resnik similarity.
///   \item lin_sim - Float, the Lin similarity.
/// }
///
/// @export
#[extendr]
fn rs_onto_similarity(
    terms: Vec<String>,
    ancestor_list: List,
    ic_list: List,
) -> extendr_api::Result<List> {
    let ancestors_map = r_list_to_hashmap_set(ancestor_list)?;
    let ic_map = ic_list_to_ic_hashmap(ic_list);

    let terms_split: Vec<(String, &[String])> = terms
        .iter()
        .enumerate()
        .map(|(i, first)| (first.clone(), &terms[i + 1..]))
        .take_while(|(_, rest)| !rest.is_empty())
        .collect();

    let onto_sim: Vec<Vec<OntoSimRes<'_>>> = terms_split
        .par_iter()
        .map(|(t1, others)| {
            let mut sim_vec: Vec<OntoSimRes<'_>> = Vec::with_capacity(others.len());
            others.iter().for_each(|t2| {
                let sim_res = onto_sim(t1, t2, &ancestors_map, &ic_map);
                sim_vec.push(sim_res);
            });
            sim_vec
        })
        .collect();

    let onto_sim = flatten_vector(onto_sim);

    let mut term1 = Vec::new();
    let mut term2 = Vec::new();
    let mut res_sim = Vec::new();
    let mut lin_sim = Vec::new();

    for sim_res in onto_sim.iter() {
        if sim_res.resnik_sim > 0.0 {
            term1.push(sim_res.t1.to_string());
            term2.push(sim_res.t2.to_string());
            res_sim.push(sim_res.resnik_sim);
            lin_sim.push(sim_res.lin_sim);
        }
    }

    Ok(list!(
        term1 = term1,
        term2 = term2,
        resnik_sim = res_sim,
        lin_sim = lin_sim
    ))
}

extendr_module! {
  mod fun_ontology;
  fn rs_onto_similarity;
}
