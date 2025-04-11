use extendr_api::prelude::*;
use std::collections::{HashMap, HashSet};

/// Structure to store the Ontology similarity results
pub struct OntoSimRes<'a> {
    pub t1: &'a str,
    pub t2: &'a str,
    pub resnik_sim: f64,
    pub lin_sim: f64,
}

/// Calculate the semantic similarity. Return Resnik and Lin similarity in one go.
pub fn onto_sim<'a>(
    t1: &'a str,
    t2: &'a str,
    ancestor_map: &HashMap<String, HashSet<String>>,
    info_content_map: &HashMap<String, f64>,
) -> OntoSimRes<'a> {
    // Default Hashmap to avoid all types of tests here...
    let default_hash: HashSet<String> =
        HashSet::from_iter(std::iter::once("I have no ancestors".to_string()));
    let ancestor_1 = ancestor_map.get(t1).unwrap_or(&default_hash);
    let ancestor_2 = ancestor_map.get(t2).unwrap_or(&default_hash);
    let mica = ancestor_1
        .intersection(ancestor_2)
        .map(|ancestor| info_content_map.get(ancestor).cloned().unwrap_or(0.0))
        .fold(0.0, f64::max);
    let t1_ic = info_content_map.get(t1).unwrap_or(&0.0);
    let t2_ic = info_content_map.get(t2).unwrap_or(&0.0);
    let lin_sim = 2.0 * mica / (t1_ic + t2_ic);

    OntoSimRes {
        t1,
        t2,
        resnik_sim: mica,
        lin_sim,
    }
}

/// Transform an R list that hopefully contains the IC into a HashMap
/// of floats
pub fn ic_list_to_ic_hashmap(r_list: List) -> HashMap<String, f64> {
    let mut hashmap = HashMap::with_capacity(r_list.len());
    for (name, x) in r_list {
        let name = name.to_string();
        let ic_val = x.as_real().unwrap_or(0.0);
        hashmap.insert(name, ic_val);
    }
    hashmap
}
