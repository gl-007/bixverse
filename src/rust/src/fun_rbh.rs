use extendr_api::prelude::*;
use rayon::prelude::*;

use crate::helpers_rbh::*;
use crate::utils_r_rust::{r_nested_list_to_rust, NestedHashMap};
use crate::utils_rust::flatten_vector;

/// Structure to store the RBH results.
pub struct RbhResult {
    pub origin: String,
    pub target: String,
    pub origin_modules: Vec<String>,
    pub target_modules: Vec<String>,
    pub similarities: Vec<f64>,
}

/// Generate reciprocal best hits based on set similarities
///
/// @description This function takes a nested list that contains gene modules/
/// sets derived from various methods and generate identifies reciprocal best
/// hits between gene modules/sets across the different origins. WARNING!
/// Incorrect use can cause kernel crashes. Wrapper around the Rust functions
/// with type checks are provided in the package.
///
/// @param module_list A nested named list. The outer list should contain the
/// origin of the gene modules, the inner list the names of the gene modules and
/// the respective genes in them.
/// @param overlap_coefficient Shall the overlap coefficient instead of the
/// Jaccard similarity be used.
/// @param min_similarity Minimum similarity that should exist between any two
/// given gene modules to actually calculate RBH pairs.
/// @param debug Boolean Boolean that activates print messages for debugging
/// purposes.
///
/// @return A list containing:
///  \itemize{
///   \item origin - The name of the origin of the gene modules.
///   \item target - The name of the target of the gene modules.
///   \item comparisons - Integer vector indicating how many RBH hits were
///   identified in this comparison
///   \item origin_modules - Names of the gene modules from the origin.
///   \item target_modules - Names of the gene modules from the target.
///   \item similarity - The similarities between the two respective gene
///   modules.
/// }
/// @export
#[extendr]
fn rs_rbh_sets(
    module_list: List,
    overlap_coefficient: bool,
    min_similarity: f64,
    debug: bool,
) -> extendr_api::Result<List> {
    let module_list: NestedHashMap = r_nested_list_to_rust(module_list)?;
    // Pull out all the keys
    let origins: Vec<String> = module_list.clone().into_keys().collect();

    let origins_split: Vec<(String, &[String])> = origins
        .iter()
        .enumerate()
        .map(|(i, first)| (first.clone(), &origins[i + 1..]))
        .take_while(|(_, rest)| !rest.is_empty())
        .collect();

    let rbh_results: Vec<Vec<RbhResult>> = origins_split
        .par_iter()
        .map(|(origin_module, target_modules)| {
            // Parallel iterator starts here
            let origin_module_data = module_list.get(origin_module).unwrap();

            // Iterate over the remaining target modules
            let res: Vec<RbhResult> = target_modules
                .iter()
                .map(|target| {
                    let target_module_data = module_list.get(target).unwrap();

                    let rbh_res: RbhTriplet = calculate_rbh_set(
                        origin_module_data,
                        target_module_data,
                        overlap_coefficient,
                        min_similarity,
                        debug,
                    );

                    let mut origin_modules = Vec::new();
                    let mut target_modules = Vec::new();
                    let mut similarities = Vec::new();

                    for (origin, target, similarity) in rbh_res {
                        origin_modules.push(origin);
                        target_modules.push(target);
                        similarities.push(similarity)
                    }

                    RbhResult {
                        origin: origin_module.to_string(),
                        target: target.to_string(),
                        origin_modules,
                        target_modules,
                        similarities,
                    }
                })
                .collect();

            res
        })
        .collect();

    // Flatten and extract relevant data.
    let rbh_results_flatten: Vec<_> = flatten_vector(rbh_results);

    let mut origin = Vec::new();
    let mut target = Vec::new();
    let mut comparisons = Vec::new();
    let mut origin_modules = Vec::new();
    let mut target_modules = Vec::new();
    let mut similarity = Vec::new();

    for module in rbh_results_flatten {
        origin.push(module.origin);
        target.push(module.target);
        origin_modules.push(module.origin_modules);
        target_modules.push(module.target_modules);
        similarity.push(module.similarities.clone());
        comparisons.push(module.similarities.len());
    }

    let origin_modules: Vec<_> = flatten_vector(origin_modules);
    let target_modules: Vec<_> = flatten_vector(target_modules);
    let similarity: Vec<_> = flatten_vector(similarity);

    Ok(list!(
        origin = origin,
        target = target,
        comparisons = comparisons,
        origin_modules = origin_modules,
        target_modules = target_modules,
        similarity = similarity
    ))
}

extendr_module! {
    mod fun_rbh;
    fn rs_rbh_sets;
}
