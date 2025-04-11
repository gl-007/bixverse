use rand::prelude::*;
use rand::seq::SliceRandom;
use rayon::prelude::*;
use statrs::distribution::{Continuous, ContinuousCDF, Normal};
use std::collections::HashSet;

/// Split a vector randomly into two chunks with one being [..x] and the other [x..]
pub fn split_vector_randomly(vec: Vec<f64>, x: usize, seed: u64) -> (Vec<f64>, Vec<f64>) {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut shuffled = vec.clone();
    shuffled.shuffle(&mut rng);

    let first_set = shuffled[..x].to_vec();
    let second_set = shuffled[x..].to_vec();

    (first_set, second_set)
}

/// Calculate the set similarity. Options are Jaccard (similarity_index = False)
/// or the similarity index calculation.
pub fn set_similarity(
    s_1: &HashSet<String>,
    s_2: &HashSet<String>,
    overlap_coefficient: bool,
) -> f64 {
    let i = s_1.intersection(s_2).count() as u64;
    let u = if overlap_coefficient {
        std::cmp::min(s_1.len(), s_2.len()) as u64
    } else {
        s_1.union(s_2).count() as u64
    };
    i as f64 / u as f64
}

/// Calculate the Hedge's g effect size and its standard error
pub fn hedge_g_effect(
    mean_a: &[f64],
    mean_b: &[f64],
    std_a: &[f64],
    std_b: &[f64],
    n_a: usize,
    n_b: usize,
    small_sample_correction: bool,
) -> (Vec<f64>, Vec<f64>) {
    let total_n = (n_a + n_b) as f64;
    let res: Vec<(f64, f64)> = mean_a
        .par_iter()
        .zip(mean_b.par_iter())
        .zip(std_a.par_iter())
        .zip(std_b.par_iter())
        .map(|(((mean_a, mean_b), std_a), std_b)| {
            let pooled_sd = (((n_a - 1) as f64 * std_a.powi(2) + (n_b - 1) as f64 * std_b.powi(2))
                / ((n_a + n_b - 2) as f64))
                .sqrt();
            let effect_size = (mean_a - mean_b) / pooled_sd;
            // Small sample correction if needed
            let effect_size = if small_sample_correction {
                let correction_factor =
                    ((total_n - 3.0) / (total_n - 2.25)) * ((total_n - 2.0) / total_n).sqrt();
                correction_factor * effect_size
            } else {
                effect_size
            };
            let standard_error = ((total_n / (n_a as f64 * n_b as f64))
                + (effect_size.powi(2) / (2.0 * total_n)))
                .sqrt();
            (effect_size, standard_error)
        })
        .collect();

    let mut effect_sizes: Vec<f64> = Vec::new();
    let mut standard_errors: Vec<f64> = Vec::new();

    for (effect_size, standard_error) in res {
        effect_sizes.push(effect_size);
        standard_errors.push(standard_error);
    }

    (effect_sizes, standard_errors)
}

/// Transform Z-scores into p-values (assuming normality).
pub fn z_scores_to_pval(z_scores: &[f64]) -> Vec<f64> {
    let normal = Normal::new(0.0, 1.0).unwrap();
    z_scores
        .iter()
        .map(|&z| {
            let abs_z = z.abs();
            if abs_z > 6.0 {
                // Deal with numeric precision problems for very large z-scores.
                let pdf = normal.pdf(abs_z);
                let p = pdf / abs_z * (1.0 - 1.0 / (abs_z * abs_z));
                2.0 * p
            } else {
                2.0 * (1.0 - normal.cdf(abs_z))
            }
        })
        .collect()
}

////////////////////////////
// Radial Basis functions //
////////////////////////////

/// Enum for the RBF function
#[derive(Debug)]
pub enum RbfType {
    Gaussian,
    Bump,
}

/// Parsing the RBF function
pub fn parse_rbf_types(s: &str) -> Option<RbfType> {
    match s.to_lowercase().as_str() {
        "gaussian" => Some(RbfType::Gaussian),
        "bump" => Some(RbfType::Bump),
        _ => None,
    }
}

/// Gaussian Radial Basis function
pub fn rbf_gaussian(dist: &[f64], epsilon: &f64) -> Vec<f64> {
    dist.par_iter()
        .map(|x| f64::exp(-((x * *epsilon).powi(2))))
        .collect()
}

/// Bump Radial Basis function
/// Will set dist >= 1 / epsilon to 0
pub fn rbf_bump(dist: &[f64], epsilon: &f64) -> Vec<f64> {
    dist.par_iter()
        .map(|x| {
            if *x < (1.0 / epsilon) {
                f64::exp(-(1_f64 / (1_f64 - (*epsilon * x).powi(2))) + 1_f64)
            } else {
                0_f64
            }
        })
        .collect()
}

// /// Calculate the Jaccard or Set similarity over a vector of HashSets.
// pub fn set_similarity_iter(
//   s_1: HashSet<String>,
//   other_s: Vec<HashSet<String>>,
//   similarity_index: bool,
// ) -> Vec<f64>{
//   let sims: Vec<f64>= other_s
//     .into_iter()
//     .map(|hash_i| {
//       set_similarity(hash_i, &s_1, similarity_index)
//     })
//     .collect();

//   sims
// }
