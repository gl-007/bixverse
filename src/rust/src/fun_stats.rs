use extendr_api::prelude::*;
use rand::prelude::*;
use rayon::prelude::*;
// use std::sync::{Arc, Mutex};

use crate::helpers_linalg::{col_means, col_sds};
use crate::utils_r_rust::r_matrix_to_faer;
use crate::utils_stats::{hedge_g_effect, split_vector_randomly};
// use std::collections::HashSet;
// use crate::utils_r_rust::r_list_to_str_vec;

/// Fast AUC calculation
///
/// @description This function calculates rapidly AUCs based on an approximation.
///
/// @param pos_scores The scores of your hits.
/// @param neg_scores The scores of your non-hits.
/// @param iters Number of iterations to run the function for.
/// Recommended size: 10000L.
/// @param seed Seed.
///
/// @return The AUC.
///
/// @export
#[extendr]
fn rs_fast_auc(pos_scores: Vec<f64>, neg_scores: Vec<f64>, iters: usize, seed: u64) -> f64 {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut count = 0;

    for _ in 0..iters {
        let pos_sample = pos_scores.choose(&mut rng).unwrap();
        let neg_sample = neg_scores.choose(&mut rng).unwrap();
        if pos_sample > neg_sample {
            count += 1;
        }
    }

    count as f64 / iters as f64
}

/// Create random AUCs
///
/// @description This function creates a random set of AUCs based on a score
/// vector and a size of the positive set. This can be used for permutation-
/// based estimation of Z-scores and subsequently p-values.
///
/// @param score_vec The overall vector of scores.
/// @param size_pos The size of the hits represented in the score_vec.
/// @param random_iters Number of random AUCs to generate.
/// @param auc_iters Number of random iterations to approximate the AUCs.
/// Recommended size: 10000L.
/// @param seed Seed.
///
/// @return A vector of random AUCs based the score vector and size of the
/// positive set.
///
/// @export
#[extendr]
fn rs_create_random_aucs(
    score_vec: Vec<f64>,
    size_pos: usize,
    random_iters: usize,
    auc_iters: usize,
    seed: u64,
) -> Vec<f64> {
    let iter_vec: Vec<usize> = (0..random_iters).collect();

    let random_aucs: Vec<_> = iter_vec
        .par_iter()
        .map(|x| {
            let scores = split_vector_randomly(score_vec.clone(), size_pos, *x as u64 + seed);

            rs_fast_auc(scores.0, scores.1, auc_iters, *x as u64 + 1 + seed)
        })
        .collect();

    random_aucs
}

/// Calculate the Hedge's G effect
///
/// @description Calculates the Hedge's G effect for two sets of matrices. The
/// function assumes that rows = samples and columns = features.
/// WARNING! Incorrect use can cause kernel crashes. Wrapper around the Rust
/// functions with type checks are provided in the package.
///
/// @param mat_a The matrix of samples and features in grp A for which to
/// calculate the Hedge's G effect.
/// @param mat_b The matrix of samples and features in grp B for which to
/// calculate the Hedge's G effect.
/// @param small_sample_correction Shall the small sample correction be applied.
///
/// @return Returns the harmonic sum according to the OT calculation.
///
/// @export
#[extendr]
fn rs_hedges_g(mat_a: RMatrix<f64>, mat_b: RMatrix<f64>, small_sample_correction: bool) -> List {
    let mat_a: faer::Mat<f64> = r_matrix_to_faer(&mat_a);
    let mat_b: faer::Mat<f64> = r_matrix_to_faer(&mat_b);

    let n_a = mat_a.nrows();
    let n_b = mat_b.nrows();

    let mean_a = col_means(&mat_a);
    let mean_b = col_means(&mat_b);

    let std_a = col_sds(&mat_a);
    let std_b = col_sds(&mat_b);

    let (es, se) = hedge_g_effect(
        &mean_a,
        &mean_b,
        &std_a,
        &std_b,
        n_a,
        n_b,
        small_sample_correction,
    );

    list!(effect_sizes = es, standard_errors = se)
}

/// Calculate a BH-based FDR
///
/// @description Rust implementation that will be faster if you have an
/// terrifying amount of p-values to adjust.
///
/// @param pvals Numeric vector. The p-values you wish to adjust.
///
/// @return The Benjamini-Hochberg adjusted p-values.
///
/// @export
#[extendr]
fn rs_fdr_adjustment(pvals: &[f64]) -> Vec<f64> {
    let n = pvals.len();
    let n_f64 = n as f64;

    let mut indexed_pval: Vec<(usize, f64)> =
        pvals.par_iter().enumerate().map(|(i, &x)| (i, x)).collect();

    indexed_pval
        .sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    let adj_pvals_tmp: Vec<f64> = indexed_pval
        .par_iter()
        .enumerate()
        .map(|(i, (_, p))| (n_f64 / (i + 1) as f64) * p)
        .collect();

    let mut current_min = adj_pvals_tmp[n - 1].min(1.0);
    let mut monotonic_adj = vec![current_min; n];

    for i in (0..n - 1).rev() {
        current_min = current_min.min(adj_pvals_tmp[i]).min(1.0);
        monotonic_adj[i] = current_min;
    }

    let mut adj_pvals = vec![0.0; n];

    for (i, &(original_idx, _)) in indexed_pval.iter().enumerate() {
        adj_pvals[original_idx] = monotonic_adj[i];
    }

    adj_pvals
}

// /// Set similarities
// ///
// /// This function calculates the Jaccard or similarity index between a given
// /// string vector and a list of other string vectors.
// ///
// /// @param string The String vector against which to calculate the setsimilarities.
// /// @param string_list The list of character vectors for which to calculate the set similarities.
// /// @param similarity_index Shall the similarity index instead of the Jaccard similarity be calculated.
// ///
// /// @export
// #[extendr]
// fn rs_set_sim_list(
//   string: Vec<String>,
//   string_list: List,
//   similarity_index: bool,
// ) -> Vec<f64> {
//     let string_vec = r_list_to_str_vec(string_list);
//     let string: HashSet<_> = string
//       .iter()
//       .collect();
//     let values: Vec<(u64, u64)> = string_vec
//       .iter()
//       .map(|s| {
//         let s_hash: HashSet<_> = s
//         .iter()
//         .collect();
//         let i = s_hash.intersection(&string).count() as u64;
//         let u = if similarity_index {
//           std::cmp::min(s_hash.len(), string.len()) as u64
//         } else {
//           s_hash.union(&string).count() as u64
//         };
//         (i, u)
//       })
//       .collect();

//     let mut sim: Vec<f64> = Vec::new();

//     for (i, u) in values {
//       let sim_i: f64 = (i as f64) / (u as f64);
//       sim.push(sim_i)
//     }

//     sim
//   }

extendr_module! {
    mod fun_stats;
    // fn rs_set_sim_list;
    fn rs_fast_auc;
    fn rs_create_random_aucs;
    fn rs_hedges_g;
    fn rs_fdr_adjustment;
}
