use extendr_api::prelude::*;

use faer::Mat;
use rayon::prelude::*;

use crate::utils_r_rust::faer_to_r_matrix;
use crate::utils_rust::array_f64_max_min;
use crate::utils_stats::*;

/// Calculate the OT harmonic sum
///
/// @param x The numeric vector (should be between 0 and 1) for which to
/// calculate the harmonic sum
///
/// @return Returns the harmonic sum according to the OT calculation.
///
/// @export
#[extendr]
fn rs_ot_harmonic_sum(mut x: Vec<f64>) -> f64 {
    x.sort_by(|a, b| b.partial_cmp(a).unwrap());

    let harmonic_sum: f64 = x
        .iter()
        .enumerate()
        .map(|(i, x)| x / (i + 1).pow(2) as f64)
        .sum();

    let max_sum: f64 = vec![1; x.len()]
        .into_iter()
        .enumerate()
        .map(|(i, x)| x as f64 / (i + 1).pow(2) as f64)
        .sum();

    harmonic_sum / max_sum
}

/// Reconstruct a matrix from a flattened upper triangle vector
///
/// @description This function takes a flattened vector of the upper triangle
/// from a symmetric matrix (think correlation matrix) and reconstructs the full
/// dense matrix for you.
///
/// @param cor_vector Numeric vector. The vector of correlation coefficients
/// that you want to use to go back to a dense matrix.
/// @param shift Integer. If you applied a shift, i.e. included the diagonal
/// values = 0; or excluded the diagonal values = 1.
/// @param n Integer. Original dimension (i.e., ncol/nrow) of the matrix to be
/// reconstructed.
///
/// @return The dense R matrix.
///
/// @export
#[extendr]
fn rs_upper_triangle_to_dense(
    cor_vector: &[f64],
    shift: usize,
    n: usize,
) -> extendr_api::RArray<f64, [usize; 2]> {
    let mut mat = Mat::<f64>::zeros(n, n);
    let mut idx = 0;
    for i in 0..n {
        for j in i..n {
            if shift == 1 && i == j {
                mat[(i, j)] = 1_f64
            } else {
                mat[(i, j)] = cor_vector[idx];
                mat[(j, i)] = cor_vector[idx];
                idx += 1;
            }
        }
    }

    faer_to_r_matrix(mat.as_ref())
}

/// Apply a Radial Basis Function
///
/// @description Applies a radial basis function (RBF) to a given distance
/// vector. Has at the moment a Gaussian version and a Bump version.
///
/// @param x Numeric vector. The distances you wish to apply the Gaussian kernel
/// onto.
/// @param epsilon Float. Epsilon parameter for the RBF.
/// @param rbf_type String. Needs to be from `c("gaussian", "bump)`.
///
/// @return The affinities after the Kernel was applied.
///
/// @export
#[extendr]
fn rs_rbf_function(x: &[f64], epsilon: f64, rbf_type: &str) -> extendr_api::Result<Vec<f64>> {
    let rbf_fun =
        parse_rbf_types(rbf_type).ok_or_else(|| format!("Invalid RBF function: {}", rbf_type))?;

    let res: Vec<f64> = match rbf_fun {
        RbfType::Gaussian => rbf_gaussian(x, &epsilon),
        RbfType::Bump => rbf_bump(x, &epsilon),
    };

    Ok(res)
}

/// Apply a range normalisation on a vector.
///
/// @description Applies a range normalisation on an R vector.
///
/// @param x Numerical vector. The data to normalise.
/// @param max_val Numeric. The upper bound value to normalise into. If set to 1,
/// the function will be equal to a min-max normalisation.
/// @param min_val Numeric. The lower bound value to normalise into. If set to 0,
/// the function will equal a min-max normalisation.
///
/// @return Normalised values
///
/// @export
#[extendr]
fn rs_range_norm(x: &[f64], max_val: f64, min_val: f64) -> Vec<f64> {
    let (x_min, x_max) = array_f64_max_min(x);
    let denom = x_max - x_min;
    let scale = (max_val - min_val) / denom;
    x.par_iter()
        .map(|x| (x - x_min) * scale + min_val)
        .collect()
}

extendr_module! {
    mod fun_helpers;
    fn rs_upper_triangle_to_dense;
    fn rs_ot_harmonic_sum;
    fn rs_rbf_function;
    fn rs_range_norm;
}
