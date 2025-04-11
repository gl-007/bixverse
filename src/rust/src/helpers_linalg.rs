use faer::Mat;
use rand::prelude::*;
use rand_distr::Normal;
use rayon::iter::*;

use crate::utils_rust::*;
use crate::utils_stats::*;

//////////////////////////////
// ENUMS, TYPES, STRUCTURES //
//////////////////////////////

/// Structure for random SVD results
pub struct RandomSvdResults {
    pub u: faer::Mat<f64>,
    pub v: faer::Mat<f64>,
    pub s: Vec<f64>,
}

/// Structure for DiffCor results
pub struct DiffCorRes {
    pub r_a: Vec<f64>,
    pub r_b: Vec<f64>,
    pub z_score: Vec<f64>,
    pub p_vals: Vec<f64>,
}

//////////////////////////////
// SCALING, COVAR, COR, PCA //
//////////////////////////////

/// Calculates the columns means of a matrix and returns it as Vec<f64>
pub fn col_means(mat: &Mat<f64>) -> Vec<f64> {
    let n_rows = mat.nrows();
    let ones = Mat::from_fn(n_rows, 1, |_, _| 1.0);
    let means = (ones.transpose() * mat) / n_rows as f64;

    means.row(0).iter().cloned().collect()
}

/// Calculates the column sums of a matrix and returns it as Vec<f64>
pub fn col_sums(mat: &Mat<f64>) -> Vec<f64> {
    let n_rows = mat.nrows();
    let ones = Mat::from_fn(n_rows, 1, |_, _| 1.0);
    let col_sums = ones.transpose() * mat;

    col_sums.row(0).iter().cloned().collect()
}

/// Calculate the column standard deviations
pub fn col_sds(mat: &Mat<f64>) -> Vec<f64> {
    let n = mat.nrows() as f64;
    let n_cols = mat.ncols();

    // Calculate means and SDs in one pass
    let (_, m2): (Vec<f64>, Vec<f64>) = (0..n_cols)
        .map(|j| {
            let mut mean = 0.0;
            let mut m2 = 0.0;
            let mut count = 0.0;

            for i in 0..mat.nrows() {
                count += 1.0;
                let delta = mat[(i, j)] - mean;
                mean += delta / count;
                let delta2 = mat[(i, j)] - mean;
                m2 += delta * delta2;
            }
            (mean, (m2 / (n - 1.0)).sqrt())
        })
        .unzip();

    m2
}

/// Scale a matrix by its mean (column wise)
pub fn scale_matrix_col(mat: &Mat<f64>, scale_sd: bool) -> Mat<f64> {
    let n_rows = mat.nrows();
    let n_cols = mat.ncols();

    let mut means = vec![0.0; n_cols];
    for j in 0..n_cols {
        for i in 0..n_rows {
            means[j] += mat[(i, j)];
        }
        means[j] /= n_rows as f64;
    }

    let mut result = mat.clone();
    for j in 0..n_cols {
        let mean = means[j];
        for i in 0..n_rows {
            result[(i, j)] -= mean;
        }
    }

    if !scale_sd {
        return result;
    }

    let mut std_devs = vec![0.0; n_cols];
    for j in 0..n_cols {
        for i in 0..n_rows {
            let val = result[(i, j)];
            std_devs[j] += val * val;
        }
        std_devs[j] = (std_devs[j] / (n_rows as f64 - 1.0)).sqrt();
        if std_devs[j] < 1e-10 {
            std_devs[j] = 1.0;
        }
    }

    for j in 0..n_cols {
        let std_dev = std_devs[j];
        for i in 0..n_rows {
            result[(i, j)] /= std_dev;
        }
    }

    result
}

/// Calculate the column-wise co-variance
pub fn column_covariance(mat: &Mat<f64>) -> Mat<f64> {
    let n_rows = mat.nrows();
    let centered = scale_matrix_col(mat, false);
    let covariance = (centered.transpose() * &centered) / (n_rows - 1) as f64;

    covariance
}

/// Calculate the column-wise correlation. Option to use spearman.
pub fn column_correlation(mat: &Mat<f64>, spearman: bool) -> Mat<f64> {
    let mat = if spearman {
        let ranked_vecs: Vec<Vec<f64>> = mat
            .par_col_iter()
            .map(|x_i| {
                let x_i: Vec<f64> = x_i.iter().copied().collect();
                rank_vector(&x_i)
            })
            .collect();

        nested_vector_to_faer_mat(ranked_vecs)
    } else {
        mat.cloned()
    };

    let scaled = scale_matrix_col(&mat, true);

    let nrow = scaled.nrows() as f64;

    let cor = scaled.transpose() * &scaled / (nrow - 1_f64);

    cor
}

/// Calculate differential correlations
pub fn calculate_diff_correlation(
    mat_a: &Mat<f64>,
    mat_b: &Mat<f64>,
    no_sample_a: usize,
    no_sample_b: usize,
    spearman: bool,
) -> DiffCorRes {
    let mut cors_a: Vec<f64> = Vec::new();
    let mut cors_b: Vec<f64> = Vec::new();

    let upper_triangle_indices = upper_triangle_indices(mat_a.ncols(), 1);

    for (&r, &c) in upper_triangle_indices
        .0
        .iter()
        .zip(upper_triangle_indices.1.iter())
    {
        cors_a.push(*mat_a.get(r, c));
        cors_b.push(*mat_b.get(r, c));
    }

    // Maybe save the original correlations... Note to myself.
    let original_cor_a = cors_a.to_vec();
    let original_cor_b = cors_b.to_vec();

    cors_a.par_iter_mut().for_each(|x| *x = x.atanh());
    cors_b.par_iter_mut().for_each(|x| *x = x.atanh());

    // Constant will depend on if Spearman or Pearson
    let constant = if spearman { 1.06 } else { 1.0 };
    let denominator =
        ((constant / (no_sample_a as f64 - 3.0)) + (constant / (no_sample_b as f64 - 3.0))).sqrt();

    let z_scores: Vec<f64> = cors_a
        .par_iter()
        .zip(cors_b.par_iter())
        .map(|(a, b)| (a - b) / denominator)
        .collect();

    let p_values = z_scores_to_pval(&z_scores);

    DiffCorRes {
        r_a: original_cor_a,
        r_b: original_cor_b,
        z_score: z_scores,
        p_vals: p_values,
    }
}

/// Get the eigenvalues and vectors from a symmetric matrix
pub fn get_top_eigenvalues(matrix: &Mat<f64>, top_n: usize) -> Vec<(f64, Vec<f64>)> {
    // Ensure the matrix is square
    assert!(matrix.nrows() == matrix.ncols(), "Matrix must be square");

    let eigendecomp = matrix.eigen_from_real().unwrap();

    let s = eigendecomp.S();
    let u = eigendecomp.U();

    // Extract the real part of the eigenvalues and vectors
    let mut eigenpairs = s
        .column_vector()
        .iter()
        .zip(u.col_iter())
        .map(|(l, v)| {
            let l_real = l.re;
            let v_real = v.iter().map(|v_i| v_i.re).collect::<Vec<f64>>();
            (l_real, v_real)
        })
        .collect::<Vec<(f64, Vec<f64>)>>();

    // Sort and return Top N
    eigenpairs.sort_by(|a, b| b.0.total_cmp(&a.0));

    let res: Vec<(f64, Vec<f64>)> = eigenpairs.into_iter().take(top_n).collect();

    res
}

/// Implementation of random Singular Value Decomposition to be faster
/// and computationally WAY more efficient.
pub fn randomised_svd(
    x: &Mat<f64>,
    rank: usize,
    seed: usize,
    oversampling: Option<usize>,
    n_power_iter: Option<usize>,
) -> RandomSvdResults {
    let ncol = x.ncols();
    let nrow = x.nrows();

    // Oversampling for better accuracy
    let os = oversampling.unwrap_or(10);
    let sample_size = (rank + os).min(ncol.min(nrow));
    let n_iter = n_power_iter.unwrap_or(2);

    // Create a random matrix
    let mut rng = StdRng::seed_from_u64(seed as u64);
    let normal = Normal::new(0.0, 1.0).unwrap();
    let omega = Mat::from_fn(ncol, sample_size, |_, _| normal.sample(&mut rng));

    // Multiply random matrix with original and use QR composition to get
    // low rank approximation of x
    let y = x * omega;

    let mut q = y.qr().compute_thin_Q();
    for _ in 0..n_iter {
        let z = x.transpose() * q;
        q = (x * z).qr().compute_thin_Q();
    }

    // Perform the SVD on the low-rank approximation
    let b = q.transpose() * x;
    let svd = b.thin_svd().unwrap();

    RandomSvdResults {
        u: q * svd.U(),
        v: svd.V().cloned(), // Use clone instead of manual copying
        s: svd.S().column_vector().iter().copied().collect(),
    }
}

/// Iterate over a distance vector with given RBF function and epsilon parameter
/// and return the column sums for the matrices based on the epsilons.
pub fn rbf_iterate_epsilons(
    dist: &[f64],
    epsilons: &[f64],
    n: usize,
    shift: usize,
    rbf_type: &str,
) -> Result<faer::Mat<f64>, String> {
    // Now specifying String as the error type
    let rbf_fun =
        parse_rbf_types(rbf_type).ok_or_else(|| format!("Invalid RBF function: {}", rbf_type))?;

    let k_res: Vec<Vec<f64>> = epsilons
        .par_iter()
        .map(|epsilon| {
            let affinity_adj = match rbf_fun {
                RbfType::Gaussian => rbf_gaussian(dist, epsilon),
                RbfType::Bump => rbf_bump(dist, epsilon),
            };
            let affinity_adj_mat = upper_triangle_to_sym_faer(&affinity_adj, shift, n);
            col_sums(&affinity_adj_mat)
        })
        .collect();

    Ok(nested_vector_to_faer_mat(k_res))
}
