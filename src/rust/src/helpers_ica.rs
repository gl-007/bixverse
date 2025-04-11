use extendr_api::prelude::*;
use faer::{
    linalg::solvers::{PartialPivLu, Solve},
    Mat,
};
use rand::prelude::*;
use rand_distr::Normal;
use rayon::iter::*;

use crate::helpers_linalg::{randomised_svd, scale_matrix_col};
use crate::utils_rust::*;

//////////////////////////////
// ENUMS, TYPES, STRUCTURES //
//////////////////////////////

/// Enum for the ICA types
#[derive(Debug)]
pub enum IcaType {
    Exp,
    LogCosh,
}

/// Type alias of the ICA results
type IcaRes = (faer::Mat<f64>, f64);

/// Structure to save ICA parameters
#[derive(Clone, Debug)]
pub struct IcaParams {
    pub maxit: usize,
    pub alpha: f64,
    pub tol: f64,
    pub verbose: bool,
}

/// Structure to save ICA CV results
#[derive(Clone, Debug)]
pub struct IcaCvData {
    pub pre_white_matrices: Vec<Mat<f64>>,
    pub k_matrices: Vec<Mat<f64>>,
}

/////////
// ICA //
/////////

/// Prepare ICA parameters
pub fn prepare_ica_params(r_list: List) -> IcaParams {
    let ica_params = r_list.into_hashmap();

    let maxit = ica_params
        .get("maxit")
        .and_then(|v| v.as_integer())
        .unwrap_or(200) as usize;
    let alpha = ica_params
        .get("alpha")
        .and_then(|v| v.as_real())
        .unwrap_or(1.0);
    let tol = ica_params
        .get("max_tol")
        .and_then(|v| v.as_real())
        .unwrap_or(1e-4);
    let verbose = ica_params
        .get("verbose")
        .and_then(|v| v.as_bool())
        .unwrap_or(false);

    IcaParams {
        maxit,
        alpha,
        tol,
        verbose,
    }
}

/// Whiten a matrix. This is needed pre-processing for ICA.
/// Has the option to use randomised SVD for faster computations.
/// Returns the scaled data and the pre-whitening matrix K.
pub fn prepare_whitening(
    x: &Mat<f64>,
    fast_svd: bool,
    seed: usize,
    rank: usize,
    oversampling: Option<usize>,
    n_power_iter: Option<usize>,
) -> (faer::Mat<f64>, faer::Mat<f64>) {
    let n = x.nrows();

    let centered = scale_matrix_col(x, false);

    let centered = centered.transpose();

    let v = centered * centered.transpose() / n as f64;

    let k = if fast_svd {
        let svd_res = randomised_svd(&v, rank, seed, oversampling, n_power_iter);
        let s: Vec<f64> = svd_res.s.iter().map(|x| 1_f64 / x.sqrt()).collect();
        let d = faer_diagonal_from_vec(s);
        d * svd_res.u.transpose()
    } else {
        let svd_res = v.thin_svd().unwrap();
        let s = svd_res
            .S()
            .column_vector()
            .iter()
            .map(|x| 1_f64 / x.sqrt())
            .collect::<Vec<_>>();
        let d = faer_diagonal_from_vec(s);
        let u = svd_res.U();
        let u_t = u.transpose();
        d * u_t
    };

    (centered.cloned(), k)
}

/// Update the mixing matrix for ICA
pub fn update_mix_mat(w: &Mat<f64>) -> faer::Mat<f64> {
    // SVD
    let svd_res = w.thin_svd().unwrap();

    let s = svd_res.S();
    let u = svd_res.U();
    let s = s
        .column_vector()
        .iter()
        .map(|x| 1_f64 / x)
        .collect::<Vec<_>>();
    let d = faer_diagonal_from_vec(s);

    u * d * u.transpose() * w
}

/// Generate a w_init matrix of size n_comp * n_comp given a random seed.
pub fn create_w_init(n_comp: usize, seed: u64) -> faer::Mat<f64> {
    let mut rng = StdRng::seed_from_u64(seed);
    let normal = Normal::new(0.0, 1.0).unwrap();
    let vec_size = n_comp.pow(2);
    let data: Vec<f64> = (0..vec_size).map(|_| normal.sample(&mut rng)).collect();

    Mat::from_fn(n_comp, n_comp, |i, j| data[i + j * n_comp])
}

/// Parsing the ICA types
pub fn parse_ica_type(s: &str) -> Option<IcaType> {
    match s.to_lowercase().as_str() {
        "exp" => Some(IcaType::Exp),
        "logcosh" => Some(IcaType::LogCosh),
        _ => None,
    }
}

/// Fast ICA implementation based on logcosh.
pub fn fast_ica_logcosh(
    x: &Mat<f64>,
    w_init: &Mat<f64>,
    tol: f64,
    alpha: f64,
    maxit: usize,
    verbose: bool,
) -> IcaRes {
    let p = x.ncols();
    let mut w = update_mix_mat(w_init);
    let mut lim = vec![1000_f64; maxit];

    let mut it = 0;

    while it < maxit && lim[it] > tol {
        let wx: Mat<f64> = &w * x;

        let gwx = Mat::from_fn(wx.nrows(), wx.ncols(), |i, j| {
            let x = wx.get(i, j);
            (alpha * x).tanh()
        });

        let v1 = &gwx * x.transpose() / p as f64;

        let gwx_2 = alpha
            * Mat::from_fn(gwx.nrows(), gwx.ncols(), |i, j| {
                let x = gwx.get(i, j);
                1_f64 - x.powi(2)
            });

        let ones = Mat::from_fn(p, 1, |_, _| 1.0);
        let row_means = (&gwx_2 * &ones) * (1.0 / p as f64);

        let row_means_vec: Vec<f64> = row_means
            .as_ref()
            .col_iter()
            .flat_map(|col| col.iter())
            .copied()
            .collect();

        let v2 = faer_diagonal_from_vec(row_means_vec) * w.clone();

        let w1 = update_mix_mat(&(v1 - v2));

        let w1_up = w1.clone() * w.transpose();

        let tol_it = w1_up
            .diagonal()
            .column_vector()
            .iter()
            .map(|x| (x.abs() - 1.0).abs())
            .fold(f64::NEG_INFINITY, f64::max);

        if it + 1 < maxit {
            lim[it + 1] = tol_it
        }

        if verbose {
            println!("Iteration: {:?}, tol: {:?}", it + 1, tol_it)
        }

        w = w1;

        it += 1;
    }

    let min_tol = array_f64_min(&lim);

    (w, min_tol)
}

/// Fast ICA implementation based on exp.
pub fn fast_ica_exp(
    x: &Mat<f64>,
    w_init: &Mat<f64>,
    tol: f64,
    maxit: usize,
    verbose: bool,
) -> IcaRes {
    let p = x.ncols();
    let mut w = update_mix_mat(w_init);
    let mut lim = vec![1000_f64; maxit];

    let mut it = 0;
    while it < maxit && lim[it] > tol {
        let wx: Mat<f64> = &w * x;

        let gwx = Mat::from_fn(wx.nrows(), wx.ncols(), |i, j| {
            let x = wx.get(i, j);
            x * (-x.powi(2) / 2.0).exp()
        });

        let v1 = &gwx * x.transpose() / p as f64;

        let gwx_2 = Mat::from_fn(wx.nrows(), wx.ncols(), |i, j| {
            let x = wx.get(i, j);
            (1.0 - x.powi(2)) * (-x.powi(2) / 2.0).exp()
        });

        let ones = Mat::from_fn(p, 1, |_, _| 1.0);
        let row_means = (&gwx_2 * &ones) * (1.0 / p as f64);

        let row_means_vec: Vec<f64> = row_means
            .as_ref()
            .col_iter()
            .flat_map(|col| col.iter())
            .copied()
            .collect();

        let v2 = faer_diagonal_from_vec(row_means_vec) * w.clone();

        let w1 = update_mix_mat(&(v1 - v2));

        let w1_up = w1.clone() * w.transpose();

        let tol_it = w1_up
            .diagonal()
            .column_vector()
            .iter()
            .map(|x| (x.abs() - 1.0).abs())
            .fold(f64::NEG_INFINITY, f64::max);

        if it + 1 < maxit {
            lim[it + 1] = tol_it
        }

        if verbose {
            println!("Iteration: {:?}, tol: {:?}", it + 1, tol_it)
        }

        w = w1;

        it += 1;
    }

    let min_tol = array_f64_min(&lim);

    (w, min_tol)
}

/// Iterate through a set of random initialisations with a given pre-whitened
/// matrix, the whitening matrix k and the respective ICA parameters.
pub fn stabilised_ica_iters(
    x_pre_whiten: &Mat<f64>,
    k: &Mat<f64>,
    no_comp: usize,
    no_iters: usize,
    ica_type: &str,
    ica_params: IcaParams,
    random_seed: usize,
) -> (Mat<f64>, Vec<bool>) {
    // Generate the random w_inits
    let w_inits: Vec<Mat<f64>> = (0..no_iters)
        .map(|iter| create_w_init(no_comp, (random_seed + iter) as u64))
        .collect();
    let k_ncol = k.ncols();
    let k_red = k.get(0..no_comp, 0..k_ncol);
    let x_whiten = k_red * x_pre_whiten;

    let ica_type = parse_ica_type(ica_type).unwrap();
    let iter_res: Vec<(Mat<f64>, f64)> = w_inits
        .par_iter()
        .map(|w_init| match ica_type {
            IcaType::Exp => fast_ica_exp(
                &x_whiten,
                w_init,
                ica_params.tol,
                ica_params.maxit,
                ica_params.verbose,
            ),
            IcaType::LogCosh => fast_ica_logcosh(
                &x_whiten,
                w_init,
                ica_params.tol,
                ica_params.alpha,
                ica_params.maxit,
                ica_params.verbose,
            ),
        })
        .collect();

    let mut convergence = Vec::new();
    let mut a_matrices = Vec::new();

    for (a, final_tol) in iter_res {
        a_matrices.push(a);
        convergence.push(final_tol < ica_params.tol);
    }

    let s_matrices: Vec<Mat<f64>> = a_matrices
        .par_iter()
        .map(|a| {
            let w = a * k_red;
            let to_solve = w.clone() * w.transpose();
            let identity = Mat::<f64>::identity(to_solve.nrows(), to_solve.ncols());
            let lu = PartialPivLu::new(to_solve.as_ref());
            let solved = lu.solve(&identity);
            w.transpose() * solved
        })
        .collect();

    let s_combined = colbind_matrices(s_matrices);

    (s_combined, convergence)
}

/// Generate cross-validation like data for ICA.
pub fn create_ica_cv_data(
    x: &Mat<f64>,
    num_folds: usize,
    seed: usize,
    rank: Option<usize>,
) -> IcaCvData {
    let no_samples = x.nrows();
    let no_features = x.ncols();
    let mut indices: Vec<usize> = (0..no_samples).collect();
    let mut rng = StdRng::seed_from_u64(seed as u64);
    indices.shuffle(&mut rng);

    let svd_rank = rank.unwrap_or(no_features);

    let fold_size = no_samples / num_folds;
    let remainder = no_samples % num_folds;

    let mut folds = Vec::with_capacity(num_folds);
    let mut start = 0;

    for idx in 0..num_folds {
        let current_fold_size = if idx < remainder {
            fold_size + 1
        } else {
            fold_size
        };

        let end = start + current_fold_size;
        folds.push(indices[start..end].to_vec());
        start = end;
    }

    let k_x_matrices: Vec<(Mat<f64>, Mat<f64>)> = folds
        .par_iter()
        .map(|test_indices| {
            let train_indices: Vec<usize> = indices
                .iter()
                .filter(|&idx| !test_indices.contains(idx))
                .cloned()
                .collect();
            let mut x_i = Mat::<f64>::zeros(train_indices.len(), no_features);

            for (new_row, old_row) in train_indices.iter().enumerate() {
                for j in 0..no_features {
                    x_i[(new_row, j)] = x[(*old_row, j)];
                }
            }

            prepare_whitening(&x_i, true, seed + 1, svd_rank, None, None)
        })
        .collect();

    let mut pre_white_matrices = Vec::with_capacity(num_folds);
    let mut k_matrices = Vec::with_capacity(num_folds);

    for (x_i, k_i) in k_x_matrices {
        pre_white_matrices.push(x_i);
        k_matrices.push(k_i);
    }

    IcaCvData {
        pre_white_matrices,
        k_matrices,
    }
}

#[allow(clippy::too_many_arguments)]
pub fn stabilised_ica_cv(
    x: Mat<f64>,
    no_comp: usize,
    num_folds: usize,
    no_iters: usize,
    ica_type: &str,
    ica_params: IcaParams,
    ica_cv_data: Option<IcaCvData>,
    seed: usize,
) -> (Mat<f64>, Vec<bool>) {
    // Generate cross-validation data if not provided
    let cv_data = match ica_cv_data {
        Some(data) => data, // Use the provided data
        None => create_ica_cv_data(&x, num_folds, seed, Some(no_comp)), // Generate new data
    };

    // Iterate through bootstrapped samples
    let cv_res: Vec<(Mat<f64>, Vec<bool>)> = cv_data
        .k_matrices
        .par_iter()
        .zip(cv_data.pre_white_matrices)
        .map(|(k_i, x_i)| {
            let (s_i, converged_i) = stabilised_ica_iters(
                &x_i,
                k_i,
                no_comp,
                no_iters,
                ica_type,
                ica_params.clone(),
                seed + 2,
            );

            (s_i, converged_i)
        })
        .collect();

    let mut s_final = Vec::new();
    let mut converged_final = Vec::new();

    for (s_i, converged_i) in cv_res {
        s_final.push(s_i);
        converged_final.push(converged_i);
    }

    let s_final = colbind_matrices(s_final);
    let converged_final = flatten_vector(converged_final);

    (s_final, converged_final)
}
