use extendr_api::prelude::*;

use crate::helpers_ica::*;
use crate::utils_r_rust::{faer_to_r_matrix, r_matrix_to_faer};

/// Prepare the data for whitening
///
/// @description Prepares the data for subsequent usag in ICA.
/// WARNING! Incorrect use can cause kernel crashes. Wrapper around the Rust
/// functions with type checks are provided in the package.
///
/// @param x The matrix to whiten. The whitening will happen over the columns.
/// @param fast_svd Boolean. Shall a randomised SVD be used. This is way faster
/// on larger data sets.
/// @param seed Integer. Only relevant with fast_svd is set to `TRUE`.
/// @param rank Integer. How many ranks to use for the fast SVD approximation.
/// If you supply `NULL`, it will default to `10L`. Only relevant with
/// fast_svd is set to `TRUE`.
/// @param oversampling Integer. Oversampling parameter to make the approximation
/// more precise. If you supply `NULL`, it will default to `10L`. Only relevant
/// with fast_svd is set to `TRUE`.
/// @param n_power_iter Integer. How much shall the QR low rank approximation be
/// powered. If you supply `NULL`, it will default to `2L`.
///
///
/// @return A list containing:
///  \itemize{
///   \item x - The preprocessed matrix.
///   \item k - The pre-whitening matrix k.
/// }
///
/// @export
#[extendr]
fn rs_prepare_whitening(
    x: RMatrix<f64>,
    fast_svd: bool,
    seed: usize,
    rank: Option<usize>,
    oversampling: Option<usize>,
    n_power_iter: Option<usize>,
) -> List {
    let x = r_matrix_to_faer(&x);
    let rank = rank.unwrap_or(10);

    let (x, k) = prepare_whitening(&x, fast_svd, seed, rank, oversampling, n_power_iter);

    list!(
        x = faer_to_r_matrix(x.as_ref()),
        k = faer_to_r_matrix(k.as_ref())
    )
}

/// Run the Rust implementation of fast ICA.
///
/// @description This function serves as a wrapper over the fast ICA
/// implementations in Rust. It assumes a whitened matrix and also an
/// intialised w_init. WARNING! Incorrect use can cause kernel crashes. Wrapper
/// around the Rust functions with type checks are provided in the package.
///
/// @param whiten Numerical matrix. The whitened matrix.
/// @param w_init Numerical matrix. The initial unmixing matrix. ncols need to
/// be equal to nrows of whiten.
/// @param ica_type String. One of 'logcosh' or 'exp'.
/// @param ica_params A list containing:
///  \itemize{
///   \item maxit - Integer. Maximum number of iterations for ICA.
///   \item alpha - Float. The alpha parameter for the logcosh version of ICA.
///   Should be between 1 to 2.
///   \item max_tol - Maximum tolerance of the algorithm
///   \item verbose - Verbosity of the function, i.e., shall individual iters
///   be shown.
/// }
/// If the list is empty or the expected elements are not found, default values
/// are used.
///
/// @return A list with the following items:
///  \itemize{
///   \item mixing - The mixing matrix for subsequent usage.
///   \item converged - Boolean if the algorithm converged.
/// }
///
/// @export
#[extendr]
fn rs_fast_ica(
    whiten: RMatrix<f64>,
    w_init: RMatrix<f64>,
    ica_type: &str,
    ica_params: List,
) -> extendr_api::Result<List> {
    // assert!(!whiten.nrows() == w_init.ncols(), "The dimensions of the provided matrices don't work");
    let ica_params = prepare_ica_params(ica_params);

    let x = r_matrix_to_faer(&whiten);
    let w_init = r_matrix_to_faer(&w_init);

    let ica_type =
        parse_ica_type(ica_type).ok_or_else(|| format!("Invalid ICA type: {}", ica_type))?;

    let a = match ica_type {
        IcaType::Exp => fast_ica_exp(
            &x,
            &w_init,
            ica_params.tol,
            ica_params.maxit,
            ica_params.verbose,
        ),
        IcaType::LogCosh => fast_ica_logcosh(
            &x,
            &w_init,
            ica_params.tol,
            ica_params.alpha,
            ica_params.maxit,
            ica_params.verbose,
        ),
    };

    Ok(list!(
        mixing = faer_to_r_matrix(a.0.as_ref()),
        converged = a.1 < ica_params.tol
    ))
}

/// Run ICA over a given no_comp with random initilisations of w_init
///
/// @description This function implements a stabilised ICA like algorithm in
/// Rust. Briefly, it generates random w_init matrices (total number being
/// no_random_init) and runs ICA given the x_processed and k data over these.
/// The function returns combined S from the different runs and a boolean
/// vector indicating if this specific run converged.
///
/// @param x1 Numerical matrix. The processed matrix (but not yet
/// whitened!)
/// @param k Numerical matrix. The whitening matrix.
/// @param no_comp Integer. Number of independent components to return.
/// @param no_random_init Integer. Number of random initialisations to test.
/// @param ica_type String. One of 'logcosh' or 'exp'.
/// @param random_seed Integer. Seed for randomisations.
/// @param ica_params A list containing:
/// \itemize{
///   \item maxit - Integer. Maximum number of iterations for ICA.
///   \item alpha - Float. The alpha parameter for the logcosh version of ICA.
///   Should be between 1 to 2.
///   \item max_tol - Maximum tolerance of the algorithm
///   \item verbose - Verbosity of the function, i.e., shall individual iters
///   be shown.
/// }
/// If the list is empty or the expected elements are not found, default values
/// are used.
///
/// @return A list containing:
/// \itemize{
///   \item s_combined - The combined matrices for S. Dimensions are nrows =
///   features; and ncols = ncomp * no_random_init.
///   \item converged - Boolean vector indicating if the respective run reached
///   convergence. Length = no_random_init
/// }
///
/// @export
#[extendr]
fn rs_ica_iters(
    x1: RMatrix<f64>,
    k: RMatrix<f64>,
    no_comp: usize,
    no_random_init: usize,
    ica_type: &str,
    random_seed: usize,
    ica_params: List,
) -> List {
    if k.nrows() < no_comp {
        panic!(
            "Number of rows in k ({}) must be at least as large as no_comp ({})",
            k.nrows(),
            no_comp
        );
    }
    let x_processed = r_matrix_to_faer(&x1);
    let k = r_matrix_to_faer(&k);

    let ica_params = prepare_ica_params(ica_params);

    let (s_combined, converged) = stabilised_ica_iters(
        &x_processed,
        &k,
        no_comp,
        no_random_init,
        ica_type,
        ica_params,
        random_seed,
    );

    list!(
        s_combined = faer_to_r_matrix(s_combined.as_ref()),
        converged = converged
    )
}

/// Run ICA with cross-validation and random initialsiation
///
/// @description This function will split the data into `no_folds` and apply
/// ICA with `no_random_inits` over that fold.
///
/// @param x Numeric matrix. The processed data (no whitening function has
/// been applied yet.)
/// @param no_comp Integer. Number of components to test for.
/// @param no_random_init Integer. Number of random initialisations.
/// @param no_folds Integer. Number of folds to use for the cross-validation.
/// @param ica_type String. Which type of ICA shall be run.
/// @param random_seed Integer. For reproducibility.
/// @param ica_params A list containing:
/// \itemize{
///   \item maxit - Integer. Maximum number of iterations for ICA.
///   \item alpha - Float. The alpha parameter for the logcosh version of ICA.
///   Should be between 1 to 2.
///   \item max_tol - Maximum tolerance of the algorithm
///   \item verbose - Verbosity of the function, i.e., shall individual iters
///   be shown.
/// }
/// If the list is empty or the expected elements are not found, default values
/// are used.
///
/// @return A list containing:
/// \itemize{
///   \item s_combined - The combined matrices for S. Dimensions are nrows =
///   features; and ncols = ncomp * no_random_init.
///   \item converged - Boolean vector indicating if the respective run reached
///   convergence. Length = no_random_init
/// }
///
/// @export
#[extendr]
fn rs_ica_iters_cv(
    x: RMatrix<f64>,
    no_comp: usize,
    no_folds: usize,
    no_random_init: usize,
    ica_type: &str,
    random_seed: usize,
    ica_params: List,
) -> List {
    let x = r_matrix_to_faer(&x);
    let ica_params = prepare_ica_params(ica_params);

    let (s_combined, converged) = stabilised_ica_cv(
        x,
        no_comp,
        no_folds,
        no_random_init,
        ica_type,
        ica_params,
        None,
        random_seed,
    );

    list!(
        s_combined = faer_to_r_matrix(s_combined.as_ref()),
        converged = converged
    )
}

extendr_module! {
  mod fun_ica;
  fn rs_prepare_whitening;
  fn rs_fast_ica;
  fn rs_ica_iters;
  fn rs_ica_iters_cv;
}
