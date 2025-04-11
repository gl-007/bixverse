use extendr_api::prelude::*;
use faer::Mat;
use std::collections::{HashMap, HashSet};

/// A double nested HashMap
pub type NestedHashMap = HashMap<String, HashMap<String, HashSet<String>>>;

/// Transforms a Robj List into a Hashmap
pub fn r_list_to_hashmap(r_list: List) -> extendr_api::Result<HashMap<String, Vec<String>>> {
    let mut result = HashMap::with_capacity(r_list.len());

    for (n, s) in r_list {
        let s_vec = s.as_string_vector().ok_or_else(|| {
            Error::Other(format!(
                "Failed to convert value for key '{}' to string vector",
                n
            ))
        })?;
        result.insert(n.to_string(), s_vec);
    }

    Ok(result)
}

/// Transforms a Robj List into a Hashmap with the values as Hashset
pub fn r_list_to_hashmap_set(
    r_list: List,
) -> extendr_api::Result<HashMap<String, HashSet<String>>> {
    let mut result = HashMap::with_capacity(r_list.len());

    for (n, s) in r_list {
        let s_vec = s.as_string_vector().ok_or_else(|| {
            Error::Other(format!(
                "Failed to convert value for key '{}' to string vector",
                n
            ))
        })?;
        let mut s_hash = HashSet::with_capacity(s_vec.len());
        for item in s_vec {
            s_hash.insert(item);
        }
        result.insert(n.to_string(), s_hash);
    }

    Ok(result)
}

/// Transforms a Robj List into an array of String arrays.
pub fn r_list_to_str_vec(r_list: List) -> extendr_api::Result<Vec<Vec<String>>> {
    let mut result = Vec::with_capacity(r_list.len());

    for (n, s) in r_list.into_iter() {
        let s_vec = s.as_string_vector().ok_or_else(|| {
            Error::Other(format!(
                "Failed to convert value to string vector at key '{}'",
                n
            ))
        })?;
        result.push(s_vec);
    }

    Ok(result)
}

/// Transforms a Robj nested list into a nested hashmap
pub fn r_nested_list_to_rust(r_nested_list: List) -> extendr_api::Result<NestedHashMap> {
    let mut result = HashMap::with_capacity(r_nested_list.len());

    for (n, obj) in r_nested_list {
        let inner_list = obj.as_list().ok_or_else(|| {
            Error::Other(format!("Failed to convert value for key '{}' to list", n))
        })?;
        let inner_hashmap = r_list_to_hashmap_set(inner_list)?;
        result.insert(n.to_string(), inner_hashmap);
    }

    Ok(result)
}

/// Transform an R matrix to a Faer one
pub fn r_matrix_to_faer(x: &RMatrix<f64>) -> faer::Mat<f64> {
    let ncol = x.ncols();
    let nrow = x.nrows();
    let data = x.data();

    Mat::from_fn(nrow, ncol, |i, j| data[i + j * nrow])
}

/// Transform a faer into an R matrix
pub fn faer_to_r_matrix(x: faer::MatRef<f64>) -> extendr_api::RArray<f64, [usize; 2]> {
    let nrow = x.nrows();
    let ncol = x.ncols();
    RArray::new_matrix(nrow, ncol, |row, column| x[(row, column)])
}
