//! Substitution matrices for sequence alignment.
//!
//! This module contains pre-defined substitution matrices for nucleotide alignments.

use std::collections::HashMap;

/// Get the EDNAFULL (also known as NUC4.4) substitution matrix for DNA/RNA.
///
/// This matrix was created by Todd Lowe (12/10/92) and uses ambiguous nucleotide codes.
/// Probabilities are rounded to the nearest integer.
/// Lowest score = -4.0, Highest score = 5.0
pub fn get_ednafull_matrix() -> HashMap<(char, char), f32> {
    let chars = [
        'A', 'T', 'G', 'C', 'S', 'W', 'R', 'Y', 'K', 'M', 'B', 'V', 'H', 'D', 'N',
    ];

    // Matrix values row by row (from the EDNAFULL file)
    #[rustfmt::skip]
    let values: [[f32; 16]; 16] = [
        [ 5.0, -4.0, -4.0, -4.0, -4.0,  1.0,  1.0, -4.0, -4.0,  1.0, -4.0, -1.0, -1.0, -1.0, -2.0, -4.0], // A
        [-4.0,  5.0, -4.0, -4.0, -4.0,  1.0, -4.0,  1.0,  1.0, -4.0, -1.0, -4.0, -1.0, -1.0, -2.0,  5.0], // T
        [-4.0, -4.0,  5.0, -4.0,  1.0, -4.0,  1.0, -4.0,  1.0, -4.0, -1.0, -1.0, -4.0, -1.0, -2.0, -4.0], // G
        [-4.0, -4.0, -4.0,  5.0,  1.0, -4.0, -4.0,  1.0, -4.0,  1.0, -1.0, -1.0, -1.0, -4.0, -2.0, -4.0], // C
        [-4.0, -4.0,  1.0,  1.0, -1.0, -4.0, -2.0, -2.0, -2.0, -2.0, -1.0, -1.0, -3.0, -3.0, -1.0, -4.0], // S
        [ 1.0,  1.0, -4.0, -4.0, -4.0, -1.0, -2.0, -2.0, -2.0, -2.0, -3.0, -3.0, -1.0, -1.0, -1.0,  1.0], // W
        [ 1.0, -4.0,  1.0, -4.0, -2.0, -2.0, -1.0, -4.0, -2.0, -2.0, -3.0, -1.0, -3.0, -1.0, -1.0, -4.0], // R
        [-4.0,  1.0, -4.0,  1.0, -2.0, -2.0, -4.0, -1.0, -2.0, -2.0, -1.0, -3.0, -1.0, -3.0, -1.0,  1.0], // Y
        [-4.0,  1.0,  1.0, -4.0, -2.0, -2.0, -2.0, -2.0, -1.0, -4.0, -1.0, -3.0, -3.0, -1.0, -1.0,  1.0], // K
        [ 1.0, -4.0, -4.0,  1.0, -2.0, -2.0, -2.0, -2.0, -4.0, -1.0, -3.0, -1.0, -1.0, -3.0, -1.0, -4.0], // M
        [-4.0, -1.0, -1.0, -1.0, -1.0, -3.0, -3.0, -1.0, -1.0, -3.0, -1.0, -2.0, -2.0, -2.0, -1.0, -1.0], // B
        [-1.0, -4.0, -1.0, -1.0, -1.0, -3.0, -1.0, -3.0, -3.0, -1.0, -2.0, -1.0, -2.0, -2.0, -1.0, -4.0], // V
        [-1.0, -1.0, -4.0, -1.0, -3.0, -1.0, -3.0, -1.0, -3.0, -1.0, -2.0, -2.0, -1.0, -2.0, -1.0, -1.0], // H
        [-1.0, -1.0, -1.0, -4.0, -3.0, -1.0, -1.0, -3.0, -1.0, -3.0, -2.0, -2.0, -2.0, -1.0, -1.0, -1.0], // D
        [-2.0, -2.0, -2.0, -2.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -2.0], // N
    ];

    let mut matrix = HashMap::new();
    for (i, &c1) in chars.iter().enumerate() {
        for (j, &c2) in chars.iter().enumerate() {
            matrix.insert((c1, c2), values[i][j]);
            // Also handle lowercase
            matrix.insert((c1.to_ascii_lowercase(), c2), values[i][j]);
            matrix.insert((c1, c2.to_ascii_lowercase()), values[i][j]);
            matrix.insert(
                (c1.to_ascii_lowercase(), c2.to_ascii_lowercase()),
                values[i][j],
            );
        }
    }
    matrix
}

/// Returns a list of available substitution matrix names.
pub fn available_matrices() -> Vec<&'static str> {
    vec!["EDNAFULL"]
}

/// Get a substitution matrix by name.
///
/// Supported matrix names:
/// - "EDNAFULL"
pub fn get_matrix_by_name(name: &str) -> Option<HashMap<(char, char), f32>> {
    match name.to_uppercase().as_str() {
        "EDNAFULL" => Some(get_ednafull_matrix()),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ednafull_basic_scores() {
        let matrix = get_ednafull_matrix();

        // Check match scores (diagonal)
        assert_eq!(matrix.get(&('A', 'A')), Some(&5.0));
        assert_eq!(matrix.get(&('T', 'T')), Some(&5.0));
        assert_eq!(matrix.get(&('G', 'G')), Some(&5.0));
        assert_eq!(matrix.get(&('C', 'C')), Some(&5.0));

        // Check mismatch scores
        assert_eq!(matrix.get(&('A', 'T')), Some(&-4.0));
        assert_eq!(matrix.get(&('A', 'G')), Some(&-4.0));
        assert_eq!(matrix.get(&('A', 'C')), Some(&-4.0));

        // Check ambiguous codes
        assert_eq!(matrix.get(&('A', 'W')), Some(&1.0)); // W = A or T
        assert_eq!(matrix.get(&('A', 'R')), Some(&1.0)); // R = A or G
        assert_eq!(matrix.get(&('N', 'N')), Some(&-1.0)); // N = any
    }

    #[test]
    fn test_ednafull_case_insensitive() {
        let matrix = get_ednafull_matrix();

        assert_eq!(matrix.get(&('a', 'a')), Some(&5.0));
        assert_eq!(matrix.get(&('A', 'a')), Some(&5.0));
        assert_eq!(matrix.get(&('a', 'A')), Some(&5.0));
        assert_eq!(matrix.get(&('a', 't')), Some(&-4.0));
    }

    #[test]
    fn test_ednafull_symmetry() {
        let matrix = get_ednafull_matrix();

        // Matrix should be symmetric
        for c1 in ['A', 'T', 'G', 'C', 'N'] {
            for c2 in ['A', 'T', 'G', 'C', 'N'] {
                assert_eq!(
                    matrix.get(&(c1, c2)),
                    matrix.get(&(c2, c1)),
                    "Matrix should be symmetric for ({}, {})",
                    c1,
                    c2
                );
            }
        }
    }

    #[test]
    fn test_get_matrix_by_name() {
        assert!(get_matrix_by_name("EDNAFULL").is_some());
        assert!(get_matrix_by_name("ednafull").is_some());
        assert!(get_matrix_by_name("NUC4.4").is_none());
        assert!(get_matrix_by_name("NUC44").is_none());
        assert!(get_matrix_by_name("UNKNOWN").is_none());
    }

    #[test]
    fn test_available_matrices() {
        let matrices = available_matrices();
        assert_eq!(matrices, vec!["EDNAFULL"]);
    }
}
