use crate::scoring::Scorer;
use crate::{Alignment, AlignmentOutput, Matrices, PathStep};

/// Records traceback directions for a cell, filtering to optimal predecessors.
fn record_optimal_directions(candidates: &[(u8, f32, f32)]) -> u8 {
    if candidates.is_empty() {
        return 0;
    }

    // Find the best transition value (maximum of all transition values)
    let best_transition = candidates.iter().map(|(_, val, _)| *val).reduce(f32::max).unwrap();

    // Filter to candidates that achieve the best transition
    let valid_candidates: Vec<&(u8, f32, f32)> = candidates
        .iter()
        .filter(|(_, val, _)| *val == best_transition)
        .collect();

    // Among valid candidates, find the maximum predecessor score
    let max_pred = valid_candidates
        .iter()
        .map(|(_, _, pred)| *pred)
        .reduce(f32::max)
        .unwrap();

    // Only record directions whose predecessors have the maximum score
    valid_candidates
        .iter()
        .filter(|(_, _, pred)| *pred == max_pred)
        .fold(0u8, |acc, (bit, _, _)| acc | bit)
}

/// Classic Needleman-Wunsch algorithm with linear gap penalty.
pub fn needleman_wunsch_linear(
    seq1: &str,
    seq2: &str,
    scorer: &dyn Scorer,
    gap_penalty: f32,
    max_paths: usize,
) -> AlignmentOutput {
    let seq1: Vec<char> = seq1.chars().collect();
    let seq2: Vec<char> = seq2.chars().collect();
    let n = seq1.len();
    let m = seq2.len();

    // Initialize score matrix
    let mut matrix = vec![vec![0.0f32; m + 1]; n + 1];
    // Direction matrix: bitmask where 1=Diagonal, 2=Up, 4=Left
    let mut directions = vec![vec![0u8; m + 1]; n + 1];

    // =========================================================================
    // Initialization: first row and column
    // =========================================================================

    for i in 1..=n {
        matrix[i][0] = (i as f32) * gap_penalty;
        directions[i][0] = 2; // Up
    }
    for j in 1..=m {
        matrix[0][j] = (j as f32) * gap_penalty;
        directions[0][j] = 4; // Left
    }

    // =========================================================================
    // Matrix Fill
    // =========================================================================

    for i in 1..=n {
        for j in 1..=m {
            let char_score = scorer.score(seq1[i - 1], seq2[j - 1]);

            let diag = matrix[i - 1][j - 1] + char_score;
            let up = matrix[i - 1][j] + gap_penalty;
            let left = matrix[i][j - 1] + gap_penalty;

            let max_val = diag.max(up).max(left);
            matrix[i][j] = max_val;

            // Use unified helper to record optimal directions
            let candidates = &[
                (1, diag, matrix[i - 1][j - 1]),
                (2, up, matrix[i - 1][j]),
                (4, left, matrix[i][j - 1]),
            ];
            directions[i][j] = record_optimal_directions(candidates);
        }
    }

    // =========================================================================
    // Traceback to find all optimal alignments
    // =========================================================================

    let mut alignments = Vec::new();
    find_all_paths_linear(
        n,
        m,
        &directions,
        &seq1,
        &seq2,
        &mut Vec::new(),
        &mut alignments,
        max_paths,
    );

    AlignmentOutput {
        score: matrix[n][m],
        alignments,
        matrices: Matrices {
            m: matrix,
        },
    }
}

/// Recursive traceback for linear gap penalty mode.
fn find_all_paths_linear(
    i: usize,
    j: usize,
    directions: &[Vec<u8>],
    seq1: &[char],
    seq2: &[char],
    current_path: &mut Vec<(usize, usize, char, char)>,
    results: &mut Vec<Alignment>,
    max_paths: usize,
) {
    if results.len() >= max_paths {
        return;
    }

    if i == 0 && j == 0 {
        // Build alignment from path
        let mut aligned1 = String::new();
        let mut aligned2 = String::new();
        let mut path = vec![PathStep { i: 0, j: 0 }];

        for &(r, c, c1, c2) in current_path.iter().rev() {
            aligned1.push(c1);
            aligned2.push(c2);
            path.push(PathStep { i: r, j: c });
        }

        results.push(Alignment {
            aligned_seq1: aligned1,
            aligned_seq2: aligned2,
            path,
        });
        return;
    }

    let dir = directions[i][j];

    // Diagonal (match/mismatch)
    if i > 0 && j > 0 && (dir & 1) != 0 {
        current_path.push((i, j, seq1[i - 1], seq2[j - 1]));
        find_all_paths_linear(i - 1, j - 1, directions, seq1, seq2, current_path, results, max_paths);
        current_path.pop();
    }

    // Up (gap in seq2)
    if i > 0 && (dir & 2) != 0 && results.len() < max_paths {
        current_path.push((i, j, seq1[i - 1], '-'));
        find_all_paths_linear(i - 1, j, directions, seq1, seq2, current_path, results, max_paths);
        current_path.pop();
    }

    // Left (gap in seq1)
    if j > 0 && (dir & 4) != 0 && results.len() < max_paths {
        current_path.push((i, j, '-', seq2[j - 1]));
        find_all_paths_linear(i, j - 1, directions, seq1, seq2, current_path, results, max_paths);
        current_path.pop();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scoring::{SimpleScorer, MatrixScorer};
    use crate::matrices::get_ednafull_matrix;

    fn simple_scorer(match_score: f32, mismatch: f32) -> SimpleScorer {
        SimpleScorer {
            match_score,
            mismatch_penalty: mismatch,
        }
    }

    #[test]
    fn test_linear_identical() {
        let scorer = simple_scorer(5.0, -4.0);
        let result = needleman_wunsch_linear("ACGT", "ACGT", &scorer, -2.0, 10);
        assert_eq!(result.score, 20.0); // 4 matches × 5
        assert_eq!(result.alignments[0].aligned_seq1, "ACGT");
        assert_eq!(result.alignments[0].aligned_seq2, "ACGT");
    }

    #[test]
    fn test_linear_with_gap() {
        let scorer = simple_scorer(5.0, -4.0);
        let result = needleman_wunsch_linear("ACGT", "ACT", &scorer, -2.0, 10);
        // ACT aligned to ACGT: AC-T vs ACGT or ACT- vs ACGT
        // Best: A C G T
        //       A C - T  → 3 matches (15.0) + 1 gap (-2.0) = 13.0
        assert_eq!(result.score, 13.0);
    }

    #[test]
    fn test_linear_with_ednafull() {
        let matrix = get_ednafull_matrix();
        let scorer = MatrixScorer {
            matrix,
            default_score: -4.0,
        };

        let result = needleman_wunsch_linear("ACGT", "ACGT", &scorer, -2.0, 10);
        assert_eq!(result.score, 20.0); // 4 matches × 5
    }
}
