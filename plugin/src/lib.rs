use serde::{Deserialize, Serialize};
use wasm_minimal_protocol::*;

pub mod algorithm;
pub mod matrices;
pub mod scoring;

pub use matrices::{available_matrices, get_ednafull_matrix, get_matrix_by_name};
use scoring::{Scorer, SimpleScorer, MatrixScorer};
use algorithm::needleman_wunsch_linear;

initiate_protocol!();

// =============================================================================
// Input Structures
// =============================================================================

/// Scoring parameters for matches and mismatches (simple mode).
#[derive(Deserialize, Clone)]
pub struct MatchScores {
    /// Score for matching characters (typically positive, e.g., 5.0)
    #[serde(rename = "match")]
    pub match_score: f32,
    /// Penalty for mismatching characters (typically negative, e.g., -4.0)
    pub mismatch: f32,
}

/// Scoring configuration - either simple match/mismatch or a substitution matrix.
#[derive(Deserialize, Clone)]
#[serde(untagged)]
pub enum ScoringConfig {
    /// Simple match/mismatch scoring
    Simple(MatchScores),
    /// Named substitution matrix (e.g., "EDNAFULL")
    Matrix { matrix: String },
}

/// Input parameters for sequence alignment.
#[derive(Deserialize)]
pub struct AlignmentInput {
    pub seq1: String,
    pub seq2: String,
    /// Scoring configuration: either { "match": 5, "mismatch": -4 } or { "matrix": "EDNAFULL" }
    pub scores: ScoringConfig,
    /// Linear gap penalty (cost per gap character)
    pub gap_penalty: f32,
    /// Maximum number of optimal alignments to return
    pub max_paths: Option<usize>,
}

// =============================================================================
// Output Structures
// =============================================================================

/// A single step in the traceback path.
#[derive(Serialize, Clone, Debug, PartialEq)]
pub struct PathStep {
    /// Row index in the matrix
    pub i: usize,
    /// Column index in the matrix
    pub j: usize,
}

/// A single alignment result with the traceback path.
#[derive(Serialize)]
pub struct Alignment {
    pub aligned_seq1: String,
    pub aligned_seq2: String,
    /// Ordered traceback path from (0,0) to (n,m)
    pub path: Vec<PathStep>,
}

/// Scoring matrices from the alignment algorithm.
#[derive(Serialize)]
pub struct Matrices {
    /// Score matrix
    pub m: Vec<Vec<f32>>,
}

/// Output from the alignment algorithm.
#[derive(Serialize)]
pub struct AlignmentOutput {
    /// The optimal alignment score
    pub score: f32,
    /// List of optimal alignments found
    pub alignments: Vec<Alignment>,
    /// The scoring matrices
    pub matrices: Matrices,
}

// =============================================================================
// WASM Entry Point
// =============================================================================

#[wasm_func]
pub fn align(input_data: &[u8]) -> Result<Vec<u8>, String> {
    let input: AlignmentInput = serde_json::from_slice(input_data)
        .map_err(|e| format!("Failed to parse input: {}", e))?;

    let result = run_alignment(input);

    serde_json::to_vec(&result)
        .map_err(|e| format!("Failed to serialize output: {}", e))
}

/// Run the Needleman-Wunsch alignment with linear gap penalty.
pub fn run_alignment(input: AlignmentInput) -> AlignmentOutput {
    // Create the appropriate scorer based on the scoring configuration
    let scorer: Box<dyn Scorer> = match &input.scores {
        ScoringConfig::Simple(match_scores) => Box::new(SimpleScorer {
            match_score: match_scores.match_score,
            mismatch_penalty: match_scores.mismatch,
        }),
        ScoringConfig::Matrix { matrix } => {
            let mat = get_matrix_by_name(matrix)
                .unwrap_or_else(|| panic!("Unknown substitution matrix: '{}'. Available matrices: {:?}", matrix, available_matrices()));
            Box::new(MatrixScorer {
                matrix: mat,
                default_score: -4.0, // Default penalty for unknown character pairs
            })
        }
    };

    needleman_wunsch_linear(
        &input.seq1,
        &input.seq2,
        scorer.as_ref(),
        input.gap_penalty,
        input.max_paths.unwrap_or(100),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_scoring_config_simple() {
        let input = AlignmentInput {
            seq1: "ACGT".to_string(),
            seq2: "ACGT".to_string(),
            scores: ScoringConfig::Simple(MatchScores {
                match_score: 5.0,
                mismatch: -4.0,
            }),
            gap_penalty: -2.0,
            max_paths: Some(1),
        };

        let result = run_alignment(input);
        assert_eq!(result.score, 20.0);
    }

    #[test]
    fn test_scoring_config_matrix() {
        let input = AlignmentInput {
            seq1: "ACGT".to_string(),
            seq2: "ACGT".to_string(),
            scores: ScoringConfig::Matrix {
                matrix: "EDNAFULL".to_string(),
            },
            gap_penalty: -2.0,
            max_paths: Some(1),
        };

        let result = run_alignment(input);
        assert_eq!(result.score, 20.0); // 4 matches × 5 (EDNAFULL)
    }
}
