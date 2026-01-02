use std::collections::HashMap;

/// Trait for scoring character pairs in sequence alignment.
pub trait Scorer {
    fn score(&self, a: char, b: char) -> f32;
}

/// Simple match/mismatch scoring.
pub struct SimpleScorer {
    pub match_score: f32,
    pub mismatch_penalty: f32,
}

impl Scorer for SimpleScorer {
    fn score(&self, a: char, b: char) -> f32 {
        if a.to_ascii_uppercase() == b.to_ascii_uppercase() {
            self.match_score
        } else {
            self.mismatch_penalty
        }
    }
}

/// Substitution matrix-based scoring.
pub struct MatrixScorer {
    pub matrix: HashMap<(char, char), f32>,
    /// Default score for pairs not in the matrix
    pub default_score: f32,
}

impl Scorer for MatrixScorer {
    fn score(&self, a: char, b: char) -> f32 {
        let a_upper = a.to_ascii_uppercase();
        let b_upper = b.to_ascii_uppercase();
        *self.matrix.get(&(a_upper, b_upper)).unwrap_or(&self.default_score)
    }
}
