#let myplugin = plugin("alignment.wasm")

// Import tiptoe for drawing arrows
#import "@preview/tiptoe:0.4.0": line as tiptoe-line, tikz

/// Perform sequence alignment using Needleman-Wunsch algorithm (linear gap penalty).
///
/// Parameters:
/// - seq1, seq2: The sequences to align
/// - scores: Either a dictionary with "match" and "mismatch" keys, or a dictionary with "matrix" key
///   - Simple scoring: (match: 5.0, mismatch: -4.0)
///   - Matrix scoring: (matrix: "EDNAFULL")
/// - gap: Linear gap penalty (default -2.0)
/// - max_paths: Maximum number of optimal alignments to return
///
/// Available matrices: EDNAFULL
#let align(seq1, seq2, scores: (match: 1.0, mismatch: -1.0), gap: -2.0, max_paths: 100) = {
  let input = json.encode((
    seq1: seq1,
    seq2: seq2,
    scores: scores,
    gap_penalty: gap,
    max_paths: max_paths,
  ))
  let result_bytes = myplugin.align(bytes(input))
  json(result_bytes)
}

/// Draw a score matrix with traceback path highlighting.
///
/// Parameters:
/// - matrix: The score matrix (2D array)
/// - s1, s2: The sequences
/// - path: The unified traceback path (list of {i, j} objects)
/// - cell-size: Fixed square cell size (default 25pt)
/// - highlight-color: Color for highlighting traceback cells (default yellow.lighten(60%))
#let draw-matrix(
  matrix,
  s1,
  s2,
  path: none,
  cell-size: 25pt,
  highlight-color: yellow.lighten(60%),
) = {
  let n = s1.len()
  let m = s2.len()

  // Calculate total grid dimensions (including header row/column)
  let total-cols = m + 2 // label + "-" + seq2 characters
  let total-rows = n + 2 // header + "-" + seq1 characters

  // Grid stroke for score cells only
  let grid-stroke = 0.5pt + luma(200)

  // Use path directly if present (no filtering needed as there is only one matrix)
  let cells-in-matrix = if path != none { path } else { () }

  // Build header row: empty + "-" + seq2 characters (no stroke)
  let header-row = ()
  header-row.push(grid.cell(stroke: none)[])
  header-row.push(grid.cell(stroke: none)[*-*])
  for c in s2.clusters() {
    header-row.push(grid.cell(stroke: none)[*#c*])
  }

  // Build data rows
  let data-rows = ()
  for i in range(n + 1) {
    // Row label: "-" for first row, otherwise sequence character (no stroke)
    let row-label = if i == 0 { [*-*] } else { [*#s1.at(i - 1)*] }
    data-rows.push(grid.cell(stroke: none)[#row-label])

    // Cell values (with stroke)
    for j in range(m + 1) {
      let val = matrix.at(i).at(j)

      // Check if this cell is in the path
      let is-highlighted = cells-in-matrix.any(step => step.i == i and step.j == j)

      let fill-color = if is-highlighted {
        highlight-color
      } else {
        white
      }

      // Format value (handle none if any)
      let display-val = if val == none {
        sym.minus + sym.infinity
      } else {
        str(val)
      }

      data-rows.push(grid.cell(fill: fill-color, stroke: grid-stroke)[
        #set text(8pt)
        #display-val
      ])
    }
  }

  // Generate arrows from consecutive pairs in the path
  // Draw arrows in traceback direction (from end to start, reversing the path order)
  let arrows = ()
  if path != none and path.len() > 1 {
    for idx in range(path.len() - 1) {
      let from = path.at(idx)
      let to = path.at(idx + 1)
      // Reverse direction: draw from 'to' to 'from' to show traceback direction
      arrows.push((to.i, to.j, from.i, from.j))
    }
  }

  // Calculate cell center position
  // Cell (i, j) in the matrix is at grid position (row: i+1, col: j+1) due to headers
  let cell-center(i, j) = {
    let x = (j + 1) * cell-size + cell-size / 2
    let y = (i + 1) * cell-size + cell-size / 2
    (x, y)
  }

  // Total box dimensions
  let box-width = total-cols * cell-size
  let box-height = total-rows * cell-size

  // Create the visualization
  box(width: box-width, height: box-height)[
    // First layer: the grid with cell values
    #grid(
      columns: (cell-size,) * total-cols,
      rows: (cell-size,) * total-rows,
      align: center + horizon,
      stroke: none, // Default no stroke, individual cells set their own
      ..header-row,
      ..data-rows,
    )

    // Second layer: traceback arrows overlaid on top using tiptoe
    #for arrow in arrows {
      let (from-i, from-j, to-i, to-j) = arrow

      // Get cell centers
      let (from-x, from-y) = cell-center(from-i, from-j)
      let (to-x, to-y) = cell-center(to-i, to-j)

      // Calculate arrow direction vector
      let dx = to-x - from-x
      let dy = to-y - from-y

      // Shorten arrows to 50% of center-to-center length, centered on midpoint
      let mid-x = (from-x + to-x) / 2
      let mid-y = (from-y + to-y) / 2
      let arrow-start-x = mid-x - dx * 0.25
      let arrow-start-y = mid-y - dy * 0.25
      let arrow-end-x = mid-x + dx * 0.25
      let arrow-end-y = mid-y + dy * 0.25

      // Draw arrow using tiptoe with tikz arrow head
      place(
        top + left,
        tiptoe-line(
          start: (arrow-start-x, arrow-start-y),
          end: (arrow-end-x, arrow-end-y),
          stroke: 0.8pt + luma(40),
          tip: tikz,
        ),
      )
    }
  ]
}

// =============================================================================
// Example: Linear Gap Penalty
// =============================================================================

#let s1 = "GATTACA"
#let s2 = "GCATGCU"
#let res_linear = align(s1, s2, scores: (match: 1.0, mismatch: -1.0), gap: -1.0)

= Needleman-Wunsch Sequence Alignment

== Linear Gap Penalty Mode

*Sequences:* #raw(s1) vs #raw(s2) \
*Scoring:* Match = +1, Mismatch = -1, Gap = -1 \
*Score:* #res_linear.score \
*Number of optimal paths:* #res_linear.alignments.len()

#for i in range(calc.min(res_linear.alignments.len(), 3)) {
  let al = res_linear.alignments.at(i)
  [=== Alignment #(i + 1)]

  raw(al.aligned_seq1 + "\n" + al.aligned_seq2)

  draw-matrix(
    res_linear.matrices.m,
    s1,
    s2,
    path: al.path,
  )

  v(10pt)
}

#pagebreak()

// =============================================================================
// Example: Using EDNAFULL Substitution Matrix
// =============================================================================

== Test alignment (EDNAFULL)

#let s1_matrix = "AGATTACA"
#let s2_matrix = "GCCCWTGCWG"
// Use EDNAFULL scoring matrix with linear gap penalty
#let res_matrix = align(s1_matrix, s2_matrix, scores: (matrix: "EDNAFULL"), gap: -2.0)

*Score:* #res_matrix.score \
*Number of optimal paths:* #res_matrix.alignments.len()

// First alignment
#let al = res_matrix.alignments.at(0)
#raw(al.aligned_seq1 + "\n" + al.aligned_seq2, block: true)

#draw-matrix(
  res_matrix.matrices.m,
  s1_matrix,
  s2_matrix,
  path: al.path,
)

// Check total score matches sum of cell values (for verification)
#let total = 0.0
#for step in al.path {
  // Skip the first step (0,0) score is always 0.0, but loop covers it.
  // Wait, the path includes (0,0), and matrix[0][0] is 0.0.
  // But standard scoring is sum of transitions.
  // The matrix cell values accumulate the score. The last cell is the total score.
  // So we don't need to sum them again. The debug check in previous code was summing cell values which is wrong for NW verification,
  // but useful if checking path integrity.
  // I'll leave the visualization as the main output.
}
