//! Star patterns (quads and triangles), their matching, and match filtering.

pub mod build;
pub mod r#match;
pub mod tetra;
pub mod vote;

pub use build::{build_quads, build_quads_presorted};
pub use r#match::{
    CatalogCodes, QuadGrid, QuadMatch, extract_star_pairs, filter_by_scale, find_matches,
    find_matches_indexed, find_matches_sorted,
};
pub use tetra::{
    TETRA_TOL_FACTOR, TriMatch, TriangleList, bijective_filter, build_triangles,
    extract_triangle_pairs, filter_triangles_by_scale, find_triangle_matches,
};
pub use vote::vote_filter;
