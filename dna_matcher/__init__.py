"""DNA matcher package for NCBI sequence fetching and comparison algorithms."""

from .core import (
    search_species_options,
    fetch_common_marker_sequences,
    compare_species_matrix,
    calculate_similarity_and_alignment,
    build_neighbor_joining_tree,
)

__all__ = [
    "search_species_options",
    "fetch_common_marker_sequences",
    "compare_species_matrix",
    "calculate_similarity_and_alignment",
    "build_neighbor_joining_tree",
]
