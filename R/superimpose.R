##
#' Superimpose structures
#'
#' Superimposes the mobile structure on the target structure
#' using Kabsch algorithm.
#'
#' @param mobile coord data of the mobile structure
#' @param target_subset,mobile_subset atom subsets to be aligned
#'     (subsets should be the same length and in the same order)
#' @return Coord data of the superimposed structure
#' @references
#' Kabsch, Wolfgang, (1976), A solution for the best rotation to relate
#' two sets of vectors. \emph{Acta Crystallographica} \strong{32}:922.
#' \url{https://dx.doi.org/10.1107/S0567739476001873}
#' @examples
#' \dontrun{
#' superimpose(pdb2$coord,
#'             selectRecords(pdb1, name = "P", resSeq = 1:10),
#'             selectRecords(pdb2, name = "P", resSeq = 2:11))
#' }
#' @export
##
superimpose = function(mobile, target_subset, mobile_subset) {
  # Extract XYZ coordinates from coord data
  target_subset_xyz = target_subset[, 9:11]
  mobile_xyz = mobile[, 9:11]
  mobile_subset_xyz = mobile_subset[, 9:11]

  # Center both atom subsets (move their centroids to [0, 0, 0])
  target_centroid = apply(as.matrix(target_subset_xyz), 2, mean)
  mobile_centroid = apply(as.matrix(mobile_subset_xyz), 2, mean)
  coord0 = t(apply(target_subset_xyz, 1, `-`, target_centroid))
  coord1 = t(apply(mobile_subset_xyz, 1, `-`, mobile_centroid))

  # Calculate rotation matrix using Kabsch method
  cov_matrix = t(coord0) %*% coord1
  svd_ = svd(cov_matrix)
  V = svd_$u
  W = svd_$v
  d = sign(det(W %*% t(V)))
  m = matrix(0, 3, 3)
  diag(m) = 1
  m[3, 3] = d
  U = W %*% m %*% t(V)

  # Center all mobile atoms based on the *mobile subset* centroid
  mobile_xyz_centered = t(apply(mobile_xyz, 1, `-`, mobile_centroid))
  # Rotate all mobile atoms
  mobile_xyz_rotated = t(apply(mobile_xyz_centered, 1, `%*%`, U))
  # Move all mobile atoms to the *target subset* centroid
  mobile_xyz_translated = t(apply(mobile_xyz_rotated, 1, `+`, target_centroid))

  mobile[, 9:11] = mobile_xyz_translated
  mobile
}

