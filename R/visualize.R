# Default color table (coloring by element name) for showStructure()
ELEMENT_COLOR_TABLE = list(
  C = "black",
  O = "red",
  H = "grey",
  N = "blue",
  P = "orange",
  S = "yellow"
)

# Default color table (coloring by chain name) for showStructure()
dummy_ = 2:37
names(dummy_) = c(LETTERS, 0:9)
CHAIN_COLOR_TABLE = as.list(dummy_)


##
#' Show structure
#'
#' Shows an interactive 3D view of the structure from given coordinates
#' in a graphics device (window, Web browser etc.)
#'
#' @param coord coordinate object (data frame)
#' @param color_by name of the column used for coloring
#' @param color_table list of colors associated with different atoms
#' @param color_other color associated with atoms not found in color_table
#' @param ... any valid parameter for \link{plot3d}
#' @export
##
showStructure = function(coord, color_by = "element", color_other = "pink",
                         color_table = if(color_by == "element")
                         ELEMENT_COLOR_TABLE else CHAIN_COLOR_TABLE,
                         type = "s", size = 1.2, box = FALSE, axes = FALSE,
                         xlab = "", ylab = "", zlab = "", ...) {
  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("`rgl` package required for visualization.")
  }
  for (value in names(color_table)) {
    coord[coord[, color_by] == value, color_by] = color_table[[value]]
  }
  coord[! coord[, color_by] %in% color_table, color_by] = color_other
  with(coord, {
    rgl::plot3d(x, y, z, col = coord[, color_by], type = type, size = size,
                box = box, axes = axes, xlab = xlab, ylab = ylab, zlab = zlab,
                ...)
  })
}
