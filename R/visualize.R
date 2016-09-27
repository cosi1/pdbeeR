# Default color table for showStructure()
DEFAULT_COLOR_TABLE = list(
  C = "black",
  O = "red",
  H = "grey",
  N = "blue",
  P = "orange",
  S = "yellow"
)


##
#' Show structure
#'
#' Shows an interactive 3D view of the structure from given coordinates
#' in a graphics device (window, Web browser etc.)
#'
#' @param coord coordinate object (data frame)
#' @param color_table list of colors associated with different atoms
#' @param color_other color associated with atoms not found in color_table
#' @param ... any valid parameter for \link{plot3d}
#' @export
##
showStructure = function(coord, color_table = DEFAULT_COLOR_TABLE,
                         color_other = "pink",
                         type = "s", size = 1.2, box = FALSE, axes = FALSE,
                         xlab = "", ylab = "", zlab = "", ...) {
  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("`rgl` package required for visualization.")
  }
  with(coord, {
    for (el in names(color_table)) {
      element[element == el] = color_table[[el]]
    }
    element[! element %in% color_table] = color_other
    rgl::plot3d(x, y, z, col = element, type = type, size = size, box = box,
                axes = axes, xlab = xlab, ylab = ylab, zlab = zlab, ...)
  })
}
