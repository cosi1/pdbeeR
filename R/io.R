##
#' Read PDB file
#'
#' Reads a .pdb file into a pdb object.
#'
#' @param filename name of the input .pdb file
#' @return pdb object (see description below)
#' @field header data frame with PDB header
#' @field crystal data frame with crystallographic data
#' @field coord data frame with coordinate data
#' @field footer data frame with PDB footer
#' @export
##
read.pdb = function(filename) {
  f = file(filename, open = "r")
  on.exit(close(f))
  pdb_data = list(header = NULL, crystal = NULL, coord = NULL, footer = NULL)

  lines = readLines(f, warn = FALSE)
  rec_names = gsub(" ", "", substring(lines, 1, 6), fixed = TRUE)
  shorter_rec_names = substring(rec_names, 1, 5)
  header_ids = rec_names %in%
    c("HEADER", "OBSLTE", "TITLE", "SPLT", "CAVEAT", "COMPND",
    "SOURCE", "KEYWDS", "EXPDTA", "NUMMDL", "MDLTYP", "AUTHOR",
    "REVDAT", "SPRSDE", "JRNL", "REMARK",
    "DBREF", "DBREF1", "DBREF2", "SEQADV", "SEQRES", "MODRES",
    "HET", "HETNAM", "HETSYN", "FORMUL", "HELIX", "SHEET",
    "SSBOND", "LINK", "CISPEP", "SITE")
  crystal_ids = shorter_rec_names %in%
    c("CRYST", "ORIGX", "SCALE", "MTRIX")
  coord_ids = rec_names %in%
    c("MODEL", "ATOM", "HETATM", "TER", "ENDMDL")
  footer_ids = rec_names %in%
    c("CONECT", "MASTER", "END")
  header_lines = cbind(lines[header_ids], rec_names[header_ids])
  crystal_lines = cbind(lines[crystal_ids], rec_names[crystal_ids])
  coord_lines = cbind(lines[coord_ids], rec_names[coord_ids])
  footer_lines = cbind(lines[footer_ids], rec_names[footer_ids])

  header = t(apply(header_lines, 1, function(l) {
    c(l[2], substring(l[1], c(7, 11), c(10, 80)))
  }))
  if (length(header) > 0) {
    pdb_data$header = as.data.frame(header, stringsAsFactors = FALSE)
    rm(header)
    names(pdb_data$header) = c("rec_name", "num", "value")
  }

  crystal = t(apply(crystal_lines, 1, function(l) {
    c(l[2], substring(l[1], 7, 80))
  }))
  if (length(crystal) > 0) {
    pdb_data$crystal = as.data.frame(crystal, stringsAsFactors = FALSE)
    rm(crystal)
    names(pdb_data$crystal) = c("rec_name", "value")
  }

  START_COORD_COLUMNS =
    c(7, 13, 17, 18, 22, 23, 27, 31, 39, 47, 55, 61, 77, 79)
  END_COORD_COLUMNS =
    c(11, 16, 17, 20, 22, 26, 27, 38, 46, 54, 60, 66, 78, 80)

  coord = as.vector(apply(coord_lines, 1, function(l) {
    c(l[2], gsub(" +", "", substring(l[1],
    START_COORD_COLUMNS, END_COORD_COLUMNS)))
  }))
  if (length(coord) > 0) {
    pdb_data$coord = as.data.frame(matrix(coord, ncol = 15, byrow = TRUE),
      stringsAsFactors = FALSE)
    rm(coord)
    names(pdb_data$coord) = c("rec_name", "serial", "name", "altLoc",
      "resName", "chainID", "resSeq", "iCode", "x", "y", "z", "occupancy",
      "tempFactor", "element", "charge")
    for (i in 9:13) {
      pdb_data$coord[, i] = as.double(pdb_data$coord[, i])
    }
    for (i in c(2, 7)) {
      pdb_data$coord[, i] = as.integer(pdb_data$coord[, i])
    }
  }

  footer = t(apply(footer_lines, 1, function(l) {
    c(l[2], substring(l[1], 7, 80))
  }))
  if (length(footer) > 0) {
    pdb_data$footer = as.data.frame(footer, stringsAsFactors = FALSE)
    rm(footer)
    names(pdb_data$footer) = c("rec_name", "value")
  }

  class(pdb_data) = "pdb"
  pdb_data
}


##
#' Write PDB file
#'
#' Writes a pdb object to a file.
#'
#' @param data pdb object
#' @param filename name of the output .pdb file
#' @export
##
write.pdb = function(data, filename) {
  if (class(data) != "pdb") {
    stop("`data` must be a pdb object.")
  }
  f = file(filename, open = "w")
  on.exit(close(f))

  if (!is.null(data$header)) {
    data$header[, 1] = format(data$header[, 1], width = 6)
    data$header[, 2] = format(data$header[, 2], width = 4)
    data$header[, 3] = format(data$header[, 3], width = 70)
    writeLines(apply(data$header, 1, paste0, collapse = ""), f)
  }

  if (!is.null(data$crystal)) {
    data$crystal[, 1] = format(data$crystal[, 1], width = 6)
    data$crystal[, 2] = format(data$crystal[, 2], width = 74)
    writeLines(apply(data$crystal, 1, paste0, collapse = ""), f)
  }

  if (!is.null(data$coord)) {
    data$coord[, 1] = format(data$coord[, 1], width = 6)
    data$coord[, 2] = format(data$coord[, 2], width = 5, justify = "right")
    data$coord[, 3] = sprintf("%5s", format(data$coord[, 3], width = 3))
    data$coord[, 4] = format(data$coord[, 4], width = 1)
    data$coord[, 5] = format(data$coord[, 5], width = 3, justify = "right")
    data$coord[, 6] = format(data$coord[, 6], width = 2, justify = "right")
    data$coord[, 7] = format(data$coord[, 7], width = 4, justify = "right")
    data$coord[, 8] = format(data$coord[, 8], width = 4)
    data$coord[, 9:11] = apply(data$coord[, 9:11], 2, function(column) {
      column = sprintf("%8.3f", column)
      column[column == "      NA"] = "        "
      column
    })
    data$coord[, 12:13] = apply(data$coord[, 12:13], 2, function(column) {
      column = sprintf("%6.2f", column)
      column[column == "    NA"] = "      "
      column
    })
    data$coord[, 14] = format(data$coord[, 14], width = 12, justify = "right")
    data$coord[, 15] = format(data$coord[, 15], width = 2)
    writeLines(apply(data$coord, 1, paste0, collapse = ""), f)
  }

  if (!is.null(data$footer)) {
    data$footer[, 1] = format(data$footer[, 1], width = 6)
    data$footer[, 2] = format(data$footer[, 2], width = 74)
    writeLines(apply(data$footer, 1, paste0, collapse = ""), f)
  }
}
