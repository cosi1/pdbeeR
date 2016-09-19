## pdbeeR - tools for manipulating PDB data
## version 0.0.0.9000
## author: Paweł Piątkowski <pawel.piatkowski@posteo.net>


##
# Reads a .pdb file into a pdb object.
# Structure:
#   $header - data frame with PDB header
#   $crystal - data frame with crystallographic data
#   $coord  - data frame with coordinate data
#   $footer - data frame with PDB footer
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

  pdb_data
}


##
# Writes a pdb object to a file.
##
write.pdb = function(data, filename) {
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


##
# Returns resSeq IDs for a coord object.
# Example:
#   getResSeq(selectRecords(pdb, chainID = "B"))
##
getResSeq = function(coord) {
  unique(coord$resSeq)
}


##
# Returns the sequence of a coord object (as a vector).
# Example:
#   getSequence(selectRecords(pdb, chainID = "A"))
##
getSequence = function(coord) {
  coord[!duplicated(coord$resSeq),]$resName
}


##
# Converts sequence vector to a sequence string in FASTA format
# Input:
#   sequence - vector containing valid residue codes
#   name - name for the output sequence
#   wrapAt - column, at which the sequence is wrapped
#     (set to 0 for unwrapped output)
# Example:
#   chainB = selectRecords(pdb, chainID = "B")
#   sequence = getSequence(chainB)
#   seqToFasta(sequence)
##
seqToFasta = function(sequence, name = "Sequence", wrapAt = 80) {
  longCodes = c("ALA", "CYS", "ASP", "GLU", "PHE",
                "GLY", "HIS", "ILE", "LYS", "LEU",
                "MET", "ASN", "PRO", "GLN", "ARG",
                "SER", "THR", "VAL", "TRP", "TYR",
                "HOH")
  shortCodes = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                 "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y",
                 "")

  mapply(function(long, short) {
    sequence <<- gsub(long, short, sequence)
  }, longCodes, shortCodes)
  sequence[nchar(sequence) > 1] = "x"
  str_seq = paste0(sequence, collapse = "")
  if (wrapAt > 0) {
    chunks = sapply(seq(1, nchar(str_seq), wrapAt),
      function(i) substr(str_seq, i, i + wrapAt - 1))
    str_seq = paste0(chunks, collapse = "\n")
  }
  sprintf(">%s\n%s", name, str_seq)
}


##
# Returns a subset of a pdb object's coord data.
# Examples:
#   selectRecords(pdb, resSeq = 1:10)
#   selectRecords(pdb, element = "C", resSeq = c(20, 22:25))
#   selectRecords(pdb, name = "P", resName = c("G", "C"), chainID = "A")
##
selectRecords = function(pdb, ...) {
  l = eval(substitute(alist(...)))
  l = lapply(l, function(x) if (is.character(x)) sprintf('"%s"', x) else x)
  e = sapply(names(l), function(name) {
    paste(c(name, "%in%", l[[name]]), collapse = " ")})
  e = paste(e, collapse = " & ")
  e = parse(text = e)
  s = substitute(subset(pdb$coord, e))
  eval(s)
}


##
# Superimposes the mobile structure on the target structure.
# Input:
#   mobile - coord data of the mobile structure
#   target_subset, mobile_subset - atom subsets to be aligned
#     (subsets should be the same length and in the same order)
# Output:
#   coord data of the superimposed structure
# Example:
#   superimpose(pdb2$coord,
#               selectRecords(pdb1, name = "P", resSeq = 1:10),
#               selectRecords(pdb2, name = "P", resSeq = 2:11))
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

