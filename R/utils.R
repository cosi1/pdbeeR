##
#' resSeq IDs
#'
#' Returns resSeq IDs for a coordinate object.
#'
#' @param coord coordinate object (data frame)
#' @return Vector of resSeq IDs
#' @examples
#' \dontrun{
#' getResSeq(selectRecords(pdb, chainID = "B"))
#' }
#' @export
##
getResSeq = function(coord) {
  if (class(coord) != "data.frame") {
    stop("`coord` must be a data frame.")
  }
  unique(coord$resSeq)
}


##
#' Get sequence
#'
#' Returns the sequence of a coordinate object (as a vector).
#'
#' @param coord coordinate object (data frame)
#' @return Vector of residue names
#' @examples
#' \dontrun{
#' getSequence(selectRecords(pdb, chainID = "A"))
#' }
#' @export
##
getSequence = function(coord) {
  if (class(coord) != "data.frame") {
    stop("`coord` must be a data frame.")
  }
  coord[!duplicated(coord$resSeq),]$resName
}


##
#' Convert sequence to FASTA
#'
#' Converts sequence vector to a sequence string in FASTA format
#'
#' @param sequence vector containing valid residue codes
#' @param name name for the output sequence
#' @param wrapAt column, at which the sequence is wrapped
#'     (set to 0 for unwrapped output)
#' @return String in FASTA format
#' @examples
#' \dontrun{
#' chainB = selectRecords(pdb, chainID = "B")
#' sequence = getSequence(chainB)
#' seqToFasta(sequence)
#' }
#' @export
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
#' Select records
#'
#' Returns a subset of a pdb object's coord data.
#'
#' @param pdb pdb object
#' @param ... selection criteria (see Examples)
#' @return Coordinate object (data frame) containing the selected subset
#' @examples
#' \dontrun{
#' selectRecords(pdb, resSeq = 1:10)
#' selectRecords(pdb, element = "C", resSeq = c(20, 22:25))
#' selectRecords(pdb, name = "P", resName = c("G", "C"), chainID = "A")
#' }
#' @export
##
selectRecords = function(pdb, ...) {
  if (class(pdb) != "pdb") {
    stop("`pdb` must be a pdb object.")
  }
  l = eval(substitute(alist(...)))
  l = lapply(l, function(x) if (is.character(x)) sprintf('"%s"', x) else x)
  e = sapply(names(l), function(name) {
    paste(c(name, "%in%", l[[name]]), collapse = " ")})
  e = paste(e, collapse = " & ")
  e = parse(text = e)
  s = substitute(subset(pdb$coord, e))
  eval(s)
}
