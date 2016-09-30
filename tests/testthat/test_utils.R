library(pdbeeR)
context("PDB data manipulation")

test_that("subsetting works", {
  pdb = read.pdb("1zci.pdb")
  chainA = selectRecords(pdb, chainID = "A")
  expect_equal(nrow(chainA), 576)
  seq = getSequence(chainA)
  expect_equal(seq[1:5], c("C", "5BU", "U", "G", "C"))
  expect_equal(seqToFasta(getSequence(chainA), "chainA"),
               ">chainA\nCxUGCUGAAGUGCACACAGCAAGK")
  cleaned_chainA = chainA[chainA$resName != "HOH", ]
  expect_equal(getResSeq(cleaned_chainA), 1:24)

  expect_error(selectRecords(chainA, invalid_name = "X"))
})
