library(pdbeeR)
context("Loading and saving .pdb files")

test_that("file is loaded correctly", {
  pdb = read.pdb("1zci.pdb")
  expect_equal(class(pdb), "pdb")
  expect_equal(nrow(pdb$header), 423)
  expect_equal(nrow(pdb$crystal), 7)
  expect_equal(nrow(pdb$coord), 2275)
  expect_equal(nrow(pdb$footer), 126)
  expect_equal(unique(pdb$coord$chainID), c("A", "B", "C", "D"))
})

test_that("file is saved correctly", {
  pdb = read.pdb("1zci.pdb")
  write.pdb(pdb, "1zci_copy.pdb")
  pdb2 = read.pdb("1zci_copy.pdb")
  expect_equal(pdb2, pdb)
  unlink("1zci_copy.pdb")

  expect_error(write.pdb(pdb$coord, "error.pdb"), "must be a pdb object")
})
