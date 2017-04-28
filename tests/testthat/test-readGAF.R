test_that("readGAF() works", {
  d<-readGAF("gene_association_head.goa_ref_human.gz")
  expect_equal(d@numberOfItems, 184)
  expect_equal(dim(d@setAnnotations)[1], 1574)
})

test_that("readGAF() with aspect works", {
  d<-readGAF("gene_association_head.goa_ref_human.gz",aspect="C")
  expect_equal(d@numberOfItems, 149)
})
