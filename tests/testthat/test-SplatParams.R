context("SplatParams")

test_that("nGroups checks work", {
  params <- newSplatParams()
  expect_error(setParam(params, "nCells", 1),
               "nCells cannot be set directly, set groupCells instead")
})
