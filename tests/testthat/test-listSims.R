context("listSims")

test_that("listSims printing works", {
    expect_output(listSims(), "Splatter currently contains")
})

test_that("listSims return works", {
    ll <- listSims(print = FALSE)
    checkmate::expect_class(ll, "data.frame")
})
