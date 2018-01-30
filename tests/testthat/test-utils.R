context("utils")

test_that("logistic function works", {
    expect_equal(logistic(0, x0 = 0, k = 1), 0.5)
})

test_that("rbindMatched works", {
    df1 <- data.frame(A = 1:3, B = 4:6, C = 7:9)
    df2 <- data.frame(D = 0)
    expect_error(rbindMatched(df1, df2),
                 "There must be at least two columns in common")
    df2 <- data.frame(A = 1:3)
    expect_error(rbindMatched(df1, df2),
                 "There must be at least two columns in common")
    df2 <- data.frame(A = 1:3, C = 7:9, D = 0)
    expect_equal(colnames(rbindMatched(df1, df2)), c("A", "C"))
})

test_that("winsorize works", {
    x <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2)
    expect_true(all(winsorize(x, q = 0.1) == 1))
})
