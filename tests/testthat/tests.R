
test_that("direct loess matches R's loess", {
  n <- 100
  x <- 1:n
  set.seed(34565)
  y <- rnorm(n)

  dc <- loess.control(surface = "direct")

  # we can only test indices 16:85 because for some reason
  # on Windows i386 only, the non-symmetric fits are off

  r1 <- stlplus:::.loess_stlplus(x = x, y = y, span = 31, degree = 0)
  r2 <- predict(loess(y ~ x, degree = 0, span = 31 / n, control = dc))
  expect_true(mean(abs(r1 - r2)[16:85]) < 1.0e-15)

  r1 <- stlplus:::.loess_stlplus(x = x, y = y, span = 31, degree = 1)
  r2 <- predict(loess(y ~ x, degree = 1, span = 31 / n, control = dc))
  expect_true(mean(abs(r1 - r2)[16:85]) < 1.0e-15)

  r1 <- stlplus:::.loess_stlplus(x = x, y = y, span = 31, degree = 2)
  r2 <- predict(loess(y ~ x, degree = 2, span = 31 / n, control = dc))
  expect_true(mean(abs(r1 - r2)[16:85]) < 1.0e-15)
})

test_that("stlplus matches stl", {
  x1 <- stlplus(co2, t = as.vector(time(co2)), n.p = 12,
    l.window = 13, t.window = 19, s.window = 11, s.degree = 1,
    s.jump = 1, t.jump = 1, l.jump = 1)
  x2 <- stl(co2, l.window = 13, t.window = 19, s.window = 11, s.degree = 1,
    s.jump = 1, t.jump = 1, l.jump = 1)
  expect_true(mean(abs(seasonal(x1) - seasonal(x2))) < 1.0e-12)
  expect_true(mean(abs(trend(x1) - trend(x2)))  < 1.0e-12)

  x1 <- stlplus(co2, t = as.vector(time(co2)), n.p = 12,
    l.window = 13, t.window = 39, s.window = 101, s.degree = 1,
    s.jump = 1, t.jump = 1, l.jump = 1)
  x2 <- stl(co2, l.window = 13, t.window = 39, s.window = 101, s.degree = 1,
    s.jump = 1, t.jump = 1, l.jump = 1)
  expect_true(mean(abs(seasonal(x1) - seasonal(x2))) < 1.0e-12)
  expect_true(mean(abs(trend(x1) - trend(x2)))  < 1.0e-12)
})

# compare to loess with missing values
# y[10] <- NA
# r1 <- .loess_stlplus(x = x, y = y, span = 31, degree = 1)
# r2 <- predict(loess(y ~ x, degree = 1, span = 31/n, control = dc), newdata = c(1:n))
# plot(r1-r2)
# stopifnot(mean(abs(r1-r2)) < 1.0e-15)

## doesn't match with loess when span > n (but matches stl())
##---------------------------------------------------------

# r1 <- .loess_stlplus(x = x, y = y, span = 131, degree = 0)
# r2 <- predict(loess(y ~ x, degree = 0, span = 131 / n, control = dc))
# stopifnot(mean(abs(r1-r2)) < 1.0e-15)
#
# r1 <- .loess_stlplus(x = x, y = y, span = 131, degree = 1)
# r2 <- predict(loess(y ~ x, degree = 1, span = 131 / n, control = dc))
# stopifnot(mean(abs(r1-r2)) < 1.0e-15)
#
# r1 <- .loess_stlplus(x = x, y = y, span = 131, degree = 2)
# r2 <- predict(loess(y ~ x, degree = 2, span = 131 / n, control = dc))
# stopifnot(mean(abs(r1-r2)) < 1.0e-15)

## compare to operator
##---------------------------------------------------------

# library(operator)
# n <- 100
# set.seed(5432)
# y <- rnorm(100)
# x1 <- stlplus(y, n.p = 12, l.window = 13, t.window = 19, s.window = 11, s.degree = 1, s.jump = 1, t.jump = 1, l.jump = 1)
# x2 <- stlOp(x1)
# stopifnot(mean(abs(seasonal(x1) - x2$seas$O %*% y)) < 1.0e-12)
# stopifnot(mean(abs(trend(x1) - x2$trend$O %*% y))  < 1.0e-12)
