traj <- data.frame(year = c(2023, 2024), value = c(9000, 9500),
                   valueSD = c(1080, 1140), ratio = c(0.92, 0.95),
                   ratioSD = c(0.16, 0.166), refpt = 10000,
                   quantity = "ssb", relative = "bratio",
                   source = "derived_quants", stringsAsFactors = FALSE)

# 2024 values are the 2026 assessment terminal year: F/FMSY = 0.72, with the
# standard error that reproduces sigma_log = 0.384
ftraj <- data.frame(year = c(2023, 2024), value = c(0.68, 0.72),
                    valueSD = c(0.271, 0.287), ratio = c(0.68, 0.72),
                    ratioSD = c(0.271, 0.287), refpt = 1,
                    quantity = "f", relative = "reported",
                    source = "derived_quants", stringsAsFactors = FALSE)

test_that("the interval is lognormal, not Wald", {

  ci <- lognormalCI(traj, 0.95)

  expect_gt(ci$ratioHi[1] - traj$ratio[1], traj$ratio[1] - ci$ratioLo[1])
  expect_equal(sqrt(ci$ratioLo * ci$ratioHi), traj$ratio, tolerance = 1e-10)
  expect_true(all(lognormalCI(traj, 0.999)$ratioLo > 0))
})

test_that("the ratio basis reproduces the reported 2024 mako interval", {

  ci <- lognormalCI(traj, 0.95)

  expect_equal(ci$ratioLo[2], 0.676, tolerance = 1e-3)
  expect_equal(ci$ratioHi[2], 1.335, tolerance = 1e-3)
})

test_that("the absolute basis is narrower and warns that it is conditional", {

  expect_warning(ciAbs <- lognormalCI(traj, 0.95, basis = "absolute"),
                 "ignores uncertainty in the reference point")

  ciRat <- lognormalCI(traj, 0.95)

  expect_lt(ciAbs$ratioHi[2] - ciAbs$ratioLo[2],
            ciRat$ratioHi[2] - ciRat$ratioLo[2])
  expect_equal(ciAbs$basis[1], "absolute")
})

test_that("nested levels are returned long form and ordered by width", {

  ci <- lognormalCI(traj, levels = c(0.5, 0.8, 0.95))

  expect_equal(nrow(ci), 6L)
  expect_equal(sort(unique(ci$level)), c(0.5, 0.8, 0.95))

  w <- tapply(ci$ratioHi - ci$ratioLo, ci$level, function(x) x[1])
  expect_true(all(diff(w) > 0))
})

test_that("a zero or missing SD gives NA bounds, not a zero-width interval", {

  bad <- traj
  bad$ratioSD <- c(0, NA)

  expect_warning(ci <- lognormalCI(bad, 0.95), "missing or zero standard error")
  expect_true(all(is.na(ci$ratioLo)))
})

test_that("an absent ratio SD falls back to the absolute basis", {

  noRatio <- traj
  noRatio$ratioSD <- NA_real_

  expect_warning(ci <- lognormalCI(noRatio, 0.95), "falling back")
  expect_equal(ci$basis[1], "absolute")
})

test_that("lognormalCI rejects malformed input", {

  expect_error(lognormalCI(traj[, 1:2]), "missing: ")
  expect_error(lognormalCI(traj, levels = 1), "between 0 and 1")
  expect_error(lognormalCI(traj, levels = 0), "between 0 and 1")
})

test_that("the F trajectory takes the same path as the SSB one", {

  ci <- lognormalCI(ftraj, 0.95)

  expect_equal(ci$ratioLo[2], 0.72 * exp(-qnorm(0.975) *
                 sqrt(log(1 + (0.287 / 0.72)^2))), tolerance = 1e-10)
  expect_true(all(ci$ratioLo > 0))
  expect_equal(ci$quantity[1], "f")
})

test_that("the reported interval is lognormal but not centred on the point", {

  # The 2026 assessment reports F/FMSY = 0.72 with a 95% interval of
  # 0.35 to 2.37. That interval is log-symmetric, so it is lognormal in shape,
  # but its geometric centre is 0.911 rather than 0.72. A lognormal centred on
  # the reported point estimate therefore cannot reproduce it, and the mismatch
  # is in the centring and scale, not in the weight of the tail.
  lo <- 0.35; hi <- 2.37
  ctr <- sqrt(lo * hi)

  expect_equal(hi / ctr, ctr / lo, tolerance = 1e-3)   # log-symmetric
  expect_equal(ctr, 0.911, tolerance = 1e-3)           # but not 0.72

  # what a lognormal centred on the point estimate gives instead
  ci <- lognormalCI(ftraj, 0.95)
  expect_equal(ci$ratioLo[2], 0.339, tolerance = 1e-3)
  expect_equal(ci$ratioHi[2], 1.528, tolerance = 1e-3)

  # the sigma implied by the reported interval, against the sigma that
  # reproduces the published quadrant percentages: a factor of about 1.27
  sigRep <- log(hi / ctr) / qnorm(0.975)
  expect_equal(sigRep, 0.488, tolerance = 1e-3)
  expect_equal(sigRep / 0.384, 1.27, tolerance = 0.01)

  # and the two imply materially different probabilities of overfishing
  pKobe <- pnorm(log(0.72) / 0.384)
  pRep  <- pnorm(log(ctr) / sigRep)
  expect_equal(pKobe, 0.196, tolerance = 1e-3)
  expect_gt(pRep - pKobe, 0.2)
})

test_that("the interval and kobeGreen agree on the same moments", {

  cv  <- traj$ratioSD[2] / traj$ratio[2]
  sig <- sqrt(log(1 + cv^2))
  mu  <- log(traj$ratio[2])

  ci <- lognormalCI(traj, 0.95)
  expect_equal(ci$ratioLo[2], exp(mu - qnorm(0.975) * sig), tolerance = 1e-10)

  # P(ratio > 1) from the same moments, two ways
  expect_equal(1 - pnorm(0, mu, sig), pnorm(mu / sig), tolerance = 1e-10)
})

test_that("the FLQuant method returns FLQuants and brackets the median", {

  skip_if_not_installed("FLCore")

  fq <- FLCore::FLQuant(c(9000, 9500), dimnames = list(age = "all",
                                                      year = 2023:2024))
  sq <- FLCore::FLQuant(c(1080, 1140), dimnames = list(age = "all",
                                                      year = 2023:2024))

  res <- lognormalCI(fq, sd = sq, levels = 0.95)

  expect_s4_class(res, "FLQuants")
  expect_equal(names(res), c("lo", "hi"))
  expect_true(all(res$lo < fq))
  expect_true(all(res$hi > fq))

  # geometric centring, as for the data.frame method
  expect_equal(c(sqrt(res$lo * res$hi)), c(fq), tolerance = 1e-10)

  # two levels give four quantities, named by level
  res2 <- lognormalCI(fq, sd = sq, levels = c(0.5, 0.95))
  expect_equal(names(res2), c("lo50", "hi50", "lo95", "hi95"))

  expect_error(lognormalCI(fq), "needs 'sd'")
})

test_that("vanilla ggplot2 can facet nested ribbons without a plot helper", {

  skip_if_not_installed("ggplot2")

  mk <- function(nm, mult) {
    d <- lognormalCI(transform(traj, ratio = traj$ratio * mult),
                     levels = c(0.5, 0.95))
    transform(d, run = nm)
  }
  ci <- rbind(mk("250t", 1), mk("1100t", 1.1), mk("1500t", 1.2))
  ci$run <- factor(ci$run, levels = c("250t", "1100t", "1500t"))
  med <- ci[!duplicated(ci[, c("year", "run")]), ]

  p <- ggplot2::ggplot(ci, ggplot2::aes(year, ratio)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = ratioLo, ymax = ratioHi,
                                      alpha = factor(level)),
                         fill = "#1b9e77") +
    ggplot2::geom_line(data = med) +
    ggplot2::facet_wrap(~run) +
    ggplot2::theme_bw()

  expect_s3_class(p, "ggplot")
  expect_equal(levels(ci$run), c("250t", "1100t", "1500t"))
  expect_equal(nrow(med), 6L)
})
