test_that("quadrant probabilities are a proper distribution", {

  m <- list(muB = log(0.95), muF = log(0.72),
            sigmaB = 0.174, sigmaF = 0.384, rho = -0.652)

  p <- vapply(c("green", "red", "yellow", "orange"),
              function(q) do.call(kobeGreen, c(m, quadrant = q)),
              numeric(1))

  expect_true(all(p >= 0 & p <= 1))
  expect_equal(sum(p), 1, tolerance = 1e-6)
})

test_that("the 2026 North Atlantic shortfin mako Kobe is reproduced", {

  # moments back-calculated from the reported 2024 marginal probabilities;
  # target values are the published point estimates
  m <- list(muB = log(0.95), muF = log(0.72),
            sigmaB = 0.174, sigmaF = 0.384, rho = -0.652)

  expect_equal(do.call(kobeGreen,  c(m, quadrant = "green")),  0.371,
               tolerance = 5e-4)
  expect_equal(do.call(kobeGreen,  c(m, quadrant = "red")),    0.183,
               tolerance = 5e-4)
  expect_equal(do.call(kobeGreen,  c(m, quadrant = "yellow")), 0.433,
               tolerance = 5e-4)
  expect_equal(do.call(kobeGreen,  c(m, quadrant = "orange")), 0.013,
               tolerance = 5e-4)
})

test_that("degenerate and infeasible moments return NA rather than 0 or 1", {

  expect_true(is.na(kobeGreen(log(0.95), log(0.72), 0, 0.384, -0.652)))
  expect_true(is.na(kobeGreen(log(0.95), log(0.72), 0.174, 0.384, -1.2)))
  expect_true(is.na(kobeGreen(NA, log(0.72), 0.174, 0.384, -0.652)))
})

test_that("kobeGreen is monotone in the biomass and F means", {

  p <- sapply(log(c(0.6, 0.8, 1.0, 1.2)),
              function(mb) kobeGreen(mb, log(0.72), 0.174, 0.384, -0.652))
  expect_true(all(diff(p) > 0))

  p <- sapply(log(c(0.5, 0.7, 0.9, 1.1)),
              function(mf) kobeGreen(log(0.95), mf, 0.174, 0.384, -0.652))
  expect_true(all(diff(p) < 0))
})

test_that("rebuilding time uses the last crossing, not the first", {

  # transient dip, as in an age-structured rebuilding projection
  prob <- data.frame(year  = c(2027, 2030, 2035, 2040, 2045, 2050, 2070),
                     green = c(0.50, 0.44, 0.55, 0.80, 0.92, 0.97, 1.00))

  r60 <- kobeRebuildTime(prob, 0.60)
  expect_equal(r60$rebuildYear, 2036.0, tolerance = 1e-6)
  expect_equal(r60$nCrossings, 1L)
  expect_equal(r60$dPdt, 0.05, tolerance = 1e-9)
  expect_true(r60$attained)

  # at 0.50 the series crosses twice; a first-crossing rule would return 2027
  r50 <- kobeRebuildTime(prob, 0.50)
  expect_equal(r50$nCrossings, 2L)
  expect_gt(r50$rebuildYear, 2032)
  expect_lt(r50$rebuildYear, 2033)
})

test_that("an unattained target is reported as such, not as the last year", {

  prob <- data.frame(year  = c(2027, 2040, 2070),
                     green = c(0.35, 0.30, 0.36))

  r <- kobeRebuildTime(prob, 0.60)
  expect_false(r$attained)
  expect_true(is.na(r$rebuildYear))
})

test_that("a series already above the target rebuilds at its first year", {

  prob <- data.frame(year = c(2027, 2040, 2070), green = c(0.7, 0.8, 0.9))

  r <- kobeRebuildTime(prob, 0.60)
  expect_true(r$attained)
  expect_equal(r$rebuildYear, 2027)
  expect_equal(r$nCrossings, 0L)
})

test_that("a flat approach to the target warns about conditioning", {

  prob <- data.frame(year  = seq(2027, 2070, by = 1),
                     green = seq(0.58, 0.62, length.out = 44))

  expect_warning(kobeRebuildTime(prob, 0.60), "poorly determined")
})

test_that("kobeRebuildTime rejects malformed input", {

  expect_error(kobeRebuildTime(data.frame(year = 2027, p = 0.5)), "columns")
  expect_error(kobeRebuildTime(data.frame(year = 2027, green = 0.5), target = 0),
               "between 0 and 1")
  expect_error(kobeRebuildTime(data.frame(year = 2027, green = 0.5), target = 1),
               "between 0 and 1")
})

test_that("kobeMomentGrid requires catch levels as names", {

  expect_error(kobeMomentGrid(list(a = NULL), 2027), "catch levels")
})

test_that("interpolating moments beats interpolating probabilities", {

  # a smooth surface for which the exact answer is known
  moments <- function(catch, year) {
    t <- year - 2025
    c(muB = -0.05 + 0.030 * t - 0.10 * exp(-((t - 7) / 5)^2) -
             catch / 1000 * 0.022 * t / (1 + t / 25),
      muF = log(0.72) - 0.010 * t + catch / 1000 * 0.45)
  }

  years <- 2025:2070
  catch <- c(0, 250, 1000, 1500, 1883)

  grid <- do.call(rbind, lapply(catch, function(c0)
    do.call(rbind, lapply(years, function(y) {
      m <- moments(c0, y)
      data.frame(catch = c0, year = y, muB = m["muB"], muF = m["muF"],
                 sigmaB = 0.174, sigmaF = 0.384, rhoLog = -0.652,
                 feasible = TRUE)
    }))))

  target <- 1100
  for (y in c(2040, 2050)) {

    m     <- moments(target, y)
    truth <- kobeGreen(m["muB"], m["muF"], 0.174, 0.384, -0.652)

    viaMoments <- kobeInterp(grid, target, y)$green

    g <- grid[grid$year == y, ]
    pr <- mapply(kobeGreen, g$muB, g$muF, g$sigmaB, g$sigmaF, g$rhoLog)
    viaProbs <- stats::splinefun(g$catch, pr, method = "monoH.FC")(target)

    expect_lt(abs(viaMoments - truth), abs(viaProbs - truth))
  }
})

test_that("the two exchange-rate estimates agree on a smooth surface", {

  years <- 2025:2070
  catch <- c(0, 250, 500, 1000, 1500, 1883)

  grid <- do.call(rbind, lapply(catch, function(c0)
    do.call(rbind, lapply(years, function(y) {
      t <- y - 2025
      data.frame(catch = c0, year = y,
                 muB = -0.05 + 0.030 * t - c0 / 1000 * 0.022 * t,
                 muF = log(0.72) - 0.010 * t + c0 / 1000 * 0.45,
                 sigmaB = 0.174, sigmaF = 0.384, rhoLog = -0.652,
                 feasible = TRUE)
    }))))

  x <- kobeExchangeRate(grid, catch = c(250, 500), target = 0.60)

  expect_true(all(is.finite(x$dTdC)))
  expect_true(all(x$dTdC > 0))            # more catch delays rebuilding
  expect_equal(x$dTdC, x$dTdCimplicit, tolerance = 0.25)
})

test_that("pAbove reproduces the 2026 mako marginals", {

  expect_equal(pAbove(0.95, 0.16656), 0.3841, tolerance = 1e-3)
  # sigma_log = 0.384 implies CV = sqrt(exp(sigma^2) - 1)
  seF <- 0.72 * sqrt(exp(0.384^2) - 1)
  expect_equal(pAbove(0.72, seF), 0.196, tolerance = 5e-3)
})

test_that("makeRebuildRefs returns RebuildAge response metrics", {

  # recovery time falls from ~2 GT to ~0.2 GT across starting biomass
  rt <- data.frame(initial = seq(0.1, 0.9, by = 0.1),
                   TRecover = 20 * (1 - seq(0.1, 0.9, by = 0.1)))
  refs <- makeRebuildRefs(rt, gt = 10)

  expect_true(all(c("Tlim", "Ttrig", "Bgt", "Bgt50", "gt") %in% names(refs)))
  expect_equal(refs$gt, 10)
  expect_true(is.finite(refs$Tlim))
  expect_true(refs$Tlim > refs$Ttrig)   # lower start -> longer rebuild
  expect_true(is.finite(refs$Bgt) && is.finite(refs$Bgt50))
  expect_true(refs$Bgt50 > refs$Bgt)    # faster recovery needs higher start
})

test_that("makeK2SM builds the three Liniers tables from a moment grid", {

  years <- c(2027, 2040, 2070)
  catch <- c(0, 250, 1000)
  grid <- do.call(rbind, lapply(catch, function(c0)
    do.call(rbind, lapply(years, function(y) {
      t <- y - 2025
      data.frame(catch = c0, year = y,
                 muB = -0.05 + 0.03 * t - c0 / 1000 * 0.02 * t,
                 muF = log(0.72) + c0 / 1000 * 0.4,
                 sigmaB = 0.174, sigmaF = 0.384, rhoLog = -0.652,
                 feasible = TRUE)
    }))))

  k2 <- makeK2SM(grid, years = years, asPercent = TRUE)

  expect_equal(names(k2), c("pNoOverfishing", "pNotOverfished", "green", "long"))
  expect_equal(nrow(k2$green), 3L)
  expect_true(all(k2$green$catch == catch))
  expect_true(all(k2$green[["2070"]] >= k2$green[["2027"]] |
                    k2$green$catch == 1000))
})
