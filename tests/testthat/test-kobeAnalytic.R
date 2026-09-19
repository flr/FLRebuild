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

test_that("degenerate moments return NA; |rho|>=1 is clamped not NA", {

  expect_true(is.na(kobeGreen(log(0.95), log(0.72), 0, 0.384, -0.652)))
  expect_true(is.na(kobeGreen(NA, log(0.72), 0.174, 0.384, -0.652)))
  # Boundary / beyond-boundary correlation: still a probability
  p <- kobeGreen(log(0.95), log(0.72), 0.174, 0.384, -1.2)
  expect_true(is.finite(p) && p >= 0 && p <= 1)
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

test_that("deterministic F = 0 does not warn and yields P(no overfishing) = 1", {

  dq <- data.frame(
    Label  = c("Bratio_2030", "F_2030"),
    Value  = c(0.95, 0),
    StdDev = c(0.16656, 0),
    stringsAsFactors = FALSE)
  covar <- data.frame(Label.i = character(0), Label.j = character(0),
                      corr = numeric(0), stringsAsFactors = FALSE)

  expect_silent(hat <- FLRebuild:::kobeHat(dq, "F_2030"))
  expect_equal(unname(hat["value"]), 0)

  expect_silent(m <- FLRebuild:::kobeMoments(dq, covar, 2030))
  expect_true(is.infinite(m$mu[2]) && m$mu[2] < 0)
  expect_equal(m$sigmaF, 0)
  expect_equal(m$fratio, 0)

  pB <- 1 - stats::pnorm(0, m$mu[1], m$sigmaB)
  expect_equal(kobeGreen(m$mu[1], m$mu[2], m$sigmaB, m$sigmaF, m$rhoLog),
               pB, tolerance = 1e-10)
  expect_equal(FLRebuild:::kobePNoOverfishing(m$mu[2], m$sigmaF), 1)

  grid <- data.frame(catch = 0, year = 2030,
                     muB = m$mu[1], muF = m$mu[2],
                     sigmaB = m$sigmaB, sigmaF = m$sigmaF,
                     rhoLog = m$rhoLog, feasible = TRUE)
  k2 <- makeK2SM(grid, years = 2030, asPercent = TRUE)
  expect_equal(as.numeric(k2$pNoOverfishing[["2030"]]), 100)
  expect_equal(as.numeric(k2$green[["2030"]]), round(100 * pB))
})

test_that("kobeMomentGrid summarises missing moments without per-year spam", {

  dq_ok <- data.frame(
    Label  = c("Bratio_2028", "F_2028"),
    Value  = c(0.9, 0.5),
    StdDev = c(0.1, 0.05),
    stringsAsFactors = FALSE)
  dq_bad <- data.frame(
    Label  = c("Bratio_2028", "F_2028"),
    Value  = c(0.9, 0.5),
    StdDev = c(0.1, 0),
    stringsAsFactors = FALSE)
  covar <- data.frame(
    Label.i = "Bratio_2028", Label.j = "F_2028", corr = -0.5,
    stringsAsFactors = FALSE)

  runs <- list(
    "250" = list(derived_quants = dq_ok, CoVar = covar),
    "500" = list(derived_quants = dq_bad, CoVar = covar))

  expect_warning(
    grid <- kobeMomentGrid(runs, years = 2028, transform = "exact"),
    "1 run x year")
  expect_true(is.finite(grid$muF[grid$catch == 250]))
  expect_true(is.na(grid$muF[grid$catch == 500]))
})

test_that("infeasible log-scale correlation is clamped and summarised once", {

  # CVs and natural-scale rho that push exact transform outside (-1, 1)
  dq <- data.frame(
    Label  = c("Bratio_2070", "F_2070"),
    Value  = c(0.5, 0.4),
    StdDev = c(0.4, 0.35),
    stringsAsFactors = FALSE)
  covar <- data.frame(
    Label.i = "Bratio_2070", Label.j = "F_2070", corr = -0.95,
    stringsAsFactors = FALSE)
  runs <- list("1000" = list(derived_quants = dq, CoVar = covar))

  expect_warning(
    grid <- kobeMomentGrid(runs, years = 2070, transform = "exact"),
    "1 run x year.*outside \\(-1, 1\\)")
  expect_false(isTRUE(grid$feasible[1]))
  expect_true(abs(grid$rhoLog[1]) < 1)
  expect_true(is.finite(kobeGreen(grid$muB[1], grid$muF[1], grid$sigmaB[1],
                                  grid$sigmaF[1], grid$rhoLog[1])))
})

test_that("kobeMomentGridJabba + makeK2SM work on synthetic kbtrj", {
  set.seed(42)
  n <- 500
  syn <- data.frame(
    run = rep(c("1-B", "1-S"), each = n),
    year = 2024,
    stock = c(exp(rnorm(n, log(0.7), 0.2)), exp(rnorm(n, log(0.9), 0.15))),
    harvest = c(exp(rnorm(n, log(1.3), 0.25)), exp(rnorm(n, log(0.8), 0.2))))
  grid <- kobeMomentGridJabba(syn, years = 2024)
  expect_equal(nrow(grid), 2L)
  expect_true(all(c("run", "muB", "muF", "sigmaB", "sigmaF", "rhoLog") %in% names(grid)))
  expect_true(all(is.finite(grid$muB)))
  k2 <- makeK2SM(grid, years = 2024)
  expect_true(all(c("run", "green", "pNoOverfishing", "pNotOverfished") %in% names(k2$long)))
  expect_true(all(k2$long$green >= 0 & k2$long$green <= 1))
})

test_that("kobeMomentGridJabba accepts catch scenario column (projections)", {
  set.seed(7)
  n <- 300
  syn <- data.frame(
    catch = rep(c(0, 250, 1200), each = n),
    year = 2070,
    stock = exp(rnorm(3 * n, log(0.85), 0.2)),
    harvest = c(rep(1e-12, n), exp(rnorm(2 * n, log(0.9), 0.25))))
  grid <- kobeMomentGridJabba(syn, years = 2070)
  expect_equal(nrow(grid), 3L)
  expect_true("catch" %in% names(grid))
  expect_true(is.finite(grid$muB[grid$catch == 0]))
  expect_true(is.infinite(grid$muF[grid$catch == 0]) && grid$muF[grid$catch == 0] < 0)
  k2 <- makeK2SM(grid, years = 2070)
  expect_true(all(c("catch", "green") %in% names(k2$long)))
  expect_true(all(is.finite(k2$long$green)))
})

test_that("kobeInterp keeps deterministic F = 0 at a fitted catch", {

  grid <- data.frame(
    catch = c(0, 250), year = 2030,
    muB = c(log(0.9), log(0.8)), muF = c(-Inf, log(0.5)),
    sigmaB = c(0.2, 0.2), sigmaF = c(0, 0.3), rhoLog = c(0, -0.5),
    feasible = TRUE)
  g0 <- kobeInterp(grid, catch = 0, years = 2030)
  expect_true(is.infinite(g0$muF) && g0$muF < 0)
  expect_equal(g0$sigmaF, 0)
  expect_equal(kobeGreen(g0$muB, g0$muF, g0$sigmaB, g0$sigmaF, g0$rhoLog),
               1 - stats::pnorm(0, g0$muB, g0$sigmaB), tolerance = 1e-10)
})

test_that("kobeRebuildCurve accepts several targets", {

  years <- 2027:2035
  grid <- do.call(rbind, lapply(c(0, 500), function(c0)
    data.frame(catch = c0, year = years,
               muB = log(0.8) + 0.04 * (years - 2027) - c0 / 5000,
               muF = if (c0 == 0) -Inf else log(0.6),
               sigmaB = 0.15, sigmaF = if (c0 == 0) 0 else 0.2,
               rhoLog = 0, feasible = TRUE)))
  crv <- suppressWarnings(
    kobeRebuildCurve(grid, catch = c(0, 500), years = years,
                     target = c(0.50, 0.60)))
  expect_equal(nrow(crv), 4L)
  expect_true(all(c(0.50, 0.60) %in% crv$target))
  expect_true(all(crv$attained[crv$catch == 0]))
})

test_that("kobeMomentCI matches lognormalCI on converted moments", {

  grid <- data.frame(catch = 250, year = 2030,
                     muB = log(0.95), muF = log(0.72),
                     sigmaB = 0.174, sigmaF = 0.384, rhoLog = -0.5)
  ci <- kobeMomentCI(grid, quant = "B", levels = 0.95)
  expect_equal(ci$ratio, exp(grid$muB), tolerance = 1e-12)
  z <- stats::qnorm(0.975)
  expect_equal(ci$ratioLo, exp(grid$muB - z * grid$sigmaB), tolerance = 1e-10)
  expect_equal(ci$ratioHi, exp(grid$muB + z * grid$sigmaB), tolerance = 1e-10)
  z0 <- kobeMomentCI(transform(grid, muF = -Inf, sigmaF = 0),
                     quant = "F", levels = 0.95)
  expect_equal(z0$ratio, 0)
  expect_equal(z0$ratioLo, 0)
})

test_that("kobeDraws samples F = 0 as harvest pinned at 0", {

  grid <- data.frame(catch = 0, year = 2070,
                     muB = log(1.1), muF = -Inf,
                     sigmaB = 0.2, sigmaF = 0, rhoLog = 0)
  set.seed(1)
  d <- kobeDraws(grid, n = 50)
  expect_equal(nrow(d), 50L)
  expect_true(all(d$harvest == 0))
  expect_true(all(d$stock > 0))
})
