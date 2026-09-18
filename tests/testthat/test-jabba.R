test_that("jabbaPriors from FLStock and FLBRP matches jabba-f inputs", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  data(ple4brp)

  pri <- jabbaPriors(ple4, eq = ple4brp)
  expect_type(pri, "list")
  expect_true(all(c("pr", "pr.sd") %in% names(pri)))
  expect_true(all(c("r", "psi", "shape") %in% dimnames(pri$pr)$params))
  expect_true(all(is.finite(c(pri$pr[c("r", "psi", "shape")]))))
  expect_gt(c(pri$pr["r"]), 0)
  expect_gt(c(pri$pr["shape"]), 0)
  expect_lt(c(pri$pr["shape"]), 1)

  pri2 <- jabbaPriors(ple4brp, stock = ple4)
  expect_equal(c(pri$pr["r"]), c(pri2$pr["r"]), tolerance = 1e-8)
})

test_that("jabbaInput from FLStock has catch and ebiomass index", {
  library(FLCore)
  data(ple4)
  inp <- jabbaInput(ple4, index = ebiomass)
  expect_true(all(c("year", "catch") %in% names(inp$catch)))
  expect_true(all(c("year", "index") %in% names(inp$index)))
  expect_equal(nrow(inp$catch), dim(ple4)[2])
  expect_equal(nrow(inp$index), dim(ple4)[2])

  inp2 <- jabbaInput(ple4, index = "ebiomass")
  expect_equal(inp$index$index, inp2$index$index, tolerance = 1e-10)
})

test_that("jabbaInput can add F or F/Fmsy (jabba-f auxiliary)", {
  library(FLCore)
  library(FLBRP)
  data(ple4)
  data(ple4brp)

  inp.f <- jabbaInput(ple4, index = ebiomass, f = fbar)
  expect_equal(inp.f$auxiliary.type, "f")
  expect_true(all(c("year", "f") %in% names(inp.f$auxiliary)))
  expect_equal(nrow(inp.f$auxiliary), dim(ple4)[2])

  inp.ff <- jabbaInput(ple4, eq = ple4brp, index = ebiomass, f = ffmsy)
  expect_equal(inp.ff$auxiliary.type, "ffmsy")
  expect_true(all(c("year", "ffmsy") %in% names(inp.ff$auxiliary)))
  expect_true(all(is.finite(inp.ff$auxiliary$ffmsy)))

  q <- ffmsy(ple4, eq = ple4brp)
  expect_s4_class(q, "FLQuant")
  expect_equal(inp.ff$auxiliary$ffmsy, as.numeric(c(q)), tolerance = 1e-10)
})

test_that("jabba2biodyn coerces a JABBA-like fit", {
  skip_if_not_installed("mpb")
  yrs <- 2000:2004
  ts <- array(1, dim = c(5, 3, 1),
              dimnames = list(year = yrs, quant = c("mu", "lci", "uci"),
                               type = "B"))
  ts[, "mu", "B"] <- 100 + seq_along(yrs)
  fit <- list(
    assessment = "ple4",
    scenario = "test",
    timeseries = ts,
    inputseries = list(
      catch = data.frame(year = yrs, catch = 10),
      cpue = data.frame(year = yrs, index = 1)
    ),
    pars = data.frame(
      Median = c(0.4, 2000, 1.2, 0.9),
      row.names = c("r", "K", "m", "psi")
    )
  )
  bd <- jabba2biodyn(fit)
  expect_s4_class(bd, "biodyn")
  expect_equal(c(bd@params["r"]), 0.4, tolerance = 1e-10)
  expect_equal(c(bd@params["k"]), 2000, tolerance = 1e-10)
  expect_equal(as.numeric(c(bd@stock)), as.numeric(ts[, "mu", "B"]))

  bd2 <- jabba2biodyn(list(fit = fit))
  expect_equal(c(bd@params["r"]), c(bd2@params["r"]))
})

test_that("jabbaTs keeps sequential years starting at 1", {
  yrs <- 1:20
  ts <- array(NA_real_, dim = c(20, 3, 2),
              dimnames = list(year = yrs, quant = c("mu", "lci", "uci"),
                               type = c("B", "procB")))
  ts[, "mu", "B"] <- seq_along(yrs)
  ts[, "lci", "B"] <- seq_along(yrs) - 0.5
  ts[, "uci", "B"] <- seq_along(yrs) + 0.5
  ts[, "mu", "procB"] <- 0.1
  ts[, "lci", "procB"] <- 0.05
  ts[, "uci", "procB"] <- 0.15
  fit <- list(timeseries = ts)
  jb <- jabbaTs(fit, vars = c(stock = "B", pe = "procB"))
  expect_equal(as.numeric(dimnames(jb$stock)$year), yrs)
  expect_equal(as.numeric(c(jb$stock)), as.numeric(yrs))
  expect_equal(as.numeric(c(jb$pe)), rep(0.1, 20))

  ci <- jabbaTsCI(fit, vars = c(stock = "B", pe = "procB"))
  expect_equal(unique(ci$year), yrs)
  expect_equal(ci$mu[ci$qname == "stock"], as.numeric(yrs))
  expect_equal(ci$lci[ci$qname == "stock"], as.numeric(yrs) - 0.5)
  expect_equal(ci$uci[ci$qname == "stock"], as.numeric(yrs) + 0.5)

  p <- plotJabbaTs(fit, qname = "stock")
  expect_s3_class(p, "ggplot")

  om1 <- FLCore::FLQuant(seq_along(yrs), dimnames = list(year = yrs))
  om2 <- FLCore::FLQuant(seq_along(yrs) / 2, dimnames = list(year = yrs))
  p2 <- plotJabbaTs(fit, qname = "stock",
                    om = FLCore::FLQuants(`OM SSB` = om1, `OM EB` = om2),
                    om.qname = "stock")
  expect_s3_class(p2, "ggplot")
})
