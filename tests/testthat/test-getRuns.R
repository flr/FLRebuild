ssStub <- function()
  list(timeseries = data.frame(Yr = 2024),
       derived_quants = data.frame(Label = "SSB_2024", Value = 1, StdDev = 0.1))

jabbaStub <- function(run = "1-S")
  list(scenario = run,
       kbtrj = data.frame(year = 2024, iter = 1L, stock = 0.9, harvest = 0.8,
                          stringsAsFactors = FALSE),
       kobe = data.frame(stock = 0.9, harvest = 0.8, level = run,
                         stringsAsFactors = FALSE))

test_that("getRuns is an S4 generic with character and list methods", {
  expect_true(isGeneric("getRuns"))
  expect_true(hasMethod("getRuns", "character"))
  expect_true(hasMethod("getRuns", "list"))
})

test_that("getRuns list method sorts SS runs by catch", {
  runs <- list("250" = ssStub(), "0" = ssStub())
  out <- getRuns(runs)
  expect_equal(names(out), c("0", "250"))
})

test_that("getRuns list method keeps JABBA run labels", {
  jb <- list("2-S" = jabbaStub("2-S"), "1-S" = jabbaStub("1-S"))
  out <- getRuns(jb)
  expect_equal(names(out), c("2-S", "1-S"))
  out1 <- getRuns(jb, runs = "1-S")
  expect_equal(names(out1), "1-S")
})

test_that("getRuns list method rejects mixed or empty input", {
  expect_error(getRuns(list()), "empty")
  mixed <- list("250" = ssStub(), "1-S" = jabbaStub())
  expect_error(getRuns(mixed), "SS_output-shaped or all JABBA")
})

test_that("getRuns character method loads SS .Rdata by catch name", {
  td <- tempfile("ssruns")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  ss <- ssStub()
  save(ss, file = file.path(td, "HW2e_0t.Rdata"))
  save(ss, file = file.path(td, "HW2e_250t.Rdata"))
  out <- getRuns(td)
  expect_equal(names(out), c("0", "250"))
  expect_true(all(c("timeseries", "derived_quants") %in% names(out[["0"]])))
})

test_that("getRuns character method auto-detects JABBA files", {
  td <- tempfile("jbrun")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  jabba <- jabbaStub("1-S")
  save(jabba, file = file.path(td, "SMA2026_1-S_jabba.rdata"))
  jabba <- jabbaStub("2-S")
  save(jabba, file = file.path(td, "SMA2026_2-S_jabba.rdata"))
  out <- getRuns(td)
  expect_equal(sort(names(out)), c("1-S", "2-S"))
  out1 <- getRuns(td, runs = "1-S")
  expect_equal(names(out1), "1-S")
})

test_that("getJabbaRuns is a source='jabba' wrapper with ensemble default", {
  td <- tempfile("jbrun")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  for (lab in c("1-B", "1-S", "2-S", "3-X")) {
    jabba <- jabbaStub(lab)
    save(jabba, file = file.path(td, paste0("SMA2026_", lab, "_jabba.rdata")))
  }
  jb <- getJabbaRuns(td)
  expect_equal(sort(names(jb)), c("1-B", "1-S", "2-S"))
  jb2 <- getRuns(td, source = "jabba", runs = c("1-B", "1-S", "2-S"))
  expect_equal(names(jb2), names(jb))
})

test_that("getRuns named character paths dispatch on name type", {
  td <- tempfile("paths")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  ss <- ssStub()
  ssFile <- file.path(td, "run.Rdata")
  save(ss, file = ssFile)
  jabba <- jabbaStub("1-S")
  jbFile <- file.path(td, "fit.rdata")
  save(jabba, file = jbFile)

  ssOut <- getRuns(c("250" = ssFile, "0" = ssFile))
  expect_equal(names(ssOut), c("0", "250"))

  jbOut <- getRuns(c("1-S" = jbFile))
  expect_equal(names(jbOut), "1-S")
  expect_true(is.data.frame(jbOut[["1-S"]]$kbtrj))
})
