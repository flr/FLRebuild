test_that("relsig is 1 at BMSY and larger at low depletion", {
  expect_equal(relsig(1), 1)
  expect_gt(relsig(0.2), relsig(1))
  expect_lt(relsig(2), relsig(1))
  expect_length(relsig(c(0.2, 1, 2)), 3)
})

test_that("pe vector sets a, c, d", {
  rs <- FLRebuild:::.jabbaRelsigPars(pe = c(a = 0.2, c = 0.35, d = 2))
  expect_equal(rs, list(a = 0.2, c = 0.35, d = 2))
  rs.pos <- FLRebuild:::.jabbaRelsigPars(pe = c(0.2, 0.35, 2))
  expect_equal(rs.pos, rs)
  rs.part <- FLRebuild:::.jabbaRelsigPars(pe = c(a = 0.2))
  expect_equal(rs.part$a, 0.2)
  expect_equal(rs.part$c, 0.223)
  expect_equal(FLRebuild:::.jabbaRelsigPars()$a, 0.128)
  expect_error(FLRebuild:::.jabbaRelsigPars(pe = c(a = 0)), "pe")
  expect_error(FLRebuild:::.jabbaRelsigPars(pe = c(b = 0.2)), "names")
})

test_that("JAGS relsig constants can be patched", {
  f <- tempfile(fileext = ".jags")
  writeLines(c("a_pe <- 0.128", "c_pe <- 0.223", "d_pe <- 1.326"), f)
  expect_true(FLRebuild:::.patchJabbaJagsRelsig(f, 0.2, 0.4, 2))
  txt <- paste(readLines(f), collapse = "\n")
  expect_match(txt, "a_pe <- 0.2", fixed = TRUE)
  expect_match(txt, "c_pe <- 0.4", fixed = TRUE)
  expect_match(txt, "d_pe <- 2", fixed = TRUE)
})

test_that("compareJabbaPE is exported", {
  expect_true(exists("compareJabbaPE", mode = "function", envir = asNamespace("FLRebuild")))
})

test_that("compareJabbaPE requires JABBA proc.dyn", {
  skip_if_not_installed("JABBA")
  has <- "proc.dyn" %in% names(formals(JABBA::build_jabba))
  if (!has) {
    expect_error(compareJabbaPE(data.frame(year = 1, catch = 1)),
                 "proc.dyn")
  } else {
    expect_true(exists("compareJabbaPE", mode = "function"))
  }
})
