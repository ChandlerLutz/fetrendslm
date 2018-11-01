## c:/Dropbox/Rpackages/fetrendslm/tests/testthat/test_deaths.R

##    Chandler Lutz
##    Questions/comments: cl.eco@cbs.dk
##    $Revisions:      1.0.0     $Date:  2018-10-04

## -- To test the deaths data -- ##

library(lfe); library(broom)

rm(list = ls())

context("deaths data")

## From mastering metrics
data(deaths)

DT <- deaths %>%
  .[, beerpercap.quintile := as.character(dplyr::ntile(beerpercap, 5))]

##For felm
deaths2 <- copy(DT) %>%
  .[, year.char := as.factor(year.char)] %>%
  .[, state.char := as.factor(state.char)] %>%
  .[, beerpercap.quintile := as.factor(beerpercap.quintile)]

DT2 <- deaths2

DT <- CLmisc::dummy_cols(DT, select_columns = "beerpercap.quintile",
                         remove_first_dummy = TRUE,
                         by_reference = TRUE)
beercap.quantile.vars <- names(DT) %>% .[grepl("beerpercap.quintile_", x = .)]


test_that("Basic Setup works against felm", {
  temp <- FeTrendsLm$new(DT = copy(DT),
                         .f = ~ state.char + year.char,
                         cluster.vars = c("state.char"),
                         chunk.size = 5)
  res <- temp$fetrendslm(y.var = "mrate", x.vars = c("legal"))

  res.felm <- felm(mrate ~ legal | state.char + year.char| 0 | state.char,
                   DT2) %>%
      broom::tidy(.) %>% setDT

  ##The coefficients
  res.coef <- res$DT.tidy[, estimate] %>% round(2)
  felm.coef <- res.felm[, estimate] %>% round(2)
  expect_true(all.equal(res.coef, felm.coef))

  ##The standard errors
  felm.se <- res.felm[, std.error] %>% round(2)
  res.se <- res$DT.tidy[, std.error] %>% round(2)
  expect_true(all.equal(felm.se, res.se, tolerance = 0.2))

})

test_that("Basic Setup works against felm with beercap.quantile", {
  temp <- FeTrendsLm$new(DT = copy(DT),
                         .f = ~ state.char + year.char,
                         cluster.vars = c("state.char"),
                         chunk.size = 5)
  res <- temp$fetrendslm(y.var = "mrate", x.vars = c("legal", beercap.quantile.vars))

  res.felm <- felm(mrate ~ legal | state.char + year.char + beerpercap.quintile | 0 |
                     state.char,
                   DT2, exactDOF = TRUE) %>%
      broom::tidy(.) %>% setDT

  ##The coefficients
  res.coef <- res$DT.tidy[, estimate] %>% round(2)
  felm.coef <- res.felm[, estimate] %>% round(2)
  expect_true(all.equal(res.coef[1], felm.coef[1]))

  ##The standard errors
  felm.se <- res.felm[, std.error] %>% round(2)
  res.se <- res$DT.tidy[, std.error] %>% round(2)
  expect_true(all.equal(felm.se[1], res.se[1], tolerance = 0.2))

})


formulas <- list(
  ~ state.char,
  ~ state.char + year.char,
  ~ year.char + state.char,
  ~ state.char + year.char + state.char:year,
  ~ year.char + state.char + state.char:year,
  ~ state.char + year.char + state.char:year + state.char:year2,
  ~ year.char + state.char + state.char:year + state.char:year2,
  ~ state.char + state.char:year,
  ~ state.char + state.char:year + state.char:year2
)

y.vars <- c("mrate")
x.vars <- list(
  c("legal"),
  c("legal", "beerpercap")
)

cluster.vars <- list(
  NA_character_,
  c("state.char"),
  c("state.char", "year.char")
)


chunk.size = c(2, 5, 10, 50)


keycolvars <- list(
  NA_character_,
  c("state.char"),
  c("state.char", "year.char")
)

weights <- c(NA_character_, "pop")

DT.test <- expand.grid(.f = formulas, y.var = y.vars, x.vars = x.vars,
                       cluster.vars = cluster.vars,
                       chunk.size = chunk.size,
                       keycolvars = keycolvars,
                       ##weights = weights,
                       weights = NA,
                       stringsAsFactors = FALSE) %>%
  setDT

for (i in 1L:nrow(DT.test)) {

  DT.temp <- DT.test[i, ]

  ##print this interation
  ##print(DT.temp); cat("\n\n")

  weights <- !is.na(DT.temp[["weights"]])
  if (weights) {
    weight.col <- DT.temp[["weights"]]
    weight.vals <- DT[[weight.col]]
  } else {
    weight.col <- NULL
    weight.vals <- NULL
  }

  use.clus.se <- !anyNA(DT.temp[["cluster.vars"]][[1]])
  if (use.clus.se) {
    cluster.vars <- DT.temp[["cluster.vars"]][[1]]
    cluster.vars.felm <- paste0(cluster.vars, collapse = " + ")
  } else {
    cluster.vars <- NULL
    cluster.vars.felm <- " 0 "
  }

  if (anyNA(DT.temp[["keycolvars"]][[1]])) {
    keycolvars <- NULL
  } else {
    keycolvars <- DT.temp[["keycolvars"]][[1]]
  }

  .f <- DT.temp[[".f"]][[1]]
  x.vars <- DT.temp[["x.vars"]][[1]]
  y.var <- DT.temp[["y.var"]][[1]]

  temp <- FeTrendsLm$new(
    DT = copy(DT),
    .f = .f,
    main.reg.vars = c(y.var, x.vars),
    chunk.size = DT.temp[["chunk.size"]][[1]],
    keycolsv = keycolvars,
    weight.var = weight.col
  )

  res <- temp$fetrendslm(y.var = y.var, x.vars = x.vars)

  f.felm <- sprintf("%s ~ %s | %s | 0 | %s",
                    y.var,
                    paste0(x.vars, collapse = "+"),
                    as.character(.f)[2],
                    cluster.vars.felm
                    )

  if (weights) {
    res.felm <- felm(as.formula(f.felm), data = DT2, weights = weight.vals)

  } else {
    res.felm <- felm(as.formula(f.felm), data = DT2)
  }
  res.felm.rse <- res.felm$rse %>% as.numeric

  res.felm <- res.felm %>% broom::tidy(.) %>% setDT %>%
    .[, std.error := res.felm.rse]

  ##The coefficients
  res.coef <- res$DT.tidy[, estimate] %>% round(2)
  felm.coef <- res.felm[, estimate] %>% round(2)
  expect_true(all.equal(res.coef, felm.coef), info = sprintf("i = %s", i))

  ##The standard errors
  felm.se <- res.felm[, std.error] %>% round(2)
  res.se <- res$DT.tidy[, std.error] %>% round(2)
  expect_true(all.equal(felm.se, res.se, tolerance = 0.2), info = sprintf("i = %s", i))

}
