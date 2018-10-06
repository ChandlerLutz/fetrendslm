## c:/Dropbox/Rpackages/fetrendslm/tests/testthat/test_petersen.R

##    Chandler Lutz
##    Questions/comments: cl.eco@cbs.dk
##    $Revisions:      1.0.0     $Date:  2018-10-05

##Clear the workspace
rm(list = ls())

context("petersen data")

data(petersen)

DT <- petersen %>%
  .[, firmid.char := as.character(firmid)] %>%
  .[, year.char := as.character(year)] %>%
  .[, year2 := year ^ 2]

DT2 <- copy(DT) %>%
  .[, firmid.char := as.factor(firmid.char)] %>%
  .[, year.char := as.factor(year.char)]

formulas <- list(
  ~ firmid.char,
  ~ year.char,
  ~ firmid.char + year.char,
  ~ firmid.char + year.char + firmid.char:year,
  ~ firmid.char + year.char + firmid.char:year + firmid.char:year2
)


y.vars <- "y"
x.vars <- "x"

chunk.size = c(2, 5, 10, 50)

keycolvars <- list(
  NA_character_,
  c("firmid"),
  c("firmid", "year.char")
)

cluster.vars <- list(
  NA_character_,
  c("firmid"),
  c("firmid", "year.char"),
  c("year.char")
)


weights <- NA_character_


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
  expect_true(all.equal(res.coef, felm.coef, tolerance = 0.1),
              info = sprintf("i = %s", i))

  ##The standard errors
  felm.se <- res.felm[, std.error] %>% round(2)
  res.se <- res$DT.tidy[, std.error] %>% round(2)
  expect_true(all.equal(felm.se, res.se, tolerance = 0.2), info = sprintf("i = %s", i))

}
