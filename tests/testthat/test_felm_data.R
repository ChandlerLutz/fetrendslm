## c:/Dropbox/Rpackages/fetrendslm/tests/testthat/test_felm_data.R

##    Chandler Lutz
##    Questions/comments: cl.eco@cbs.dk
##    $Revisions:      1.0.0     $Date:  2018-10-04

library(lfe);

rm(list = ls())

context("felm generated data data")

set.seed(12345)
## create covariates
x <- rnorm(1000)
x2 <- rnorm(length(x))


## individual and firm
id <- factor(sample(20,length(x),replace=TRUE))
firm <- factor(sample(13,length(x),replace=TRUE))

## effects for them
id.eff <- rnorm(nlevels(id))
firm.eff <- rnorm(nlevels(firm))

## left hand side
u <- rnorm(length(x))
y <- x + 0.5*x2 + id.eff[id] + firm.eff[firm] + u

DT <- data.table(y=y,x=x,x2=x2,id=as.character(id),firm=as.character(firm)) %>%
  .[, t := 1:.N, by = firm] %>%
  .[, t.char := as.character(t)] %>%
  .[, t2 := t ^ 2] %>%
  .[, t2.char := as.character(t)]

DT2 <- copy(DT) %>%
  .[, firm := as.factor(firm)] %>%
  .[, id := as.factor(id)] %>%
  .[, t.char := as.factor(t.char)] %>%
  .[, t := as.numeric(t)] %>%
  .[, t2.char := as.factor(t.char)] %>%
  .[, t2 := as.numeric(t2)]

formulas <- list(
  ~ id,
  ~ firm,
  ~ firm + id,
  ~ t.char,
  ~ id + t.char,
  ~ firm + t.char,
  ~ id + firm + t.char,
  ~ id + firm + t.char + id:t,
  ~ id + firm + t.char + firm:t,
  ~ id + firm + t.char + id:t + id:t2,
  ~ id + firm + t.char + firm:t + firm:t2,
  ~ id + firm + t.char + id:t + id:t2 + firm:t + firm:t2
)

y.vars <- "y"
x.vars <- list(c("x"), c("x", "x2"))

chunk.size = c(2, 5, 10, 50)

keycolvars <- list(
  NA_character_,
  c("id"),
  c("id", "firm"),
  c("id", "firm", "t.char")
)

cluster.vars <- list(
  NA_character_
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
