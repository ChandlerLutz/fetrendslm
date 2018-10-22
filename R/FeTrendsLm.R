# c:/Dropbox/Rpackages/fetrendslm/R/FeTrensLm.R

##    Chandler Lutz
##    Questions/comments: cl.eco@cbs.dk
##    $Revisions:      1.0.0     $Date:  2018-10-04


## ----------- ##
## -- UTILS -- ##
## ----------- ##


## The initialization function.
initialize <- function(DT, .f, main.reg.vars = NULL, cluster.vars = NULL,
                       chunk.size = 100, keycolsv = NULL,
                       weight.var = NULL) {

  ##Just for you, R CMD check
  cluster <- char.var <- temp.index <- temp.char. <- NULL;


  ## -- Assign inputs -- ##
  private$..N <- nrow(DT)

  private$..f <- .f; rm(.f)
  if (length(private$..f) != 2)
    stop("Error: .f must be a RHS formula only")

  if (!is.null(main.reg.vars))
    private$..main.reg.vars <- main.reg.vars; rm(main.reg.vars)




  ##Chunk size -- save for later
  private$..chunk.size <- chunk.size; rm(chunk.size)

  ##key the data
  if (!is.null(keycolsv)) {
    setkeyv(DT, keycolsv)
    private$..keycolsv <- keycolsv; rm(keycolsv)
  }

  ##The clustered standard error information
  if (!is.null(cluster.vars)) {
    private$..cluster.vars <- cluster.vars
    p1 <- function(...) paste(..., sep = "_")
    private$..DT.cluster <- DT[, .SD, .SDcols = cluster.vars] %>%
      .[, cluster := Reduce(p1, .SD)] %>%
      ##delete the cluster var columns, so we only have one
      ##cluster column
      .[, (cluster.vars) := NULL]

    ##Also get the number of clusters
    private$..clusters.num <- uniqueN(private$..DT.cluster[["cluster"]])
    rm(p1, cluster.vars)
  }


  if (!is.null(weight.var)) {
    stop("Error: Weights not currently implemented")
    if (!is.character(weight.var) | length(weight.var) > 1)
      stop("Error: weight.var must be a character string with column name of weights")
    private$..weight.var <- weight.var; rm(weight.var)
    private$..weight.vals <- DT[[private$..weight.var]]
    private$..weight.vals <- private$..weight.vals / sum(private$..weight.vals)
  }

  ##Parse the formula
  rhs <- as.character(private$..f)[2]

  ##the variables to partial out
  po.vars <- strsplit(rhs,"\\+") %>% unlist %>% stringr::str_trim(.) %>%
    ##split for multiplication
    strsplit(., split = ":") %>%
    ##Only unique vars
    unique(.)

  po.vars.unique <- unlist(po.vars) %>% unique

  ##The variable classes
  DT.var.classes <- sapply(DT, class)


  ##Get a data.table with the classes
  DT.po.vars <- data.table(po.vars = po.vars) %>%
    .[, class := lapply(po.vars, function(x) DT.var.classes[x])]

  ##Helper function to get the the character vectors
  ##given a vector of classes
  f_get_char_var <- function(vars, classes) {
    vars[classes == "character"]
  }

  ##Make sure that there are no character missing classes
  lapply(DT.po.vars[["class"]], function(x) {
    if (anyNA(x)) {
      print(DT.po.vars)
      stop("Error: Some variables in .f have missing classes. Are you sure they are in the dataset?")
      }
    return(invisible(NULL))
  })

  ##Get the character vector for the
  DT.po.vars <- DT.po.vars[, char.var := Map(f_get_char_var, po.vars, class)]

  ##Check to make sure that each entry only has one character vector
  if (any(DT.po.vars[, sapply(char.var, length)] > 1))
    stop("Error: Cannot interact two character columns in .f")

  DT.po.vars <- DT.po.vars %>%
    ##Convert char.var in DT.po.vars a character column
    .[, char.var := sapply(char.var, function(x) x[1])] %>%
    ##A dummy if the po.vars are an interaction
    ##(there is more than variable in po.vars8)
    .[, interaction := sapply(po.vars, function(x) as.integer(length(x) > 1))]

  ##Order by the character vector
  DT.po.vars <- DT.po.vars[order(char.var)]
  private$..DT.po.vars <- DT.po.vars

  ##The variables to keep
  if (is.null(private$..main.reg.vars)) {
    private$..main.reg.vars <- names(DT) %>% .[!(. %in% po.vars.unique)]
  }

  ##Only include numeric vars in main.reg.vars
  numeric.vars <- DT.var.classes %>% .[. %in% c("numeric", "integer")] %>%
    names(.)
  private$..main.reg.vars <- private$..main.reg.vars %>%
    .[. %in% c(numeric.vars)]

  ##Leave Y as a sparse (dense) matrix
  Y <- as(as.matrix(DT[, .SD, .SDcols = private$..main.reg.vars]), "sparseMatrix")

  ##Remove main.reg.vars from the data.table
  delete.vars <- private$..main.reg.vars %>% .[!(. %in% po.vars.unique)] %>%
    unique(.)
  if (length(delete.vars)) {
    DT[, c(delete.vars) := NULL]; rm(delete.vars)
  }

  ##Only keep columns that are going to be used in DT
  delete.vars <- names(DT) %>% .[!(. %in% po.vars.unique)] %>% unique(.)
  if (length(delete.vars) > 0) {
    DT <- DT[, c(delete.vars) := NULL]; rm(delete.vars)
  }

  ## --- Partial out all char vars and their interactions --- ##

  ##Add the index to the data to help with creating sparse matrices
  DT <- DT[, temp.index := 1L:.N]

  char.vars.unique <- unique(DT.po.vars[["char.var"]])

  ##Loop through the character vectors and partial out all relevent columns
  for (i in 1:length(char.vars.unique)) {

    ##the character vector
    char.temp <- char.vars.unique[i]

    ##change its name by reference for easy access later
    setnames(DT, char.temp, "temp.char.")

    char.temp.unique <- unique(DT[["temp.char."]])

    counter <- 0
    while(length(char.temp.unique) > 0) {
      counter <- counter + 1

      ##Get the unique values of the character vector for this chunk
      char.temp.vals <- char.temp.unique[1:private$..chunk.size] %>% .[!is.na(.)]

      DT.i.j.vals <- DT[temp.char. %chin% c(char.temp.vals),
                        .(temp.index, temp.vals = 1L, temp.grp = .GRP),
                        by = temp.char.]

      ## The Sparse X matrix
      X <- sparseMatrix(
        i = DT.i.j.vals[["temp.index"]], j = DT.i.j.vals[["temp.grp"]],
        x = DT.i.j.vals[["temp.vals"]],
        dims = c(nrow(DT), length(char.temp.vals))
      )

      ## Loop through the different potential interactions for
      ## char.temp and add any interactions
      interact.vars <- DT.po.vars[char.var == c(char.temp) & interaction == 1L,
                                  po.vars]

      ##If there are interaction terms
      if (length(interact.vars) > 0) {

        ##Loop through the interaction terms
        for (j in 1:length(interact.vars)) {
          temp.vars <- interact.vars[[j]]

          ##There are interaction terms, partial out those as well
          other.vars <- temp.vars %>% .[. != char.temp]

          if (length(other.vars) == 1) {
            ##only one other var -- a vector
            X.numeric.vec <- DT[[other.vars]]

          } else {
            ##Multiple other variables
            X.numeric.vec <- DT[, .(temp.col = Reduce("*", .SD)), .SDcols = other.vars]
          }

          ##Multiply X times X.1 for the interactions
          if (j == 1) {
            X.1 <- cbind(X * X.numeric.vec)
          } else {
            X.1 <- cbind(X.1, X * X.numeric.vec)
          }
        }

        if (DT.po.vars[char.var == (char.temp), sum(interaction) / .N == 1]) {
          ##If only interaction terms, set X.1 to X and ignore the non-interaction
          ##terms
          X <- X.1
        } else {
          ##Use both the raw dummies and the interaction terms
          X <- cbind(X, X.1)
        }

        ##clean up
        rm(X.1)

      } ## end of length(interact.vars) > 0 IF

      ##Get the number of columns and add to private$..po.ncol
      private$..po.ncol = private$..po.ncol + ncol(X)

      ##Partial out X times the interaction terms.
      if (!is.null(private$..weight.vals)) {
        ##weighthed least squares

        ##See https://stackoverflow.com/a/24196694
        ## See the wikipedia article on how you can
        ## just multiply X and Y by the square root of the diagonal
        ##https://en.wikipedia.org/wiki/Weighted_least_squares
        ## W <- Diagonal(x = private$..weight.vals)
        ## Y <- Y - X %*% solve(t(X) %*% W %*% X, t(X) %*% W %*% Y)
        Y <- Y - X %*% (solve(crossprod(X * sqrt(private$..weight.vals)),
                              crossprod(X * sqrt(private$..weight.vals),
                                        Y * sqrt(private$..weight.vals))))
      } else {
        ##non-weighted least squares
        beta <- solve(crossprod(X), crossprod(X, Y))
        Yhat <- X %*% beta; rm(beta)
        Y <- Y - Yhat; rm(Yhat)
      }

      ##Remove the values from char.temp.unique for this
      char.temp.unique <- char.temp.unique[-c(1:private$..chunk.size)]
      ##clean up
      rm(X, DT.i.j.vals, interact.vars)

    } ##End of char.var.unique private$..chunk.size while loop

    ##Change the name back
    setnames(DT, "temp.char.", char.temp)

  }  ## end of character vars for


  rm(DT)

  self$Y = Y

}


## The regression function
fetrendslm <- function(y.var, x.vars) {

  ##Just for you, R CMD check
  X.mat <- X.mat <- U.mat <- statistic <- std.error <- p.value <- cluster <- NULL

  Y <- self$Y[, y.var, drop = FALSE]
  X <- self$Y[, x.vars, drop = FALSE]

  if (!is.null(private$..weight.vals)) {
    ##weighted least squares
    ## See the wikipedia article on how you can
    ## just multiply X and Y by the square root of the diagonal
    ##https://en.wikipedia.org/wiki/Weighted_least_squares
    ## W <- Diagonal(x = private$..weight.vals)
    ## beta <- solve(t(X) %*% W %*% X, t(X) %*% W %*% Y)
    beta <- solve(crossprod(X * sqrt(private$..weight.vals)),
                  crossprod(X * sqrt(private$..weight.vals),
                            Y * sqrt(private$..weight.vals)))

  } else {
    beta <- solve(crossprod(X), crossprod(X, Y))
  }

  ##The residuals
  resid <- Y - X %*% beta
  colnames(resid) <- "resid"



  if (is.null(private$..cluster.vars)) {
    ##No clustering, regular robust standard errors
    ## (X'X)^{-1} X' Omega X (X'X)^{-1}

    ##The degrees of freedom correction -- HC1 from sandwich
    dfc <- private$..N / (private$..N - ncol(X) - private$..po.ncol)

    omega <- tcrossprod(resid) %>% diag %>% Diagonal(x = .)
    ##the meat --  X' Omega X
    meat <- crossprod(X, omega) %*% X
    rm(omega)

    ##The bread
    XXinv <- solve(crossprod(X))
    se <- (XXinv %*% meat %*% XXinv * dfc) %>% diag %>% sqrt
    rm(XXinv)

  } else {
    ##clustering
    m <- private$..clusters.num
    N <- private$..N

    ##Degrees of Freedom correction from Cameron and Miller (2015) eq. 12
    dfc <- (m / (m - 1)) * ((N / (N - ncol(X) - private$..po.ncol)))
    X.DT <- X %>% as.matrix %>% as.data.table
    resid.DT <- resid %>% as.matrix %>% as.data.table

    ##helper function to get the meat of the var-cov matrix
    f_meat <- function(X, U) {
      XU <- crossprod(X, U)
      UX <- crossprod(U, X)
      out <- XU %*% UX
    }

    ## Formula for meet from Cameron and Miller (2015) eq. 11
    meat <- cbind(X.DT, resid.DT, private$..DT.cluster) %>%
      ##Collapse by the cluster variables
      .[, .(X.mat = list(as(as.matrix(.SD), "sparseMatrix")),
            U.mat = list(as(as.matrix(resid), "dgeMatrix"))),
        .SDcols = x.vars,
        by = cluster] %>%
      ##get the meet
      .[, .(meat = Map(f_meat, X.mat, U.mat))] %>%
      ##sum over all clusters
      .[, Reduce("+", meat)]

    ##The bread
    XXinv <- solve(crossprod(X))
    se <- (XXinv %*% meat %*% XXinv * dfc) %>% diag %>% sqrt

  }

  ## -- get ready for the return -- ##
  beta <- beta %>% as.numeric %>% setNames(x.vars)
  se <- se %>% as.numeric %>% setNames(x.vars)

  ##A tidy data.table matching broom::tidy(.)
  DT.tidy <- data.table(term = x.vars, estimate = beta, std.error = se) %>%
    .[, statistic := beta / std.error]
  ##the degrees of freedom for the t-test
  t.df <- private$..N - length(x.vars) - private$..po.ncol
  ##Add in the p-value from a two-sided t-test
  DT.tidy <- DT.tidy[, p.value := (1 - pt(abs(statistic), df = c(t.df))) * 2]




  ##get a small lm model to trick stargazer

  DT.temp <- cbind(Y[1:50, ], X[1:50, ]) %>% as.matrix %>% as.data.table %>%
    setnames(names(.), c(y.var, x.vars))

  f <- sprintf("%s ~ 0 + .", y.var)
  mod.lm <- lm(as.formula(f), DT.temp)
  mod.lm$coefficients <- beta

  return(list(coef = beta, se = se, N = private$..N, DT.tidy = DT.tidy,
              lm.mod = mod.lm))

}


## ------------------------- ##
## -- FeTrendsLm R6 Class -- ##
## ------------------------- ##

##R6 documenation from https://stackoverflow.com/a/45603005

#' FeTrendsLm
#'
#'
#' Large Data Estimation linear fixed effect models with trends
#'
#' The FeTrendsLm class provides an R6 class to estimate linear models
#' fixed effects and/or trends with large data. Consider the following
#' model from Angrist and Pischke "Mastering Metrics" (2015), equation
#' (5.6) that is going to estimate a difference-in-differences model
#' that measures the impact of legalized drinking on deaths (see
#' examples below for implentation):
#'
#' \eqn{Y_{st} = \alpha + \gamma_{rDD} \textit{LEGAL}_{st} + \gamma_i + \gamma_t + \gamma_{st} + \varepsilon_{it}}
#'
#' \eqn{\gamma_i} and \eqn{\gamma_t} are state and time fixed
#' effects. \eqn{\gamma_{st}} is a state time trend (\eqn{\gamma_i \times
#' t}). The state and time fixed effects along with the state time
#' trends, yield a large number of parameters to be estimated.
#'
#' See also equation (5.2.7) from Angrist and Pischke (2009) employs
#' an individual level panel, state fixed effects, and state time trends:
#'
#' \eqn{Y_{ist} = \gamma_{0s} + \gamma_{1st} + \lambda_t + \beta D_{st} + X^\prime_{ist} \delta + \varepsilon_{ist}}
#'
#' @section FeTrendsLm$new() Method:
#' The `FeTrendsLm$new()` method initializes FeTrendsLm class and partial out variables in
#' preperation for estimation. `$new()` takes the following parameters
#' (see examples below):
#'
#' **Usage**:
#'
#' `out <- FeTrendsLm$new(DT, .f, main.reg.vars = NULL, cluster.vars = NULL,
#'                       chunk.size = 100, keycolsv = NULL,
#'                       weight.var = NULL)`
#'
#' **Parameters**:
#' \describe{
#'    \item{`DT`}{A `data.table` with the data. `DT` should have no
#'                `NA` values for variables to be used in estimation.
#'                Note that the `FeTrendsLm$new()` function
#'                **updates DT by reference**. So, `DT` will be different
#'                in the global environment after running `FeTrendsLm$new()`.
#'                Setting `data.table` `keys` based on the dimension of
#'                the fixed effects to be partialed out, will speed up
#'                estiamtion. See the data.table documenation for more.}
#'    \item{`.f`}{A right-hand-side formula with variables to be
#'                partialed out. Fixed effects must be character columns.
#'                Character vectors can be interacted with a numeric column
#'                (e.g. a time trend)}
#'    \item{`main.reg.vars`}{Optional character vector of the variables to be
#'                           used in the final regression. If `NULL`, `main.reg.vars`
#'                           with defualt to all variables that are not in `~.f`}
#'    \item{`cluster.vars`}{Optional character vector of variables to be used
#'                          in clustering standard errors. If `NULL`, robust `HC1`
#'                          standard errors will be used. If not `NULL` clustered standard
#'                          errors will be computed using the using the `HC1` correction}
#'    \item{`chunk.size`}{Integer of length 1. The size of the the chunks with which to
#'                        partial out columns. The defualt is `100`, meaning that columns
#'                        will be partialed out in batches of 100 at a time. Larger values
#'                        may be faster, but will use more memory.}
#'    \item{`keycolsv`}{A character vector of the columns which will be used as `data.table`
#'                      keys. If making multiple calls to the same dataset, `DT`, it will
#'                      faster to set the keys beforehand.}
#'    \item{`weight.var`}{A string with the variable to weight. **NOT CURRENTLY IMPLEMENTED**}
#' }
#'
#' **Return Value**:
#' An object of class `FeTrendsLm`
#'
#'
#' @section $fetrendslm():
#' The `$fetrendslm()` method estimates a regression after partialling out variables
#'
#' **Usage**:
#'
#' See examples below
#'
#' **Parameters**:
#' \describe{
#'    \item{`y.var`}{A string with the dependent, left-hand-side variable}
#'    \item{`x.vars`}{A character vector with the independent, right-hand-side variables}
#' }
#'
#' **Return Value**:
#' A list with the following items
#' \describe{
#'    \item{coef}{A numeric vector with the regression coefficients}
#'    \item{se}{A numeric vector with the regression standard errors}
#'    \item{N}{The number of observations used in the regression}
#'    \item{DT.tidy}{A `data.table` with parameter estimation output that
#'                   matches `broom::tidy(.)`}
#'    \item{lm.mod}{A small linear model, that can be used with stargazer.
#'                  Note that the standard errors will have to be applied manually
#'                  to stargazer}
#' }
#'
#' @return Object of `R6Class` with methods and data related to estiamtion.
#' @name FeTrendsLm
#' @examples
#'
#' ## From Angrist and Pischke (2015) eqn. 5.6
#' ## Difference-in-differences.
#'
#' library(lfe) ##To compare results to lfe:felm()
#' library(stargazer) ##for pretty printing of output
#' data(deaths)
#'
#' ##Dataset for felm()
#' DT.felm <- copy(deaths)
#' DT.felm <- DT.felm[, year.char := as.factor(year.char)]
#' DT.felm <- DT.felm[, state.char := as.factor(state.char)]
#'
#' ##so the original deaths dataset isn't updated by reference
#' DT1 <- copy(deaths)
#'
#' ## -- State and time fixed effects with clustering --##
#'
#' est1 <- FeTrendsLm$new(DT1, .f = ~ state.char + year.char,
#'                        cluster.vars = "state.char")
#' mod1 <- est1$fetrendslm(y.var = "mrate", x.vars = c("legal"))
#'
#' ##The parameters
#' print(mod1$coef)
#' #'The standard errors
#' print(mod1$se)
#' ##Print the number of observations
#' print(mod1$N)
#' ##print the data.table with the same outpub as broom::tidy()
#' print(mod1$DT.tidy)
#'
#' ##compare to lfe::felm() and print to stargazer
#' mod1.felm <- felm(mrate ~ legal | state.char + year.char | 0 | state.char,
#'                   DT.felm)
#'
#' ##Stargazer comparison
#' ##Notice how the number of obersvations is wrong
#' ##for mod1$lm.mod -- you must fix that by hand
#' stargazer(mod1$lm.mod, mod1.felm, type = "text",
#'           se = list(mod1$se, mod1.felm$cse),
#'           keep.stat = "n")
#'
#' ## -- State FE, time FE, State time tremds; with clustering --##
#'
#' DT2 <- copy(deaths)
#'
#' est2 <- FeTrendsLm$new(DT2, .f = ~ state.char + year.char + state.char:year,
#'                        cluster.vars = "state.char")
#' mod2 <- est2$fetrendslm(y.var = "mrate", x.vars = c("legal"))
#'
#' ##compare to lfe::felm() and print to stargazer
#' mod2.felm <- felm(mrate ~ legal | state.char + year.char + state.char:year | 0 | state.char,
#'                   DT.felm)
#'
#' ##Stargazer comparison
#' ##Notice how the number of obersvations is wrong
#' ##for mod1$lm.mod -- you must fix that by hand
#' stargazer(mod2$lm.mod, mod2.felm, type = "text",
#'           se = list(mod2$se, mod2.felm$cse),
#'           keep.stat = "n")
NULL

#' @importFrom R6 R6Class
#' @importFrom methods as
#' @importFrom stats as.formula lm pt setNames
#' @export
FeTrendsLm <- R6Class(
  classname = "FeTrendsLm",

  ## ---------------------------------- ##
  ## -- Public methods - Class Utils -- ##
  ## ---------------------------------- ##
  public = list(
    initialize = initialize,
    fetrendslm = fetrendslm,
    Y = NULL ##The matrix Y with partialed out data
  ),


  private = list(

    ## --------------------- ##
    ## -- Private Members -- ##
    ## --------------------- ##
    ..main.reg.vars = NULL, ##The regression variables to keep
    ..cluster.vars = NULL, ##The variables to using in clustering
    ..DT.cluster = NULL, ##data.table with the cluster information
    ..clusters.num = NULL, ##the number of clusters
    ..f = NULL, ##The regression formula
    ..keycolsv = NULL, ##The keys to set
    ..chunk.size = NULL, ##The Size of the Chunk
    ..weight.var = NULL, ##The weight variable as a character string
    ..weight.vals = NULL, ##The sparse weight matrix

    ..DT.po.vars = NULL, ##data.table with infomration on the partialed out vars
    ..N = NULL, ##The number of observations in the dataset DT
    ..po.ncol = 0 ##the number of partialed out columns -- initialize to 0
  ),
  active = list(
    main.reg.vars = function() private$..main.reg.vars,
    cluster.vars = function() private$..cluster.vars,
    DT.cluster = function() private$..DT.cluster,
    clusters.num = function() private$..clusters.num,
    .f = function() private$..f,
    keycolsv = function() private$..keycolsv,
    chunk.size = function() private$..chunk.size,
    weight.var = function() private$..weight.var,
    weight.mat = function() private$..weight.mat,
    DT.po.vars = function() private$..DT.po.vars,
    N = function() private$..N,
    po.ncol = function() private$..po.ncol
    )
)
