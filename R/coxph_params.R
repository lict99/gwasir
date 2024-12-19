#' Create Cox Proportional Hazards Model Parameters
#'
#' This function prepares the parameters required to fit a Cox proportional
#' hazards model using the [survival::coxph()] function.
#'
#' @param data A data frame containing the ids, covariates and survival data.
#' @param id_col A string specifying the column name for the subject IDs.
#' @param time_col A string specifying the column name for the survival time.
#' @param event_col A string specifying the column name for the event indicator
#' (1 if event occurred, 0 otherwise).
#' @param covariates_col A vector of strings specifying the column names for the
#' covariates.
#' @param ... Additional arguments to be passed to the [survival::coxph()]
#' function except for the `formula` and `data` arguments.
#'
#' @return A `model_params` object.
#'
#' @examples
#' df <- data.frame(
#'   sample_id = 1:10,
#'   time = rexp(10),
#'   event = sample(0:1, 10, replace = TRUE),
#'   age = rnorm(10, 50, 10),
#'   treatment = sample(0:1, 10, replace = TRUE)
#' )
#' mk_coxph_params(
#'   data = df,
#'   id_col = "sample_id",
#'   time_col = "time",
#'   event_col = "event",
#'   covariates_col = c("age", "treatment")
#' )
#' @export
mk_coxph_params <- function(
    data,
    id_col,
    time_col,
    event_col,
    covariates_col,
    ...) {
  .chk_model_params(c("formula", "data"), ...)
  extract_coxph <- function(formula, data, ...) {
    fit <- survival::coxph(formula = formula, data = data, ...)
    fit_smr <- summary(fit)
    return(
      data.frame(
        coef = fit_smr[["coefficients"]]["genotype", "coef"],
        se = fit_smr[["coefficients"]]["genotype", "se(coef)"],
        p_value = fit_smr[["coefficients"]]["genotype", "Pr(>|z|)"],
        hr = fit_smr[["conf.int"]]["genotype", "exp(coef)"],
        hr_l95 = fit_smr[["conf.int"]]["genotype", "lower .95"],
        hr_u95 = fit_smr[["conf.int"]]["genotype", "upper .95"],
        n = fit_smr[["n"]],
        n_event = fit_smr[["nevent"]]
      )
    )
  }

  fml <- stats::as.formula(
    paste(
      sprintf("survival::Surv(time = %s, event = %s)", time_col, event_col),
      paste(c("genotype", covariates_col), collapse = " + "),
      sep = " ~ "
    ),
    env = new.env()
  )

  return(
    mk_model_params(
      data = data,
      id_col = id_col,
      data_col = c(time_col, event_col, covariates_col),
      fun = extract_coxph,
      fun_params = c(list(formula = fml, data = NULL), list(...)),
      header = c(
        "coef", "se", "p_value", "hr", "hr_l95", "hr_u95", "n", "n_event"
      ),
      gt = "genotype"
    )
  )
}
