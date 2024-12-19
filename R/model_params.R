#' Create Model Parameters
#'
#' This function creates a `model_params` object.
#'
#' @param data A data frame containing the data.
#' @param id_col A string specifying the column name in `data` that contains
#' the IDs.
#' @param data_col A string specifying the column name in `data` that contains
#' the data.
#' @param fun A function to be used for modeling.
#' @param fun_params A list of parameters to be passed to `fun`.
#' @param header A character vector specifying the header names.
#' @param gt A string specifying the genotype. Default is "genotype".
#'
#' @return A list with class "model_params" containing the following elements:
#' \describe{
#'   \item{id}{The IDs extracted from `data` based on `id_col`.}
#'   \item{gt}{The genotype string.}
#'   \item{data}{The data extracted from `data` based on `data_col`.}
#'   \item{fun}{The function to be used for modeling.}
#'   \item{fun_params}{The list of parameters to be passed to `fun`.}
#'   \item{header}{The header names.}
#' }
#'
#' @details todo!
#'
#' @examples
#' data <- data.frame(id = letters[1:5], x = rnorm(5), y = rnorm(5))
#' fun <- function(formula, data, ...) {
#'   res <- extract_lm(formula, data, ...)
#'   return(data.frame(coef, se, p))
#' }
#' fun_params <- list(formula = y ~ genotype + x, data = NULL)
#' header <- c("coef", "se", "p")
#' mk_model_params(data, "id", c("x", "y"), fun, fun_params, header)
#'
#' @export
mk_model_params <- function(
    data,
    id_col,
    data_col,
    fun,
    fun_params,
    header,
    gt = "genotype") {
  if (!is.character(gt) || length(gt) != 1L) {
    stop("Parameter `gt` should be a character", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("Parameter `data` should be a data frame", call. = FALSE)
  }
  if (!is.function(fun)) {
    stop("Parameter `fun` should be a function", call. = FALSE)
  }
  if (!is.list(fun_params)) {
    stop("Parameter `fun_params` should be a list", call. = FALSE)
  }
  if (!is.character(header) || length(header) == 0L) {
    stop("Parameter `header` should be a character vector", call. = FALSE)
  }
  return(
    structure(
      list(
        id = data[[id_col]],
        gt = gt,
        data = data[data_col],
        fun = fun,
        fun_params = fun_params,
        header = header
      ),
      class = "model_params"
    )
  )
}
