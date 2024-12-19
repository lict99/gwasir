#' Check Model Parameters
#'
#' This function checks the validity of model parameters passed through `...`.
#' It ensures that no positional arguments are used and that certain parameters
#' specified in the `ignore` list are not included.
#'
#' @param ignore A character vector of parameter names that should be ignored.
#' @param ... Additional parameters to be checked.
#'
#' @return Returns `TRUE` if all checks pass.
#' @noRd
.chk_model_params <- function(ignore, ...) {
  input_params <- list(...)
  if (length(input_params) > 0 && is.null(names(input_params))) {
    stop("Do not support positional arguments in `...`", call. = FALSE)
  }
  if (length(input_params) > 0 && any("" %in% names(input_params))) {
    stop("Do not support positional arguments in `...`", call. = FALSE)
  }
  for (i in ignore) {
    if (i %in% names(input_params)) {
      stop("Parameter `", i, "` should be ignored", call. = FALSE)
    }
  }
  return(TRUE)
}
