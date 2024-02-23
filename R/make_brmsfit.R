#' @title Make \code{brmsfit} Object
#' 
#' @description This is primarily used internally, wherein the \code{cmdstanr} 
#' object is converted into a \code{brmsfit} object.
#' 
#' @param object An object of class \code{chkpt_brms}
#' 
#' @param formula An object of class \code{\link[stats]{formula}}, 
#'                \code{\link[brms]{brmsformula}}, or \code{\link{brms}{mvbrmsformula}}.
#'                Further information can be found in  \code{\link[brms]{brmsformula}}.
#'                
#' 
#' @param data  An object of class \code{data.frame} (or one that can be coerced to that class) 
#'             containing data of all variables used in the model.
#' 
#' @param path Character string. The path to the folder, that is used for 
#'             saving the checkpoints.
#'             
#' @param ... Additional arguments to be passed to \code{brm}.
#' 
#' @importFrom brms brm 
#' 
#' @return An object of class \code{brmsfit}
#' 
#' @note This is primarily an internal function that constructs
#' a \code{brmsfit} object.
#' 
#' @export
make_brmsfit <- function(object, formula = NULL, data = NULL, path, ...) {
  
  if(is.null(formula)){
    formula <- object$args$formula
  }
  
  if(is.null(data)){
    data <- object$args$data
  }

  fit <- brms::brm(formula = formula,
                   data = data, 
                   empty = TRUE,
                   ...)

  file_names <- list.files(paste0(path, "/cmd_fit/"))
  
  checkpoints <- length(file_names)
  
  ordered_file_names <- paste0(path, "/cmd_fit/cmd_fit_",
                               sort(as.numeric(
                                 gsub(".*fit_(.+).rds.*", "\\1", file_names)
                               )), ".rds")
  
  stanfit_to_brms <- readRDS(ordered_file_names[1])
  
  chains <- stanfit_to_brms@sim$chains
  
  for (i in seq_len(chains)) {
    stanfit_to_brms@sim$samples[[i]] <-
      do.call(rbind, lapply(seq_len(checkpoints), function(x) {
        readRDS(ordered_file_names[x])@sim$samples[[i]]
      }))
    
  }
  
  stanfit_to_brms@sim$n_save <- stanfit_to_brms@sim$n_save * checkpoints
  
  fit$fit <- stanfit_to_brms
  
  fit <- brms::rename_pars(fit)
  
  return(fit)
}