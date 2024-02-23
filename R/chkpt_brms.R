#' @title Checkpoint Sampling: brms
#' 
#' @description Fit Bayesian generalized (non-)linear multivariate multilevel models using brms
#'              with checkpointing.
#' 
#' @param formula An object of class \code{\link[stats]{formula}}, 
#'                \code{\link[brms]{brmsformula}}, or \code{\link{brms}{mvbrmsformula}}.
#'                Further information can be found in  \code{\link[brms]{brmsformula}}.
#'                
#' @param data An object of class \code{data.frame} (or one that can be coerced to that class) 
#'             containing data of all variables used in the model.
#' 
#' @param iter_warmup (positive integer) The number of warmup iterations to run 
#'                    per chain (defaults to 1000).
#' 
#' @param iter_sampling (positive integer) The number of post-warmup iterations 
#'                      to run per chain (defaults to 1000).
#' 
#' @param iter_per_chkpt (positive integer). The number of iterations per 
#'                        checkpoint. Note that \code{iter_sampling} is divided
#'                        by \code{iter_per_chkpt} to determine the number of
#'                        checkpoints. This must result in an integer 
#'                        (if not, there will be an error).
#' 
#' @param iter_typical (positive integer) The number of iterations in the 
#'                     initial warmup, which finds the so-called typical set.
#'                     This is an initial phase, and not included in 
#'                     \code{iter_warmup}. Note that a large enough value
#'                     is required to ensure convergence (defaults to 150).
#'                     
#' 
#' @param parallel_chains (positive integer) The \emph{maximum number} of MCMC 
#'                        chains to run in parallel. If parallel_chains is not 
#'                        specified then the default is to look for the option 
#'                        \code{mc.cores}, which can be set for an entire R session by 
#'                        \code{options(mc.cores=value)}. If the \code{mc.cores} 
#'                        option has not been set then the default is \code{1}.
#' 
#' @param threads_per (positive integer) Number of threads to use in within-chain 
#'                     parallelization (defaults to \code{1}).
#' 
#' @param chkpt_progress logical. Should the \code{chkptstanr} progress 
#'                       be printed (defaults to \code{TRUE}) ? If set to 
#'                       \code{FALSE}, the standard \code{cmdstanr} progress
#'                       bar is printed for each checkpoint 
#'                       (which does not actually keep track of 
#'                       checkpointing progress)
#' 
#' @param control A named list of parameters to control the sampler's behavior. 
#'                It defaults to NULL so all the default values are used.
#'                For a comprehensive overview see \code{\link[rstan]{stan}}.
#' 
#' @param brmsfit Logical. Should a \code{brmsfit} object be returned 
#'                (defaults to \code{TRUE}). 
#' 
#' @param seed (positive integer). The seed for random number generation to 
#'             make results reproducible.
#'             
#' @param stop_after (positive integer). The number of checkpoints to sample
#'                  before stopping. If \code{NULL}, then all checkpoints are
#'                  sampled (defaults to \code{NULL}). Used only for testing
#' 
#' @param path Character string. The path to the folder, that is used for 
#'             saving the checkpoints.
#' 
#' @param ... Additional arguments based to \code{\link[brms]{make_stancode}},
#'            including, for example, user-defined prior distributions and the 
#'            \code{\link[brms]{brmsfamily}} (e.g., \code{family = poisson())}.
#'
#' @return An object of class \code{brmsfit} (with \code{brmsfit = TRUE}) 
#'         or \code{chkpt_brms} (with \code{brmsfit = FALSE}).
#'
#' @importFrom brms make_standata make_stancode threading bf
#' @importFrom rstan read_stan_csv
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(brms)
#' library(cmdstanr)
#' 
#' # path for storing checkpoint info
#' path <- create_folder(folder_name  = "chkpt_folder_fit1")
#' 
#' # "random" intercept
#' fit1 <- chkpt_brms(bf(formula = count ~ zAge + zBase * Trt + (1|patient),
#'                       family = poisson()), 
#'                    data = epilepsy, , 
#'                    iter_warmup = 1000, 
#'                    iter_sampling = 1000, 
#'                    iter_per_chkpt = 250, 
#'                    path = path)
#'                    
#' # brmsfit output
#' fit1
#' 
#' # path for storing checkpoint info
#'  path <- create_folder(folder_name  = "chkpt_folder_fit2")
#' 
#' # remove "random" intercept (for model comparison)
#' fit2 <- chkpt_brms(bf(formula = count ~ zAge + zBase * Trt, 
#'                       family = poisson()), 
#'                    data = epilepsy, , 
#'                    iter_warmup = 1000, 
#'                    iter_sampling = 1000, 
#'                    iter_per_chkpt = 250, 
#'                    path = path)
#'                    
#' # brmsfit output
#' fit2
#' 
#' # compare models
#' loo(fit1, fit2)
#' 
#' 
#' # using custom priors
#' path <- create_folder(folder_name = "chkpt_folder_fit3")
#' 
#' # priors
#' bprior <- prior(constant(1), class = "b") +
#'   prior(constant(2), class = "b", coef = "zBase") +
#'   prior(constant(0.5), class = "sd")
#' 
#' # fit model
#' fit3 <-
#'   chkpt_brms(
#'     bf(
#'       formula = count ~ zAge + zBase + (1 | patient),
#'       family = poisson()
#'     ),
#'     data = epilepsy,
#'     path  = path,
#'     prior = bprior,
#'     iter_warmup = 1000,
#'     iter_sampling = 1000,
#'     iter_per_chkpt = 250, 
#'   )
#' 
#' 
#' # check priors
#' prior_summary(fit3)
#' 
#' }
chkpt_brms <- function(formula,
                       data,
                       iter_warmup = 1000,
                       iter_sampling = 1000,
                       iter_per_chkpt = 100,
                       iter_typical = 150,
                       parallel_chains = 2,
                       threads_per = 1,
                       chkpt_progress = TRUE,
                       control = NULL,
                       brmsfit = TRUE,
                       seed = 1,
                       stop_after = NULL,
                       path, ...){
  
  
  args <- c(as.list(environment()), list(...)) 
  
  if(!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("Please install the '", "cmdstanr", "' package.")
  }
  
  if(isFALSE(requireNamespace("brms", quietly = TRUE))) {
    stop("Please install the '", "brms", "' package.")
  }

  if (isFALSE(is(formula, "formula") |
              is(formula, "brmsformula"))) {
    stop("formula must be of class formula or brmsformula")
  }

  if (isFALSE(is.character(path))) {
    stop("path must be a character string to parent_folder")
  }
  
  if (threads_per == 1) {
    
    stan_data <- brms::make_standata(formula = formula, data = data, ...)
    
    stan_code <- brms::make_stancode(formula = formula, data = data, ...)
    
  } else {
  
    stan_data <- brms::make_standata(formula = formula, 
                                     data = data,  
                                     threads = brms::threading(threads_per),
                                     ...)
    
    stan_code <- brms::make_stancode(
      formula = formula,
      data = data,
      threads = brms::threading(threads_per),
      ...
    )
  }
  
  stan_code_path <- paste0(path, "/stan_model/model.stan")
  
  if (isFALSE(check_for_model("model.stan", path))) {

    stan_code_path <- cmdstanr::write_stan_file(
      code = stan_code,
      dir = paste0(path, "/stan_model"),
      basename = "model"
    )

  }
  
  model_threads_name <- ifelse(.Platform$OS.type == "unix", 
                               "model_threads", 
                               "model_threads.exe")
  
  if (isFALSE(check_for_model(eval(model_threads_name), path))) {
    stan_m3 <- cmdstanr::cmdstan_model(stan_file = stan_code_path,
                                       cpp_options = list(stan_threads = TRUE))
    saveRDS(stan_m3, file = paste0(path, "/stan_model/comp.rds"))
    saveRDS(args, paste0(path, "/stan_model/args.rds"))
  } else {

    initial_args <- readRDS(paste0(path, "/stan_model/args.rds"))

    if (isFALSE(check_restart(args, initial_args))) {
      stop("invalid restart (arguments have been changed)",
           call. = FALSE)
    }

    stan_m3 <- readRDS(paste0(path, "/stan_model/comp.rds"))

    message("Sampling next checkpoint")
  }

  chkpt_set_up <- chkpt_setup(iter_warmup = iter_warmup,
                              iter_sampling = iter_sampling,
                              iter_per_chkpt = iter_per_chkpt)

  warmup_chkpts <- chkpt_set_up$warmup_chkpts

  total_chkpts <- chkpt_set_up$total_chkpts

  cmd_args <- list(data = stan_data,
                   chains = parallel_chains,
                   parallel_chains = parallel_chains,
                   threads_per_chain = threads_per)

  if (is_zero(list.files(paste0(path, "/cp_info")))) {

    last_chkpt <- 0

    message("Initial Warmup (Typical Set)")

    sample_chunk <- chkpt_typical(
      seed = seed,
      control = control,
      model = stan_m3$sample,
      iter_typical = iter_typical,
      cmd_args = cmd_args,
      progress = chkpt_progress
    )

    stan_state <- extract_stan_state(sample_chunk, "warmup")

  } else {

    cp_fn_list <- list.files(paste0(path, "/cp_info"))

    checkpoints <- as.numeric(gsub(".*info_(.+).rds.*", "\\1", cp_fn_list))

    last_chkpt <- max(checkpoints)

    stan_state <-  readRDS(file = paste0(path, "/cp_info/",
                           cp_fn_list[which.max(checkpoints)]))

  }

  if (last_chkpt == total_chkpts) {
   
     message("Checkpointing complete")
    
    
     if (brmsfit) {
      
      if (is.null(args$prior)) {
        
        returned_object <- make_brmsfit(formula = formula,
                                        data = data,
                                        path = path)
        
      } else {
        
        returned_object <- make_brmsfit(
          formula = formula,
          prior = args$prior,
          data = data,
          path = path
        )
      }
      
    } else {
      
      returned_object <- list(args = args)
      class(returned_object) <- "chkpt_brms"
      
    }
  
  return(returned_object)
  
  } else {

    cp_seq <- seq(last_chkpt + 1, total_chkpts)

  }

  for (i in cp_seq) {
    
    if (!is.null(stop_after) && i > stop_after) {
      message("Stopping after ", stop_after, " checkpoints")
      stop_quietly()
    }

    if (i <= warmup_chkpts) {

      stan_phase <- "warmup"

      sample_chunk <- chkpt_sample(
        model = stan_m3$sample,
        cmd_args = cmd_args,
        control = control,
        progress = chkpt_progress,
        cp_cmd_args =  cp_cmd_args(
          seed = seed + i,
          phase = stan_phase,
          stan_state = stan_state,
          iter_per_chkpt = iter_per_chkpt
        )
      )

    } else {

      stan_phase <- "sample"

      sample_chunk <- chkpt_sample(
        model = stan_m3$sample,
        cmd_args = cmd_args,
        control = control,
        progress = chkpt_progress,
        cp_cmd_args =  cp_cmd_args(
          seed = seed + i,
          phase = stan_phase,
          stan_state = stan_state,
          iter_per_chkpt = iter_per_chkpt
        )
      )

      sample_chunk$save_object(paste0(path, "/cp_samples/",
                                      "samples_", i, ".rds"))

      # this is optional
      if (brmsfit) {

        saveRDS(
          object = rstan::read_stan_csv(sample_chunk$output_files()),
          file = paste0(path, "/cmd_fit/cmd_fit_", i, ".rds")
        )
      }

    }

    stan_state <- extract_stan_state(sample_chunk, stan_phase)

    saveRDS(object = stan_state,
            file = paste0(path, "/cp_info/cp_info_", i, ".rds"))

    cat(progress_bar(chkpt_set_up, phase = stan_phase, i = i))

    if (i == total_chkpts) {

      message("Checkpointing complete")

    }

  }

  if (brmsfit) {
    
    if (is.null(args$prior)) {
      
      returned_object <- make_brmsfit(formula = formula,
                                      data = data,
                                      path = path)
      
    } else {
      returned_object <- make_brmsfit(
        formula = formula,
        prior = args$prior,
        data = data,
        path = path
      )
    }
    
  } else {
    
    returned_object <- list(args = args)
    class(returned_object) <- "chkpt_brms"
    
  }
  
  returned_object$path <- path
  return(returned_object)

}


#' @title Print \code{chkpt_brms} Objects 
#' 
#' @param x Object of class \code{chkpt_brms} 
#' @param ... Currently ignored
#' 
#' @note  
#' 
#' This function mainly avoids printing out a list,
#' and it is only used when \code{brmsfit = "FALSE"} in \code{\link[chkptstanr]{chkpt_brms}}.
#' 
#' Typically, after fitting, the posterior draws should be summarized with
#' \code{\link[chkptstanr]{combine_chkpt_draws}} (assuming \code{brmsfit = "FALSE"}).
#' 
#' @return No return value, and used to print the \code{chkpt_brms} object.
#' 
#' @export
print.chkpt_brms <- function(x, ...){
  cat("chkptstanr \n")
  cat("----- \n")
  cat("Date:", date(), "\n")
  cat("(see 'combine_chkpt_draws') \n")
}
