#' @title Checkpoint Sampling: Stan
#' 
#' @description Fit Bayesian  models using Stan with checkpointing.
#' 
#' @param model_code Character string corresponding to the Stan model.
#' 
#' @param data A named list of R objects (like for RStan). 
#'            Further details can be found in \code{\link[cmdstanr]{sample}}.
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
#'                     is required to ensure converge (defaults to 150).
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
#' @param seed (positive integer). The seed for random number generation to 
#'             make results reproducible.
#' 
#' @param path Character string. The path to the folder, that is used for 
#'             saving the checkpoints.
#' 
#' @param ... Currently ignored.
#'
#' @return An objet of class \code{chkpt_stan}
#' 
#' @export
#' 
#' @examples 
#' \dontrun{
#' 
#' # path for storing checkpoint info
#' path <- create_folder(folder_name = "chkpt_folder_fit1")
#' 
#' stan_code <- make_stancode(bf(formula = count ~ zAge + zBase * Trt + (1|patient),
#'                               family = poisson()),
#'                            data = epilepsy)
#' stan_data <- make_standata(bf(formula = count ~ zAge + zBase * Trt + (1|patient),
#'                               family = poisson()),
#'                            data = epilepsy)
#' 
#' # "random" intercept
#' fit1 <- chkpt_stan(model_code = stan_code, 
#'                    data = stan_data,
#'                    iter_warmup = 1000,
#'                    iter_sampling = 1000,
#'                    iter_per_chkpt = 250,
#'                    path = path)
#' 
#' draws <- combine_chkpt_draws(object = fit1)
#' 
#' posterior::summarise_draws(draws)
#' 
#' 
#' # eight schools example
#'  
#' # path for storing checkpoint info
#' path <- create_folder(folder_name = "chkpt_folder_fit2")
#' 
#' stan_code <- "
#' data {
#'  int<lower=0> n;
#'   real y[n]; 
#'   real<lower=0> sigma[n]; 
#' }
#' parameters {
#'   real mu;
#'   real<lower=0> tau; 
#'   vector[n] eta; 
#' }
#' transformed parameters {
#'   vector[n] theta; 
#'   theta = mu + tau * eta; 
#' }
#' model {
#'   target += normal_lpdf(eta | 0, 1); 
#'   target += normal_lpdf(y | theta, sigma);  
#' }
#' "
#' stan_data <- schools.data <- list(
#'   n = 8,
#'   y = c(28,  8, -3,  7, -1,  1, 18, 12),
#'   sigma = c(15, 10, 16, 11,  9, 11, 10, 18)
#' )
#'
#' fit2 <- chkpt_stan(model_code = stan_code, 
#'                    data = stan_data,
#'                    iter_warmup = 1000,
#'                    iter_sampling = 1000,
#'                    iter_per_chkpt = 250,
#'                    path = path)
#' 
#' draws <- combine_chkpt_draws(object = fit2)
#' 
#' posterior::summarise_draws(draws)
#' }
chkpt_stan <- function(model_code,
                       data,
                       iter_warmup = 1000,
                       iter_sampling = 1000,
                       iter_per_chkpt = 100,
                       iter_typical = 150,
                       parallel_chains = 2,
                       threads_per = 1,
                       chkpt_progress = TRUE,
                       control = NULL,
                       seed = 1,
                       path, ...){
  
  
  
  if(!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("Please install the '", "cmdstanr", "' package.")
  }
  
  args <- c(as.list(environment()), list(...)) 
  
  if (isFALSE(is.character(path))) {
    stop("path must be a character string to parent_folder")
  }
  
  if(isFALSE(is.list(data))){
    stop("data must be a list. see examples.")
  }
  
  stan_data <- data
  
  if (isFALSE(check_for_model("model.stan", path))) {
    
    stan_code_path <- cmdstanr::write_stan_file(
      code = model_code,
      dir = paste0(path, "/stan_model"),
      basename = "model"
    )

  } else {
    stan_code_path <- paste0(path, "/stan_model/model.stan")
  }
  
  model_threads_name <- ifelse(.Platform$OS.type == "unix", 
                               "model", 
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
  
  # checkpoint test
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
    returned_object <- list(args = args)
    class(returned_object) <- "chkpt_stan"
    return(returned_object)
    
    
  } else {
    
    cp_seq <- seq(last_chkpt + 1, total_chkpts)
    
  }
  
  for (i in cp_seq) {
    
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
    }
    
    stan_state <- extract_stan_state(sample_chunk, stan_phase)
    
    saveRDS(object = stan_state,
            file = paste0(path, "/cp_info/cp_info_", i, ".rds"))
    
    cat(progress_bar(chkpt_set_up, phase = stan_phase, i = i))
    
    if (i == total_chkpts) {
      
      message("Checkpointing complete")
      
    }
    
  }
  
  returned_object <- list(args = args)
  class(returned_object) <- "chkpt_stan"
  return(returned_object)
  
}


#' @title Print \code{chkpt_stan} Objects 
#' 
#' @param x Object of class \code{chkpt_stan} 
#' @param ... Currently ignored
#' 
#' @note  
#' 
#' This function mainly avoids printing out a list.
#' 
#' Typically, after fitting, the posterior draws should be summarized with
#' \code{\link[chkptstanr]{combine_chkpt_draws}}.
#' 
#' @return No return value, and used to print the \code{chkpt_stan} object.
#' 
#' @export
print.chkpt_stan <- function(x, ...){
  cat("chkptstanr \n")
  cat("----- \n")
  cat("Date:", date(), "\n")
  cat("(see 'combine_chkpt_draws') \n")
}