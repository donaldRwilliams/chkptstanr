#' @importFrom utils capture.output packageVersion

# find so-called `typical set`
chkpt_typical <- function(model,
                          cmd_args,
                          control,
                          iter_typical,
                          seed,
                          progress) {
  
  cmd_args$iter_sampling <- 0
  cmd_args$seed <- seed
  cmd_args$save_warmup <- TRUE
  cmd_args$iter_warmup <- iter_typical
  
  if (progress) {
    capture_output <- capture.output({
      suppressMessages({
        sample_chunk <- do.call(what = model,
                                args = c(cmd_args, control))
        
      })
      
    })
    
  } else {
    sample_chunk <- do.call(what = model,
                            args = c(cmd_args, control))
    
    
  }
  
  return(sample_chunk)
}


# check point sampling
chkpt_sample <- function(model, 
                         cmd_args, 
                         cp_cmd_args, 
                         control, 
                         progress){
  
  if(progress){
    
    capture_output <- capture.output({
      suppressMessages({
        sample_chunk <-
          do.call(model, args = c(cmd_args, cp_cmd_args, control))
      })
    })
    
  } else {
    
    sample_chunk <-
      do.call(model, args = c(cmd_args, cp_cmd_args, control))
  }
  
  return(sample_chunk)
}


# sampling cmdstanr arguments
cp_cmd_args <- function(seed,
                        phase,
                        stan_state,
                        iter_per_chkpt){
  
  if (phase == "warmup") {
    
    returned_object <- list(
      iter_sampling = 0,
      seed = seed,
      save_warmup = TRUE,
      iter_warmup = iter_per_chkpt,
      adapt_engaged = FALSE,
      init = stan_state$inits,
      inv_metric = stan_state$inv_metric,
      step_size = stan_state$step_size_adapt
    )
    
  } else if (phase == "sample") {
    
    returned_object <- list(
      iter_sampling = iter_per_chkpt,
      seed = seed,
      save_warmup = FALSE,
      iter_warmup = 0,
      adapt_engaged = FALSE,
      init = stan_state$inits,
      inv_metric = stan_state$inv_metric,
      step_size = stan_state$step_size_adapt
    )
    
  } else {
    
    stop("phase must be warmup or sample")
    
  }
  
  return(returned_object)
  
}

# matrix initial values
get_init <- function(draws, max_iter, chain){
  
  last_draw <- draws[max_iter, chain, -1]
  
  param_names <- dimnames(last_draw)$variable
  
  flat_draws <- as.numeric(last_draw)
  
  unique_names <- unique(gsub(pattern = "\\[.*?\\]",
                              replacement =  "",
                              x = param_names))
  
  returned_object <-
    lapply(seq_along(unique_names), function(x) {
      ids <- grepl(pattern = paste0("\\b", unique_names[x] , "\\b"),
                   x = param_names)
      
      bracket_check <-
        regmatches(x = param_names[ids],
                   m = gregexpr(pattern = "\\[.*?\\]",
                                text = param_names[ids]))
      
      if (length(unlist(bracket_check)) == 0) {
        
        return(flat_draws[ids])
        
      } else {
        
        indices <- sapply(strsplit(
          gsub(
            pattern = "\\[|\\]",
            replacement = "",
            x = bracket_check
          ),
          split = ","
        ), as.numeric)
        
        if (is(indices, "matrix")) {
          
          return(matrix(
            data = flat_draws[ids],
            nrow = max(indices[1, ]),
            ncol = max(indices[2, ])
          ))
          
        } else {
          
          return(flat_draws[ids])
          
        }
        
      }
      
    })
  
  names(returned_object) <- unique_names
  
  return(returned_object)
  
}

check_for_model <- function(x, path) {
  return(any(x == list.files(paste0(
    path, "/stan_model"
  ))))
}

is_zero <- function(x){
  return(length(x) == 0)
}

check_restart <- function(initial_args, restart_args) {
  return(isTRUE(suppressWarnings({
    all.equal(initial_args, restart_args)
  })))
  
}


stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}
