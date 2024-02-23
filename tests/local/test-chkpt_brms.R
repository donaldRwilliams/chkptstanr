test_that('chkpt_brms picks up after stopping', {
  path <- create_folder(folder_name  = "tests/local/chkpt_folder_m1", 
                        path = here::here())
  # clean up
  on.exit(unlink(path, recursive = TRUE))
  
  # simplified example from vignete for faster execution
  bf_m1 <- brms::bf(formula = count ~ zAge + zBase,
              family = poisson())
  
  cat("\n\nRunning for 1 checkpoint then stopping\n\n")
  
  # run for 1 checkpoints then stop
  fit_m1 <- try(chkpt_brms(
    formula = bf_m1,
    data = brms::epilepsy,
    path  = path,
    iter_warmup = 100,
    iter_sampling = 200,
    iter_per_chkpt = 100,
    stop_after = 1
  ), silent = T)
  
  cat("\n\nTrying to pick up where we stopped\n\n")
  
  res <- try(chkpt_brms(
    formula = bf_m1,
    data = brms::epilepsy,
    path  = path,
    iter_warmup = 100,
    iter_sampling = 200,
    iter_per_chkpt = 100,
  ), silent=T)
  
  expect_false(is(res, 'try-error'))
})


test_that('chkpt_brms works with data2', {
  # generate data 
  # from: https://cran.r-project.org/web/packages/brms/vignettes/brms_phylogenetics.html
  phylo <- ape::read.nexus("https://paul-buerkner.github.io/data/phylo.nex")
  data_simple <- read.table(
    "https://paul-buerkner.github.io/data/data_simple.txt",
    header = TRUE
  )
  A <- ape::vcv.phylo(phylo)
  
  path <- create_folder(folder_name  = "tests/local/chkpt_folder_m2",
                        path = here::here())
  
  # clean up
  on.exit(unlink(path, recursive = TRUE))
  
  prior <- c(
    brms::prior(normal(0, 10), "b"),
    brms::prior(normal(0, 50), "Intercept"),
    brms::prior(student_t(3, 0, 20), "sd"),
    brms::prior(student_t(3, 0, 20), "sigma")
  )
  
  checkpoint_model <- try(chkpt_brms(
    phen ~ cofactor + (1|gr(phylo, cov = A)),
    data = data_simple,
    family = gaussian(),
    data2 = list(A = A),
    prior = prior,
    path  = path,
    iter_warmup = 100,
    iter_sampling = 200,
    iter_per_chkpt = 100
  ), silent = T)
  
  expect_false(is(checkpoint_model, 'try-error'))
})
