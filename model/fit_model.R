#!/usr/bin/env Rscript

library(optparse)

here::i_am("model/fit_model.R")

# Define command line arguments
option_list <- list(
  make_option(c("-m", "--mode"), type = "character", default = "fixed",
              help = "Mode for generating the model object: 'random', 'fixed', or 'download'", 
              metavar = "MODE"),
  make_option(c("-n", "--nruns"), type = "integer", default = 1),
  make_option(c("-d", "--diagnostics"), type = "logical", default = FALSE)
)

source(here::here("model/thermaltimemodel.R"))

fixed_model_params <- tibble(
  sex = c("FEMALE", "FEMALE", "MALE", "MALE"),
  event = c("begin", "end", "begin", "end"),
  seed = c(475902884, 614687001, 1960616843, 1938696768)
)

opt <- parse_args(OptionParser(option_list = option_list))

valid_modes <- c("random", "fixed", "download")
if (!(opt$mode %in% valid_modes)) {
  stop("Error: Invalid mode. Choose from 'random', 'fixed', or 'download'.")
}

if (opt$mode == "download") {
  options(timeout = 6000)
  # This URL points to the file stored in DOI 10.5281/zenodo.14597582e
  zenodo_url <- "https://zenodo.org/records/14597582/files/modells.rds?download=1"
  output_path <- here::here("model/modells.rds")
  download.file(zenodo_url, output_path, mode = "wb")
  quit(save = "no", status = 0)
}

for (i in 1:opt$nruns) {
  message("Model run ", i, " of ", opt$nruns, " in mode: ", opt$mode)
  if (opt$mode == "random") {
    # copy fixed_model_params to a temporary and then
    # modify each of the seed values to be a randomly
    # generated integer
    random_model_params <- fixed_model_params %>%
      mutate(seed = sample(1:2^31, nrow(.), replace = TRUE))
    fit_thermaltimemodel(random_model_params)
  } else if (opt$mode == "fixed") {
    fits <- fit_thermaltimemodel(fixed_model_params)
    # write out package versions and seeds
    pkgversions <- fits[[1]]$version
    modelseeds <- list(fb = rstan::get_seed(fits[[1]]$fit),
                        fe = rstan::get_seed(fits[[2]]$fit),
                        mb = rstan::get_seed(fits[[3]]$fit),
                        me = rstan::get_seed(fits[[4]]$fit))
    saveRDS(list(vers = pkgversions, seeds = modelseeds), "output/model_meta.rds")
  }
  # if diagnostics are requested, run the diagnostics
  if (opt$diagnostics) {
    message("Running model diagnostics")
    source(here::here("analysis/1.diagnostics.R"))
  }
}
