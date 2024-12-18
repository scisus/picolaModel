#!/usr/bin/env Rscript

# When running this script in Rstudio or on windows, the 4 model fits will be
# sequential. To fit them in parallel, run this script from the command line.

here::i_am("run_all.R")

# Find all .rds files, recursively, and delete them.
file.remove(list.files(path = ".", pattern = ".rds$", recursive = TRUE))

# Choose which mode to fit the model in here if running interactively.
cmdargs = c("--mode=fixed")
# cmdargs = c("--mode=random")
# cmdargs = c("--mode=download")

# If running from the command line, you can pass --mode=(fixed|random|download)
# to customize how the model is fit. This is ignored if running interactively.
if (!interactive()) {
    # Parse command-line arguments if running as script.
    library(optparse)
    option_list <- list(make_option(c("-m", "--mode"), type = "character", default = "fixed",
                    help = "Mode for generating the model object: 'random', 'fixed', or 'download'",
                    metavar = "MODE")
    )
    opt <- parse_args(OptionParser(option_list=option_list))
    valid_modes <- c("random", "fixed", "download")
    if (!(opt$mode %in% valid_modes)) {
        stop("Error: Invalid mode. Choose from 'random', 'fixed', or 'download'.")
    }
    cmdargs <- c(paste0("--mode=", opt$mode))
}

# Prepare inputs.
callr::rscript(here::here("inputs/calc_forcing.R"))

# Fit (or download) the model.
callr::rscript(script=here::here("model/fit_model.R"), cmdargs = cmdargs)

# Run the analysis.
callr::rscript(here::here("analysis/run_all_analysis.R"))

# Produce outputs for use elsewhere.
callr::rscript(here::here("output/build_all_outputs.R"))
