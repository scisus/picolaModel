library(assertthat)
library(dplyr)
# functions
#
# calculate censoring status of each observation. left, right, and interval. Assumes no uncensored data - uncensored data will be treated as left or right censored - and that all trees flowered. Returns a dataframe with identifying information for each observation and a censorship code. takes wide form of data
add_censor_indicator <- function(phenevent) {

  #calculate indexes of trees with before flowering or past flowering observations

  have1 <- unique(phenevent[which(phenevent$State == 1), "Index"]) # trees with an event 1/state 1 obs
  have4 <- unique(phenevent[which(phenevent$State == 3), "Index"]) # trees with an event 4/state 3 obs

  # add censoring labels
  censorind <- phenevent %>%
    mutate(censored = case_when(State %in% c(1,3) ~ "interval", # because all trees have a flowering observation, all state 1 and 3 observations are interval censored
                                # if no event/state 1 was observed, then event 2 obs are left censored, otherwise they are interval censored
                                Event_Obs == 2 & !Index %in% have1 ~ "left",
                                Event_Obs == 2 & Index %in% have1 ~ "interval",
                                # if no event 4/state 3 was observed, then event 3 obs are left censored, otherwise they are interval censored
                                Event_Obs == 3 & !Index %in% have4 ~ "right",
                                Event_Obs == 3 & Index %in% have4 ~ "interval"))

  assertthat::are_equal(which(is.na(censorind$censored)), integer(0))
  assertthat::are_equal(nrow(censorind), nrow(phenevent))

  return(censorind)
}


# munge phenology data. Remove duplicate observations, add censoring information, combine with forcing data, and standardize sum_forcing. phendat should be flowers::lodgepole_phenology_event or structured similarly. clim is a dataframe of daily climate with Site + Year + DoY columns that can be matched against the ones in phendat. spu is a dataframe of information on spus with a provenance column SPU_Name and Orchard column that can be matched against the one in phendat
prepare_data <- function(phendat, clim, spu) {
  # 4 trees were observed by both Wagner and Walsh at PGTIS in 2006 - drop 1 copy of them (16 duplicate observations).
  rmidx <- phendat %>%
    group_by(Index) %>%
    summarize(count = n()) %>%
    filter(count > 4)

  spu <- select(spu, SPU_Name, Orchard) # drop unnecessary columns

  #add information about censoring
  phen <- phendat %>%
    filter(! (Index %in% rmidx$Index & Source == "Rita Wagner")) %>% # Drop duplicates
    add_censor_indicator() %>% # add censoring type
    mutate(censored_lronly = case_when(censored == "interval" ~ "none",
                                       censored %in% c("left", "right") ~ censored)) %>% # exclude interval censoring
    # add bound labels for interval censoring models
    mutate(bound = case_when(Event_Obs == 1 | Event_Obs == 2 & censored == "left" ~ "lower",
                             Event_Obs == 2 & censored == "interval" ~ "upper",
                             Event_Obs == 3 ~ "lower",
                             Event_Obs == 4 ~ "upper"))

  phenf <- phen %>%
    dplyr::left_join(clim) %>%
    dplyr::left_join(spu) %>%
    dplyr::mutate(Year = as.character(Year), Genotype = as.character(Genotype)) %>%
    dplyr::rename(Provenance = SPU_Name) %>%
    distinct() %>%

  return(phenf)
}

# filter phenology data by sex and event in order to build models for male and female, begin and end. sex is one of "FEMALE" or "MALE" event is one of "begin" or "end, and dat is output from `prepare_data` or similarly structured
filter_sex_event <- function(sex, event, dat = phenf) {
  assertthat::assert_that(sex %in% c("MALE", "FEMALE"))
  assertthat::assert_that(event %in% c("begin", "end"))

  if (event == "begin") {
    event_obs <- c(1,2)
  } else {
    event_obs <- c(3,4)
  }

  dat_prepped <- dat %>%
    dplyr::filter(Sex == sex & Event_Obs %in% event_obs) %>%
    #dplyr::select(-DoY, -Date, -State, -contains("Event"), -mean_temp, -forcing, -sum_forcing) %>%
    dplyr::select(-DoY, -Date, -State, -contains("Event"), -mean_temp, -forcing) %>%
    distinct() %>%
   # tidyr::pivot_wider(names_from = bound, values_from = sum_forcing_centered, values_fill = 0) %>%
   tidyr::pivot_wider(names_from = bound, values_from = sum_forcing, values_fill = 0) %>%
    dplyr::rename(sum_forcing = lower) %>%
    dplyr::mutate(event = event) %>%
    filter(!is.na(MAT))

  return(dat_prepped)
}

# convenience function for gathering population mean draws from the model `mod`
gather_intercept_draws <- function(mod) {
  draws <- mod %>% tidybayes::gather_draws(b_Intercept, ndraws = nsamp, seed = seed)
  return(draws)
}

gather_slope_draws <- function(mod) {
  draws <- mod %>% tidybayes::gather_draws(b_MAT, ndraws = nsamp, seed = seed)
  return(draws)
}

# convenience function for gathering generation parameters draws from the model `mod`
gather_gen_draws <- function(mod) {
  draws <- mod %>% spread_draws(`.*moGen.*`, regex = TRUE, ndraws = nsamp, seed = seed) %>%
    select(-contains("prior"))
  return(draws)
}

# convenience function for gathering sd (population and offset) draws from the model `mod`
gather_var_draws <- function(mod) {
  draws <- mod %>% gather_draws(`sd_.*`, `sigma`, regex = TRUE, ndraws = nsamp, seed = seed)
  return(draws)
}

# convenience function for gathering offset delta draws from the model `mod`
gather_offset_draws <- function(mod) {
  draws <- mod %>% gather_draws(`r_.*`, regex = TRUE, ndraws = nsamp, seed = seed)
  return(draws)
}

# prepare dataframes a and b for interval finding by splitting into lists and ensuring that the climate list (a) and the phenology list (b) contain information for the same sites and years. a and b must both have Site and Year columns
split_df_to_lists <- function(a, b) {
  lista <- split(a, f = list(a$Site, a$Year), drop = TRUE)
  listb <- split(b, f = list(b$Site, b$Year), drop = TRUE)

  # check that all sites and years in the phenology frame are in the climate frame
  assertthat::assert_that(all(names(listb) %in% names(lista))) # all entries in B must be in A

  # subset lista so it only contains Site x Year that also occur in B
  ainb <- lista[names(listb)]

  assertthat::are_equal(names(ainb), names(listb))

  return(list(listainb = ainb, listb = listb))
}

# a = a,  b= forcing for translation, aforce = name of forcing column in a, bforce = name of forcing columb in b
forcing_to_doy <- function(a, b, aforce, bforce, newdoycolname) {
  # prepare dataframes for interval finding by splitting into lists
  splitdfs <- split_df_to_lists(a, b)

  df <- purrr::map2(splitdfs$listainb, splitdfs$listb, find_day_of_forcing, aforce = aforce, bforce = bforce) %>% # find DoY in A corresponding to each sum_forcing in B
    purrr::map_dfr(bind_rows) # combine into a single dataframe

  names(df)[which(names(df) == "newdoycol")] <- newdoycolname

  return(df)
}

# given dataframes adf (climate) and bdf (phenology) identify the day of year in adf corresponding to reaching each sum_forcing bcol in bdf. adf must have a sum_forcing column identified by name aforce and a DoY column and b must have a sum_forcing column identified with name bforce.
find_day_of_forcing <- function(adf, bdf, aforce, bforce) {

  # what row in a contains the interval for entries in b. Add 1 to the index because phenological events require the threshold to be reached. this introduces error, but is unavoidable in some direction.
  a_index <- findInterval(bdf[[bforce]], adf[[aforce]]) + 1


  # add a column to b for Day of Year and extract the correct Day of year from a using the index
  bdf$newdoycol <- adf$DoY[a_index]

  # when sum_forcing in b is exactly identical to sum_forcing in b, a_index will (incorrectly) be +1 day. Re-write those to be the correct day - the day the sum_forcing is reached
  identical_forcing_index <- which(bdf[[bforce]] %in% adf[[aforce]])
  bdf$newdoycol[identical_forcing_index] <- bdf$newdoycol[identical_forcing_index] - 1

  return(bdf)
}

# useful when you have to do multiple layers of indexing
find_day_of_forcing_mapper <- function(alist, bdf, bforce = ".epred") {
  message("running find_day_of_forcing_mapper for ", alist$Site[1])
  res <- map_dfr(split(alist, f = list(alist$index), drop = TRUE),
          find_day_of_forcing, .id = ".id",
          bdf = bdf, aforce = "sum_forcing", bforce = bforce,
          .progress = TRUE)
  gc()
  return(res)
}
