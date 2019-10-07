require(gdata)
require(tidyverse)


compute_derived_quantities_by_sample <- function(control_lookup_table, smr_data) {
  
  
  derived_quantities_by_sample <- control_lookup_table %>% rowwise() # initialize
  temp_hellinger_metrics_df <- data.frame()
  for (row_num in 1:nrow(derived_quantities_by_sample)) {
    current_row <- derived_quantities_by_sample[row_num,]
    hellinger_metrics <- compute_hellinger_metrics(data = smr_data, 
                                                   id1 = current_row$expt_id,
                                                   id2 = current_row$control_id,
                                                   n_bootstrap = 100)
    
    current_row_updated <- bind_cols(current_row, as.data.frame(hellinger_metrics))
    temp_hellinger_metrics_df <- bind_rows(temp_hellinger_metrics_df, current_row_updated)
  }
  
  derived_quantities_by_sample <- temp_hellinger_metrics_df # this will be the function output
  
  derived_quantities_by_sample <- derived_quantities_by_sample %>%
    ungroup()
  
  return(derived_quantities_by_sample)
  
}
 
resample_vector <- function(vector) {
  require(gdata)
  resampled <- gdata::resample(vector, 
                               size = length(vector),
                               replace = TRUE)
  return(resampled)
}

downsample_vector <- function(vector, n) {
  require(gdata)
  if (n >= length(vector)) {
    return(vector)
  } else if (n < length(vector)) {
    downsampled <- gdata::resample(vector,
                                   size = n,
                                   replace = FALSE)
    return(downsampled)
  }
}

compute_hellinger_distribution <- function(data, id1, id2, N_ITER = 500, downsample_size = NA) {
  control_masses <- get_masses(data, id1)
  #print(paste0("length(control_masses):", length(control_masses)))
  tmz_masses <- get_masses(data, id2)
  #print(paste0("length(tmz_masses):", length(tmz_masses)))
  safe_hellinger <- possibly(custom_hellinger, NA)
  
  hellinger_distances_resampled <- c() 
  for (ii in 1:N_ITER) {
    if (is.na(downsample_size)) {
      resampled_control <- resample_vector(control_masses)
      resampled_tmz <- resample_vector(tmz_masses)
    } else if (!is.na(downsample_size)) {
      #print(paste0("downsampling to ", downsample_size))
      resampled_control <- downsample_vector(control_masses, n = downsample_size)
      resampled_tmz <- downsample_vector(tmz_masses, n = downsample_size)
    }
    
    hellinger_distance <- safe_hellinger(x = resampled_control,
                                         y = resampled_tmz)
    hellinger_distances_resampled[ii] <- hellinger_distance
  }
  
  return(hellinger_distances_resampled)
}



compute_hellinger_statistics <- function(hellinger_distances_resampled) {
  hellinger_metrics <- list(
    hellinger_dist = mean(hellinger_distances_resampled, na.rm = T),
    hellinger_sd = sd(hellinger_distances_resampled, na.rm = T),
    hellinger_05 = as.numeric(quantile(hellinger_distances_resampled, probs = 0.05, na.rm = T)),
    hellinger_95 = as.numeric(quantile(hellinger_distances_resampled, probs = 0.95, na.rm = T))
  )
  
  return(hellinger_metrics)
  
}

compute_hellinger_metrics <- function(data, id1, id2, n_bootstrap = 500) {
  print(paste0("Calculating Hellinger statistics for ", id1, "/", id2))
  hellinger_distances_resampled <- compute_hellinger_distribution(data = data,
                                                                  id1 = id1,
                                                                  id2 = id2,
                                                                  N_ITER = n_bootstrap)
  hellinger_metrics <- compute_hellinger_statistics(hellinger_distances_resampled)
  
  return(hellinger_metrics)
}

#
find_matching_control <- function(metadata, EXPT_ID) {
  
  CURRENT_SAMPLE_METADATA <- metadata %>% filter(expt_id == EXPT_ID)
  
  EXPT_INFO <- list(
    "line" = CURRENT_SAMPLE_METADATA$line,
    "treatment" = CURRENT_SAMPLE_METADATA$treatment
  )
  
  if (EXPT_INFO$treatment == "DMSO") {
    print(sprintf('%s is already a control sample', EXPT_ID))
    matching_control_id <- NA
    return(matching_control_id)
  }
  
  CURRENT_CONTROL_METADATA <- metadata %>%
    filter(line == EXPT_INFO$line,
           treatment == "DMSO")
  
  if (nrow(CURRENT_CONTROL_METADATA) == 0) {
    print(sprintf('No control found matching %s', EXPT_ID))
    matching_control_id <- NA
  } else if (nrow(CURRENT_CONTROL_METADATA) >= 2) {
    print(sprintf('Found two control experiments matching %s', EXPT_ID))
    matching_control_id <- NA
  } else if (nrow(CURRENT_CONTROL_METADATA) == 1) {
    matching_control_id <- CURRENT_CONTROL_METADATA$expt_id
  }
  
  return(matching_control_id)
}

#
match_control_to_treatment <- function(smr_metadata) {
  drug_metadata <- smr_metadata %>% 
    filter(treatment != "DMSO")
  
  control_lookup_table <- drug_metadata %>%
    rowwise() %>%
    mutate(control_id = find_matching_control(metadata = smr_metadata, EXPT_ID = expt_id)) %>%
    ungroup() %>%
    select(line, treatment, expt_id, control_id)
  
  return(control_lookup_table)
}

#
get_masses <- function(data, expt_id) {
  masses <- data[data$expt_id == expt_id, "mass"]
  return(masses)
}

#
custom_hellinger <- function (x, y, lower = 0, upper = 1000, method = 1, ...) 
{
  require(statip)
  
  fx <- densityfun(x, ...)
  fy <- densityfun(y, ...)
  if (method == 1) {
    g <- function(z) (fx(z)^0.5 - fy(z)^0.5)^2
    stats::integrate(g, lower, upper, subdivisions=500L)$value/2
  }
  else if (method == 2) {
    g <- function(z) (fx(z) * fy(z))^0.5
    1 - stats::integrate(g, lower, upper, subdivisions=500L)$value
  }
  else {
    stop("incorrect 'method' argument")
  }
}