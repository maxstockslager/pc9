###########################################################################################################
read_google_sheet <- function(filename, sheetname) {
  
  require(googlesheets)
  require(tidyverse)
  
  data <- gs_ls() %>% 
    filter(sheet_title == filename) %>% 
    .$sheet_key %>%
    gs_key() %>%
    gs_gs() %>%
    gs_read(ws = sheetname)
  
  
  return(data)
}

###########################################################################################################
add_file_paths <- function(raw_metadata, data_root) {
  metadata <- raw_metadata %>%
    rowwise() %>%
    mutate(full_file = paste0(data_root,
                              parent_folder,
                              '/', 
                              experiment_folder,
                              '/masses.csv')) %>%
    ungroup()
}

###########################################################################################################
clean_metadata <- function(metadata) {

  # check if all the files described in the spreadsheet exist, return a warning if not
  metadata$data_exists <- file.exists(metadata$full_file)
  if (any(metadata$data_exists == FALSE)) {
    n_missing_files = sum(metadata$data_exists == FALSE)
    if (n_missing_files == 1) {
      print(paste0("Warning: no data found for 1 file:"))
    } else {
      print(paste0("Warning: no data found for ", as.character(n_missing_files), " files:"))
    }
  
    missing_expt_ids <- metadata$expt_id[metadata$data_exists == FALSE]
    for (ii in seq_along(missing_expt_ids)) {
      print(paste0(missing_expt_ids[ii]))
    }
  }
  
  # remove files that don't exist from the metadata dataframe
  metadata <- metadata %>%
    filter((exclude_from_analysis != TRUE) | is.na(metadata$exclude_from_analysis)) %>%
    filter(data_exists == TRUE)
  
  return(metadata)
}

###########################################################################################################
load_metadata <- function(sheetname) {
  
  require(assertthat)
  
  raw_metadata_mb1 <- suppressMessages(read_google_sheet(filename = 'mass-blaster-1',
                                                         sheetname = sheetname))
  metadata_mb1 <- add_file_paths(raw_metadata_mb1, 
                                 data_root = "Z:/maxstock/222 systems/Data - Mass blaster/")
  
  
  raw_metadata_mb2 <- suppressMessages(read_google_sheet(filename = 'mass-blaster-2',
                                                         sheetname = sheetname))
  metadata_mb2 <- add_file_paths(raw_metadata_mb2,
                                 data_root = 'Z:/maxstock/222 systems/Data - Blue system/')
  
  metadata <- bind_rows(metadata_mb1, metadata_mb2)
  
  if (!all(names(metadata_mb1) == names(metadata_mb2))) {
    stop("Column names don''t match for the two metadata spreadsheets.")
  }
  
  
  return(metadata)
}

###########################################################################################################
load_smr_data <- function(metadata) {
  data <- map(metadata$full_file, ~suppressMessages(read_csv(., 
                                                             col_names = c("mass", "transit_time", "sensor_number", "peak_height"), 
                                                             skip = 1)))
  temp_data <- data.frame()
  for (ii in seq_along(data)) {
    expt_data <- bind_cols(data.frame(expt_id = rep(x = metadata$expt_id[ii], 
                                                    times = nrow(data[[ii]]))),
                           data[[ii]])
    temp_data <- suppressWarnings(bind_rows(temp_data, expt_data))
  }
  data <- temp_data
  rm(temp_data)
  
  return(data)
}


###########################################################################################################
load_genomics <- function() {
  genomics <- suppressMessages(read_csv("raw/genomics.csv"))
  genomics$pt_mgmt <- as.factor(genomics$pt_mgmt)
  genomics$pdcl_mgmt = as.factor(genomics$pdcl_mgmt)
  genomics$mmr = as.factor(genomics$mmr)
  return(genomics)
}

###########################################################################################################
load_clinical_data <- function() {
  clinical <- suppressMessages(read_csv("raw/overall_survival.csv"))
  clinical$overall_survival <- as.numeric(clinical$overall_survival)
  return(clinical)
}

#############################################################
exclude_lines <- function(biomarkers, lines_to_exclude) {
  for (row_number in 1:nrow(lines_to_exclude)) {
    current_row <- raw_data$lines_to_exclude[row_number,]
    
    if (current_row$batch == "any") {
      biomarkers <- filter(biomarkers, !(line == current_row$line))
      print(paste0("Excluded ", current_row$line, ". Reason: ", current_row$reason))
    } else if (current_row$batch != "any" & current_row$line != "any") {
      biomarkers <- filter(biomarkers, !(line == current_row$line & batch == current_row$batch))
      print(paste0("Excluded ", current_row$line, " from Batch ", current_row$batch, ". Reason: ", current_row$reason))
    } else if (current_row$batch != "any" & current_row$line == "any") {
      biomarkers <- filter(biomarkers, !(batch == current_row$batch))
      print(paste0("Excluded all lines from Batch ", current_row$batch, ". Reason: ", current_row$reason))
    }
  }
  
  return(biomarkers)
}