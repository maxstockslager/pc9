# load R packages required for the analysis
library(tidyverse)
library(R.utils)
sourceDirectory("R/", modifiedOnly=TRUE)

CELL_COUNT_THRESHOLD <- 500 # minimum number of cells measured for us to process a sample. samples with less are not included

smr_metadata_raw <- load_metadata(sheetname = '20191002-PC9') # read from Google Sheets
smr_metadata <- clean_metadata(smr_metadata_raw)               # make sure the data exists
smr_data <- load_smr_data(smr_metadata)                        # load the SMR data from Rowley

control_lookup_table <- match_control_to_treatment(smr_metadata) # describes which experiment IDs correspond to paired treatment/controls
derived_quantities_by_sample <- compute_derived_quantities_by_sample(control_lookup_table, smr_data)

######
DOSE <- tribble(
  ~name, ~dose,
  "Gefitinib-Dose1", 10000/1000, 
  "Gefitinib-Dose2", 4000/1000,
  "Gefitinib-Dose3", 1600/1000,
  "Gefitinib-Dose4", 640/1000,
  "Gefitinib-Dose5", 256/1000,
  "Gefitinib-Dose6", 102/1000,
  "Gefitinib-Dose7", 41/1000,
  "Gefitinib-Dose8", 16/1000
)

#######
smr_data %>%
  left_join(smr_metadata, by = "expt_id") %>%
  ggplot(mapping = aes(x = treatment, y = mass, fill = line)) + 
  geom_violin() + 
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  theme_bw() + 
  facet_grid(line ~ .) + 
  theme(panel.grid.major.x = element_blank())

########
smr_metadata %>% 
  left_join(derived_quantities_by_sample) %>%
  #filter(treatment %in% DOSE$name) %>%
  left_join(DOSE, by = c("treatment" = "name")) %>%
  ggplot(mapping = aes(x = treatment, y = hellinger_dist, color = line)) +
  facet_grid(line ~ .) +
  geom_point() + 
  geom_errorbar(mapping = aes(ymin= hellinger_05, ymax = hellinger_95)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        panel.grid.major.x = element_blank()) + 
 # scale_x_log10() + 
  scale_y_log10()

#########
smr_data %>%
  left_join(smr_metadata, by = "expt_id") %>%
  group_by(treatment, line) %>%
  summarize(mean_mass = mean(mass),
            n_cells = n(),
            sd_mass = sd(mass),
            mean_mass_se = sd_mass/sqrt(n_cells)) %>%
  left_join(DOSE, by = c("treatment" = "name")) %>%
  ggplot(mapping = aes(x = treatment, y = mean_mass)) + 
  geom_point() + 
  geom_errorbar(mapping = aes(ymin = mean_mass - mean_mass_se, ymax = mean_mass + mean_mass_se)) +
  facet_grid(line ~ .) + 
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())


###############
counts <- smr_data %>%
  left_join(smr_metadata, by = "expt_id") %>%
  group_by(line, treatment) %>% 
  summarize(n_cells = n()) 

###############
counts %>%
  ggplot(mapping = aes(x = treatment, y = n_cells)) + 
  geom_col() + 
  facet_grid(line ~.) + 
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

####################################
# make plots (one example included here)


smr_metadata %>%
  left_join(derived_quantities_by_sample) %>%
  filter(status == "live",
         treatment != "DMSO",
         treatment != "Ponatinib") %>%
  left_join(DOSE, by = c("treatment" = "name")) %>%
  ggplot(mapping = aes(x = dose, y = hellinger_dist, color = line)) + 
  geom_point() + 
  geom_errorbar(mapping = aes(ymin = hellinger_dist - hellinger_sd, ymax = hellinger_dist + hellinger_sd)) +
  theme_bw() + 
  xlab("Treatment") + 
  ylab("SMR response score (Hellinger distance)") + 
  ggtitle("BaF3 WT/T315I imatinib dose response") + 
  coord_cartesian(ylim = c(0, 0.075)) + 
  scale_x_log10() + 
  facet_grid(. ~ fct_rev(line))

SAVE_PLOT <- F
if (SAVE_PLOT) {
  ggsave("./output/figs-raw/imatinib-dose-response.svg",
         width = 5, 
         height = 4)
}
