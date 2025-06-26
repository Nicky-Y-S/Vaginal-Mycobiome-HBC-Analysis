#Nicky Sun 

library(phyloseq)
library(dplyr)
library(stringr)
library(writexl)
library(decontam)
library(phyloseq)
library(vegan)
library(pheatmap)
library(tidyverse)
library(Matrix)
library(readxl)
library(scales)
library(performance)
library(ggforce)
library(lme4)
library(lmerTest)
library(patchwork)
library(tidyr)
library(ggrepel)          

getwd()
setwd("/Users/yiningsun/Desktop/Tetel Lab")

data <- readRDS("/Users/yiningsun/Desktop/Tetel Lab/Walther-Antonio_Project_022_ITS2.rds")
sample_label <- read_excel("/Users/yiningsun/Desktop/Tetel Lab/cleaned_samplesv2.xlsx")

#change column name of sample_label for merging later
#View(sample_label)
sample_label$SampleID <- sample_label$qr
#View(sample_label)

###################################################################################################################

otuT <- otu_table(data)
taxaT <- tax_table(data)

#View(otuT)
#View(taxaT)

###################################################################################################################
#DATA PROCESSING

#step 1: create new phyloseq object (new_phyloseq) with desired sample data, OTU, and taxonomy table

#convert OTU table to a plain data.frame and transpose the table (taxa x samples -> samples x taxa)
otu_table.data <- as.data.frame(otu_table(data))
otu_table.data <- as.data.frame(t(otu_table.data))
#View(otu_table.data) #should be transposed now

#process the otu_table.data object to:
#extract row names (assumed to be sample IDs) into a new column called SampleID
#create a logical column (is_blank) that indicates whether the sample ID contains the word "BLANK"
#clean up the SampleID column by removing everything after the first dot
#retain only the SampleID and is_blank columns in the final output, which is stored in metadata
metadata <- otu_table.data %>% 
  mutate(SampleID=rownames(otu_table.data),
         is_blank=as.logical(ifelse(str_detect(SampleID, "BLANK"), "TRUE", "FALSE")),
         SampleID= sub("\\..*", "", SampleID)) %>%
  select(SampleID, is_blank)
#View(metadata)

#create new OTU table object
otu_table_obj <- otu_table(otu_table.data, taxa_are_rows = FALSE)
#View(otu_table_obj) #original OTU table should be transposed now

#create new sample-data object
sample_data_obj <- sample_data(metadata)
#View(sample_data_obj)

#assemble the phyloseq object with the metadata, transposed OTU table, and original taxonomy table
new_phyloseq <- phyloseq(otu_table_obj, sample_data_obj, tax_table(data))

###################################################################################################################

#step 2: filter contaminant from new_phyloseq

#identify contaminants (non-biological samples) based on prevalence
contam_prev <- isContaminant(new_phyloseq, method = "prevalence", neg = "is_blank")
contaminants <- contam_prev$contaminant

###################################################################################################################

#step 3: filter the OTU and taxa table in the phyloseq object

new_phyloseq_no_contam <- prune_taxa(!contaminants, new_phyloseq)

#filters OTU and taxa table to only include samples classified as Fungi and has non‐blank Phylum
fungal_phyloseq_subset <- subset_taxa(new_phyloseq_no_contam, Kingdom == "Fungi" & !is.na(Phylum) & Phylum != "")

fungal_otu_tab <- otu_table(fungal_phyloseq_subset)
fungal_taxa_tab <- tax_table(fungal_phyloseq_subset)

#View(fungal_otu_table)
#View(fungal_taxa_table)
#View(sample_data(fungal_phyloseq_subset))

###################################################################################################################

#step 4: merge metadata with sample_label

#convert metadata to df
fungal_metadata_dataframe <- as.data.frame(as.matrix(sample_data(fungal_phyloseq_subset)))
#View(fungal_metadata_df)

#filter samples
fungal_metadata_dataframe_IDs <- fungal_metadata_dataframe %>% 
  left_join(sample_label, by="SampleID") %>% 
  filter(!is.na(biome_id))
#View(fungal_metadata_dataframe_IDs)

###################################################################################################################

#step 5: update OTU and taxa table and build filtered phyloseq object

#strips off everything from the first dot onward, leaving "D101", "F103", etc., so that they line up exactly with the SampleID values in metadata
rownames(fungal_otu_tab) <- sub("\\..*", "", rownames(fungal_otu_tab))
rownames(fungal_taxa_tab) <- sub("\\..*", "", rownames(fungal_taxa_tab))

#And set the row names of metadata as sampleID
fungal_metadata_dataframe_IDs_rownames <- fungal_metadata_dataframe_IDs$SampleID
rownames(fungal_metadata_dataframe_IDs) <- fungal_metadata_dataframe_IDs$SampleID
#View(fungal_metadata_dataframe_IDs)

otu_tab_filtered <- fungal_otu_tab[fungal_metadata_dataframe_IDs_rownames, , drop = FALSE]
taxa_tab_filtered <- fungal_taxa_tab[colnames(otu_tab_filtered), , drop = FALSE]

fungal_phyloseq_filtered <- phyloseq(otu_table(otu_tab_filtered), sample_data(fungal_metadata_dataframe_IDs), tax_table(taxa_tab_filtered))

###################################################################################################################

#step 6: get OTU table, taxa table, and sample data from filtered phyloseq object

fungal_otu_tab.data <- (otu_table(fungal_phyloseq_filtered))
fungal_taxa_tab.data <- tax_table(fungal_phyloseq_filtered)
fungal_sample_data_phyloseq <- as.data.frame(as.matrix(sample_data(fungal_phyloseq_filtered)))

#View(fungal_otu_tab.data)
#View(fungal_taxa_tab.data)
#View(fungal_sample_data_phyloseq)

###################################################################################################################

#step 7: filter for just vaginal

vaginal_sample_data_phyloseq <- fungal_sample_data_phyloseq %>% 
  filter(sampleType=="vaginal")
vaginal_fungal_otu <- fungal_otu_tab.data[rownames(vaginal_sample_data_phyloseq), , drop = FALSE]
vaginal_fungal_taxa <- fungal_taxa_tab.data[colnames(vaginal_fungal_otu), , drop = FALSE]

#vaginal phyloseq obj
vaginal_phyloseq <- phyloseq(otu_table(vaginal_fungal_otu, taxa_are_rows=FALSE), 
                             sample_data(vaginal_sample_data_phyloseq), tax_table(vaginal_fungal_taxa))

vaginal_tax_tab <- as.data.frame(tax_table(vaginal_phyloseq))
#View(vaginal_tax_tab) #bad species format, need cleaning

###################################################################################################################

#step 8: clean tax table 
vaginal_tax_tab <- vaginal_tax_tab %>% 
  mutate(
    Species_full = case_when(
      !is.na(Genus)   & Genus   != "" &
        !is.na(Species) & Species != "" &
        str_starts(Species, paste0(Genus, "_")) &        # begins with the genus
        !str_ends(Species, "_sp")                        # and does NOT end in _sp
      ~ Species,
            !is.na(Genus)   & Genus   != "" &
        !is.na(Species) & Species != "" &
        !str_ends(Species, "_sp")                        # exclude Genus_sp
      ~ paste0(Genus, "_", Species),
      TRUE ~ "Unknown"
    )
  )

###################################################################################################################

#vaginal dominant species
species_lookup <- vaginal_tax_tab$Species_full   # already just Genus_species or "Unknown"
vaginal_dominant_species <- apply(otu_table(vaginal_phyloseq), 1, function(x) {
  if (all(x == 0, na.rm = TRUE)) return("Unknown")
  lab <- species_lookup[ which.max(x) ]
  lab
})

sample_data(vaginal_phyloseq)$DominantSpecies <- vaginal_dominant_species

###################################################################################################################

#vaginal alpha diversity
vaginal_alpha_diversity <- estimate_richness(vaginal_phyloseq, measures = c("Shannon"))$Shannon
sample_data(vaginal_phyloseq)$Shannon <- vaginal_alpha_diversity

###################################################################################################################

#total reads per sample
taxa_are_rows(vaginal_phyloseq)
total_reads <- rowSums(otu_table(vaginal_phyloseq))
sample_data(vaginal_phyloseq)$total_reads <- total_reads
detection_limit <- 1 / total_reads

###################################################################################################################

#transform to rel. abundances
vaginal_phyloseq_rel <- transform_sample_counts(vaginal_phyloseq, function(x) x / sum(x))

#add C. albicans rel. abundance to sample data
ca_phy <- subset_taxa(vaginal_phyloseq_rel, Species == "Candida_albicans")
ca_abund <- rowSums( otu_table(ca_phy)[ , , drop = FALSE ] )
sample_data(vaginal_phyloseq_rel)$CA_abund <- ca_abund[sample_names(vaginal_phyloseq_rel)]

ca_abund_limited <- ifelse(ca_abund == 0, detection_limit[names(ca_abund)], ca_abund)
sample_data(vaginal_phyloseq_rel)$CA_abund_limited <- ca_abund_limited

###################################################################################################################

#convert sample_data to data frame
vaginal_rel_metadata_df <- as(sample_data(vaginal_phyloseq_rel), "data.frame")
vaginal_rel_metadata_df$biome_id <- as.integer(vaginal_rel_metadata_df$biome_id)
vaginal_rel_metadata_df$qr <- NULL
vaginal_rel_metadata_df$is_blank <- NULL
vaginal_rel_metadata_df$status <- NULL
#View(vaginal_rel_metadata_df) #1140 entries

###################################################################################################################

#DATA CLEANING

#medication cleaning
med_df <- read.csv("cleaned_Report 9-Volunteer Medical History.csv",
                             header = TRUE, stringsAsFactors = FALSE) %>%
  filter(birthControl != "Orilissa (Elagolix)") %>%
  select(biome_id, birthControl, taken_antibiotics) %>%
  mutate(biome_id = as.integer(biome_id),
         birthControl = factor(birthControl,
                               levels = c("None", "Local P", "Systemic P only", "Systemic Combined (E&P)"))) #fewer number of participants for medical history survey

vaginal_rel_df_matched <- vaginal_rel_metadata_df %>%
  left_join(
    med_df %>% select(biome_id, birthControl, taken_antibiotics),
    by = "biome_id"
  )

vaginal_rel_df_matched <- vaginal_rel_df_matched %>%
  filter(!is.na(birthControl))
#View(vaginal_rel_df_matched) #1084 entries

#DASS cleaning
dass_df <- read.csv("/Users/yiningsun/Desktop/Tetel Lab/DASS_0503_2024-final_df.csv")
dass_df$mood_date <- as.Date(dass_df$Timestamp)
dass_df$biome_id <- dass_df$study_id

dass_df <- dass_df %>%
  mutate(
    stress_severity = case_when(
      stress_score >= 0  & stress_score <= 14 ~ 0,
      stress_score >= 15 & stress_score <= 18 ~ 1,
      stress_score >= 19 & stress_score <= 25 ~ 2,
      stress_score >= 26 & stress_score <= 33 ~ 3,
      stress_score >= 34                      ~ 4,
      TRUE ~ NA_real_
    )
  )

vaginal_rel_df_matched$logDate <- as.Date(vaginal_rel_df_matched$logDate)

dass_clean <- dass_df %>%
  select(biome_id, mood_date, stress_score, stress_severity) %>%
  distinct()

closest_matches <- vaginal_rel_df_matched %>%
  inner_join(dass_clean, by = "biome_id") %>%
  mutate(date_diff = abs(as.numeric(difftime(logDate, mood_date, units = "days")))) %>%
  group_by(biome_id, logDate) %>%
  slice_min(order_by = date_diff, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    stress_score = ifelse(date_diff <= 7, stress_score, NA)  # NA if >7 days
    # do NOT change stress_severity – keep it even if score is NA
  )

vaginal_rel_df_matched <- vaginal_rel_df_matched %>%
  left_join(
    closest_matches %>%
      select(biome_id, logDate, mood_date, stress_score, stress_severity),
    by = c("biome_id", "logDate")
  )

#sport, athlete, field hockey cleaning
athlete_df <- dass_df %>%
  select(biome_id, sport, athlete) %>% 
  distinct(biome_id, .keep_all = TRUE)

vaginal_rel_df_matched <- vaginal_rel_df_matched %>% 
  left_join(athlete_df, by = "biome_id")

vaginal_rel_df_matched <- vaginal_rel_df_matched %>% 
  mutate(
    field_hockey = case_when(
      sport == "In-Season" ~ "Field Hockey",  # samples taken during FH season
      sport == "None"      ~ "None",          # non-athletes
      TRUE                 ~ "Other"          # Club, Off-Season, Nordic Ski
    ),
    field_hockey = factor(field_hockey, levels = c("None", "Other", "Field Hockey"))
  )

#View(vaginal_rel_df_matched) #1084 entries

#NA check
na_samples <- vaginal_rel_df_matched %>% 
  filter(is.na(sport))            
na_samples #participant 31 has NA info

###################################################################################################################
#DATA ANALYSIS -> C. ALBICANS RELATIVE ABUNDANCE AND LIFESTYLE FACTORS

#C. albicans rel. abundance by HBC sina plot (w/ detection limit)
hbc_cols <- c("None" = "#D73027", "Local P" = "#4575B4", "Systemic P only" = "#91CF60", "Systemic Combined (E&P)" = "#8073AC")

box_width <- 0.75
jitter_width <- box_width/2  

p1 <- ggplot(vaginal_rel_df_matched,
       aes(x = birthControl, y = CA_abund_limited)) +
  geom_sina(aes(color = birthControl),        
            size = 1.5, alpha = 0.7,
            maxwidth = box_width * 0.4,
            show.legend = FALSE) +                
  geom_violin(aes(fill = birthControl),
              color = "black", size = 0.3,
              alpha = 0.4, width = box_width) +
  geom_boxplot(width = 0.1, outlier.shape = NA,
               color = "black", size = 0.5) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(name = "Birth Control",
                    values = hbc_cols) +
  scale_color_manual(values = hbc_cols, guide = "none") +
  theme_minimal(base_size = 14) +
  labs(title = "Vaginal Mycobiome: C. albicans Relative Abundance by Birth Control Type (Detection Limit Applied)",
       x = "Birth Control", y = "Relative Abundance of C. albicans") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", size = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.ticks.length = unit(3, "pt"),
        legend.position = "right")

###################################################################################################################

#C. albicans rel. abundance over time by HBC (w/ detection limit)
all_days <- seq.Date(as.Date("2022-10-13"), as.Date("2022-12-16"), by = "day")
vaginal_rel_df_matched$study_day <- match(as.Date(vaginal_rel_df_matched$logDate), all_days) - 1

p_hbc <- ggplot(vaginal_rel_df_matched, aes(x = study_day, y = CA_abund_limited, color = birthControl)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(se = FALSE, method = "loess", size = 0.5) +
  scale_color_manual(values = hbc_cols) +
  labs(title = "Vaginal Mycobiome: C. albicans Relative Abundance Over Time by Birth Control Type (Detection Limit Applied)",
       x = "Study Day",
       y = "Relative Abundance of C. albicans",
       color = "Birth Control") +
  theme_minimal()

###################################################################################################################

#C. albicans rel. abundance overtime by athlete status (w/ detection limit)
athlete_cols <- c("0" = "#D36FA4", "1" = "#1A9850")

p_ath <- ggplot(
  vaginal_rel_df_matched,
  aes(study_day, CA_abund_limited,
      colour = factor(athlete, exclude = NULL))
) +
  geom_point(alpha = 0.40, size = 1) +
  geom_smooth(method = "loess", se = FALSE, size = 0.5) +
  scale_colour_manual(
    values  = athlete_cols,
    breaks  = c("0", "1"),
    labels  = c("Non-athlete", "Athlete"),
    name    = "Athlete Status",
    na.translate = TRUE,
    na.value     = na_col
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Vaginal Mycobiome: C. albicans Relative Abundance Over Time by Athlete Status (Detection Limit Applied)",
    x     = "Study Day",
    y     = "Relative Abundance of C. albicans"
  ) +
  theme_minimal()

print(p_ath)

#C. albicans rel. abundance and athlete mixed effect quadratic model
vaginal_rel_df_matched <- vaginal_rel_df_matched %>%
  mutate(day_c = scale(study_day, center = TRUE, scale = FALSE))

m_rel_athlete_quad <- lmer(CA_abund_limited ~ athlete + day_c + I(day_c^2) + (1 | biome_id), data = vaginal_rel_df_matched)
summary(m_rel_athlete_quad)
r2(m_rel_athlete_quad)

###################################################################################################################

#C. albicans rel. abundance overtime by sports category (w/ detection limit)
vaginal_rel_df_matched <- vaginal_rel_df_matched %>% 
  mutate(sport = factor(sport, levels = c(NA, "None", "Club", "In-Season", "Off-Season", "Nordic Ski")))

sports_cols <- c("None" = "#B0B0B0", "Club" = "#795548", "Off-Season" = "#FB8072", "In-Season" = "#8DD3C7", "Nordic Ski"  = "#80B1D3")

ggplot(vaginal_rel_df_matched,
       aes(x = study_day, y = CA_abund_limited, colour = sport)) +
  geom_point(size = 1, alpha = 0.55) +
  geom_smooth(method = "loess", se = FALSE, size = 0.5) +
  scale_color_manual(values = sports_cols,
                     name = "Sport Category") +
  scale_y_continuous(limits = c(0, 1)) +         
  labs(
    title = "Vaginal Mycobiome: C. albicans Relative Abundance Over Time by Sport Category (Detection Limit Applied)",
    x = "Study Day",
    y = "Relative Abundance of C. albicans"        
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90", size = 0.5),
    legend.position    = "right"
  )

###################################################################################################################

#C. albicans rel. abundance overtime by field hockey category (w/ detection limit)
fh_cols <- c("None" = "#B0B0B0", "Other" = "#FDD95E", "Field Hockey" = "#8DD3C7")

p_fh <- ggplot(
  vaginal_rel_df_matched,
  aes(study_day, CA_abund_limited, colour = field_hockey)
) +
  geom_point(alpha = 0.40, size = 1) +
  geom_smooth(method = "loess", se = FALSE, size = 0.5) +
  scale_colour_manual(
    values = fh_cols,
    name   = "Sport Group"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Vaginal Mycobiome: C. albicans Relative Abundance Over Time\nNone vs. Other Sports vs. Field Hockey (Detection Limit Applied)",
    x     = "Study Day",
    y     = "Relative Abundance of C. albicans"
  ) +
  theme_minimal()

print(p_fh)

#C. albicans rel. abundance and field hockey mixed effect quadratic model
m_rel_FH_quad <- lmer(CA_abund_limited ~ field_hockey + day_c + I(day_c^2) + (1 | biome_id), data = vaginal_rel_df_matched)
summary(m_rel_FH_quad)
r2(m_rel_FH_quad)

###################################################################################################################

#comparison
combined_panel <- (p_hbc / p_ath / p_fh) +
  plot_layout(guides = "keep")  # each keeps its own legend
combined_panel

###################################################################################################################

#C. albicans rel. abundance and stress correlation
ggplot(vaginal_rel_df_matched, aes(stress_score, CA_abund_limited)) +
  geom_point(alpha = .35, size = 1) +
  geom_smooth(method = "loess", se = TRUE, colour = "black") +
  labs(x = "Stress score", y = "C. albicans relative abundance") +
  theme_minimal(base_size = 14)

cor.test(vaginal_rel_df_matched$stress_score, vaginal_rel_df_matched$CA_abund_limited,  method = "spearman")
cor.test(vaginal_rel_df_matched$stress_score, vaginal_rel_df_matched$CA_abund_limited,  method = "pearson")

###################################################################################################################

#individual C. albicans rel. abundance over time
ggplot(vaginal_rel_df_matched, aes(x = study_day, y = CA_abund_limited)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_line(aes(group = biome_id), alpha = 0.4) +
  facet_wrap(~ biome_id) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Vaginal Mycobiome: C. albicans Relative Abundance Over Time by Participant",
    x = "Study Day",
    y = "Relative Abundance of C. albicans"
  ) +
  theme(
    strip.text = element_text(size = 8),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.ticks.length = unit(3, "pt")
  )

ggplot(vaginal_rel_df_matched, aes(x = study_day, y = CA_abund_limited)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "loess", se = FALSE, size = 0.5, color = "black") +
  facet_wrap(~ biome_id) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Vaginal Mycobiome: Smoothed C. albicans Relative Abundance Over Time by Participant",
    x = "Study Day",
    y = "Relative Abundance of C. albicans"
  ) +
  theme(
    strip.text = element_text(size = 8),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.ticks.length = unit(3, "pt")
  )

###################################################################################################################

#troubleshoot

#uneven sampling overtime, missing data for certain days 
vaginal_rel_df_matched %>%
  count(study_day) %>%
  ggplot(aes(x = study_day, y = n)) +
  geom_col() +
  labs(title = "Sample Count per Study Day (Unfiltered Data)", x = "Study Day", y = "Number of Samples")

table(vaginal_rel_df_matched$study_day)

#test C. albicans between dropouts and completers
vaginal_rel_df_matched <- vaginal_rel_df_matched %>%
  mutate(study_week = floor(study_day / 7))

weeks_per_person <- vaginal_rel_df_matched %>%
  distinct(biome_id, study_week) %>%          
  count(biome_id, name = "weeks_with_sample")

total_weeks <- 10

completer_flag <- weeks_per_person %>%
  mutate(completer = weeks_with_sample >= 8) #80% coverage

median_CA <- vaginal_rel_df_matched %>%      
  group_by(biome_id) %>%
  summarise(baseline_CA = median(CA_abund, na.rm = TRUE),
            .groups = "drop") %>%
  left_join(completer_flag, by = "biome_id")

attr_p <- wilcox.test(baseline_CA ~ completer, data = median_CA)$p.value
attr_p

#re-graph plot with completer cohort
completer_ids <- completer_flag %>% filter(completer) %>% pull(biome_id)

comp_data <- vaginal_rel_df_matched %>%
  filter(biome_id %in% completer_ids)

ggplot(comp_data, aes(x = study_day, y = CA_abund, color = birthControl)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(se = FALSE, method = "loess", size = 0.5) +
  scale_color_manual(values = hbc_cols) +
  labs(
    title = paste0("1 Sample Each Week For At Least 8 Weeks"),
    x = "Study Day", y = "Relative Abundance of C. albicans",
    color = "Birth Control"
  ) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "right")

###################################################################################################################
#DATA ANALYSIS -> SHANNON DIVERSITY

#Shannon and HBC violin plot
ggplot(vaginal_rel_df_matched, aes(x = birthControl, y = Shannon)) +
  geom_sina(aes(color = birthControl),
            size = 1.5, alpha = 0.7,
            maxwidth = box_width * 0.4,
            show.legend = FALSE) +
  geom_violin(aes(fill = birthControl),
              color = "black", size = 0.3,
              alpha = 0.4, width = box_width) +
  geom_boxplot(width = 0.1,
               outlier.shape = NA,
               color = "black", size = 0.5) +
  scale_y_continuous(name = "Shannon Diversity") +
  scale_fill_manual(name = "Birth Control", values = hbc_cols) +
  scale_color_manual(values = hbc_cols, guide = "none") +
  theme_minimal(base_size = 14) +
  labs(title = "Vaginal Mycobiome: Shannon Diversity by Birth Control Type",
       x = "Birth Control") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", size = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.ticks.length = unit(3, "pt"),
        legend.position = "right")

###################################################################################################################

#Shannon over time by HBC
ggplot(vaginal_rel_df_matched, aes(x = study_day, y = Shannon, color = birthControl)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(se = FALSE, method = "loess", size = 0.5) +
  scale_color_manual(values = hbc_cols) +
  labs(title = "Vaginal Mycobiome: Shannon Diversity Over Time by Birth Control Type",
       x = "Study Day",
       y = "Shannon Diversity",
       color = "Birth Control") +
  theme_minimal() #spline reveals a linear mixed effect model seems appropriate

#Shannon and HBC mixed effect model 
m_shannon_HBC_time <- lmer(Shannon ~ birthControl + study_day + (1 | biome_id), data = vaginal_rel_df_matched)
summary(m_shannon_HBC_time)
r2(m_shannon_HBC_time)

###################################################################################################################

#average Shannon by HBC boxplot
vaginal_avg_shannon_df <- vaginal_rel_df_matched %>%
  group_by(biome_id, birthControl) %>%
  summarise(
    n_samples    = n(),
    Mean_Shannon = mean(Shannon, na.rm = TRUE),
    SD_Shannon   = sd(Shannon,   na.rm = TRUE),
    Min_Shannon  = min(Shannon,  na.rm = TRUE),
    Max_Shannon  = max(Shannon,  na.rm = TRUE)
  ) %>%
  ungroup()

ggplot(vaginal_avg_shannon_df, aes(x = birthControl, y = Mean_Shannon)) +
  geom_sina(aes(color = factor(biome_id)), size = 1.5, alpha = 0.7, maxwidth = box_width * 0.4, show.legend = FALSE) +
  geom_violin(aes(fill = birthControl), color = "black", size = 0.3, alpha = 0.4, width = box_width) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", size = 0.5) +
  scale_fill_manual(name = "Birth Control", values = hbc_cols) +
  scale_color_viridis_d(option = "turbo", guide = FALSE) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Vaginal Mycobiome: Average Shannon Diversity by Birth Control Type",
    x = "Birth Control",
    y = "Average Shannon Diversity"
  ) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", size = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.ticks.length = unit(3, "pt"),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

#avg Shannon diversity and HBC linear model
m_avg_shannon_hbc <- lm(Mean_Shannon ~ birthControl, data = vaginal_avg_shannon_df)
summary(m_avg_shannon_hbc)

###################################################################################################################
#DATA ANALYSIS -> SHANNON, DOMINANT SPECIES, AND HBC

#filter top 3 most frequent dominant species (at sample-level)
top3 <- vaginal_rel_df_matched %>%
  count(DominantSpecies, sort = TRUE) %>%
  filter(DominantSpecies != "Unknown") %>%
  slice_head(n = 3) %>%
  pull(DominantSpecies)

vaginal_rel_df_matched <- vaginal_rel_df_matched %>%
  mutate(
    DomSpec_grp = case_when(
      DominantSpecies %in% top3 ~ DominantSpecies,
      DominantSpecies == "Unknown" ~ "Unknown",
      TRUE ~ "Other"
    ),
    DomSpec_grp = factor(
      DomSpec_grp,
      levels = c(top3, "Other", "Unknown")
    )
  )

species_cols <- c(setNames(c("#D95F02", "#1B9E77", "#E7298A"), top3),  "Other" = "#7570B3", "Unknown" = "#666666")

###################################################################################################################

#Shannon by dominant species boxplot
ggplot(vaginal_rel_df_matched,              
       aes(x = DomSpec_grp, y = Shannon, fill = DomSpec_grp)) +
  geom_boxplot(outlier.shape = NA, width = 0.55,
               colour = "black", alpha = 0.5) +
  geom_jitter(width = 0.15, shape = 21,
              size = 1.2, alpha = 0.5, colour = "black") +
  scale_fill_manual(values = species_cols,      
                    name   = "Dominant Species") +
  labs(
    title = "Vaginal Mycobiome: Shannon Diversity by Dominant Species",
    x     = "Dominant Species",
    y     = "Shannon Diversity"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x        = element_text(angle = 0, hjust = 0.5, size = 9),
    panel.grid.major.x = element_blank(),
    legend.title       = element_text(size = 10),
    legend.text        = element_text(size = 9),
    legend.key.size    = unit(0.4, "cm")
  )

###################################################################################################################

#Shannon over time by dominant species
ggplot(vaginal_rel_df_matched,
       aes(x = study_day, y = Shannon, colour = DomSpec_grp)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "loess", se = FALSE, size = 0.5) +
  scale_colour_manual(values = species_cols, name = "Dominant Species") +
  labs(title = "Vaginal Mycobiome: Shannon Diversity Over Time",
       x = "Study Day", y = "Shannon Diversity") +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.x = element_blank())

#Shannon diversity and dominant species mixed effect model
m_shannon_dom_time <- lmer(Shannon ~ DomSpec_grp + study_day + (1 | biome_id), data = vaginal_rel_df_matched)
summary(m_shannon_dom_time)
r2(m_shannon_dom_time)

###################################################################################################################

#dominant species and HBC bar chart
sample_counts <- vaginal_rel_df_matched %>%
  count(birthControl, DomSpec_grp, name = "n_samples") %>%
  group_by(birthControl) %>%
  mutate(prop = n_samples / sum(n_samples)) %>%
  ungroup()

hbc_totals <- vaginal_rel_df_matched %>%
  count(birthControl, name = "n_total")

ggplot(sample_counts,
       aes(x = birthControl, y = prop, fill = DomSpec_grp)) +
  geom_col(width = 0.8, colour = "black") +
  geom_text(
    data = hbc_totals,
    aes(x = birthControl, y = 1.02, label = paste0("n = ", n_total)),
    vjust = 0, size = 3.5, inherit.aes = FALSE
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1.08), expand = c(0, 0)) +
  scale_fill_manual(values = species_cols, name = "Dominant Species") +
  labs(title = "Vaginal Mycobiome: Sample-Level Dominant Species by Birth Control",
       x = "Birth Control", y = "Samples (%)") +
  theme_minimal(base_size = 13) +
  coord_cartesian(clip = "off")

###################################################################################################################

#Shannon by dominant species and HBC boxplot
ggplot(vaginal_rel_df_matched,
       aes(x = DomSpec_grp, y = Shannon, fill = birthControl)) +
  geom_boxplot(outlier.shape = NA, width = 0.55,
               colour = "black", alpha = 0.50) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15,
                                              dodge.width  = 0.55),
              shape = 21, size = 1.2, alpha = 0.50, colour = "black") +
  scale_fill_manual(values = hbc_cols, name = "Birth Control") +
  labs(title = "Vaginal Mycobiome: Shannon Diversity by Dominant Species and Birth Control",
       x = "Dominant Species", y = "Shannon Diversity") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9),
        panel.grid.major.x = element_blank())

#Shannon diversity, dominant species, and HBC mixed effect model
m_shannon_dom_hbc_time <- lmer(Shannon ~ DomSpec_grp + birthControl + study_day + (1 | biome_id), data = vaginal_rel_df_matched)
summary(m_shannon_dom_hbc_time)
r2(m_shannon_dom_hbc_time)

#Shannon diversity, dominant species and HBC interaction model
m_shannon_dom_hbc_int_time <- lmer(Shannon ~ DomSpec_grp * birthControl + study_day + (1 | biome_id), data = vaginal_rel_df_matched)
summary(m_shannon_dom_hbc_int_time)
r2(m_shannon_dom_hbc_int_time)

###################################################################################################################



