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
library(MuMIn)

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
        str_starts(Species, paste0(Genus, "_")) &       
        !str_ends(Species, "_sp")                        
      ~ Species,
            !is.na(Genus)   & Genus   != "" &
        !is.na(Species) & Species != "" &
        !str_ends(Species, "_sp")           
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

#DATA CLEANING & MERGING

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

#DASS cleaning
dass_df <- read.csv("/Users/yiningsun/Desktop/Tetel Lab/DASS_0503_2024-final_df.csv") %>%
  mutate(
    mood_date = as.Date(Timestamp),
    biome_id = study_id,
    stress_severity = case_when(
      stress_score >= 0  & stress_score <= 14 ~ 0,
      stress_score >= 15 & stress_score <= 18 ~ 1,
      stress_score >= 19 & stress_score <= 25 ~ 2,
      stress_score >= 26 & stress_score <= 33 ~ 3,
      stress_score >= 34                      ~ 4,
      TRUE ~ NA_real_
    )
  )

dass_clean <- dass_df %>%
  select(biome_id, mood_date, stress_score, stress_severity) %>%
  distinct()

vaginal_rel_df_matched$logDate <- as.Date(vaginal_rel_df_matched$logDate)

closest_matches <- vaginal_rel_df_matched %>%
  inner_join(dass_clean, by = "biome_id") %>%
  mutate(date_diff = abs(as.numeric(difftime(logDate, mood_date, units = "days")))) %>%
  group_by(biome_id, logDate) %>%
  slice_min(date_diff, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    stress_score = ifelse(date_diff <= 7, stress_score, NA_real_),
    stress_severity = ifelse(is.na(stress_score), NA_real_, stress_severity),
    stress_severity_char = case_when(
      is.na(stress_severity) ~ NA_character_,
      stress_severity == 0 ~ "Normal",
      stress_severity == 1 ~ "Mild",
      stress_severity == 2 ~ "Moderate",
      stress_severity == 3 ~ "Severe",
      stress_severity == 4 ~ "Extremely Severe"
    ),
    stress_severity_char = factor(
      stress_severity_char,
      levels = c("Normal", "Mild", "Moderate", "Severe", "Extremely Severe")
    )
  )

vaginal_rel_df_matched <- vaginal_rel_df_matched %>%
  left_join(
    closest_matches %>%
      select(biome_id, logDate, stress_score, stress_severity, stress_severity_char),
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
      is.na(sport)           ~ NA_character_,     # Preserve NA
      sport == "In-Season"   ~ "Field Hockey",
      sport == "None"        ~ "None",
      TRUE                   ~ "Other"            # Club, Off-Season, Nordic Ski
    ),
    field_hockey = factor(field_hockey, levels = c("None", "Other", "Field Hockey"))
  )

#climate cleaning
weather_df <- read.csv("/Users/yiningsun/Desktop/Tetel Lab/LCD_USW00094746_2022.csv") %>%
  select(Date = DATE, 
         DailyAverageDryBulbTemperature,
         DailyAverageWetBulbTemperature,
         DailyAverageDewPointTemperature,
         DailyAverageRelativeHumidity,
         DailyPrecipitation,
         ) %>%
  mutate(Date = as.Date(Date)) %>%
  filter(Date >= as.Date("2022-10-13") & Date <= as.Date("2022-12-15")) %>%
  distinct(Date, .keep_all = TRUE)

weather_df <- weather_df %>%
  rename(logDate = Date)

vaginal_rel_df_matched <- vaginal_rel_df_matched %>%
  left_join(weather_df, by = "logDate")

#menses cleaning 
menses_df <- read.csv("/Users/yiningsun/Desktop/Tetel Lab/imputed_menstruation_data_3_11.csv")
menses_df <- menses_df %>% 
  rename_with(~gsub("X2022.", "2022.", .), starts_with("X2022.")) %>% 
  rename_with(~gsub("\\.", "-", .))

menses_df_long <- menses_df %>% 
  pivot_longer(cols=starts_with("2022-"), names_to="logDate", values_to="menses_status")
menses_df_long$logDate <- as.Date(menses_df_long$logDate)

vaginal_rel_df_matched <- vaginal_rel_df_matched %>% 
  left_join(menses_df_long, by=c("biome_id", "logDate"))

vaginal_rel_df_matched <- vaginal_rel_df_matched %>% 
  mutate(menses_day = ifelse(menses_status %in% c(1,2,3,7,9,78), "menses", 
                             ifelse(menses_status %in% c(4,5,6,10), "not_menses", NA)))

#study day imputing
all_days <- seq.Date(as.Date("2022-10-13"), as.Date("2022-12-16"), by = "day")
vaginal_rel_df_matched$study_day <- match(as.Date(vaginal_rel_df_matched$logDate), all_days) - 1

#quadratic day centering
vaginal_rel_df_matched <- vaginal_rel_df_matched %>%
  mutate(day_c = scale(study_day, center = TRUE, scale = FALSE))

###################################################################################################################
#DATA ANALYSIS -> C. ALBICANS RELATIVE ABUNDANCE AND LIFESTYLE FACTORS

#C. albicans rel. abundance by HBC sina plot (w/ detection limit)
candida_hbc_filtered <- vaginal_rel_df_matched %>%
  filter(!is.na(birthControl))

hbc_cols <- c("None" = "#D73027", "Local P" = "#4575B4", "Systemic P only" = "#91CF60", "Systemic Combined (E&P)" = "#8073AC")

box_width <- 0.75
jitter_width <- box_width/2  

b_hbc <- ggplot(candida_hbc_filtered,
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

print(b_hbc)

###################################################################################################################

#C. albicans rel. abundance over time by HBC (w/ detection limit)
x_breaks <- seq(0, max(vaginal_rel_df_matched$study_day, na.rm = TRUE), by = 5)

p_hbc <- ggplot(candida_hbc_filtered, aes(study_day, CA_abund_limited, color = birthControl)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(se = FALSE, method = "loess", size = 0.5) +
  scale_color_manual(values = hbc_cols, name = "Birth Control") +
  scale_y_continuous(name = "Relative Abundance of C. albicans", limits = c(0, 1)) +
  scale_x_continuous(name = "Study Day", breaks = x_breaks) +
  labs(
    title = "Vaginal Mycobiome: C. albicans Relative Abundance Over Time by Birth Control Type (Detection Limit Applied)"
  ) +
  theme_minimal()

print(p_hbc)

###################################################################################################################

#C. albicans rel. abundance overtime by athlete status (w/ detection limit)
candida_ath_filtered <- vaginal_rel_df_matched %>%
  filter(!is.na(athlete))

athlete_cols <- c("0" = "#D36FA4", "1" = "#1A9850")

p_ath <- ggplot(candida_ath_filtered, aes(study_day, CA_abund_limited, color = factor(athlete))) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "loess", se = FALSE, size = 0.5) +
  scale_color_manual(values = athlete_cols, breaks = c("0", "1"),
                     labels = c("Non-athlete", "Athlete"),
                     name = "Athlete Status") +
  scale_y_continuous(name = "Relative Abundance of C. albicans", limits = c(0, 1)) +
  scale_x_continuous(name = "Study Day", breaks = x_breaks) +
  labs(
    title = "Vaginal Mycobiome: C. albicans Relative Abundance Over Time by Athlete Status (Detection Limit Applied)"
  ) +
  theme_minimal()

print(p_ath)

#C. albicans rel. abundance and athlete mixed effect quadratic model
m_rel_athlete_quad <- lmer(CA_abund_limited ~ athlete + day_c + I(day_c^2) + (1 | biome_id), data = vaginal_rel_df_matched)
summary(m_rel_athlete_quad)
r2(m_rel_athlete_quad)

###################################################################################################################

#C. albicans rel. abundance overtime by sports category (w/ detection limit)
candida_sport_filtered <- vaginal_rel_df_matched %>%
  filter(!is.na(sport))

sports_cols <- c("None" = "#B0B0B0", "Club" = "#795548", "Off-Season" = "#FB8072", "In-Season" = "#8DD3C7", "Nordic Ski"  = "#80B1D3")

p_sport <- ggplot(candida_sport_filtered, aes(study_day, CA_abund_limited, color = sport)) +
  geom_point(alpha = 0.55, size = 1) +
  geom_smooth(method = "loess", se = FALSE, size = 0.5) +
  scale_color_manual(values = sports_cols, name = "Sport Category") +
  scale_y_continuous(name = "Relative Abundance of C. albicans", limits = c(0, 1)) +
  scale_x_continuous(name = "Study Day", breaks = x_breaks) +
  labs(
    title = "Vaginal Mycobiome: C. albicans Relative Abundance Over Time by Sport Category (Detection Limit Applied)"
  ) +
  theme_minimal()

print(p_sport)

###################################################################################################################

#C. albicans rel. abundance overtime by field hockey category (w/ detection limit)
candida_fh_filtered <- vaginal_rel_df_matched %>%
  filter(!is.na(field_hockey))

fh_cols <- c("None" = "#B0B0B0", "Other" = "#FDD95E", "Field Hockey" = "#8DD3C7")

p_fh <- ggplot(candida_fh_filtered, aes(study_day, CA_abund_limited, color = field_hockey)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "loess", se = FALSE, size = 0.5) +
  scale_color_manual(values = fh_cols, name = "Sport Group") +
  scale_y_continuous(name = "Relative Abundance of C. albicans", limits = c(0, 1)) +
  scale_x_continuous(name = "Study Day", breaks = x_breaks) +
  labs(
    title = "Vaginal Mycobiome: C. albicans Relative Abundance Over Time\nNone vs. Other Sports vs. Field Hockey (Detection Limit Applied)"
  ) +
  theme_minimal()

print(p_fh)

#C. albicans rel. abundance and field hockey mixed effect quadratic model
m_rel_fh_quad <- lmer(CA_abund_limited ~ field_hockey + day_c + I(day_c^2) + (1 | biome_id), data = vaginal_rel_df_matched)
summary(m_rel_fh_quad)
r2(m_rel_fh_quad)

###################################################################################################################

#field hockey players analysis
games_df <- data.frame(
  game_date = as.Date(c("2022-10-15", "2022-10-18", "2022-10-22", "2022-10-29", "2022-11-01")), 
  W = c(0, 0, 1, 0, 0),
  home = c(1, 0, 0, 1, 0),
  score_diff = c(-4, -1, +1, -1, -2)
)

games_df <- games_df %>%
  arrange(game_date) %>%
  mutate(
    game_number = row_number(),
    result = ifelse(W == 1, "W", "L")
    )

fh_df <- vaginal_rel_df_matched %>%
  filter(field_hockey == "Field Hockey") %>%
  mutate(sample_date = as.Date(logDate)) %>%
  left_join(games_df, by = c("sample_date" = "game_date")) %>%
  mutate(
    is_game_day = !is.na(W),
    label_text = ifelse(is_game_day, paste0("Day ", study_day), NA),
    biome_id = factor(biome_id)
  )

View(fh_df)

ggplot(fh_df, aes(x = study_day, y = CA_abund_limited, group = biome_id, color = biome_id)) +
  geom_smooth(se = FALSE, method = "loess", size = 0.5) +
  geom_point(size = 1.5, alpha = 0.8) +
  geom_text_repel(
    data = fh_df %>% filter(is_game_day),
    aes(label = label_text),
    size = 3,
    show.legend = FALSE,
    max.overlaps = 50,
    segment.color = "grey60",
    box.padding = 0.4,
    point.padding = 0.3
  ) +
  scale_y_continuous(name = "C. albicans Relative Abundance", limits = c(0, 1)) +
  scale_x_continuous(name = "Study Day", breaks = x_breaks) +
  labs(
    title = "Vaginal Mycobiome: C. albicans Relative Abundance in Field Hockey Players Over Time (Detection Limit Applied)",
    subtitle = "Game Days Labeled by Study Day",
    color = "Participant ID"
  ) +
  theme_minimal()

###################################################################################################################

#C. albicans rel. abundance over time by menstruation status (w/ detection limit)
candida_menses_filtered <- vaginal_rel_df_matched %>%
  filter(!is.na(menses_day))

menses_cols <- c("menses" = "orange", "not_menses" = "brown")

p_menses <- ggplot(candida_menses_filtered, aes(study_day, CA_abund_limited, color = menses_day)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = FALSE, size = 0.5) +
  scale_color_manual(values = menses_cols,
                     labels = c("menses" = "Menstruating", "not_menses" = "Not Menstruating"),
                     name = "Menstruation Status") +
  scale_y_continuous(name = "Relative Abundance of C. albicans", limits = c(0, 1)) +
  scale_x_continuous(name = "Study Day", breaks = x_breaks) +
  labs(
    title = "Vaginal Mycobiome: C. albicans Relative Abundance Over Time by Menstruation Status (Detection Limit Applied)"
  ) +
  theme_minimal()

print(p_menses)

#C. albicans rel. abundance and menstruation mixed effect quadratic model
m_rel_men_quad <- lmer(CA_abund_limited ~ menses_day + day_c + I(day_c^2) + (1 | biome_id), data = vaginal_rel_df_matched)
summary(m_rel_men_quad)
r2(m_rel_men_quad)

###################################################################################################################

#C. albicans rel. abundance over time by stress severity (w/ detection limit)
candida_stress_filtered <- vaginal_rel_df_matched %>%
  filter(!is.na(stress_severity_char))

stress_severity_cols <- c("Normal" = "#F3E5F5", "Mild" = "#CE93D8", "Moderate" = "#AB47BC", "Severe" = "#6A1B9A", "Extremely Severe" = "#4A0072")

p_stress <- ggplot(candida_stress_filtered, aes(study_day, CA_abund_limited, color = stress_severity_char)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = FALSE, size = 0.5) +
  scale_color_manual(values = stress_severity_cols, 
                     name = "Stress Severity") +
  scale_y_continuous(name = "Relative Abundance", limits = c(0, 1)) +
  scale_x_continuous(name = "Study Day", breaks = x_breaks) +
  labs(
    title = "Vaginal Mycobiome: C. albicans Relative Abundance Over Time by Stress Severity (Detection Limit Applied)"
  ) +
  theme_minimal()

print(p_stress)

###################################################################################################################

#comparison
p_hbc_clean <- p_hbc + scale_y_continuous(name = "Relative Abundance") + labs(title = NULL) + ggtitle("Hormonal Birth Control")
p_ath_clean <- p_ath + scale_y_continuous(name = "Relative Abundance") + labs(title = NULL) + ggtitle("Athlete Status")
p_fh_clean <- p_fh + scale_y_continuous(name = "Relative Abundance") + labs(title = NULL) + ggtitle("Type of Sports")
p_menses_clean <- p_menses  + scale_y_continuous(name = "Relative Abundance") + labs(title = NULL) + ggtitle("Menstruation Status")

combined_panel <- (p_hbc_clean / p_ath_clean / p_fh_clean / p_menses_clean) +
  plot_layout(guides = "keep") +
  plot_annotation(
    title = "Vaginal Mycobiome: C. albicans Relative Abundance Over Time by Different Lifestyle Factors (Detection Limit Applied)",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  )

combined_panel

###################################################################################################################

#C. albicans rel. abundance and climate conditions (dry bulb, wet bulb, dew point, and relative humidity)

cor.test(vaginal_rel_df_matched$CA_abund_limited, vaginal_rel_df_matched$DailyAverageDryBulbTemperature, method = "spearman")
cor.test(vaginal_rel_df_matched$CA_abund_limited, vaginal_rel_df_matched$DailyAverageWetBulbTemperature, method = "spearman")
cor.test(vaginal_rel_df_matched$CA_abund_limited, vaginal_rel_df_matched$DailyAverageDewPointTemperature, method = "spearman")
cor.test(vaginal_rel_df_matched$CA_abund_limited, vaginal_rel_df_matched$DailyAverageRelativeHumidity, method = "spearman")

common_theme <- theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank())

smooth_layer <- geom_smooth(method = "loess", color = "blue", se = FALSE, size = 0.5)
point_layer <- geom_point(alpha = 0.6, size = 1)

c1 <- ggplot(vaginal_rel_df_matched, aes(x = DailyAverageDryBulbTemperature, y = CA_abund_limited)) +
  point_layer +
  smooth_layer +
  labs(
    title = "Dry Bulb Temperature   (Spearman = -0.144, p < 0.001)",
    x = "Daily Average Temperature (°C)",
    y = "C. albicans Relative Abundance"
  ) + common_theme

c2 <- ggplot(vaginal_rel_df_matched, aes(x = DailyAverageWetBulbTemperature, y = CA_abund_limited)) +
  point_layer +
  smooth_layer +
  labs(
    title = "Wet Bulb Temperature   (Spearman = -0.174, p < 0.001)",
    x = "Daily Average Temperature (°C)",
    y = "C. albicans Relative Abundance"
  ) + common_theme

c3 <- ggplot(vaginal_rel_df_matched, aes(x = DailyAverageDewPointTemperature, y = CA_abund_limited)) +
  point_layer +
  smooth_layer +
  labs(
    title = "Dew Point Temperature   (Spearman = -0.173, p < 0.001)",
    x = "Daily Average Temperature (°C)",
    y = "C. albicans Relative Abundance"
  ) + common_theme

c4 <- ggplot(vaginal_rel_df_matched, aes(x = DailyAverageRelativeHumidity, y = CA_abund_limited)) +
  point_layer +
  smooth_layer +
  labs(
    title = "Relative Humidity   (Spearman = -0.145, p = 0.001)",
    x = "Daily Average Relative Humidity",
    y = "C. albicans Relative Abundance"
  ) + common_theme

(c1 | c2) / (c3 | c4)

p2 <- ggplot(vaginal_rel_df_matched, aes(x = logDate, y = DailyAverageRelativeHumidity)) +
  geom_line(color = "grey") +
  labs(y = "Relative Humidity") +
  theme_minimal()

p3 <- ggplot(vaginal_rel_df_matched, aes(x = logDate, y = DailyAverageRelativeHumidity)) +
  geom_smooth(se = FALSE, color = "firebrick") +
  labs(y = "Relative Humidity") +
  theme_minimal()

p4 <- ggplot(vaginal_rel_df_matched, aes(x = logDate, y = CA_abund_limited)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE, color = "steelblue") +
  labs(y = "Relative abundance of C. albicans") +
  theme_minimal()

p2 / p3 / p4

#variable selection
m_cli_var <- lmer(CA_abund_limited ~ DailyAverageRelativeHumidity + DailyAverageDryBulbTemperature + DailyAverageWetBulbTemperature + DailyAverageDewPointTemperature + day_c + I(day_c^2) + (1 | biome_id), data = vaginal_rel_df_matched)
options(na.action = "na.fail")  
best_cli_var <- dredge(m_cli_var)
get.models(best_cli_var, 1)[[1]] #dredge() says dry bulb and wet bulb are the best predictors

###################################################################################################################

#overall model (selected by dredge())
complete_data <- vaginal_rel_df_matched %>%
  select(CA_abund_limited, DailyAverageRelativeHumidity, DailyAverageDryBulbTemperature,
         DailyAverageWetBulbTemperature, DailyAverageDewPointTemperature, birthControl,
         field_hockey, stress_score, menses_day, day_c, biome_id) %>%
  drop_na()

m_all <- lmer(CA_abund_limited ~ DailyAverageRelativeHumidity + DailyAverageDryBulbTemperature + DailyAverageWetBulbTemperature + DailyAverageDewPointTemperature + birthControl + field_hockey + stress_score + menses_day + day_c + I(day_c^2) + (1 | biome_id), data = complete_data)
summary(m_all)
r2(m_all)

options(na.action = "na.fail")  
model_set <- dredge(m_all, trace = TRUE)
get.models(model_set, 1)[[1]]

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
#DATA ANALYSIS -> SHANNON DIVERSITY AND LIFESTYLE FACTORS

#Shannon and HBC violin plot
b_shannon_hbc <- ggplot(candida_hbc_filtered, aes(x = birthControl, y = Shannon)) +
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

print(b_shannon_hbc)

###################################################################################################################

#Shannon over time by HBC
p_shannon_hbc <- ggplot(candida_hbc_filtered, aes(x = study_day, y = Shannon, color = birthControl)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(se = FALSE, method = "loess", size = 0.5) +
  scale_color_manual(values = hbc_cols) +
  labs(title = "Vaginal Mycobiome: Shannon Diversity Over Time by Birth Control Type",
       x = "Study Day",
       y = "Shannon Diversity",
       color = "Birth Control") +
  theme_minimal() 

print(p_shannon_hbc)

#Shannon and HBC mixed effect model 
m_shannon_HBC_time <- lmer(Shannon ~ birthControl + study_day + (1 | biome_id), na.action = na.omit, data = vaginal_rel_df_matched)
summary(m_shannon_HBC_time)
r2(m_shannon_HBC_time)

###################################################################################################################

#average Shannon by HBC boxplot
vaginal_avg_shannon_hbc_df <- candida_hbc_filtered %>%
  group_by(biome_id, birthControl) %>%
  summarise(
    n_samples = n(),
    Mean_Shannon = mean(Shannon, na.rm = TRUE),
    SD_Shannon = sd(Shannon, na.rm = TRUE),
    Min_Shannon = min(Shannon, na.rm = TRUE),
    Max_Shannon = max(Shannon, na.rm = TRUE)
  ) %>%
  ungroup()

b_avg_shannon_hbc <- ggplot(vaginal_avg_shannon_hbc_df, aes(x = birthControl, y = Mean_Shannon)) +
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

print(b_avg_shannon_hbc)

#avg Shannon and HBC linear model
m_avg_shannon_hbc <- lm(Mean_Shannon ~ birthControl, data = vaginal_avg_shannon_hbc_df)
summary(m_avg_shannon_hbc)

###################################################################################################################

#Shannon by menstruation status violin plot
b_shannon_menses <- ggplot(candida_menses_filtered, aes(x = menses_day, y = Shannon)) +
  geom_sina(aes(color = menses_day), size = 1.5, alpha = 0.7,
            maxwidth = box_width * 0.4, show.legend = FALSE) +
  geom_violin(aes(fill = menses_day), color = "black", size = 0.3,
              alpha = 0.3, width = box_width) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", size = 0.5) +
  scale_fill_manual(
    values = menses_cols,
    labels = c("menses" = "Menstruating", "not_menses" = "Not Menstruating"),
    name = "Menstruation Status"
  ) +
  scale_color_manual(
    values = menses_cols,
    guide = "none"
  ) +
  scale_x_discrete(
    name = "Menstruation Status",
    labels = c("menses" = "Menstruating", "not_menses" = "Not Menstruating")
  ) +
  scale_y_continuous(name = "Shannon Diversity") +
  labs(title = "Vaginal Mycobiome: Shannon Diversity by Menstruation Status") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", size = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.ticks.length = unit(3, "pt")
  )

print(b_shannon_menses)

###################################################################################################################

#Shannon over time by menstruation status
p_shannon_menses <- ggplot(candida_menses_filtered, aes(study_day, Shannon, color = menses_day)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = FALSE, size = 0.5) +
  scale_color_manual(values = menses_cols,
                     labels = c("menses" = "Menstruating", "not_menses" = "Not Menstruating"),
                     name = "Menstruation Status") +
  scale_y_continuous(name = "Shannon Diversity") +
  scale_x_continuous(name = "Study Day", breaks = x_breaks) +
  labs(
    title = "Vaginal Mycobiome: Shannon Diversity Over Time by Menstruation Status"
  ) +
  theme_minimal()

print(p_shannon_menses)

###################################################################################################################

#avg Shannon by menstruation status test
vaginal_avg_shannon_menses_df <- candida_menses_filtered %>%
  filter(!is.na(Shannon)) %>%
  group_by(biome_id, menses_day) %>%
  summarise(
    Mean_Shannon = mean(Shannon),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = menses_day,
    values_from = Mean_Shannon,
    names_prefix = "Shannon_"
  ) %>%
  filter(!is.na(Shannon_menses), !is.na(Shannon_not_menses))

wilcox.test(vaginal_avg_shannon_menses_df$Shannon_menses, vaginal_avg_shannon_menses_df$Shannon_not_menses, paired = TRUE)

###################################################################################################################

#Shannon by HBC and menses boxplot
clean_df <- vaginal_rel_df_matched %>%
  filter(!is.na(Shannon), !is.na(birthControl), !is.na(menses_day)) %>%
  mutate(
    menses_day = factor(menses_day, levels = c("not_menses", "menses"),
                        labels = c("Not Menstruating", "Menstruating"))
  )

ggplot(clean_df, aes(x = birthControl, y = Shannon, fill = menses_day)) +
  geom_violin(position = position_dodge(width = 0.8),
              color = "black", alpha = 0.3, width = 0.7) +
  geom_boxplot(position = position_dodge(width = 0.8),
               width = 0.2, outlier.shape = NA, color = "black", alpha = 0.6) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    size = 1.2, alpha = 0.5, shape = 21, color = "black"
  ) +
  scale_fill_manual(values = menses_cols, name = "Menstruation Status") +
  labs(
    title = "Vaginal Mycobiome: Shannon Diversity by Birth Control and Menstruation",
    x = "Birth Control",
    y = "Shannon Diversity"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.ticks.length = unit(3, "pt")
  )
###################################################################################################################

#Shannon by HBC and menses lienar mixed effect model
m_shannon_hbc_menses <- lmer(Shannon ~ menses_day + birthControl + study_day + (1 | biome_id), na.action = na.omit, data = vaginal_rel_df_matched)
summary(m_shannon_hbc_menses)
r2(m_shannon_hbc_menses)

###################################################################################################################

#DATA ANALYSIS -> SHANNON & DOMINANT SPECIES

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

###################################################################################################################

#Shannon by dominant species boxplot
species_cols <- c(setNames(c("#D95F02", "#1B9E77", "#E7298A"), top3),  "Other" = "#7570B3", "Unknown" = "#666666")

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
  theme_minimal()

###################################################################################################################

#Shannon over time by dominant species
ggplot(vaginal_rel_df_matched,
       aes(x = study_day, y = Shannon, colour = DomSpec_grp)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "loess", se = FALSE, size = 0.5) +
  scale_colour_manual(values = species_cols, name = "Dominant Species") +
  labs(title = "Vaginal Mycobiome: Shannon Diversity Over Time",
       x = "Study Day", y = "Shannon Diversity") +
  theme_minimal()

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

#transitions
dominant_species_transitions <- vaginal_rel_df_matched %>%
  group_by(biome_id, study_day) %>%
  summarise(dom = first(DominantSpecies), .groups = "drop") %>%
  arrange(biome_id, study_day) %>%
  group_by(biome_id) %>%
  mutate(switch = dom != lag(dom, default = first(dom))) %>%
  summarise(
    n_switches = sum(switch),
    n_timepoints = n(),
    first_dom = first(dom),
    final_dom = last(dom),
    .groups = "drop"
  )

View(dominant_species_transitions)

