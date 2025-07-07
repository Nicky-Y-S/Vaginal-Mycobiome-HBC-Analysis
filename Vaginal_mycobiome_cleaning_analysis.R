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
library(readxl)

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

###################################################################################################################

#step 8: amplification bias check

otu_df <- as.data.frame(otu_table(vaginal_phyloseq))
if (taxa_are_rows(vaginal_phyloseq)) {
  otu_df <- t(otu_df)
}
otu_df <- as.data.frame(otu_df)

asv_ids <- colnames(otu_df)
asv_sequences <- as.character(taxa_names(vaginal_phyloseq))  # DNA strings
names(asv_sequences) <- asv_ids

repeat_info <- apply(otu_df, 1, function(counts) {
  total <- sum(counts)
  max_idx <- which.max(counts)
  max_count <- counts[max_idx]
  max_seq <- asv_sequences[max_idx]
  list(max_seq = max_seq, max_count = max_count, total = total, prop = max_count / total)
})

repeat_df <- tibble::tibble(
  SampleID = rownames(otu_df),
  MaxSeq = sapply(repeat_info, `[[`, "max_seq"),
  MaxCount = sapply(repeat_info, `[[`, "max_count"),
  TotalReads = sapply(repeat_info, `[[`, "total"),
  MaxProp = sapply(repeat_info, `[[`, "prop")
)

View(repeat_df)

ggplot(repeat_df, aes(x = TotalReads, y = MaxProp)) +
  geom_point(alpha = 0.6) +
  
  # Threshold lines
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 200, linetype = "dashed", color = "orange") +
  geom_vline(xintercept = 1000, linetype = "dashed", color = "orange") +
  
  # Vertical labels
  annotate("text", x = 200, y = 1.08, label = "200 reads", angle = 90, vjust = 0, hjust = 0, color = "orange") +
  annotate("text", x = 1000, y = 1.08, label = "1000 reads", angle = 90, vjust = 0, hjust = 0, color = "orange") +
  
  # Horizontal labels
  annotate("text", x = min(repeat_df$TotalReads), y = 1, label = "MaxProp = 1.0", hjust = -0.1, vjust = -0.5, color = "red") +
  annotate("text", x = min(repeat_df$TotalReads), y = 0.7, label = "MaxProp = 0.7", hjust = -0.1, vjust = -0.5, color = "red") +
  
  scale_x_log10() +
  labs(
    title = "Sample Quality: Max ASV Proportion vs Total Reads",
    x = "Total Reads (log10 scale)",
    y = "Proportion of Top ASV"
  ) +
  coord_cartesian(ylim = c(0, 1.1), clip = "off") +
  theme_minimal() +
  theme(
    plot.margin = margin(t = 20, r = 10, b = 10, l = 30)
  )

hist_data_log <- ggplot_build(
  ggplot(repeat_df, aes(x = TotalReads)) +
    geom_histogram(bins = 50) +
    scale_x_log10()
)$data[[1]]
y_max_log <- max(hist_data_log$count)

ggplot(repeat_df, aes(x = TotalReads)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 200, linetype = "dashed", color = "orange") +
  geom_vline(xintercept = 1000, linetype = "dotted", color = "orange") +
  annotate("text", x = 200, y = y_max_log + 2, label = "200 reads", angle = 90, vjust = 0, color = "orange") +
  annotate("text", x = 1000, y = y_max_log + 2, label = "1000 reads", angle = 90, vjust = 0, color = "orange") +
  scale_x_log10() +
  labs(
    title = "Distribution of Total Reads per Sample",
    x = "Total Reads (log10 scale)",
    y = "Sample Count"
  ) +
  coord_cartesian(ylim = c(0, y_max_log + 5)) +
  theme_minimal()

repeat_df <- repeat_df %>%
  mutate(Flagged = (MaxProp == 1 & TotalReads < 200))

vaginal_phyloseq <- prune_samples(!(sample_names(vaginal_phyloseq) %in% repeat_df$SampleID[repeat_df$Flagged]), vaginal_phyloseq)

###################################################################################################################

#step 9: clean tax table 

vaginal_tax_tab <- as.data.frame(tax_table(vaginal_phyloseq))
#View(vaginal_tax_tab) #bad species format, need cleaning

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
View(vaginal_rel_metadata_df) #1092 after pruning badly sequenced samples

###################################################################################################################

#DATA CLEANING & MERGING

#HBC cleaning
hbc_df <- read.csv("/Users/yiningsun/Desktop/Tetel Lab/cleaned_Report 9-Volunteer Medical History.csv",
                             header = TRUE, stringsAsFactors = FALSE) %>%
  filter(birthControl != "Orilissa (Elagolix)") %>%
  select(biome_id, birthControl) %>%
  mutate(biome_id = as.integer(biome_id),
         birthControl = factor(birthControl,
                               levels = c("None", "Local P", "Systemic P only", "Systemic Combined (E&P)"))) 

vaginal_rel_df_matched <- vaginal_rel_metadata_df %>%
  left_join(
    hbc_df %>% select(biome_id, birthControl, taken_antibiotics),
    by = "biome_id"
  )

#medication cleaning
med_df <- read_xlsx("/Users/yiningsun/Desktop/Tetel Lab/medications.xlsx")

ssri_summary <- med_df %>%
  group_by(study_id) %>%
  summarise(took_ssri = any(SSRI == 1), .groups = "drop") %>%
  mutate(biome_id = as.integer(study_id)) %>%
  select(biome_id, took_ssri)

vaginal_rel_df_matched <- vaginal_rel_df_matched %>%
  left_join(ssri_summary, by = "biome_id")

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
      is.na(sport)           ~ NA_character_,     
      sport == "In-Season"   ~ "Field Hockey",
      sport == "None"        ~ "None",
      TRUE                   ~ "Other"            
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

#View(vaginal_rel_df_matched)

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

#C. albicans rel. abundance over time by menstruation status (by day)
candida_menses_filtered <- vaginal_rel_df_matched %>%
  filter(!is.na(menses_day))

menses_day_cols <- c("menses" = "brown", "not_menses" = "orange")

p_menses_day <- ggplot(candida_menses_filtered, aes(study_day, CA_abund_limited, color = menses_day)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = FALSE, size = 0.5) +
  scale_color_manual(values = menses_day_cols,
                     labels = c("menses" = "Menstruating", "not_menses" = "Not Menstruating"),
                     name = "Menstruation Status (By Day)") +
  scale_y_continuous(name = "Relative Abundance of C. albicans", limits = c(0, 1)) +
  scale_x_continuous(name = "Study Day", breaks = x_breaks) +
  labs(
    title = "Vaginal Mycobiome: C. albicans Relative Abundance Over Time by Menstruation Status (By Day)"
  ) +
  theme_minimal()

print(p_menses_day)

#C. albicans rel. abundance and menstruation mixed effect quadratic model
m_rel_men_quad <- lmer(CA_abund_limited ~ menses_day + day_c + I(day_c^2) + (1 | biome_id), data = vaginal_rel_df_matched)
summary(m_rel_men_quad)
r2(m_rel_men_quad)

#align each participant
candida_menses_aligned <- candida_menses_filtered %>%
  group_by(biome_id) %>%
  # Only include those with at least one "menses" day
  filter(any(menses_day == "menses", na.rm = TRUE)) %>%
  mutate(
    first_menses_day = min(study_day[menses_day == "menses"], na.rm = TRUE),
    days_since_menses = study_day - first_menses_day
  ) %>%
  ungroup()

ggplot(candida_menses_aligned, 
       aes(x = days_since_menses, y = CA_abund_limited, 
           color = menses_day, shape = menses_day)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(
    values = menses_day_cols,
    labels = c("menses" = "Menstruating", "not_menses" = "Not Menstruating"),
    name = "Menstruation Status (By Day)"
  ) +
  scale_shape_manual(
    values = c("menses" = 17, "not_menses" = 16),
    labels = c("menses" = "Menstruating", "not_menses" = "Not Menstruating"),
    name = "Menstruation Status (By Day)"
  ) +
  facet_wrap(~ biome_id, scales = "free_x") +
  scale_x_continuous(
    name = "Days Since First Menstruation",
    breaks = scales::pretty_breaks(n = 10)(seq(-100, 100, by = 5))
  ) +
  scale_y_continuous(
    name = "Relative Abundance of C. albicans",
    limits = c(0, 1)
  ) +
  labs(title = "Vaginal Mycobiome: Participant-Specific Trends in C. albicans Relative Abundance Aligned to First Menstruation Day") +
  theme_minimal()

###################################################################################################################

#summarize mean C. albicans rel. abundance per participant by menstruation status color coded by HBC
abund_means <- candida_menses_filtered %>%
  filter(!is.na(menses_day)) %>%
  group_by(biome_id, menses_day) %>%
  summarise(mean_abund = mean(CA_abund_limited, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = menses_day,
    values_from = mean_abund,
    names_prefix = "mean_"
  ) %>%
  filter(!is.na(mean_menses) & !is.na(mean_not_menses)) %>%
  left_join(
    candida_menses_filtered %>%
      select(biome_id, birthControl) %>%
      distinct(),
    by = "biome_id"
  )

ggplot(abund_means, aes(x = mean_not_menses, y = mean_menses, color = birthControl)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_color_manual(values = hbc_cols) +
  labs(
    x = "Not Menses (Mean C. albicans Relative Abundance)",
    y = "Menses (Mean C. albicans Relative Abundance)",
    color = "Birth Control",
    title = "Vaginal Mycobiome: Mean C. albicans Relative Abundance by Menstruation Status"
  ) +
  theme_minimal()

###################################################################################################################

#summarize mean C. albicans rel. abundance per participant by menstruation status color coded by SSRI
candida_ssri_filtered <- vaginal_rel_df_matched %>%
  filter(!is.na(took_ssri))

abund_means_ssri <- candida_ssri_filtered %>%
  filter(!is.na(menses_day)) %>%
  group_by(biome_id, menses_day) %>%
  summarise(mean_abund = mean(CA_abund_limited, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = menses_day,
    values_from = mean_abund,
    names_prefix = "mean_"
  ) %>%
  filter(!is.na(mean_menses) & !is.na(mean_not_menses)) %>%
  left_join(
    candida_ssri_filtered %>%
      select(biome_id, took_ssri) %>%
      distinct(),
    by = "biome_id"
  ) %>%
  mutate(ssri_label = ifelse(took_ssri, "SSRI", "No SSRI"))

ggplot(abund_means_ssri, aes(x = mean_not_menses, y = mean_menses, color = ssri_label)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_color_manual(values = c("SSRI" = "tomato", "No SSRI" = "steelblue")) +
  labs(
    x = "Not Menses (Mean C. albicans Relative Abundance)",
    y = "Menses (Mean C. albicans Relative Abundance)",
    color = "SSRI Use",
    title = "Vaginal Mycobiome: Mean C. albicans Relative Abundance by Menstruation Status and SSRI"
  ) +
  theme_minimal()

###################################################################################################################

#C. albicans rel. abundance over time by stress severity
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

#C. albicans rel. abundance and climate conditions (dry bulb, wet bulb, dew point, and relative humidity)

#dry bulb temp.
c1 <- ggplot(vaginal_rel_df_matched, aes(x = logDate, y = CA_abund_limited)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE, color = "steelblue") +
  labs(x = "Date",
       y = "Relative Abundance of C. albicans") +
  theme_minimal()

d1 <- ggplot(weather_df, aes(x = logDate, y = DailyAverageDryBulbTemperature)) +
  geom_line(color = "firebrick") +
  labs(title = "Side-by-side Comparison of Daily Average Dry Bulb Temperature (°C) and C. albicans Over Time",
       x = "Date",
       y = "Daily Average Temperature (°C)") +
  theme_minimal()

d2 <- ggplot(weather_df, aes(x = logDate, y = DailyAverageDryBulbTemperature)) +
  geom_smooth(se = FALSE, color = "firebrick") +
  labs(x = "Date",
       y = "Daily Average Temperature (°C)") +
  theme_minimal()

d1 / d2 / c1

#wet bulb temp.
w1 <- ggplot(weather_df, aes(x = logDate, y = DailyAverageWetBulbTemperature)) +
  geom_line(color = "firebrick") +
  labs(title = "Side-by-side Comparison of Daily Average Wet Bulb Temperature (°C) and C. albicans Over Time",
       x = "Date",
       y = "Daily Average Temperature (°C)") +
  theme_minimal()

w2 <- ggplot(weather_df, aes(x = logDate, y = DailyAverageWetBulbTemperature)) +
  geom_smooth(se = FALSE, color = "firebrick") +
  labs(x = "Date",
       y = "Daily Average Temperature (°C)") +
  theme_minimal()

w1 / w2 / c1

#dew point temp.
dp1 <- ggplot(weather_df, aes(x = logDate, y = DailyAverageDewPointTemperature)) +
  geom_line(color = "firebrick") +
  labs(title = "Side-by-side Comparison of Daily Average Dew Point Temperature (°C) and C. albicans Over Time",
       x = "Date",
       y = "Daily Average Temperature (°C)") +
  theme_minimal()

dp2 <- ggplot(weather_df, aes(x = logDate, y = DailyAverageDewPointTemperature)) +
  geom_smooth(se = FALSE, color = "firebrick") +
  labs(x = "Date",
       y = "Daily Average Temperature (°C)") +
  theme_minimal()

dp1 / dp2 / c1

#rel. humidity
r1 <- ggplot(weather_df, aes(x = logDate, y = DailyAverageRelativeHumidity)) +
  geom_line(color = "firebrick") +
  labs(title = "Side-by-side Comparison of Daily Average Relative Humidity (%) and C. albicans Over Time",
       x = "Date",
       y = "Relative Humidity (%)") +
  theme_minimal()

r2 <- ggplot(weather_df, aes(x = logDate, y = DailyAverageRelativeHumidity)) +
  geom_smooth(se = FALSE, color = "firebrick") +
  labs(x = "Date",
       y = "Relative Humidity (%)") +
  theme_minimal()

r1 / r2 / c1

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
m_shannon_HBC_time <- lmer(Shannon ~ birthControl + day_c + (1 | biome_id), data = vaginal_rel_df_matched)
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
  geom_sina(aes(color = birthControl), size = 1.5, alpha = 0.7, maxwidth = box_width * 0.4, show.legend = FALSE) +
  geom_violin(aes(fill = birthControl), color = "black", size = 0.3, alpha = 0.4, width = box_width) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", size = 0.5) +
  scale_fill_manual(name = "Birth Control", values = hbc_cols) +
  scale_color_manual(values = hbc_cols) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Vaginal Mycobiome: Average Shannon Diversity by Birth Control Type",
    x = "Birth Control",
    y = "Average Shannon Diversity"
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
    values = menses_day_cols,
    labels = c("menses" = "Menstruating", "not_menses" = "Not Menstruating"),
    name = "Menstruation Status (By Day)"
  ) +
  scale_color_manual(
    values = menses_day_cols,
    guide = "none"
  ) +
  scale_x_discrete(
    name = "Menstruation Status",
    labels = c("menses" = "Menstruating", "not_menses" = "Not Menstruating")
  ) +
  scale_y_continuous(name = "Shannon Diversity") +
  labs(title = "Vaginal Mycobiome: Shannon Diversity by Menstruation Status (By Day)") +
  theme_minimal()

print(b_shannon_menses)

###################################################################################################################

#Shannon over time by menstruation status (by day)
p_shannon_menses <- ggplot(candida_menses_filtered, aes(study_day, Shannon, color = menses_day)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = FALSE, size = 0.5) +
  scale_color_manual(values = menses_day_cols,
                     labels = c("menses" = "Menstruating", "not_menses" = "Not Menstruating"),
                     name = "Menstruation Status (By Day)") +
  scale_y_continuous(name = "Shannon Diversity") +
  scale_x_continuous(name = "Study Day", breaks = x_breaks) +
  labs(
    title = "Vaginal Mycobiome: Shannon Diversity Over Time by Menstruation Status (By Day)"
  ) +
  theme_minimal()

print(p_shannon_menses)

###################################################################################################################

#avg Shannon by menstruation status colored by HBC
shannon_means <- vaginal_rel_df_matched %>%
  filter(!is.na(menses_day), !is.na(Shannon)) %>%
  group_by(biome_id, menses_day) %>%
  summarise(mean_shannon = mean(Shannon, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = menses_day,
    values_from = mean_shannon,
    names_prefix = "mean_"
  ) %>%
  filter(!is.na(mean_menses) & !is.na(mean_not_menses)) %>%
  left_join(
    vaginal_rel_df_matched %>% select(biome_id, birthControl) %>% distinct(),
    by = "biome_id"
  )

ggplot(shannon_means, aes(x = mean_not_menses, y = mean_menses, color = birthControl)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  coord_equal(xlim = c(0, max(shannon_means$mean_not_menses, shannon_means$mean_menses, na.rm = TRUE)),
              ylim = c(0, max(shannon_means$mean_not_menses, shannon_means$mean_menses, na.rm = TRUE))) +
  scale_color_manual(values = hbc_cols, name = "Birth Control") +
  labs(
    x = "Not Menses (Mean Shannon Diversity)",
    y = "Menses (Mean Shannon Diversity)",
    title = "Vaginal Mycobiome: Mean Shannon Diversity by Menstruation Status"
  ) +
  theme_minimal()

###################################################################################################################

#avg Shannon by menstruation status colored by SSRI
shannon_means_ssri <- vaginal_rel_df_matched %>%
  filter(!is.na(menses_day), !is.na(Shannon), !is.na(took_ssri)) %>%
  group_by(biome_id, menses_day) %>%
  summarise(mean_shannon = mean(Shannon, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = menses_day,
    values_from = mean_shannon,
    names_prefix = "mean_"
  ) %>%
  filter(!is.na(mean_menses) & !is.na(mean_not_menses)) %>%
  left_join(
    vaginal_rel_df_matched %>%
      select(biome_id, took_ssri) %>%
      distinct(),
    by = "biome_id"
  ) %>%
  mutate(ssri_label = factor(
    ifelse(took_ssri, "SSRI", "No SSRI"),
    levels = c("No SSRI", "SSRI")
  ))

ggplot(shannon_means_ssri, aes(x = mean_not_menses, y = mean_menses, color = ssri_label)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  coord_equal(xlim = c(0, max(shannon_means_ssri$mean_not_menses, shannon_means_ssri$mean_menses, na.rm = TRUE)),
              ylim = c(0, max(shannon_means_ssri$mean_not_menses, shannon_means_ssri$mean_menses, na.rm = TRUE))) +
  scale_color_manual(values = c("No SSRI" = "steelblue", "SSRI" = "tomato"), name = "SSRI Use") +
  labs(
    x = "Not Menses (Mean Shannon Diversity)",
    y = "Menses (Mean Shannon Diversity)",
    title = "Vaginal Mycobiome: Mean Shannon Diversity by Menstruation Status"
  ) +
  theme_minimal()

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
