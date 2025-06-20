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

###################################################################################################################

#vaginal dominant species
vaginal_dominant_species <- apply(otu_table(vaginal_phyloseq), 1, function(x) {
  spec <- vaginal_tax_tab$Species[which.max(x)]
  spec <- ifelse(spec == "", NA, spec) 
  ifelse(is.na(spec), "Unknown", spec) #if dominant species is empty, turn it into NA and then replace with Unknown
})
sample_data(vaginal_phyloseq)$DominantSpecies <- vaginal_dominant_species

###################################################################################################################

#vaginal alpha diversity
vaginal_alpha_diversity <- estimate_richness(vaginal_phyloseq, measures = c("Shannon"))$Shannon
sample_data(vaginal_phyloseq)$Shannon <- vaginal_alpha_diversity

###################################################################################################################

#transform to rel. abundances
vaginal_phyloseq_rel <- transform_sample_counts(vaginal_phyloseq, function(x) x / sum(x))

#add C. albicans rel. abundance to sample data
ca_phy <- subset_taxa(vaginal_phyloseq_rel, Species == "Candida_albicans")
ca_abund <- rowSums( otu_table(ca_phy)[ , , drop = FALSE ] )
sample_data(vaginal_phyloseq_rel)$CA_abund <- ca_abund[ sample_names(vaginal_phyloseq_rel) ]

#add globosa rel. abundance to sample data
cg_phy <- subset_taxa(vaginal_phyloseq_rel, Species == "globosa")
cg_abund <- rowSums(otu_table(cg_phy)[, , drop = FALSE])
sample_data(vaginal_phyloseq_rel)$CG_abund <- cg_abund[ sample_names(vaginal_phyloseq_rel) ]

###################################################################################################################

#check all other species' rel abundance is 0 for samples with C. albicans rel abundance of 1
otu_rel_mat <- as.data.frame(otu_table(vaginal_phyloseq_rel))
samples_with_CA_1 <- sample_names(vaginal_phyloseq_rel)[sample_data(vaginal_phyloseq_rel)$CA_abund == 1]
otu_rel_CA1 <- otu_rel_mat[samples_with_CA_1, , drop = FALSE]
ca_taxa <- taxa_names(subset_taxa(vaginal_phyloseq_rel, Species == "Candida_albicans"))
otu_rel_CA1_non_CA <- otu_rel_CA1[, !colnames(otu_rel_CA1) %in% ca_taxa, drop = FALSE]
sum_nonzero_other_species <- sum(otu_rel_CA1_non_CA != 0)
cat("Non-zero counts of other species in samples where C. albicans = 1:", sum_nonzero_other_species, "\n")

###################################################################################################################

#convert sample_data to csv
vaginal_rel_metadata_df <- as(sample_data(vaginal_phyloseq_rel), "data.frame")
vaginal_rel_metadata_df$biome_id <- as.integer(vaginal_rel_metadata_df$biome_id)
#View(vaginal_rel_metadata_df) #1140 entries

###################################################################################################################

#read in participant data
HBC_data <- read.csv("cleaned_Report 9-Volunteer Medical History.csv",
                             header = TRUE, stringsAsFactors = FALSE) %>%
  filter(birthControl != "Orilissa (Elagolix)") %>%
  select(biome_id, birthControl) %>%
  mutate(biome_id = as.integer(biome_id),
         birthControl = factor(birthControl,
                               levels = c("None", "Local P", "Systemic P only", "Systemic Combined (E&P)"))) #fewer number of participants for medical history survey

#merge HBC with rel. abundance data
vaginal_rel_metadata_hbc_df <- vaginal_rel_metadata_df %>%
  left_join(
    HBC_data %>% select(biome_id, birthControl),
    by = "biome_id"
  )

#drop anyone with NA HBC info
vaginal_rel_metadata_hbc_df <- vaginal_rel_metadata_hbc_df %>%
  filter(!is.na(birthControl))
View(vaginal_rel_metadata_hbc_df) #1084 entries

###################################################################################################################

#C. albicans rel. abundance by HBC boxplot (%)
hbc_cols <- c("None" = "#D73027", "Local P" = "#4575B4", "Systemic P only" = "#91CF60", "Systemic Combined (E&P)" = "#8073AC")

box_width <- 0.75
jitter_width <- box_width/2  

ggplot(vaginal_rel_metadata_hbc_df, aes(x = birthControl, y = CA_abund)) +
  geom_jitter(aes(color = factor(biome_id)),
              shape = 16, size = 1, alpha = 0.6, width = jitter_width, show.legend = FALSE) +
  geom_boxplot(aes(fill = birthControl),
               width = box_width, colour = "black", size = 0.3, alpha = 0.5) +
  scale_y_continuous(labels = percent_format(1), limits = c(0, 1)) +
  scale_fill_manual(name = "Birth Control", values = hbc_cols) +
  scale_color_viridis_d(guide = FALSE, option = "turbo") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", size = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.ticks.length = unit(3, "pt"),
        legend.position = "right") +
  labs(title = "C. albicans Relative Abundance by Birth Control Type",
       x = "Birth Control",
       y = "Relative Abundance (%)")

###################################################################################################################

#C. albicans rel. abundance by HBC boxplot (log and constant)
vaginal_rel_metadata_hbc_df$CA_log <- log(vaginal_rel_metadata_hbc_df$CA_abund + 1e-6)
ggplot(vaginal_rel_metadata_hbc_df, aes(x = birthControl, y = CA_abund)) +
  geom_jitter(aes(color = factor(biome_id)),
              shape = 16, size = 1, alpha = 0.6, width = jitter_width, show.legend = FALSE) +
  geom_boxplot(aes(fill = birthControl),
               width = box_width, colour = "black", size = 0.3, alpha = 0.5,
               outlier.shape = NA) +
  scale_y_log10(
    breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 0.5),
    limits = c(1e-4, 0.5)
  ) +
  scale_fill_manual(name = "Birth Control", values = hbc_cols) +
  scale_color_viridis_d(guide = FALSE, option = "turbo") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", size = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.ticks.length = unit(3, "pt"),
        legend.position = "right") +
  labs(title = "C. albicans (Log Relative Abundance) by Birth Control Type",
       x = "Birth Control",
       y = "Log Relative Abundance")

###################################################################################################################

#C. albicans rel. abundance by HBC boxplot (logit; drop samples with rel. abundance = 1 or 0)
vaginal_rel_metadata_hbc_logit <- vaginal_rel_metadata_hbc_df %>%
  filter(CA_abund > 0 & CA_abund < 1) %>%
  mutate(CA_logit = log(CA_abund / (1 - CA_abund)))

ggplot(vaginal_rel_metadata_hbc_logit, aes(x = birthControl, y = CA_logit)) +
  geom_jitter(aes(color = factor(biome_id)),
              shape = 16, size = 1, alpha = 0.6, width = jitter_width, show.legend = FALSE) +
  geom_boxplot(aes(fill = birthControl),
               width = box_width, colour = "black", size = 0.3, alpha = 0.5) +
  scale_fill_manual(name = "Birth Control", values = hbc_cols) +
  scale_color_viridis_d(guide = FALSE, option = "turbo") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", size = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.ticks.length = unit(3, "pt"),
        legend.position = "right") +
  labs(title = "C. albicans (Logit Relative Abundance) by Birth Control Type",
       x = "Birth Control",
       y = "Logit Relative Abundance")

###################################################################################################################

#C. albicans rel. abundance by HBC boxplot (square root)
vaginal_rel_metadata_hbc_df$CA_sqrt <- sqrt(vaginal_rel_metadata_hbc_df$CA_abund)
ggplot(vaginal_rel_metadata_hbc_df, aes(x = birthControl, y = CA_sqrt)) +
  geom_jitter(aes(color = factor(biome_id)),
              shape = 16, size = 1, alpha = 0.6, width = jitter_width, show.legend = FALSE) +
  geom_boxplot(aes(fill = birthControl),
               width = box_width, colour = "black", size = 0.3, alpha = 0.5) +
  scale_fill_manual(name = "Birth Control", values = hbc_cols) +
  scale_color_viridis_d(guide = FALSE, option = "turbo") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", size = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.ticks.length = unit(3, "pt"),
        legend.position = "right") +
  labs(title = "C. albicans (Square Root Relative Abundance) by Birth Control Type",
       x = "Birth Control",
       y = "Square Root Relative Abundance")

###################################################################################################################

#C. albicans rel. abundance by HBC violin and sina plot (%)
ggplot(vaginal_rel_metadata_hbc_df, aes(x = birthControl, y = CA_abund)) +
  geom_sina(aes(color = factor(biome_id)), size = 1.5, alpha = 0.7, maxwidth = box_width * 0.4, show.legend = FALSE) +
  geom_violin(aes(fill = birthControl), color = "black", size = 0.3, alpha = 0.4, width = box_width) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", size = 0.5) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_fill_manual(name = "Birth Control", values = hbc_cols) +
  scale_color_viridis_d(option = "turbo", guide = FALSE) +
  theme_minimal(base_size = 14) +
  labs(title = "C. albicans Relative Abundance by Birth Control Type", x = "Birth Control", y = "Relative Abundance (%)") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", size = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.ticks.length = unit(3, "pt"),
        legend.position = "right")

###################################################################################################################

#C. albicans rel. abundance over time by HBC spaghetti plot
all_days <- seq.Date(as.Date("2022-10-13"), as.Date("2022-12-16"), by = "day")
vaginal_rel_metadata_hbc_df$study_day <- match(as.Date(vaginal_rel_metadata_hbc_df$logDate), all_days) - 1

ggplot(vaginal_rel_metadata_hbc_df, aes(x = study_day, y = CA_abund, color = birthControl)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(se = FALSE, method = "loess", size = 0.5) +
  scale_color_manual(values = hbc_cols) +
  labs(title = "C. albicans Relative Abundance Over Time by Birth Control Type",
       x = "Study Day",
       y = "Relative Abundance of C. albicans",
       color = "Birth Control") +
  theme_minimal() #splines shows a M shape, linear mixed effect model may not be appropriate

###################################################################################################################

#Shannon diversity over time by HBC spaghetti plot
ggplot(vaginal_rel_metadata_hbc_df, aes(x = study_day, y = Shannon, color = birthControl)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(se = FALSE, method = "loess", size = 0.5) +
  scale_color_manual(values = hbc_cols) +
  labs(title = "Shannon Diversity Over Time by Birth Control Type",
       x = "Study Day",
       y = "Shannon Diversity",
       color = "Birth Control") +
  theme_minimal() #spline reveals a linear mixed effect model seems appropriate

###################################################################################################################

#Shannon diversity and HBC (sample-level) mixed effect model 
m_shannon_HBC_time <- lmer(Shannon ~ birthControl + study_day + (1 | biome_id), data = vaginal_rel_metadata_hbc_df)
summary(m_shannon_HBC_time)
r2(m_shannon_HBC_time)

###################################################################################################################

#average Shannon diversity by HBC boxplot
vaginal_avg_shannon_df <- vaginal_rel_metadata_hbc_df %>%
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
    title = "Average Shannon Diversity by Birth Control Type",
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
###################################################################################################################

#avg Shannon diversity and HBC linear model
m_avg_shannon_hbc <- lm(Mean_Shannon ~ birthControl, data = vaginal_avg_shannon_df)
summary(m_avg_shannon_hbc)

###################################################################################################################

#Shannon Diversity vs. C. albicans relative abundance by HBC
ggplot(vaginal_rel_metadata_hbc_df, aes(x = CA_abund, y = Shannon, color = birthControl)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(se = FALSE, method = "loess", size = 0.5) +
  scale_color_manual(values = hbc_cols) +
  labs(
    title = "Shannon Diversity vs. C. albicans Abundance by Birth Control Type",
    x = "Relative Abundance of C. albicans",
    y = "Shannon Diversity",
    color = "Birth Control") +
  theme_minimal()

###################################################################################################################

#Other species vs. C. albicans relative abundance by HBC
ggplot(vaginal_rel_metadata_hbc_df, aes(x = CA_abund, y = CG_abund, color = birthControl)) +
  geom_point(alpha = 0.4, size = 1) +
  #geom_smooth(se = FALSE, method = "loess", size = 0.5) +
  scale_color_manual(values = hbc_cols) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Globosa vs. C. albicans Abundance by Birth Control Type",
    x = "Relative Abundance of C. albicans",
    y = "Relative Abundance of Globosa",
    color = "Birth Control") +
  theme_minimal()

###################################################################################################################

#DASS importing

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

vaginal_rel_metadata_hbc_df$logDate <- as.Date(vaginal_rel_metadata_hbc_df$logDate)

dass_clean <- dass_df %>%
  select(biome_id, mood_date, stress_score, stress_severity) %>%
  distinct()

closest_matches <- vaginal_rel_metadata_hbc_df %>%
  inner_join(dass_clean, by = "biome_id") %>%
  mutate(date_diff = abs(as.numeric(difftime(logDate, mood_date, units = "days")))) %>%
  group_by(biome_id, logDate) %>%
  slice_min(order_by = date_diff, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    stress_score = ifelse(date_diff <= 7, stress_score, NA)  # NA if >7 days
    # do NOT change stress_severity – keep it even if score is NA
  )

vaginal_rel_metadata_hbc_df_matched <- vaginal_rel_metadata_hbc_df %>%
  left_join(
    closest_matches %>%
      select(biome_id, logDate, mood_date, stress_score, stress_severity),
    by = c("biome_id", "logDate")
  )

#View(vaginal_rel_metadata_hbc_df_matched) #1084 entries

#this dataset contains stress score of participants WITH valid HBC info (1084 entries vs. 1140 entries for entire mycobiome data)
#even if stress_score is NA because no mood log is within ±7 days of the swab, stress_severity is still retained from the nearest mood_date, even if it's outside the ±7 day window
#this is Alice's approach, is that reasonable?
#only 1 sample using systemic P with stress severity of 4
sum(vaginal_rel_metadata_hbc_df_matched$stress_severity == 4 & vaginal_rel_metadata_hbc_df_matched$birthControl == "Systemic P only")

###################################################################################################################

#stress score over time by HBC spaghetti plot
ggplot(vaginal_rel_metadata_hbc_df_matched, aes(x = study_day, y = stress_score, color = birthControl)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(se = FALSE, method = "loess", size = 0.5) +
  scale_color_manual(values = hbc_cols) +
  labs(
    title = "Stress Score Over Time by Birth Control Type",
    x = "Study Day",
    y = "Stress Score",
    color = "Birth Control"
  ) +
  theme_minimal()

###################################################################################################################

#C. albicans abundance and stress score over time by HBC
p1 <- ggplot(vaginal_rel_metadata_hbc_df_matched, 
             aes(x = study_day, y = stress_score, color = birthControl)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(se = FALSE, method = "loess", size = 0.5) +
  scale_color_manual(values = hbc_cols) +
  labs(title = "Stress Score Over Time", x = "Study Day", y = "Stress Score", color = "Birth Control") +
  theme_minimal()

p2 <- ggplot(vaginal_rel_metadata_hbc_df_matched, 
             aes(x = study_day, y = CA_abund, color = birthControl)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(se = FALSE, method = "loess", size = 0.5) +
  scale_color_manual(values = hbc_cols) +
  labs(title = "C. albicans Relative Abundance Over Time", x = "Study Day", y = "C. albicans Abundance", color = "Birth Control") +
  theme_minimal()

p1 / p2 
###################################################################################################################

#Shannon diversity and stress score over time by HBC

p3 <- ggplot(vaginal_rel_metadata_hbc_df_matched, 
             aes(x = study_day, y = stress_score, color = birthControl)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(se = FALSE, method = "loess", size = 0.5) +
  scale_color_manual(values = hbc_cols) +
  labs(title = "Stress Score Over Time", x = "Study Day", y = "Stress Score", color = "Birth Control") +
  theme_minimal()

p4 <- ggplot(vaginal_rel_metadata_hbc_df_matched, 
             aes(x = study_day, y = Shannon, color = birthControl)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(se = FALSE, method = "loess", size = 0.5) +
  scale_color_manual(values = hbc_cols) +
  labs(title = "Shannon Diversity Over Time", x = "Study Day", y = "Shannon Diversity", color = "Birth Control") +
  theme_minimal()

p3 / p4
###################################################################################################################

#
m_stress_hbc_day <- lmer(stress_score ~ birthControl + study_day + (1 | biome_id),
                         data = vaginal_rel_metadata_hbc_df_matched)
summary(m_stress_hbc_day)

aov_stress <- aov(stress_score ~ birthControl, data = vaginal_rel_metadata_hbc_df_matched)
summary(aov_stress)

###################################################################################################################
#C. albicans rel. abundance by stress severity and HBC
vaginal_rel_metadata_hbc_df_matched$stress_severity <- factor(
  vaginal_rel_metadata_hbc_df_matched$stress_severity,
  levels = 0:4,
  labels = c("Normal", "Mild", "Moderate", "Severe", "Extremely Severe")
)

# Filter out NA severity (Participant 31 never filled out DASS survey)
df_plot <- vaginal_rel_metadata_hbc_df_matched %>%
  filter(!is.na(stress_severity))

ggplot(df_plot, aes(x = stress_severity, y = CA_abund, fill = birthControl)) +
  geom_jitter(aes(color = factor(biome_id)), shape = 16, size = 1, alpha = 0.6,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              show.legend = FALSE) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.5,
               colour = "black", size = 0.3, alpha = 0.5, outlier.shape = NA) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_fill_manual(name = "Birth Control", values = hbc_cols) +
  scale_color_viridis_d(option = "turbo", guide = FALSE) +
  theme_minimal(base_size = 14) +
  labs(title = "C. albicans Relative Abundance by Stress Severity and Birth Control Type",
       x = "Stress Severity", y = "Relative Abundance (%)") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", size = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.ticks.length = unit(3, "pt"),
        legend.position = "right")

###################################################################################################################

#Shannon diversity by stress severity and HBC
ggplot(df_plot, aes(x = stress_severity, y = Shannon, fill = birthControl)) +
  geom_jitter(aes(color = factor(biome_id)), shape = 16, size = 1, alpha = 0.6,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              show.legend = FALSE) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.5,
               colour = "black", size = 0.3, alpha = 0.5, outlier.shape = NA) +
  scale_fill_manual(name = "Birth Control", values = hbc_cols) +
  scale_color_viridis_d(option = "turbo", guide = FALSE) +
  theme_minimal(base_size = 14) +
  labs(title = "Shannon Diversity by Stress Severity and Birth Control Type",
       x = "Stress Severity", y = "Shannon Diversity") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", size = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.ticks.length = unit(3, "pt"),
        legend.position = "right")

###################################################################################################################

ggplot(vaginal_rel_metadata_hbc_df_matched, aes(x = stress_score, y = Shannon, color = birthControl)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE) +
  facet_wrap(~ birthControl) +
  scale_color_manual(name = "Birth Control", values = hbc_cols) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Shannon Diversity vs. Stress Score by Birth Control Type",
    x = "Stress Score",
    y = "Shannon Diversity"
  ) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "right"
  )

ggplot(vaginal_rel_metadata_hbc_df_matched, aes(x = stress_score, y = CA_abund, color = birthControl)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE) +
  facet_wrap(~ birthControl) +
  scale_color_manual(name = "Birth Control", values = hbc_cols) +
  theme_minimal(base_size = 14) +
  labs(
    title = "C. albicans Relative Abundance vs. Stress Score by Birth Control Type",
    x = "Stress Score",
    y = "C. albicans Relative Abundance"
  ) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "right"
  )

sum(is.na(vaginal_rel_metadata_hbc_df_matched$stress_score))
table(vaginal_rel_metadata_hbc_df_matched$birthControl)
