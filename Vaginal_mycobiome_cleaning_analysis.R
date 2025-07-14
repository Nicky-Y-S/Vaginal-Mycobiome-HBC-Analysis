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

summary(vaginal_rel_df_matched$total_reads)
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

sample_metadata <- as.data.frame(as.matrix(sample_data(vaginal_phyloseq)))
sample_metadata$SampleID <- rownames(sample_metadata)

repeat_df <- repeat_df %>%
  left_join(sample_metadata %>% select(SampleID, biome_id), by = "SampleID")

biome_ids <- repeat_df %>%
  filter(!is.na(biome_id)) %>%
  distinct(biome_id) %>%
  arrange(biome_id) %>%
  pull(biome_id)

biome_chunks <- split(biome_ids, ceiling(seq_along(biome_ids) / 10))

plot_list <- list()

for (i in seq_along(biome_chunks)) {
  this_chunk <- biome_chunks[[i]]
  this_df <- repeat_df %>% filter(biome_id %in% this_chunk)
  
  p <- ggplot(this_df, aes(x = TotalReads, y = MaxProp)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 0.7, linetype = "dashed", color = "red") +
    geom_vline(xintercept = 200, linetype = "dashed", color = "orange") +
    geom_vline(xintercept = 1000, linetype = "dashed", color = "orange") +
    annotate("text", x = 200, y = 1.08, label = "200", angle = 90, vjust = 0, hjust = 0, color = "orange") +
    annotate("text", x = 1000, y = 1.08, label = "1000", angle = 90, vjust = 0, hjust = 0, color = "orange") +
    annotate("text", x = min(this_df$TotalReads), y = 1, label = "1.0", hjust = -0.1, vjust = -0.5, color = "red") +
    annotate("text", x = min(this_df$TotalReads), y = 0.7, label = "0.7", hjust = -0.1, vjust = -0.5, color = "red") +
    scale_x_log10() +
    labs(
      x = "Total Reads (log10 scale)",
      y = "Proportion of Top ASV"
    ) +
    scale_x_log10(limits = c(10, 1e5), breaks = c(100, 1000, 10000, 100000)) +
    coord_cartesian(ylim = c(0, 1.1), clip = "off") +
    facet_wrap(~ biome_id, ncol = 5) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 8),
      plot.margin = margin(t = 20, r = 10, b = 10, l = 30)
    )
  
  plot_list[[i]] <- p
}

plot_list[[1]]
plot_list[[2]]
plot_list[[3]]
plot_list[[4]]
plot_list[[5]]
plot_list[[6]]
plot_list[[7]]

ggplot(repeat_df, aes(x = TotalReads, y = MaxProp)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 200, linetype = "dashed", color = "orange") +
  geom_vline(xintercept = 1000, linetype = "dashed", color = "orange") +
  annotate("text", x = 200, y = 1.08, label = "200 reads", angle = 90, vjust = 0, hjust = 0, color = "orange") +
  annotate("text", x = 1000, y = 1.08, label = "1000 reads", angle = 90, vjust = 0, hjust = 0, color = "orange") +
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

#HBC, sexually active, sexuality, activity level, athlete, sport
volunteer_report_df <- read.csv("cleaned_Report 9-Volunteer Medical History.csv",
                                header = TRUE, stringsAsFactors = FALSE) %>%
  filter(birthControl != "Orilissa (Elagolix)") %>%
  select(biome_id, birthControl, sexuallyActive, sexuality, activity_level, sport) %>%
  mutate(
    biome_id = as.integer(biome_id),
    birthControl = factor(birthControl,
                          levels = c("None", "Local P", "Systemic P only", "Systemic Combined (E&P)")),
    athlete = if_else(sport == "None", 0, 1)
  )

#participant 7 was sexually active on 11/29, None, oral, female
volunteer_report_df <- volunteer_report_df %>%
  mutate(sexuallyActive = if_else(biome_id == 7, 1, sexuallyActive))

vaginal_rel_df_matched <- vaginal_rel_metadata_df %>%
  left_join(volunteer_report_df %>% select(biome_id, birthControl, sexuallyActive, sexuality, activity_level, sport, athlete), by = "biome_id")

#sexual activity
sex_df <- read_csv("cleaned_Report 3-Sexual Activity.csv") %>%
  mutate(logDate = as.Date(logDate)) %>%
  select(biome_id, logDate, type_of_intercourse, gender_of_partner) %>%
  distinct() %>%
  mutate(had_sex = 1)

vaginal_rel_df_matched <- vaginal_rel_df_matched %>%
  mutate(logDate = as.Date(logDate)) %>%
  left_join(sex_df, by = c("biome_id", "logDate")) %>%
  mutate(
    had_sex = if_else(is.na(had_sex), NA_real_, had_sex),  # preserve NAs
    type_of_intercourse = if_else(had_sex == 1, type_of_intercourse, NA_character_),
    gender_of_partner = if_else(had_sex == 1, gender_of_partner, NA_character_)
  )

#field hockey
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

#ssri
med_df <- read_xlsx("medications.xlsx")

ssri_summary_df <- med_df %>%
  group_by(study_id) %>%
  summarise(took_ssri = any(SSRI == 1), .groups = "drop") %>%
  mutate(biome_id = as.integer(study_id)) %>%
  select(biome_id, took_ssri)

vaginal_rel_df_matched <- vaginal_rel_df_matched %>%
  left_join(ssri_summary_df, by = "biome_id")

#DASS
dass_df <- read.csv("cleaned_dass.csv") %>%
  mutate(
    mood_date = as.Date(Timestamp),
    biome_id = as.integer(biome_id)
  )

dass_clean <- dass_df %>%
  select(biome_id, mood_date,
         depression_score, anxiety_score, stress_score,
         depressionseverity, anxietyseverity, stressseverity) %>%
  distinct()

vaginal_rel_df_matched$logDate <- as.Date(vaginal_rel_df_matched$logDate)

closest_dass <- vaginal_rel_df_matched %>%
  inner_join(dass_clean, by = "biome_id") %>%
  mutate(date_diff = abs(as.numeric(difftime(logDate, mood_date, units = "days")))) %>%
  group_by(biome_id, logDate) %>%
  slice_min(date_diff, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(across(c(depression_score, anxiety_score, stress_score,
                  depressionseverity, anxietyseverity, stressseverity),
                ~ ifelse(date_diff <= 7, ., NA)))

vaginal_rel_df_matched <- vaginal_rel_df_matched %>%
  left_join(
    closest_dass %>%
      select(biome_id, logDate,
             depression_score, anxiety_score, stress_score,
             depressionseverity, anxietyseverity, stressseverity),
    by = c("biome_id", "logDate")
  )

#menses 
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

#study day
all_days <- seq.Date(as.Date("2022-10-13"), as.Date("2022-12-16"), by = "day")
vaginal_rel_df_matched$study_day <- match(as.Date(vaginal_rel_df_matched$logDate), all_days) - 1

#quadratic day
vaginal_rel_df_matched <- vaginal_rel_df_matched %>%
  mutate(day_c = scale(study_day, center = TRUE, scale = FALSE))

View(vaginal_rel_df_matched)

###################################################################################################################

#summary data set for C. albicans rel. abundance and Shannon
mycobiome_summary_df <- vaginal_rel_df_matched %>%
  group_by(biome_id) %>%
  summarise(
    mean_CA = mean(CA_abund_limited, na.rm = TRUE),
    mean_shannon = mean(Shannon, na.rm = TRUE)
  )

###################################################################################################################

#DATA ANALYSIS -> C. ALBICANS RELATIVE ABUNDANCE AND LIFESTYLE FACTORS

#C. albicans rel. abundance by HBC sina plot
hbc_cols <- c("None" = "#D73027", "Local P" = "#4575B4", "Systemic P only" = "#91CF60", "Systemic Combined (E&P)" = "#8073AC")

box_width <- 0.75
jitter_width <- box_width/2  

x_breaks <- seq(0, max(vaginal_rel_df_matched$study_day, na.rm = TRUE), by = 5)

ggplot(vaginal_rel_df_matched %>% filter(!is.na(birthControl)),
       aes(birthControl, CA_abund_limited)) +
  geom_violin(aes(fill = birthControl), color = "black", size = 0.3, alpha = 0.4) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", size = 0.5) +
  geom_sina(aes(color = birthControl), size = 1.5, alpha = 0.7, show.legend = FALSE) +
  scale_y_continuous(name = "Relative Abundance of C. albicans", limits = c(0, 1)) +
  scale_fill_manual(name = "Birth Control", values = hbc_cols) +
  scale_color_manual(values = hbc_cols, guide = "none") +
  labs(title = "Vaginal Mycobiome: C. albicans Relative Abundance by Birth Control",
       x = "Birth Control") +
  theme_minimal()

#C. albicans rel. abundance over time by HBC
vaginal_rel_df_matched_avg <- vaginal_rel_df_matched %>%
  group_by(biome_id, study_day) %>%
  mutate(
    CA_abund_limited = mean(CA_abund_limited, na.rm = TRUE),
    Shannon = mean(Shannon, na.rm = TRUE)
  ) %>%
  slice(1) %>%
  ungroup()

ggplot(vaginal_rel_df_matched_avg %>% filter(!is.na(birthControl)), 
       aes(study_day, CA_abund_limited, color = birthControl)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(se = FALSE, method = "loess", size = 0.5) +
  scale_color_manual(values = hbc_cols, name = "Birth Control") +
  scale_y_continuous(name = "Relative Abundance of C. albicans", limits = c(0, 1)) +
  scale_x_continuous(name = "Study Day", breaks = x_breaks) +
  labs(
    title = "Vaginal Mycobiome: C. albicans Relative Abundance Over Time by Birth Control"
  ) +
  theme_minimal()

#C. albicans rel. abundance by sexual activity sina plot
ggplot(vaginal_rel_df_matched %>% filter(!is.na(sexuallyActive), !is.na(CA_abund_limited)),
       aes(factor(sexuallyActive), CA_abund_limited)) +
  geom_violin(aes(fill = factor(sexuallyActive)), color = "black", size = 0.3, alpha = 0.4) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", size = 0.5) +
  geom_sina(aes(color = factor(sexuallyActive)), size = 1.5, alpha = 0.7, show.legend = FALSE) +
  scale_x_discrete(labels = c("0" = "No", "1" = "Yes")) +
  scale_y_continuous(name = "Relative Abundance of C. albicans", limits = c(0, 1)) +
  scale_fill_manual(name = "Sexually Active", values = c("0" = "#F8766D", "1" = "#00BFC4"),
                    labels = c("0" = "No", "1" = "Yes")) +
  scale_color_manual(values = c("0" = "#F8766D", "1" = "#00BFC4"), guide = "none") +
  labs(title = "Vaginal Mycobiome: C. albicans Relative Abundance by Sexual Activity",
       x = "Sexually Active") +
  theme_minimal()

#C. albicans rel. abundance over time by sexual activity
ggplot(vaginal_rel_df_matched_avg %>% filter(!is.na(sexuallyActive)),
       aes(study_day, CA_abund_limited, color = factor(sexuallyActive))) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "loess", se = FALSE, size = 0.5) +
  scale_color_manual(name = "Sexually Active",
                     values = c("0" = "#F8766D", "1" = "#00BFC4"),
                     labels = c("0" = "No", "1" = "Yes")) +
  scale_y_continuous(name = "Relative Abundance of C. albicans", limits = c(0, 1)) +
  scale_x_continuous(name = "Study Day", breaks = x_breaks) +
  labs(title = "Vaginal Mycobiome: C. albicans Relative Abundance Over Time by Sexual Activity") +
  theme_minimal()

#HBC interacted with sexually active
ggplot(vaginal_rel_df_matched %>% filter(!is.na(CA_abund_limited), !is.na(birthControl), !is.na(sexuallyActive)),
       aes(birthControl, CA_abund_limited)) +
  geom_violin(aes(fill = factor(sexuallyActive)), color = "black", size = 0.3, alpha = 0.4,
              position = position_dodge(width = 0.8)) +
  geom_boxplot(aes(group = interaction(birthControl, sexuallyActive)),
               width = 0.1, outlier.shape = NA, color = "black", size = 0.5,
               position = position_dodge(width = 0.8)) +
  geom_sina(aes(color = factor(sexuallyActive), group = interaction(birthControl, sexuallyActive)),
            size = 1.5, alpha = 0.7, show.legend = FALSE, position = position_dodge(width = 0.8)) +
  scale_y_continuous(name = "Relative Abundance of C. albicans", limits = c(0, 1)) +
  scale_fill_manual(name = "Sexually Active",
                    values = c("0" = "#F8766D", "1" = "#00BFC4"),
                    labels = c("0" = "No", "1" = "Yes")) +
  scale_color_manual(values = c("0" = "#F8766D", "1" = "#00BFC4"), guide = "none") +
  labs(title = "Vaginal Mycobiome: C. albicans Relative Abundance by Birth Control and Sexual Activity",
       x = "Birth Control") +
  theme_minimal()

#avg C. albicans by sexual activity sina plot
sex_activity_summary <- vaginal_rel_df_matched %>%
  select(biome_id, sexuallyActive, birthControl) %>%
  distinct() %>%
  left_join(sex_df %>%
              mutate(vaginal_with_male = grepl("vaginal", tolower(type_of_intercourse)) & tolower(gender_of_partner) == "male",
                     sex_activity_class = case_when(
                       vaginal_with_male ~ "Vaginal Sex with Male Partner",
                       TRUE ~ "Ambiguous Sexual Activity")) %>%
              select(biome_id, sex_activity_class) %>%
              distinct(), by = "biome_id") %>%
  mutate(sex_activity_class = case_when(
    sexuallyActive == 0 ~ "Not Sexually Active",
    sexuallyActive == 1 & is.na(sex_activity_class) ~ "Ambiguous Sexual Activity",
    TRUE ~ sex_activity_class)) %>%
  distinct(biome_id, .keep_all = TRUE)

#43 mislabeled
sex_activity_summary <- sex_activity_summary %>%
  mutate(sex_activity_class = if_else(biome_id == 43, "Vaginal Sex with Male Partner", sex_activity_class))

avg_candida <- vaginal_rel_df_matched %>%
  filter(!is.na(CA_abund_limited)) %>%
  group_by(biome_id) %>%
  summarise(avg_CA = mean(CA_abund_limited, na.rm = TRUE), .groups = "drop")

avg_candida_df <- avg_candida %>%
  left_join(sex_activity_summary, by = "biome_id") %>%
  mutate(
    sexuallyActive = factor(sexuallyActive, levels = c(0, 1), labels = c("No", "Yes")),
    sex_activity_class = factor(
      sex_activity_class,
      levels = c("Not Sexually Active", "Vaginal Sex with Male Partner", "Ambiguous Sexual Activity")))

sex_class_colors <- c(
  "Not Sexually Active" = "#F8766D",
  "Vaginal Sex with Male Partner" = "#1B9E77",
  "Ambiguous Sexual Activity" = "#D95F02"
)

ggplot(data = filter(avg_candida_df, !is.na(avg_CA), !is.na(sex_activity_class)),
       aes(x = sexuallyActive, y = avg_CA)) +
  geom_violin(aes(fill = sexuallyActive), color = "black", size = 0.3, alpha = 0.4,
              position = position_dodge(width = 0.8)) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", size = 0.5,
               position = position_dodge(width = 0.8)) +
  geom_sina(aes(color = sex_activity_class, group = sexuallyActive),
            size = 1.5, alpha = 0.7, position = position_dodge(width = 0.8)) +
  scale_color_manual(name = "Sexual Activity Type", values = sex_class_colors) +
  scale_fill_manual(name = "Sexually Active", values = c("No" = "#F8766D", "Yes" = "#00BFC4")) +
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  scale_y_continuous(name = "Mean Relative Abundance of C. albicans", limits = c(0, 1.1)) +
  labs(
    title = "Vaginal Mycobiome: Mean C. albicans Relative Abundance by Sexual Activity",
    x = "Sexually Active"
  ) +
  theme_minimal()

#avg C. albicans by sexual activity + hbc sina plot
ggplot(data = filter(avg_candida_df, !is.na(avg_CA), !is.na(sex_activity_class)),
       aes(x = birthControl, y = avg_CA)) +
  geom_boxplot(
    aes(fill = sexuallyActive, group = interaction(birthControl, sexuallyActive)),
    width = 0.3, color = "black", outlier.shape = NA,
    position = position_dodge(width = 0.6)
  ) +
  geom_point(
    aes(color = sex_activity_class, group = interaction(birthControl, sexuallyActive)),
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
    size = 1.5, alpha = 0.7
  ) +
  scale_fill_manual(name = "Sexually Active",
                    values = c("No" = "#F8766D80", "Yes" = "#00BFC480")) +
  scale_color_manual(name = "Sexual Activity Type",
                     values = sex_class_colors) +
  scale_y_continuous(name = "Mean Relative Abundance of C. albicans", limits = c(0, 1)) +
  labs(
    title = "Vaginal Mycobiome: Mean C. albicans Relative Abundance by Sexual Activity and Hormonal Birth Control",
    x = "Hormonal Birth Control"
  ) +
  theme_minimal()

#participant 44 & 64
plot_participant_CA <- function(participant_id, vaginal_rel_df_matched, sex_df) {
  study_start <- as.Date("2022-10-13")
  
  # Filter for this participant's sex activity
  sex_df_sub <- sex_df %>%
    filter(biome_id == participant_id, had_sex == 1, !is.na(logDate)) %>%
    mutate(study_day = as.integer(logDate - study_start))
  
  sex_study_days <- sex_df_sub$study_day
  
  # Get all sample days with valid C. albicans abundance
  sample_days <- vaginal_rel_df_matched %>%
    filter(biome_id == participant_id, !is.na(CA_abund_limited)) %>%
    pull(study_day)
  
  sex_days_no_sample <- setdiff(sex_study_days, sample_days)
  
  # Data for this participant
  df <- vaginal_rel_df_matched %>%
    filter(biome_id == participant_id, !is.na(CA_abund_limited)) %>%
    mutate(
      sex_activity_label = if_else(study_day %in% sex_study_days,
                                   "Sexual Activity Recorded", "No Record")
    )
  
  # Menstruation shading
  menses_windows <- df %>%
    filter(menses_day == "menses") %>%
    distinct(study_day) %>%
    mutate(start = study_day - 0.5, end = study_day + 0.5)
  
  # 1-day post-sex shading
  yellow_windows <- data.frame(
    start = sex_study_days + 1 - 0.5,
    end = sex_study_days + 1 + 0.5
  )
  
  # All study days
  full_study_days <- data.frame(study_day = 0:63)
  
  # Tick mark formatting
  x_labels <- full_study_days %>%
    mutate(label = as.character(study_day)) %>%
    mutate(
      label = if_else(
        study_day %in% sex_days_no_sample,
        paste0(study_day, "*"),
        label
      ),
      label_color = if_else(
        study_day %in% sex_days_no_sample,
        "orange", "black"
      ),
      label_face = if_else(
        study_day %in% sex_days_no_sample,
        "bold", "plain"
      )
    )
  
  # Plot
  p <- ggplot(df, aes(x = study_day, y = CA_abund_limited)) +
    geom_rect(data = yellow_windows,
              aes(xmin = start, xmax = end),
              ymin = -Inf, ymax = Inf,
              fill = "#FFEB99", alpha = 0.4, inherit.aes = FALSE) +
    geom_rect(data = menses_windows,
              aes(xmin = start, xmax = end),
              ymin = -Inf, ymax = Inf,
              fill = "#FADADD", alpha = 0.4, inherit.aes = FALSE) +
    geom_line(color = "gray30", size = 0.4) +
    geom_point(aes(color = sex_activity_label), size = 2.5, alpha = 0.85) +
    scale_color_manual(values = c(
      "Sexual Activity Recorded" = "#F8766D",
      "No Record" = "#00BFC4"
    )) +
    scale_y_continuous(name = "Relative Abundance of C. albicans", limits = c(0, 1)) +
    scale_x_continuous(
      name = "Study Day",
      breaks = x_labels$study_day,
      labels = x_labels$label
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(
        angle = 0,
        hjust = 0.5,
        color = x_labels$label_color,
        face = x_labels$label_face
      ),
      legend.position = "right"
    ) +
    labs(
      title = paste("Participant", participant_id, ": C. albicans Relative Abundance in the Vaginal Mycobiome Over Time"),
      subtitle = "Shaded pink = menstruation; shaded yellow = 1-day post-sex window; orange star = sex without sample",
      color = "Sexual Activity"
    )
  
  return(p)
}

plot_participant_CA(44, vaginal_rel_df_matched, sex_df)
plot_participant_CA(64, vaginal_rel_df_matched, sex_df)
plot_participant_CA(4, vaginal_rel_df_matched, sex_df)
plot_participant_CA(33, vaginal_rel_df_matched, sex_df)
plot_participant_CA(10, vaginal_rel_df_matched, sex_df)
plot_participant_CA(43, vaginal_rel_df_matched, sex_df)

###################################################################################################################

#C. albicans rel. abundance by athlete status sina plot
athlete_cols <- c("0" = "#F8766D", "1" = "#00BFC4")

ggplot(vaginal_rel_df_matched %>% filter(!is.na(athlete), !is.na(CA_abund_limited)),
       aes(factor(athlete), CA_abund_limited)) +
  geom_violin(aes(fill = factor(athlete)), color = "black", size = 0.3, alpha = 0.4) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", size = 0.5) +
  geom_sina(aes(color = factor(athlete)), size = 1.5, alpha = 0.7, show.legend = FALSE) +
  scale_x_discrete(labels = c("0" = "Non-athlete", "1" = "Athlete")) +
  scale_y_continuous(name = "Relative Abundance of C. albicans", limits = c(0, 1)) +
  scale_fill_manual(name = "Athlete Status",
                    values = c("0" = "#F8766D", "1" = "#00BFC4"),
                    labels = c("0" = "Non-athlete", "1" = "Athlete")) +
  scale_color_manual(values = c("0" = "#F8766D", "1" = "#00BFC4"), guide = "none") +
  labs(title = "Vaginal Mycobiome: C. albicans Relative Abundance by Athlete Status",
       x = "Athlete Status") +
  theme_minimal()

#C. albicans rel. abundance over time by athlete status
ggplot(vaginal_rel_df_matched_avg %>% filter(!is.na(athlete), !is.na(CA_abund_limited)),
       aes(study_day, CA_abund_limited, color = factor(athlete))) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "loess", se = FALSE, size = 0.5) +
  scale_color_manual(name = "Athlete Status",
                     values = athlete_cols,
                     labels = c("0" = "Non-athlete", "1" = "Athlete")) +
  scale_y_continuous(name = "Relative Abundance of C. albicans", limits = c(0, 1)) +
  scale_x_continuous(name = "Study Day", breaks = x_breaks) +
  labs(title = "Vaginal Mycobiome: C. albicans Relative Abundance Over Time by Athlete Status") +
  theme_minimal()

#C. albicans rel. abundance and athlete mixed effect quadratic model
m_rel_athlete_quad <- lmer(CA_abund_limited ~ athlete + day_c + I(day_c^2) + (1 | biome_id), data = vaginal_rel_df_matched)
summary(m_rel_athlete_quad)
r2(m_rel_athlete_quad)

###################################################################################################################

#C. albicans rel. abundance over time by menstruation status (by day)
menses_day_cols <- c("menses" = "brown", "not_menses" = "orange")

ggplot(vaginal_rel_df_matched_avg %>% filter(!is.na(menses_day), !is.na(CA_abund_limited)), 
       aes(study_day, CA_abund_limited, color = menses_day)) +
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

#C. albicans rel. abundance and menstruation mixed effect quadratic model
m_rel_men_quad <- lmer(CA_abund_limited ~ menses_day + day_c + I(day_c^2) + (1 | biome_id), data = vaginal_rel_df_matched)
summary(m_rel_men_quad)
r2(m_rel_men_quad)

#avg C. albicans rel. abundance per participant by menstruation status color coded by HBC
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

#avg C. albicans rel. abundance per participant by menstruation status color coded by SSRI
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

#DASS and C. albicans rel. abundance 
vaginal_rel_df_matched <- vaginal_rel_df_matched %>%
  mutate(
    stress_severity_label = factor(
      case_when(
        stressseverity == 0 ~ "Normal",
        stressseverity == 1 ~ "Mild",
        stressseverity == 2 ~ "Moderate",
        stressseverity == 3 ~ "Severe",
        stressseverity == 4 ~ "Extremely Severe"
      ),
      levels = c("Normal", "Mild", "Moderate", "Severe", "Extremely Severe")
    ),
    depression_severity_label = factor(
      case_when(
        depressionseverity == 0 ~ "Normal",
        depressionseverity == 1 ~ "Mild",
        depressionseverity == 2 ~ "Moderate",
        depressionseverity == 3 ~ "Severe",
        depressionseverity == 4 ~ "Extremely Severe"
      ),
      levels = c("Normal", "Mild", "Moderate", "Severe", "Extremely Severe")
    ),
    anxiety_severity_label = factor(
      case_when(
        anxietyseverity == 0 ~ "Normal",
        anxietyseverity == 1 ~ "Mild",
        anxietyseverity == 2 ~ "Moderate",
        anxietyseverity == 3 ~ "Severe",
        anxietyseverity == 4 ~ "Extremely Severe"
      ),
      levels = c("Normal", "Mild", "Moderate", "Severe", "Extremely Severe")
    )
  )

#C. albicans rel. abundance by stress severity box plot
candida_stress_filtered <- vaginal_rel_df_matched %>%
  filter(!is.na(stressseverity))

stress_severity_colors <- c(
  "Normal" = "#A6CEE3",
  "Mild" = "#7FB8D1",
  "Moderate" = "#56A0C6",
  "Severe" = "#2C7BB6",
  "Extremely Severe" = "#1F78B4"
)

ggplot(candida_stress_filtered,
       aes(x = stress_severity_label, y = CA_abund_limited)) +
  geom_violin(aes(fill = stress_severity_label),
              color = "black", size = 0.3,
              alpha = 0.4, width = box_width) +
  geom_boxplot(width = 0.1, outlier.shape = NA,
               color = "black", size = 0.5) +
  geom_sina(aes(color = stress_severity_label),
            size = 1.5, alpha = 0.7,
            maxwidth = box_width * 0.4,
            show.legend = FALSE) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(name = "Stress Severity", values = stress_severity_colors) +
  scale_color_manual(values = stress_severity_colors, guide = "none") +
  theme_minimal() +
  labs(title = "Vaginal Mycobiome: C. albicans Relative Abundance by Stress Severity",
       x = "Stress Severity", y = "Relative Abundance of C. albicans")

#interact with HBC
ggplot(vaginal_rel_df_matched %>%
         filter(!is.na(stress_severity_label),
                !is.na(CA_abund_limited),
                !is.na(birthControl)),
       aes(x = stress_severity_label, y = CA_abund_limited, fill = birthControl)) +
  geom_boxplot(outlier.shape = NA, width = 0.55,
               colour = "black", alpha = 0.50) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15,
                                              dodge.width = 0.55),
              shape = 21, size = 1.2, alpha = 0.50, colour = "black") +
  scale_fill_manual(values = hbc_cols, name = "Birth Control") +
  labs(title = "Vaginal Mycobiome: C. albicans Relative Abundance by Stress Severity and Birth Control",
       x = "Stress Severity", y = "Relative Abundance of C. albicans") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal()

#stress score over time by HBC
ggplot(vaginal_rel_df_matched_avg %>%
         filter(!is.na(stress_score), !is.na(birthControl)),
       aes(x = study_day, y = stress_score, color = birthControl)
) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "loess", se = FALSE, size = 0.5) +
  scale_color_manual(values = hbc_cols, name = "Birth Control") +
  scale_y_continuous(name = "Stress Score", limits = c(0, NA)) +
  scale_x_continuous(name = "Study Day", breaks = x_breaks) +
  labs(title = "Stress Score Over Time by Birth Control") +
  theme_minimal()

#C. albicans rel. abundance by depression severity box plot
candida_depression_filtered <- vaginal_rel_df_matched %>%
  filter(!is.na(depressionseverity))

depression_severity_colors <- c(
  "Normal" = "#FED976",
  "Mild" = "#FEB24C",
  "Moderate" = "#FD8D3C",
  "Severe" = "#FC4E2A",
  "Extremely Severe" = "#E31A1C"
)

ggplot(candida_depression_filtered,
       aes(x = depression_severity_label, y = CA_abund_limited)) +
  geom_violin(aes(fill = depression_severity_label),
              color = "black", size = 0.3,
              alpha = 0.4, width = box_width) +
  geom_boxplot(width = 0.1, outlier.shape = NA,
               color = "black", size = 0.5) +
  geom_sina(aes(color = depression_severity_label),
            size = 1.5, alpha = 0.7,
            maxwidth = box_width * 0.4,
            show.legend = FALSE) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(name = "Depression Severity", values = depression_severity_colors) +
  scale_color_manual(values = depression_severity_colors, guide = "none") +
  theme_minimal() +
  labs(title = "Vaginal Mycobiome: C. albicans Relative Abundance by Depression Severity",
       x = "Depression Severity", y = "Relative Abundance of C. albicans")

#interact with HBC
ggplot(vaginal_rel_df_matched %>%
         filter(!is.na(depression_severity_label),
                !is.na(CA_abund_limited),
                !is.na(birthControl)),
       aes(x = depression_severity_label, y = CA_abund_limited, fill = birthControl)) +
  geom_boxplot(outlier.shape = NA, width = 0.55,
               colour = "black", alpha = 0.50) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15,
                                              dodge.width = 0.55),
              shape = 21, size = 1.2, alpha = 0.50, colour = "black") +
  scale_fill_manual(values = hbc_cols, name = "Birth Control") +
  labs(title = "Vaginal Mycobiome: C. albicans Relative Abundance by Depression Severity and Birth Control",
       x = "Depression Severity", y = "Relative Abundance of C. albicans") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal()

#depression score over time by HBC
ggplot(vaginal_rel_df_matched_avg %>%
         filter(!is.na(depression_score), !is.na(birthControl)),
       aes(x = study_day, y = depression_score, color = birthControl)
) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "loess", se = FALSE, size = 0.5) +
  scale_color_manual(values = hbc_cols, name = "Birth Control") +
  scale_y_continuous(name = "Depression Score", limits = c(0, NA)) +
  scale_x_continuous(name = "Study Day", breaks = x_breaks) +
  labs(title = "Depression Score Over Time by Birth Control") +
  theme_minimal()

#C. albicans rel. abundance by anxiety severity box plot 
candida_anxiety_filtered <- vaginal_rel_df_matched %>%
  filter(!is.na(anxietyseverity))

anxiety_severity_colors <- c(
  "Normal" = "#E5F5E0",
  "Mild" = "#A1D99B",
  "Moderate" = "#74C476",
  "Severe" = "#31A354",
  "Extremely Severe" = "#006D2C"
)

ggplot(candida_anxiety_filtered,
       aes(x = anxiety_severity_label, y = CA_abund_limited)) +
  geom_violin(aes(fill = anxiety_severity_label),
              color = "black", size = 0.3,
              alpha = 0.4, width = box_width) +
  geom_boxplot(width = 0.1, outlier.shape = NA,
               color = "black", size = 0.5) +
  geom_sina(aes(color = anxiety_severity_label),
            size = 1.5, alpha = 0.7,
            maxwidth = box_width * 0.4,
            show.legend = FALSE) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(name = "Anxiety Severity", values = anxiety_severity_colors) +
  scale_color_manual(values = anxiety_severity_colors, guide = "none") +
  theme_minimal() +
  labs(title = "Vaginal Mycobiome: C. albicans Relative Abundance by Anxiety Severity",
       x = "Anxiety Severity", y = "Relative Abundance of C. albicans")

#interact with HBC
ggplot(vaginal_rel_df_matched %>%
         filter(!is.na(anxiety_severity_label),
                !is.na(CA_abund_limited),
                !is.na(birthControl)),
       aes(x = anxiety_severity_label, y = CA_abund_limited, fill = birthControl)) +
  geom_boxplot(outlier.shape = NA, width = 0.55,
               colour = "black", alpha = 0.50) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15,
                                              dodge.width = 0.55),
              shape = 21, size = 1.2, alpha = 0.50, colour = "black") +
  scale_fill_manual(values = hbc_cols, name = "Birth Control") +
  labs(title = "Vaginal Mycobiome: C. albicans Relative Abundance by Anxiety Severity and Birth Control",
       x = "Anxiety Severity", y = "Relative Abundance of C. albicans") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal()

#anxiety score over time by HBC
ggplot(vaginal_rel_df_matched_avg %>%
         filter(!is.na(anxiety_score), !is.na(birthControl)),
       aes(x = study_day, y = anxiety_score, color = birthControl)
) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "loess", se = FALSE, size = 0.5) +
  scale_color_manual(values = hbc_cols, name = "Birth Control") +
  scale_y_continuous(name = "Anxiety Score", limits = c(0, NA)) +
  scale_x_continuous(name = "Study Day", breaks = x_breaks) +
  labs(title = "Anxiety Score Over Time by Birth Control") +
  theme_minimal()

#C. albicans rel. abundace by anxiety score
ggplot(vaginal_rel_df_matched %>%
         filter(!is.na(anxiety_score), !is.na(CA_abund_limited), !is.na(birthControl)),
       aes(x = anxiety_score, y = CA_abund_limited, color = birthControl)
) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "loess", se = FALSE, size = 0.5) +
  scale_color_manual(values = hbc_cols, name = "Birth Control") +
  scale_x_continuous(name = "Anxiety Score", limits = c(0, NA)) +
  scale_y_continuous(name = "Relative Abundance of C. albicans", limits = c(0, 1)) +
  labs(title = "Vaginal Mycobiome: C. albicans Relative Abundance by Anxiety Score") +
  theme_minimal()

###################################################################################################################
#DATA ANALYSIS -> SHANNON DIVERSITY AND LIFESTYLE FACTORS

#Shannon and HBC sina plot
ggplot(candida_hbc_filtered, aes(x = birthControl, y = Shannon)) +
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

#sexually active
ggplot(vaginal_rel_df_matched %>%
         filter(!is.na(Shannon),
                !is.na(sexuallyActive)),
       aes(x = factor(sexuallyActive), y = Shannon, fill = factor(sexuallyActive))) +
  geom_boxplot(outlier.shape = NA, width = 0.55,
               colour = "black", alpha = 0.5) +
  geom_jitter(position = position_jitter(width = 0.15),
              shape = 21, size = 1.2, alpha = 0.5, colour = "black") +
  scale_fill_manual(
    name = "Sexually Active",
    values = c("0" = "#F8766D", "1" = "#00BFC4"),
    labels = c("0" = "No", "1" = "Yes")
  ) +
  scale_x_discrete(
    name = "Sexual Activity",
    labels = c("0" = "No", "1" = "Yes")
  ) +
  labs(
    title = "Vaginal Mycobiome: Shannon Diversity by Sexual Activity",
    y = "Shannon Diversity"
  ) +
  theme_minimal()

#HBC interact with sexually active
ggplot(vaginal_rel_df_matched %>%
         filter(!is.na(Shannon),
                !is.na(birthControl),
                !is.na(sexuallyActive)),
       aes(x = birthControl, y = Shannon, fill = factor(sexuallyActive))) +
  geom_boxplot(aes(group = interaction(birthControl, sexuallyActive)),
               outlier.shape = NA, width = 0.55,
               colour = "black", alpha = 0.5,
               position = position_dodge(width = 0.55)) +
  geom_jitter(aes(group = interaction(birthControl, sexuallyActive)),
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.55),
              shape = 21, size = 1.2, alpha = 0.5, colour = "black") +
  scale_fill_manual(
    name = "Sexually Active",
    values = c("0" = "#F8766D", "1" = "#00BFC4"),
    labels = c("0" = "No", "1" = "Yes")
  ) +
  labs(
    title = "Vaginal Mycobiome: Shannon Diversity by Birth Control and Sexual Activity",
    x = "Birth Control", y = "Shannon Diversity"
  ) +
  theme_minimal()

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

#PCoA by dominant species
ordination_dom <- ordinate(vaginal_phyloseq_rel, method = "PCoA", distance = "bray")
ordination_dom <- as.data.frame(ordination_dom$vectors) %>%
  tibble::rownames_to_column("SampleID") %>%
  left_join(
    vaginal_rel_df_matched %>% select(SampleID, DominantSpecies),
    by = "SampleID"
  ) %>%
  filter(!is.na(DominantSpecies))

ggplot(ordination_dom, aes(x = Axis.1, y = Axis.2, color = DominantSpecies)) +
  geom_point(alpha = 0.6, size = 1.8) +
  stat_ellipse(type = "norm", size = 0.6, linetype = "dashed") +
  scale_color_manual(values = species_cols) +
  labs(
    title = "PCoA (Bray-Curtis) of Vaginal Mycobiome by Dominant Species",
    x = paste0("PCoA1 (", round(ordination_dom$values$Relative_eig[1] * 100, 1), "%)"),
    y = paste0("PCoA2 (", round(ordination_dom$values$Relative_eig[2] * 100, 1), "%)"),
    color = "Dominant Species"
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
  theme_minimal()

#Shannon diversity, dominant species, and HBC mixed effect model
m_shannon_dom_hbc_time <- lmer(Shannon ~ DomSpec_grp + birthControl + study_day + (1 | biome_id), data = vaginal_rel_df_matched)
summary(m_shannon_dom_hbc_time)
r2(m_shannon_dom_hbc_time)

###################################################################################################################
#DATA ANALYSIS -> INDIVIDUALS WITH MULTIPLE SAMPLES PER DAY

abundant_daily_samples <- vaginal_rel_df_matched %>%
  filter(!is.na(CA_abund_limited), !is.na(timestamp)) %>%
  group_by(biome_id, study_day) %>%
  filter(n() >= 3) %>%
  ungroup() %>%
  mutate(time_of_day = format(ymd_hms(timestamp), "%H:%M:%S"))

unique_ids <- unique(abundant_daily_samples$biome_id)

p_within_day <- lapply(unique_ids, function(id) {
  ggplot(
    abundant_daily_samples %>% filter(biome_id == id),
    aes(x = hms::as_hms(time_of_day), y = CA_abund_limited)
  ) +
    geom_point(size = 1.5, alpha = 0.7, color = "#D95F02") +
    geom_line(aes(group = logDate), color = "#D95F02", alpha = 0.5, size = 0.4) +
    facet_wrap(~ logDate, scales = "free_x", ncol = 4) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title = paste("Participant", id, "– C. albicans by Time of Day"),
      x = "Time of Day",
      y = "Relative Abundance of C. albicans"
    ) +
    theme_minimal()
})

print(p_within_day[[1]])
print(p_within_day[[2]])
print(p_within_day[[3]])
print(p_within_day[[4]])


