# prep_data.R
# One-time script to build the pre-joined data file for the DIA Shiny app.
# Run from the psup_shiny directory: Rscript dia/prep_data.R
# Mirrors the data loading pattern from dia-analyses.Rmd lines 180-214.

library(tidyverse)

dia_path <- "../calculate-psup-dia"

# 1. Load pSup data (skip comment line starting with #)
raw_data <- read_tsv(
  file.path(dia_path, "results/psups-three-species.txt"),
  comment = "#"
) %>%
  select(orf, scer_ortholog, sample_id, pSup) %>%
  mutate(sample_id = as.character(sample_id))

# 2. Load experimental conditions
conditions <- read_tsv(
  file.path(dia_path, "data/conditions.txt")
) %>%
  mutate(sample_id = as.character(identifier)) %>%
  select(-identifier)

# 3. Join psups + conditions
data <- raw_data %>%
  left_join(conditions, by = "sample_id")

# 4. Add species column (full species names)
species_map <- c(
  "BY4743"     = "S. cerevisiae",
  "FM1071"     = "S. kudriavzevii",
  "DMKU3-1042" = "K. marxianus"
)
data <- data %>%
  mutate(species = species_map[strain])

# 5. Load SGD features for gene name lookup (~6600 genes)
features <- read_tsv(
  file.path(dia_path, "data/200325-scer-features.txt"),
  comment = "#"
) %>%
  select(ORF, gene)

data <- data %>%
  left_join(features, by = c("scer_ortholog" = "ORF"))

# Gene name fallback chain: SGD gene name -> scer_ortholog -> native orf
data <- data %>%
  mutate(gene = case_when(
    !is.na(gene) & gene != "" ~ gene,
    !is.na(scer_ortholog) & scer_ortholog != "" ~ scer_ortholog,
    TRUE ~ orf
  ))

# 6. Load functional annotations (optional, ~895 proteins)
# Deduplicate: keep one annotation per ORF (some ORFs appear in multiple categories)
anno <- read_tsv(
  file.path(dia_path, "data/scer-anno-proteins.txt")
) %>%
  select(orf, annotation) %>%
  distinct(orf, .keep_all = TRUE) %>%
  rename(anno_orf = orf)

data <- data %>%
  left_join(anno, by = c("scer_ortholog" = "anno_orf"))

# 7. Select final columns (keep rows with NA scer_ortholog for native ID queries)
data <- data %>%
  select(orf, gene, scer_ortholog, species, timepoint, start_temp, end_temp, biorep, pSup)

# 8. Write output
write_tsv(data, "dia/data/dia_psup_data.tsv")

cat("Wrote", nrow(data), "rows to dia/data/dia_psup_data.tsv\n")
cat("Columns:", paste(names(data), collapse = ", "), "\n")
cat("Unique genes:", n_distinct(data$gene), "\n")
cat("Unique species:", paste(unique(data$species), collapse = ", "), "\n")
cat("Rows with NA scer_ortholog:", sum(is.na(data$scer_ortholog)), "\n")
