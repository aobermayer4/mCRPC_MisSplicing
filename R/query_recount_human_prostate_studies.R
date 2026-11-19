# ============================================================
# recount3 metadata screening: prostate cancer studies
# Goal: identify datasets for AR mis-splicing burden analysis
# ============================================================

# 1. Packages
# install.packages(c("readr", "dplyr", "stringr"))

library(readr)
library(dplyr)
library(stringr)

# 2. Load your recount3 metadata export
# Change this to your file path.
file <- "data/metadata/recount3_human_prostate_study_search.csv"

metadata <- if (grepl("\\.csv$", file, ignore.case = TRUE)) {
  read_csv(file, show_col_types = FALSE)
} else {
  read_tsv(file, show_col_types = FALSE)
}

# Inspect the available columns
print(names(metadata))
cat("Total rows:", nrow(metadata), "\n")


# 3. Identify the study/project column
# Edit this manually if the automatic choice is incorrect.
possible_study_cols <- c(
  "study_accession",
  "study",
  "project",
  "project_id",
  "project_accession",
  "project_name",
  "study_title",
  "srp",
  "gse"
)

study_col <- intersect(possible_study_cols, names(metadata))[1]

if (is.na(study_col)) {
  stop(
    "No study column found. Check names(metadata) and set study_col manually."
  )
}

cat("Grouping studies by:", study_col, "\n")

# 4. Combine searchable text fields
text_cols <- names(metadata)[vapply(metadata, is.character, logical(1))]

metadata$search_text <- do.call(
  paste,
  c(metadata[text_cols], sep = " | ")
)

# 5. Flag relevant keywords
metadata <- metadata %>%
  mutate(
    study_id = as.character(.data[[study_col]]),
    prostate_flag = str_detect(
      search_text,
      regex("prostat(e|ic)|\\bPRAD\\b", ignore_case = TRUE)
    ),
    cancer_flag = str_detect(
      search_text,
      regex(
        "cancer|carcinoma|adenocarcinoma|tumou?r|neoplas|malignan|CRPC|metasta",
        ignore_case = TRUE
      )
    ),
    cellline_flag = str_detect(
      search_text,
      regex(
        "cell.line|LNCaP|22Rv1|VCaP|DU145|PC-3|PC3",
        ignore_case = TRUE
      )
    )
  )


# 6. Summarize at the STUDY level
study_summary <- metadata %>%
  filter(prostate_flag) %>%
  mutate(
    priority = case_when(
      cancer_flag & !cellline_flag ~ "1 - Review first",
      cancer_flag & cellline_flag ~ "2 - Check sample types",
      TRUE ~ "3 - Prostate, cancer unclear"
    )
  ) %>%
  arrange(priority, desc(n_samples))


write_tsv(study_summary, "data/metadata/recount3_prostate_study_summary.txt")


study_summary_t1 <- study_summary[
  which(
    study_summary$priority == "1 - Review first"
  ),
]
