# Command line arguments --------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]
fastq_dir <- args[2]

if (any(grepl(x = args[1], pattern = "help"), length(args) != 2)) {
  message("\nUsage: Rscript renamer.R [TABLE] [DIRECTORY]")
  message(
    "Uses sample names from a TABLE to rename fastq files in DIRECTORY, ",
    "matching with the index sequence."
  )
  message(
    "Provided TABLE should be the spreadsheet (.xlsx) downloaded from the GSC ",
    "submission portal, containing all the information for the overall ",
    "submission and list of samples/indices in each pool."
  )
  q(save = "no")
}


# Setup -------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readxl)
  library(purrr)
  library(stringr)
  library(dplyr)
})

load_excel <- function(excel_file) {
  suppressMessages(
    skip_lines <- read_xlsx(excel_file) %>%
      pull(1) %>%
      grep(x = ., pattern = "Original Source Name")
  )
  read_xlsx(path = excel_file, skip = skip_lines) %>%
    janitor::clean_names()
}

rev_comp <- function(x) {
  stringi::stri_reverse(chartr(old = "ATGC", new = "TACG", x = x))
}

list_files <- function(path_, pattern_) {
  found_file <- list.files(path = path_, pattern = pattern_)
  if (length(found_file) == 0) {
    return(NA)
  } else {
    found_file
  }
}


# Load and clean GSC sheet ------------------------------------------------

summary_sheet <- load_excel(excel_file = input_file)

remaining_sheets_1 <- lapply(
  seq(2, nrow(summary_sheet) + 1),
  function(x) {
    read_xlsx(path = input_file, sheet = x) %>%
      janitor::clean_names()
  }
) %>% setNames(summary_sheet$pool_id)

remaining_sheets_2 <- bind_rows(remaining_sheets_1, .id = "pool_id")

samples_01 <- left_join(
  summary_sheet,
  remaining_sheets_2,
  by = "pool_id",
  multiple = "all"
) %>%
  select(sub_library_id, pool_id, tube_label, index, everything())

output_file <- gsub(
  x = input_file,
  pattern = "\\.xlsx",
  replacement = "_clean.tsv"
)

write.table(
  x = samples_01,
  file = output_file,
  sep = "\t",
  row.names = FALSE
)


# Fix up indices for renaming ---------------------------------------------

samples_02 <- samples_01 %>%
  rename("index_fwdi7_fwdi5" = index) %>%
  mutate(index_fwdi7_revi5 = paste0(
    str_extract(index_fwdi7_fwdi5, '^[A-Z]{8}(?=-)'),
    "-",
    map_chr(str_extract(index_fwdi7_fwdi5, '(?<=-)[A-Z]{8}$'), rev_comp)
  )) %>%
  select(sub_library_id, index_fwdi7_fwdi5, index_fwdi7_revi5)


# Create old names (file searching) and new names (table data) ------------

if (!grepl("/", fastq_dir)) {
  fastq_dir <- paste0(fastq_dir, "/")
}

samples_03 <- samples_02 %>%
  mutate(
    old_R1 = map_chr(
      index_fwdi7_revi5,
      ~list_files(
        path_ = fastq_dir,
        pattern_ = paste0(.x, "_1_")
      )
    ),
    old_R2 = map_chr(
      index_fwdi7_revi5,
      ~list_files(
        path_ = fastq_dir,
        pattern_ = paste0(.x, "_2_")
      )
    ),
    new_R1 = paste0(sub_library_id, "_R1.fastq.gz"),
    new_R2 = paste0(sub_library_id, "_R2.fastq.gz")
  )

# Check that we found any files before proceeding
if (any(is.na(samples_03$old_R1))) {
  message(
    "Error: Couldn't find any R1 files with provided information. Printing ",
    "offending samples and quitting."
  )
  print(filter(samples_03, is.na(old_R1)))
  q(save = "no")
}


# Write the rename commands -----------------------------------------------

# Rename with "mv -v" so we can see the output, and save to a log file
pwalk(select(samples_03, old_R1, new_R1), ~system(paste0(
  "mv -v ", fastq_dir, ..1, " ", fastq_dir, ..2
)))

pwalk(select(samples_03, old_R2, new_R2), ~system(paste0(
  "mv -v ", fastq_dir, ..1, " ", fastq_dir, ..2
)))
