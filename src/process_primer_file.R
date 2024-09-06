# Load required libraries
library(dplyr)
library(tidyr)
library(Biostrings)
library(data.table)

# Read in the primers data
primers <- fread("Downloads/primers_output.csv")

# Select the necessary columns
primer_bed <- primers %>%
  select(Forward_Primer, Reverse_Primer, Forward_Start, Forward_End, Reverse_Start, Reverse_End, Amplicon)
# Function to find primer positions in the sequence
find_primer_positions <- function(sequence, primer_seq) {
  match <- matchPattern(primer_seq, sequence)
  
  if (length(match) > 0) {
    start_pos <- start(match)
    end_pos <- end(match)
    return(list(start = start_pos, end = end_pos))
  } else {
    return(list(start = NA, end = NA))
  }
}

# Read the fasta file (replace 'your_fasta_file.fasta' with the actual path)
fasta_file <- "mnt/mpox_seqs/reference.fasta"
fasta_seq <- readDNAStringSet(fasta_file)

# Assume we are working with the first sequence in the fasta file
sequence <- fasta_seq[[1]]

# Read the primers data
primers <- fread("Downloads/primers_output.csv")
updated_primers <- primers %>%
  rowwise() %>%
  mutate(
    # Get the positions of the forward primer
    Forward_Positions = list(find_primer_positions(sequence, Forward_Primer)),
    Forward_Start = Forward_Positions$start,
    Forward_End = Forward_Positions$end,
    
    # Get the reverse complement of the reverse primer
    Reverse_Complement = as.character(reverseComplement(DNAString(Reverse_Primer))),
    
    # Get the positions of the reverse primer (searching the reverse complement)
    Reverse_Positions = list(find_primer_positions(sequence, Reverse_Complement)),
    Reverse_Start = Reverse_Positions$start,
    Reverse_End = Reverse_Positions$end
  ) %>%
  # Drop the helper columns
  select(-Forward_Positions, -Reverse_Complement, -Reverse_Positions)


# Collapse data and handle both forward and reverse primers
collapsed_data <- updated_primers %>%
  # Create rows for both forward and reverse primers
  tidyr::pivot_longer(
    cols = c(Forward_Primer, Reverse_Primer),
    names_to = "type",
    values_to = "primer_seq"
  ) %>%
  # Add start and end positions based on the primer type (forward or reverse)
  mutate(
    start = ifelse(type == "Forward_Primer", Forward_Start, Reverse_Start),
    end = ifelse(type == "Forward_Primer", Forward_End, Reverse_End),
    direction = ifelse(type == "Forward_Primer", "Forward", "Reverse")
  ) %>%
  # Select only the necessary columns, including Amplicon
  select(primer_seq, start, end, direction, Amplicon)

# Apply transformations to collapsed data
collapsed_data <- collapsed_data %>%
  mutate(
    pool = 1,
    handedness = ifelse(grepl("Forward", direction), "+", "-"),
    name = paste0("value_", Amplicon),
    # Append "LEFT_" for handedness "+" and "RIGHT_" for handedness "-"
    name = ifelse(handedness == "+", paste0(name, "_LEFT"), paste0(name, "_RIGHT")),
    seq = "MPXV-M5312_HM12_Rivers"
  )

collapsed_data <- collapsed_data %>%
  # Compute start1 and end1
  mutate(start1 = pmin(start, end),
         end1 = pmax(start, end))%>%
  # Select columns, removing original start and end, and renaming start1 and end1
  select(-start, -end, start = start1, end = end1)

# Reorder the columns as requested
collapsed_data <- collapsed_data %>%
  select(seq, start, end, name, pool, handedness, primer_seq)


# Define output file path
output_file_path <- "mnt/mpox_seqs/primer.bed"

# Write the output to a file
write.table(collapsed_data, file = output_file_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


