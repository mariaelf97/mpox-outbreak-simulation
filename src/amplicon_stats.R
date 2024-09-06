library(data.table)

# Function to read amplicon statistics files for a given isolate
read_amp_stat_files <- function(isolate){
  # Remove any unwanted suffixes (e.g., ".1")
  isolate_clean <- sub("\\.\\d+$", "", isolate)
  
  # Construct file path
  file_path <- paste0("mnt/mpox_seqs/simulations/pure_samples/",isolate,"/", isolate_clean, "/amplicons/amplicon_stats.csv")
  
  # Read the file
  df <- fread(file_path)
  
  # Add the isolate name to the dataframe
  df$isolate <- isolate_clean
  
  return(df)
}

lineages_B2 <- fread("mnt/mpox_seqs/lineages_B1.tsv",header = FALSE)
# Read the list of isolates
isolates <- fread("mnt/mpox_seqs/isolates_B1.txt", header = FALSE)

# Apply the function to each isolate and combine results into a single data frame
df_all <- rbindlist(mapply(read_amp_stat_files, isolates$V1, SIMPLIFY = FALSE))

# df_all now contains the combined data

df_all %>%
  filter(amplicon_length == 0) %>%
  count(amplicon_number, name = "zero_length_count")%>% View()
