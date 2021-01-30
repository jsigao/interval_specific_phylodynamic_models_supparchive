# this script reads in a FASTA file (containing SARS-CoV-2 sequences) downloaded from GISAID
# and then perform the first round of filteration to get ready for alignment
library(ape)
submission_dateby <- "041920"
collection_dateby <- "030820"
reference_genome_name_string <- "EPI_ISL_402125"
nuc_unambig <- c("a", "c", "g", "t")

#####################
# read in alignment #
#####################
# for applying this script, modify the sequence data path here
sequences <- read.dna(paste0("./data/sequence_data/sequences_raw/combined_seqs/sequences_raw_", submission_dateby, ".fasta"), 
                      format = "fasta", as.character = T)
sequence_accession <- as.character(sapply(names(sequences), function(x) unlist(strsplit(x, "\\|"))[2]))

#####################
# read in meta data #
#####################
# for applying this script, modify the metadata path here
metadata <- read.csv(paste0("./data/meta_data/metainfo_", submission_dateby, ".csv"), stringsAsFactors = F, check.names = F)
nometainfo_exclude <- names(sequences)[!(sequence_accession %in% metadata$accession_id)]

# exclude sequences without a specific date
metadata <- metadata[metadata$accession_id %in% sequence_accession, ]
collection_date <- as.Date(metadata$collection_date, format = "%Y-%m-%d", tz = "UTC")
nospecificdate_ordatetoorecent_exclude <- names(sequences)[match(metadata$accession_id[is.na(collection_date) | 
                                                                                         collection_date > as.Date(collection_dateby, format = "%m%d%y", tz = "UTC")], 
                                                                 sequence_accession)]

#####################################
# information from reference genome #
#####################################
reference_genome <- sequences[[grep(reference_genome_name_string, names(sequences))]]
reference_genome_length <- length(reference_genome)
if (reference_genome_length != "29903") stop("something wrong")

reference_genome_53primeUTRtrimmed <- reference_genome[266:29674]
reference_genome_53primeUTRtrimmed_length <- length(reference_genome_53primeUTRtrimmed)

##############################
# first round of fliteration #
##############################
# exclude nonhuman sequences
nonhuman_exclude <- grep("/pangolin/|/bat/|/tiger/|/cat/|/canine/|IVDC-HB-env", names(sequences), value = T)

# sequence with too much missing and/or too short
sequence_nomissing_length <- sapply(sequences, function(x) length(x) - sum(x %in% c("-", "?", "n")))
toomuchmissing_andor_tooshort_exclude <- names(sequences)[sequence_nomissing_length < 29000]

# duplicate sequences
# for applying this script, modify the metadata path here (downloaded from https://github.com/nextstrain/ncov/blob/master/defaults/exclude.txt)
duplicate_exclude_string <- scan("./data/meta_data/exclude_seqs.txt", what = character(), sep = "\n", comment.char = "#")
duplicate_exclude <- unlist(sapply(duplicate_exclude_string, function(x) grep(x, names(sequences), value = T)))

# combining all the exclusion criteria
sequence_name_exclude <- unique(c(nonhuman_exclude, toomuchmissing_andor_tooshort_exclude, duplicate_exclude,
                                  nospecificdate_ordatetoorecent_exclude, nometainfo_exclude))
sequence_processed <- sequences[!(names(sequences) %in% sequence_name_exclude)]

# exclude sequences with either too much missing, too short, or too different from the reference genome
# first round is to quickly exclude sequences that can be split by the start and end strings of the reference genome (and too much missing or too short)
reference_genome_53primeUTRtrimmed_startstring <- paste(reference_genome_53primeUTRtrimmed[1:15], collapse = "")
reference_genome_53primeUTRtrimmed_endstring <- paste(reference_genome_53primeUTRtrimmed[(reference_genome_53primeUTRtrimmed_length - 14):reference_genome_53primeUTRtrimmed_length], collapse = "")

sequence_string <- sapply(sequence_processed, function(x) paste(x, collapse = ""))
sequence_53primeUTRtrimmed_stringsplit <- lapply(sequence_string, function(x) unlist(strsplit(x, paste(c(reference_genome_53primeUTRtrimmed_startstring, reference_genome_53primeUTRtrimmed_endstring), collapse = "|"))))

sequence_53primeUTRtrimmed_string <- sapply(sequence_53primeUTRtrimmed_stringsplit[lengths(sequence_53primeUTRtrimmed_stringsplit) == 3], "[[", 2)
sequence_53primeUTRtrimmed <- lapply(sequence_53primeUTRtrimmed_string, function(x) unlist(strsplit(x, "")))

sequence_53primeUTRtrimmed_nomissing_length <- sapply(sequence_53primeUTRtrimmed, function(x) sum(x %in% nuc_unambig))
toomuch_missing_exclude <- names(sequence_53primeUTRtrimmed_nomissing_length)[sequence_53primeUTRtrimmed_nomissing_length < 29000]

sequence_processed <- sequence_processed[!(names(sequence_processed) %in% toomuch_missing_exclude)]

# write processed sequences out to file
# for applying this script, modify the output path here
sequence_combined_path <- paste0("./data/sequence_data/sequences_raw/combined_seqs/sequences_quickqualitycontroled_", submission_dateby, "_", collection_dateby, ".fasta")
write.dna(sequence_processed, file = sequence_combined_path, format = "fasta", nbcol = -1, colsep = "", colw = 3.2e4)
