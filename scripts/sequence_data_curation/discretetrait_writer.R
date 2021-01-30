# this script read in a sequence alignment (that will be used in the inference)
# as well as the associated metadata table to generate a table containing
# the assigned geographic area and the age (in unit of number of days) of each sequence
library(ape)
library(stringr)
source('./main_config.R')

submission_dateby <- "092220"
collection_dateby <- "030820"

####################
# data preparation #
####################
# read in sequences
# for applying this script, modify the path here
alignment_codon <- read.dna(paste0("./data/sequence_data/for_analysis/", submission_dateby, "_", collection_dateby, 
                                   "/hcov19_", submission_dateby, "_", collection_dateby, ".fasta"), format = "fasta", as.character = T)
sequence_accession <- paste0("EPI_ISL_", as.character(sapply(rownames(alignment_codon), function(x) unlist(strsplit(x, "_"))[1])))

# read in meta data
# here involves some complication as the two (early and late) datasets use different metadata tables
# so that they require different post-processing here
# for applying this script, modify the path here
if (submission_dateby != "041920") {
  metadata <- read.table(paste0("./data/gisaid_directly/metadata_", submission_dateby, ".tsv"), header = T, sep = "\t", 
                         quote = "\"", stringsAsFactors = F, check.names = F, comment.char = "")
  colnames(metadata)[colnames(metadata) == "gisaid_epi_isl"] <- "accession_id"
  colnames(metadata)[colnames(metadata) == "date"] <- "collection_date"
} else if (submission_dateby == "041920") {
  metadata <- read.table(paste0("./data/meta_data/metainfo_", submission_dateby, ".csv"), header = T, sep = ",", 
                         quote = "\"", stringsAsFactors = F, check.names = F, comment.char = "")
}

metadata <- metadata[metadata$accession_id %in% sequence_accession, ]
if (any(!(sequence_accession %in% metadata$accession_id))) {
  stop("included sequence cannot be found in the metadata table.\n")
}
metadata <- metadata[match(sequence_accession, metadata$accession_id), ]

if (submission_dateby != "041920") {
  metadata <- metadata[, c("strain", "accession_id", "collection_date", "region", "country", "division")]
} else if (submission_dateby == "041920") {
  
  # for applying this script, modify the path here
  metadata_gisaid <- read.table(paste0("./data/gisaid_directly/metadata_092220.tsv"), header = T, sep = "\t", 
                                quote = "\"", stringsAsFactors = F, check.names = F, comment.char = "")
  colnames(metadata_gisaid)[colnames(metadata_gisaid) == "gisaid_epi_isl"] <- "accession_id"
  colnames(metadata_gisaid)[colnames(metadata_gisaid) == "date"] <- "collection_date"
  metadata_gisaid <- metadata_gisaid[, c("strain", "accession_id", "collection_date", "region", "country", "division")]
  
  metadata_gisaid <- metadata_gisaid[metadata_gisaid$accession_id %in% sequence_accession, ]
  metadata_gisaid <- metadata_gisaid[!is.na(as.Date(metadata_gisaid$collection_date, format = "%Y-%m-%d")), ]
  metadata_gisaid <- metadata_gisaid[metadata$collection_date[match(metadata_gisaid$accession_id, metadata$accession_id)] == metadata_gisaid$collection_date, ]
  
  if (any(!sequence_accession %in% metadata_gisaid$accession_id)) {
    
    metadata_amend <- metadata[which(!sequence_accession %in% metadata_gisaid$accession_id), ]
    region_amend <- gsub("^ | $", "", str_match(metadata_amend$location, "^(.*?)/")[, 2])
    metadata_amend$location <- gsub("^ | $", "", gsub("^(.*?)/", "", metadata_amend$location))
    
    country_amend <- sapply(strsplit(metadata_amend$location, " /|/"), "[[", 1)
    metadata_amend$location[!grepl("/", metadata_amend$location)] <- ""
    metadata_amend$location <- gsub("^ | $", "", gsub("^(.*?)/", "", metadata_amend$location))
    
    division_amend <- gsub("^ | $", "", gsub("/(.*?)$", "", metadata_amend$location))
    metadata_amend <- cbind(metadata_amend$virus_name, metadata_amend$accession_id, metadata_amend$collection_date, region_amend, country_amend, division_amend)
    colnames(metadata_amend) <- colnames(metadata_gisaid)
    
    metadata_gisaid <- rbind(metadata_gisaid, metadata_amend)
  }
  
  metadata <- metadata_gisaid
  metadata <- metadata[match(sequence_accession, metadata$accession_id), ]
}
metadata$collection_date <- as.Date(metadata$collection_date, format = "%Y-%m-%d")
metadata[metadata$accession_id == reference_genome_name_string, "division"] <- "Hubei"
metadata$division[which(metadata$division == "Hangzhou")] <- "Zhejiang"

####################
# region partition #
####################
regions_all <- list(regions_coarse, regions_fine)
names(regions_all) <- c("coarse", "fine")
sequence_region_all <- vector("list", length(regions_all))
names(sequence_region_all) <- c("coarse", "fine")

for (k in 1:length(regions_all)) {
  
  regions <- regions_all[[k]]
  sequence_region <- character(nrow(metadata))
  
  for (i in 1:nrow(metadata)) {
    area_id <- metadata$country[i]
    seq_region <- names(regions)[grep(paste0("\"", area_id, "\""), regions)]
    
    if (length(seq_region) == 0) {
      area_id <- paste0(metadata$country[i], "_", metadata$division[i])
      seq_region <- names(regions)[grep(paste0("\"", area_id, "\""), regions)]
      
      if (length(seq_region) == 0) {
        stop(paste(area_id, "is/are not in the region list."))
      } 
    }
    
    if (length(seq_region) == 1) {
      sequence_region[i] <- seq_region
    }
  }
  
  sequence_region_all[[k]] <- sequence_region
}

############################
# write out the info table #
############################
for (k in 1:length(regions_all)) {
  regions <- regions_all[[k]]
  sequence_region <- sequence_region_all[[k]]
  
  sequence_name <- rownames(alignment_codon)
  sequence_region_code <- LETTERS[match(sequence_region, names(regions))]

  dates <- as.Date(sapply(strsplit(sequence_name, "_"), "[[", 3), format = "%m%d%y", tz = "UTC")
  dates_backwards <- as.integer(max(dates) - dates) + 1L

  discretetrait_df <- cbind(sequence_name, sequence_region_code, dates_backwards)
  colnames(discretetrait_df) <- c("taxon_id", "geography", "date")
  # for applying this script, modify the path here
  write.table(discretetrait_df, file = paste0("./data/sequence_data/for_analysis/", submission_dateby, "_", collection_dateby, 
                                              "/geography_", names(regions_all)[k], ".txt"), quote = F, sep = ",", row.names = F)
}
