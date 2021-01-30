# this script reads in a flitered alignment, 
# renames them to make their names shorter and human readable
# and then partition them in various forms (e.g., by gene and/or by codon positions)
library(ape)

#####################
# read in alignment #
#####################
submission_dateby <- "092220"
collection_dateby <- "030820"
reference_genome_name_string <- "EPI_ISL_402125"

if (submission_dateby == "041920") {
  alignment_path <- paste0("./data/sequence_data/alignments/not_partitioned/ncov_", submission_dateby, "_", collection_dateby, "_aligned_codon.fasta")
} else {
  alignment_path <- paste0("./data/gisaid_directly/", submission_dateby, "_", collection_dateby, "_aligned_codon.fasta")
}
alignment_codon <- read.dna(alignment_path, format = "fasta", as.character = T)
sequence_accession <- as.character(sapply(rownames(alignment_codon), function(x) unlist(strsplit(x, "\\|"))[2]))

#####################################
# information from reference genome #
#####################################
reference_inalignment_codon <- alignment_codon[grep(reference_genome_name_string, rownames(alignment_codon)), ]
reference_genome_annotation_aa <- matrix(c(1, 180, 181, 818, 819, 2763, 2764, 3263, 3264, 3569, 3570, 3859, 
                                           3860, 3942, 3943, 4140, 4141, 4253, 4254, 4392,
                                           4393, 5324, 5325, 5925, 5926, 6452, 6453, 6798, 6799, 7096, 
                                           7097, 7781, 7782, 8369, 8370, 8644, 8645, 8769, 8770, 8941, 8942, 9002, 
                                           9003, 9123, 9124, 9166, 9167, 9287, 9288, 9706, 9707, 9744), ncol = 2, byrow = T)
gene_names <- c("nsp1", "nsp2", "nsp3", "nsp4", "nsp5", "nsp6", "nsp7", "nsp8", "nsp9", "nsp10", 
                "nsp12", "nsp13", "nsp14", "nsp15", "nsp16", "S1", "S2", "orf3a", "E", "M", 
                "orf6", "orf7a", "orf7b", "orf8", "N", "orf10")
rownames(reference_genome_annotation_aa) <- gene_names
colnames(reference_genome_annotation_aa) <- c("start", "end")

reference_genome_annotation_codon <- reference_genome_annotation_aa
reference_genome_annotation_codon[, 1] <- (reference_genome_annotation_codon[, 1] - 1) * 3 + 1 
reference_genome_annotation_codon[, 2] <- reference_genome_annotation_codon[, 2] * 3

##########################
# modify alignment names #
##########################
# read in meta data
# for applying this script, modify the path here
if (submission_dateby == "041920") {
  metadata <- read.csv(paste0("./data/meta_data/metainfo_", submission_dateby, ".csv"), stringsAsFactors = F, check.names = F)
} else {
  metadata <- read.table(paste0("./data/gisaid_directly/metadata_", submission_dateby, ".tsv"), header = T, sep = "\t", 
                         quote = "\"", stringsAsFactors = F, check.names = F, comment.char = "")
  colnames(metadata)[colnames(metadata) == "gisaid_epi_isl"] <- "accession_id"
  colnames(metadata)[colnames(metadata) == "date"] <- "collection_date"
}

metadata <- metadata[metadata$accession_id %in% sequence_accession, ]
if (any(!(sequence_accession %in% metadata$accession_id))) {
  stop("cannot find some sequences that are supposed to be in the metadata table")
}
metadata <- metadata[match(sequence_accession, metadata$accession_id), ]

collection_date <- as.Date(metadata$collection_date, format = "%Y-%m-%d")
if (any(is.na(collection_date))) {
  stop("cannot find date")
}

alignment_codon <- alignment_codon[collection_date <= as.Date(collection_dateby, format = "%m%d%y", tz = "UTC"), ]
sequence_accession <- as.character(sapply(rownames(alignment_codon), function(x) unlist(strsplit(x, "\\|"))[2]))

# match sequences in the alignment with the metadata table again after excluding sequences that were sampled after our specified boundary
metadata <- metadata[metadata$accession_id %in% sequence_accession, ]
metadata <- metadata[match(sequence_accession, metadata$accession_id), ]
collection_date <- as.Date(metadata$collection_date, format = "%Y-%m-%d", tz = "UTC")

# each sequence name comprises three parts: 6-digit GISAID sequence accession,
# sampling country, and sampling date
sequence_name_firststring <- sapply(strsplit(sequence_accession, "_"), "[[", 3)

# modify some country names to make them shorter
if ("country" %in% colnames(metadata)) {
  sequence_name_secondstring <- metadata$country
} else {
  sequence_name_secondstring <- sapply(strsplit(metadata$location, " / |/| /|/ "), "[[", 2)
}
sequence_name_secondstring[sequence_name_secondstring == "Korea"] <- "South Korea"
sequence_name_secondstring[sequence_name_secondstring == "Democratic Republic of the Congo"] <- "DRC"
sequence_name_secondstring[sequence_name_secondstring == "Czech Republic"] <- "Czech"
sequence_name_secondstring[sequence_name_secondstring == "England" | 
                             sequence_name_secondstring == "Scotland" | 
                             sequence_name_secondstring == "United Kingdom"] <- "UK"
sequence_name_secondstring[sequence_name_secondstring == "United Arab Emirates"] <- "UAE"
sequence_name_secondstring[sequence_name_secondstring == "Saint Barthélemy"] <- "StBarts"

sequence_name_thirdstring <- format(collection_date, "%m%d%y")

sequence_name <- character(length(sequence_name_firststring))
for (i in 1:length(sequence_name)) {
  sequence_name[i] <- paste(c(sequence_name_firststring[i], sequence_name_secondstring[i], sequence_name_thirdstring[i]), collapse = "_")
}
sequence_name <- gsub(" ", "", sequence_name)
rownames(alignment_codon) <- sequence_name

###########################################
# partition alignments and write them out #
###########################################
# for applying this script, modify the path here
outputdir_path <- paste0("./data/sequence_data/for_analysis/", submission_dateby, "_", collection_dateby)
if (!dir.exists(outputdir_path)) dir.create(outputdir_path)

# full alignment
# for applying this script, modify the path here
write.dna(alignment_codon, 
          file = paste0(outputdir_path, "/hcov19_", submission_dateby, "_", collection_dateby, ".fasta"), 
          format = "fasta", nbcol = -1, colsep = "", colw = 3.2e4)

# alignments partitioned by codon positions
outputdir_path_current <- paste0(outputdir_path, "/cp")
if (!dir.exists(outputdir_path_current)) dir.create(paste0(outputdir_path_current))
num_codon <- ncol(alignment_codon) / 3
# for applying this script, modify the path here
for (j in 1:3) {
  write.dna(alignment_codon[, (1:num_codon - 1) * 3 + j], 
            file = paste0(outputdir_path_current, "/hcov19_cp", j, "_", submission_dateby, "_", collection_dateby, ".fasta"), 
            format = "fasta", nbcol = -1, colsep = "", colw = 3.2e4)
}

# alignments partitioned by gene and by codon position
outputdir_path_current <- paste0(outputdir_path, "/genecp")
if (!dir.exists(outputdir_path_current)) dir.create(outputdir_path_current)
for (i in 1:nrow(reference_genome_annotation_codon)) {
  alignment_codon_current <- alignment_codon[, reference_genome_annotation_codon[i, 1]:reference_genome_annotation_codon[i, 2]]
  num_codon <- ncol(alignment_codon_current) / 3
  # for applying this script, modify the path here
  for (j in 1:3) {
    write.dna(alignment_codon_current[, (1:num_codon - 1) * 3 + j], 
              file = paste0(outputdir_path_current, "/hcov19_", rownames(reference_genome_annotation_codon)[i], "_cp", j, "_", submission_dateby, "_", collection_dateby, ".fasta"), 
              format = "fasta", nbcol = -1, colsep = "", colw = 3.2e4)
  }
}

# alignments partitioned by orf1/remaining and by codon position
outputdir_path_current <- paste0(outputdir_path, "/orf1restcp")
if (!dir.exists(outputdir_path_current)) dir.create(paste0(outputdir_path_current))
gene_regions <- list(orf1 = grep("nsp", rownames(reference_genome_annotation_codon), value = T),
                     SEMNorf3to10 = grep("nsp", rownames(reference_genome_annotation_codon), invert = T, value = T))

for (i in 1:length(gene_regions)) {
  reference_genome_annotation_codon_current <- reference_genome_annotation_codon[rownames(reference_genome_annotation_codon) %in% gene_regions[[i]], ]
  if (!is.matrix(reference_genome_annotation_codon_current)) {
    reference_genome_annotation_codon_current <- matrix(reference_genome_annotation_codon_current, ncol = 2, byrow = T)
    rownames(reference_genome_annotation_codon_current) <- gene_regions[[i]]
    colnames(reference_genome_annotation_codon_current) <- colnames(reference_genome_annotation_codon)
  }
  
  nsites <- sum(apply(reference_genome_annotation_codon_current, 1, function(x) diff(x) + 1))
  alignment_codon_current <- matrix(ncol = nsites, nrow = nrow(alignment_codon), byrow = T)
  
  k <- 0
  for (j in 1:nrow(reference_genome_annotation_codon_current)) {
    alignment_codon_current[, (0:diff(reference_genome_annotation_codon_current[j, ]) + 1) + k] <- alignment_codon[, reference_genome_annotation_codon_current[j, 1]:reference_genome_annotation_codon_current[j, 2]]
    k <- k + diff(reference_genome_annotation_codon_current[j, ]) + 1
  }
  rownames(alignment_codon_current) <- rownames(alignment_codon)
  
  num_codon <- ncol(alignment_codon_current) / 3
  # for applying this script, modify the path here
  for (j in 1:3) {
    write.dna(alignment_codon_current[, (1:num_codon - 1) * 3 + j], 
              file = paste0(outputdir_path_current, "/hcov19_", names(gene_regions)[i], "_cp", j, "_", submission_dateby, "_", collection_dateby, ".fasta"), 
              format = "fasta", nbcol = -1, colsep = "", colw = 3.2e4)
  }
}

# alignments partitioned by orf1a, orf1b, S1, S2, orf3a, EM, orf6to8, Norf10
outputdir_path_current <- paste0(outputdir_path, "/mainregions1cp")
if (!dir.exists(outputdir_path_current)) dir.create(paste0(outputdir_path_current))
gene_regions <- list(orf1a = paste0("nsp", 1:10), orf1b = paste0("nsp", 12:16), S1 = c("S1"), S2 = c("S2"), orf3a = c("orf3a"),
                     EM = c("E", "M"), orf6to8 = c("orf6", "orf7a", "orf7b", "orf8"), Norf10 = c("N", "orf10"))

for (i in 1:length(gene_regions)) {
  
  reference_genome_annotation_codon_current <- reference_genome_annotation_codon[rownames(reference_genome_annotation_codon) %in% gene_regions[[i]], ]
  if (!is.matrix(reference_genome_annotation_codon_current)) {
    reference_genome_annotation_codon_current <- matrix(reference_genome_annotation_codon_current, ncol = 2, byrow = T)
    rownames(reference_genome_annotation_codon_current) <- gene_regions[[i]]
    colnames(reference_genome_annotation_codon_current) <- colnames(reference_genome_annotation_codon)
  }
  
  nsites <- sum(apply(reference_genome_annotation_codon_current, 1, function(x) diff(x) + 1))
  alignment_codon_current <- matrix(ncol = nsites, nrow = nrow(alignment_codon), byrow = T)
  
  k <- 0
  for (j in 1:nrow(reference_genome_annotation_codon_current)) {
    alignment_codon_current[, (0:diff(reference_genome_annotation_codon_current[j, ]) + 1) + k] <- alignment_codon[, reference_genome_annotation_codon_current[j, 1]:reference_genome_annotation_codon_current[j, 2]]
    k <- k + diff(reference_genome_annotation_codon_current[j, ]) + 1
  }
  rownames(alignment_codon_current) <- rownames(alignment_codon)
  
  num_codon <- ncol(alignment_codon_current) / 3
  # for applying this script, modify the path here
  for (j in 1:3) {
    write.dna(alignment_codon_current[, (1:num_codon - 1) * 3 + j], 
              file = paste0(outputdir_path_current, "/hcov19_", names(gene_regions)[i], "_cp", j, "_", submission_dateby, "_", collection_dateby, ".fasta"), 
              format = "fasta", nbcol = -1, colsep = "", colw = 3.2e4)
  }
}
