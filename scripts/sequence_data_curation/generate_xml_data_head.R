# this script read in partitioned sequence alignments and
# the the associated table that contains the assigned geographic area and the age (in unit of number of days) of each sequence
# to generate head of a BEAST XML script that contains the taxa chunk
library(ape)

submission_dateby <- "092220"
collection_dateby <- "030820"

####################
# data preparation #
####################
# read in sequences
# for applying this script, modify the path here
outputdir_path <- paste0("./data/sequence_data/for_analysis/", submission_dateby, "_", collection_dateby)
# to use other partition schemes, change the subdirectory path here
outputdir_path_current <- paste0(outputdir_path, "/orf1restcp")
alignment_paths <- list.files(outputdir_path_current, full.names = T, pattern = ".fasta")
alignment_all <- lapply(alignment_paths, function(x) read.dna(file = x, format = "fasta", as.character = T))

# sanity check
if (!(is.list(alignment_all) && is.matrix(alignment_all[[1]]))) {
  if (is.matrix(alignment_all)) {
    alignment_all <- list(alignment_all)
  } else {
    stop("something_wrong")
  }
}

# fetch the sampling date to generate the age of each sequence
sequence_name <- rownames(alignment_all[[1]])
dates <- as.Date(sapply(strsplit(sequence_name, "_"), "[[", 3), format = "%m%d%y")
dates_backwards <- as.integer(max(dates) - dates) + 1

x <- c("<?xml version=\"1.0\" standalone=\"yes\"?>", "<beast>", "</beast>")
x <- append(x, c("\t<taxa id=\"taxa\">", "\t</taxa>\n"), after = 2)

# read in the sampling geographic araa information
discrete_trait_path <- paste0("./data/sequence_data/for_analysis/", submission_dateby, "_", collection_dateby, "/geography_fine.txt")
states_dat <- read.table(discrete_trait_path, header = T, sep = ",", stringsAsFactors = F)

# assume there is only one discrete trait for now, could be easily relaxed
# also assumes the name of that given discrete trait is either geography or host; could be easily relaxed as well
discrete_trait_name <- colnames(states_dat)[colnames(states_dat) == "geography" | colnames(states_dat) == "host"]
colnames(states_dat)[colnames(states_dat) == "geography" | colnames(states_dat) == "host"] <- "trait"
colnames(states_dat) <- gsub("taxon_id", "taxon", colnames(states_dat))

states <- sort(as.vector(unique(states_dat$trait)))
if (any(states == "?")) {
  states <- states[-which(states == "?")]
}
states_num <- length(states)

##################
# XML generation #
##################
# insert geographic-area data
for (i in 1:length(sequence_name)) {
  taxon <- paste0("\t\t<taxon id=\"", sequence_name[i], "\">\n", 
                  "\t\t\t<date value=\"", dates_backwards[i], "\" direction=\"backwards\" units=\"days\"/>\n")
  taxon <- paste0(taxon, "\t\t\t<attr name=\"", discrete_trait_name, "\">", states_dat$trait[i], "</attr>\n")
  x <- append(x, paste0(taxon, "\t\t</taxon>"), after = length(x) - 2)
}

# instert list of geographic areas
x <- append(x, c(paste0("\t<generalDataType id=\"", discrete_trait_name, ".dataType\">"), "\t</generalDataType>"), after = length(x) - 1)
for (i in 1:states_num) {
  state_line <- paste0("\t\t<state code=\"", states[i], "\"/>")
  x <- append(x, state_line, after = length(x) - 2)
}

# insert the sequence alignment by looping over the partitions
for (i in 1:length(alignment_all)) {
  x <- append(x, c(paste0("\n\t<alignment id=\"alignment", i, "\" dataType=\"nucleotide\">"), "\t</alignment>\n"), after = length(x) - 1)
  
  for (j in 1:length(sequence_name)) {
    taxon_stateseq <- paste0("\t\t<sequence>\n",
                             "\t\t\t<taxon idref=\"", sequence_name[j], "\"/>\t\n",
                             "\t\t\t", toupper(paste(alignment_all[[i]][j, ], collapse = "")), "\n",
                             "\t\t</sequence>")
    x <- append(x, taxon_stateseq, after = length(x) - 2)
  }
}

# write the XML data head to file
cat(x, file = paste0(outputdir_path_current, "/posterior_data.xml"), sep = "\n")
