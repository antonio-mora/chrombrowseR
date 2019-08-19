################
# Get Rao Hi-C data
# AUTHOR: Antonio Mora
# Created: Mar.2016
# Last update: 28.04.2016
# PAPER: http://www.cell.com/cell/fulltext/S0092-8674%2814%2901497-4
# GEO: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525
# METHOD: In-situ Hi-C
# CELL TYPE: GM12878, HMEC, HUVEC, HeLa, IMR90, K562, KBM7, NHEK, CH12-LX 
# GENOME ASSEMBLY: b37(hg19) for human cells, and mm9 for mouse cells
# YEAR: 2014
################
get_rao <- function(cell="GM12878") {
	source("get_tables.R")

	file_location1 <- paste(getwd(), "/", cell, "_domainlist.txt.gz", sep="")
	file_location2 <- paste(getwd(), "/", cell, "_looplist.txt.gz", sep="")
	file_location3 <- paste(getwd(), "/", cell, "_looplist_with_motifs.txt.gz", sep="")

	txt_file1 <- paste(getwd(), "/", cell, "_domainlist.txt", sep="")
	txt_file2 <- paste(getwd(), "/", cell, "_looplist.txt", sep="")
	txt_file3 <- paste(getwd(), "/", cell, "_looplist_with_motifs.txt", sep="")

	if (cell == "GM12878") {
		url1 <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt.gz"
		url2 <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_GM12878_primary+replicate_HiCCUPS_looplist.txt.gz"
		url3 <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_GM12878_primary+replicate_HiCCUPS_looplist_with_motifs.txt.gz"
	} else {
		url1 <- paste("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_", cell, "_Arrowhead_domainlist.txt.gz", sep="")
		url2 <- paste("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_", cell, "_HiCCUPS_looplist.txt.gz", sep="")
		url3 <- paste("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_", cell, "_HiCCUPS_looplist_with_motifs.txt.gz", sep="")
	}

	if (cell == "KBM7") {
		domains <- get_tables(file_location1, txt_file1, url1)
		loops <- get_tables(file_location2, txt_file2, url2)
		result <- list(domains=domains, loops=loops)
	} else {
		domains <- get_tables(file_location1, txt_file1, url1)
		loops <- get_tables(file_location2, txt_file2, url2)
		loops_motifs <- get_tables(file_location3, txt_file3, url3)
		result <- list(domains=domains, loops=loops, loops_motifs=loops_motifs)
	}

	result
}
