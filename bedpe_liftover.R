bedpe_liftover <- function(bedpe_file, from="mm8", to="mm9") {
	library(rtracklayer)

	# Get the liftOver chain from UCSC:
	to <- paste(toupper(substring(to, 1,1)), substring(to, 2), sep="")
	conversion <- paste("/", from, "To", to, ".over.chain", sep="")

	chain_file <- paste(getwd(), conversion, sep="")
	chain_zip <- paste(getwd(), conversion, ".gz", sep="")
	url <- paste("http://hgdownload.cse.ucsc.edu/gbdb/", from, "/liftOver", conversion, ".gz", sep="")

	if (file.exists(chain_file) == TRUE) {
		cat("Reading available chain file...\n")
	} else {
		cat("Downloading liftOver chain to working directory...\n")
		download.file(url, destfile=chain_zip)
		gunzip(chain_zip, exdir=getwd())
	}

	# Use rtracklayer to practice liftOver to both ends of each pair:
	chain = import.chain(paste(from, "To", to, ".over.chain", sep=""))

	df1 <- data.frame(chr=bedpe_file[,1], start=bedpe_file[,2], end=bedpe_file[,3], strand=bedpe_file[,9], name=bedpe_file[,7], score=bedpe_file[,8], samplenumber=bedpe_file[,11], ids=1:dim(bedpe_file)[1]) 
	df2 <- data.frame(chr=bedpe_file[,4], start=bedpe_file[,5], end=bedpe_file[,6], strand=bedpe_file[,10], name=bedpe_file[,7], score=bedpe_file[,8], samplenumber=bedpe_file[,11], ids=1:dim(bedpe_file)[1]) 

	gr1 <- makeGRangesFromDataFrame(df1, TRUE)
	new_df1 <- as.data.frame(unlist(liftOver(gr1, chain)))
	gr2 <- makeGRangesFromDataFrame(df2, TRUE)
	new_df2 <- as.data.frame(unlist(liftOver(gr2, chain)))

	# Deal with multiplied rows (single rows in the original file that can be located in more than one place in the resulting file):
	if (length(which(duplicated(new_df1$ids)==T)) > 0) {
		new_df1 <- new_df1[-which(duplicated(new_df1$ids)==T),] # remove only the second one of a duplicated pair
	}
	if (length(which(duplicated(new_df2$ids)==T)) > 0) {
		new_df2 <- new_df2[-which(duplicated(new_df2$ids)==T),]
	}

	# Combine results into new bedpe:
	new_df_file <- merge(new_df1, new_df2, by="ids", all=TRUE)

	# Deal with lost regions:
	completely_lost_regions <- setdiff(df1$ids, new_df_file$ids)	# ids that are already removed for not being in any of the sets
	a <- is.na(new_df_file$seqnames.x)
	b <- is.na(new_df_file$seqnames.y)
	partially_lost_regions <- which(a == TRUE | b == TRUE)	# count rows with NA in either x or y

	clean_new_df_file <- new_df_file[-partially_lost_regions, c(2,3,4,10,11,12,7,8,6,14,9)]	# remove NA rows, and organize bedpe
	colnames(clean_new_df_file) <- c("chr1", "start1", "end1", "chr2", "start2", "end2", "name", "score", "strand1", "strand2", "samplenumber")

	result <- list(df=clean_new_df_file, index_lost_regions=c(completely_lost_regions, partially_lost_regions))
}
