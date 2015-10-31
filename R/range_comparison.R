range_comparison <- function(df_a, df_b, operation="intersection") {

	df_a_gr <- with(df_a, GRanges(chr, IRanges(start, end, names=1:dim(df_a)[1])))
	df_b_gr <- with(df_b, GRanges(chr, IRanges(start, end, names=1:dim(df_b)[1])))

	if (operation=="intersection") {
		cat("Ranges of the first dataframe that overlap with ranges of the second dataframe:")

		overlaps <- findOverlaps(df_a_gr, df_b_gr)
		match_hit <- data.frame(seqnames(df_a_gr)[queryHits(overlaps)],
		ranges(df_a_gr)[queryHits(overlaps)],
       	        names(df_b_gr)[subjectHits(overlaps)],
       	        stringsAsFactors=F)
		intersect_index <- as.numeric(unique(match_hit[,5]))

		result <- df_a[intersect_index,]
	}

	if (operation=="difference") {
		cat("Ranges of the first dataframe that do not overlap with ranges of the second dataframe:")

		overlaps <- findOverlaps(df_a_gr, df_b_gr)
		match_hit <- data.frame(seqnames(df_a_gr)[queryHits(overlaps)],
		ranges(df_a_gr)[queryHits(overlaps)],
       	        names(df_b_gr)[subjectHits(overlaps)],
       	        stringsAsFactors=F)
		intersect_index <- as.numeric(unique(match_hit[,5]))

		result <- df_a[-intersect_index,]
	}

	if (operation=="union") {
		cat("Merged ranges:")

		merged_gr <- reduce(c(df_a_gr, df_b_gr))
		merged_df <- data.frame(seqnames(merged_gr), ranges(merged_gr))
		rownames(merged_df) <- NULL
		colnames(merged_df) <- c("chr", "start", "end")

		result <- merged_df[,c(1,2,3)]
	}

	result
}

