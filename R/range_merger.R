range_merger <- function(df_a, df_b) {
	df_a_gr <- with(df_a, GRanges(chr, IRanges(start, end, names=1:dim(df_a)[1])))
	df_b_gr <- with(df_b, GRanges(chr, IRanges(start, end, names=1:dim(df_b)[1])))

	merged_gr <- reduce(c(df_a_gr, df_b_gr))
	merged_df <- data.frame(seqnames(merged_gr), ranges(merged_gr))
	rownames(merged_df) <- NULL
	colnames(merged_df) <- c("chr", "start", "end")

	result <- merged_df[,c(1,2,3)]
}

