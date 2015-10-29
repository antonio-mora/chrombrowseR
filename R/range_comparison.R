range_comparison <- function(df_base, df_subt, operation="intersection") {
	df_base_gr <- with(df_base, GRanges(chr, IRanges(start, end, names=1:dim(df_base)[1])))
	df_subt_gr <- with(df_subt, GRanges(chr, IRanges(start, end, names=1:dim(df_subt)[1])))

	overlaps <- findOverlaps(df_base_gr, df_subt_gr)
	match_hit <- data.frame(seqnames(df_base_gr)[queryHits(overlaps)],
		ranges(df_base_gr)[queryHits(overlaps)],
                names(df_subt_gr)[subjectHits(overlaps)],
                stringsAsFactors=F)

	intersect_index <- as.numeric(unique(match_hit[,5]))
	if (operation=="intersection") {
		result <- df_base[intersect_index,]
	} else {
		result <- df_base[-intersect_index,]
	}
	result
}

