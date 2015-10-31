range_annotator <- function(df_base, df_anno) {
	df_base_gr <- with(df_base, GRanges(chr, IRanges(start, end, names=1:nrow(df_base))))
	df_anno_gr <- with(df_anno, GRanges(chr, IRanges(start, end, names=annotation)))

	overlaps <- findOverlaps(df_base_gr, df_anno_gr)
	match_hit <- data.frame(seqnames(df_base_gr)[queryHits(overlaps)],
		ranges(df_base_gr)[queryHits(overlaps)],
                names(df_anno_gr)[subjectHits(overlaps)],
                stringsAsFactors=F
                )

	df_base[,"annotation"] <- rep("NA", nrow(df_base))
	tmp <- aggregate(match_hit, by=list(match_hit[,5]), paste)
	df_base[tmp[,1],ncol(df_base)] <- unlist(lapply(tmp[,7], paste, collapse=","))
	df_base
}

