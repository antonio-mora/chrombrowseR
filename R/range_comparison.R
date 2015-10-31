range_comparison <- function(df_a, df_b, operation="intersection", info="on") {

	df_a_gr <- with(df_a, GRanges(chr, IRanges(start, end, names=1:dim(df_a)[1])))
	df_b_gr <- with(df_b, GRanges(chr, IRanges(start, end, names=1:dim(df_b)[1])))

	if (operation=="intersection") {
		if (info=="on") { cat("Ranges of the first dataframe that overlap with ranges of the second dataframe:\n") }

		overlaps <- findOverlaps(df_a_gr, df_b_gr)
		match_hit <- data.frame(seqnames(df_a_gr)[queryHits(overlaps)],
		ranges(df_a_gr)[queryHits(overlaps)],
       	        names(df_b_gr)[subjectHits(overlaps)],
       	        stringsAsFactors=F)
		intersect_index <- as.numeric(unique(match_hit[,5]))

		result <- df_a[intersect_index,]
	}

	if (operation=="difference") {
		if (info=="on") { cat("Ranges of the first dataframe that do not overlap with ranges of the second dataframe:\n") }

		overlaps <- findOverlaps(df_a_gr, df_b_gr)
		match_hit <- data.frame(seqnames(df_a_gr)[queryHits(overlaps)],
		ranges(df_a_gr)[queryHits(overlaps)],
       	        names(df_b_gr)[subjectHits(overlaps)],
       	        stringsAsFactors=F)
		intersect_index <- as.numeric(unique(match_hit[,5]))

		result <- df_a[-intersect_index,]
	}

	if (operation=="union") {
		if (info=="on") { cat("Merged ranges:\n") }

		merged_gr <- reduce(c(df_a_gr, df_b_gr))
		merged_df <- data.frame(seqnames(merged_gr), ranges(merged_gr))
		rownames(merged_df) <- NULL
		colnames(merged_df) <- c("chr", "start", "end")

		result <- merged_df[,c(1,2,3)]
	}

#	if (operation=="plot") {
#		draw_venn <- function(bed1, bed2, name1, name2) {
#			library(VennDiagram)
#			bed1_GR <- GRanges( seqnames=bed1$chr, range=IRanges(start=bed1$start, end=bed1$end, names=paste(bed1$chr,bed1$start,sep=":")),strand="*" )
#			bed2_GR <- GRanges( seqnames=bed2$chr, range=IRanges(start=bed2$start, end=bed2$end, names=paste(bed2$chr,bed2$start,sep=":")),strand="*" )
#			inters_indexes <- bed1_GR %over% bed2_GR
#			inters <- length(bed1_GR[inters_indexes])
#			draw.pairwise.venn( length(bed1_GR), length(bed2_GR), inters, category=c(name1,name2))
#		}
#	}

	result
}

