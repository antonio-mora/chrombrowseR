##############################
# Format Rao Hi-C data
# AUTHOR: Antonio Mora
# Created: Apr.2016
# Last update: 06.05.2016
# DEPENDS ON: get_rao()
##############################
format_rao <- function(rao, file_type="loops", cell_type="GM12878") {
	source("bedpe_liftover.R")

	# Order each cell-type file:
	this_cell <- rao[[file_type]]
	if (is.null(this_cell) == FALSE) {
		rao_order <- this_cell[with(this_cell, order(chr1, x1)), ]

		# Build bedpe format:
		empty_columns <- rep(NA, dim(rao_order)[1])
		fixed_chrom1 <- do.call(paste, list(rep("chr",dim(rao_order)[1]), rao_order[,1], sep=""))
		fixed_chrom2 <- do.call(paste, list(rep("chr",dim(rao_order)[1]), rao_order[,4], sep=""))
		if (file_type=="loops") {
			rao_pref <- paste("Rao-loops", cell_type, sep="_")
			rao_names <- do.call(paste, list(rep(rao_pref,dim(rao_order)[1]), 1:dim(rao_order)[1], sep="_"))
			rao_bedpe <- data.frame(fixed_chrom1, rao_order[,2:3], fixed_chrom2, rao_order[,5:6], rao_names, empty_columns, empty_columns, empty_columns, empty_columns)
			colnames(rao_bedpe) <- c("chr1", "start1", "end1", "chr2", "start2", "end2", "name", "score", "strand1", "strand2", "samplenumber")
		}
		if (file_type=="domains") {
			rao_pref <- paste("Rao-domains", cell_type, sep="_")
			rao_names <- do.call(paste, list(rep(rao_pref,dim(rao_order)[1]), 1:dim(rao_order)[1], sep="_"))
			rao_bedpe <- data.frame(fixed_chrom1, rao_order[,2:3], rao_names, empty_columns, empty_columns)
			colnames(rao_bedpe) <- c("chr", "start", "end", "name", "score", "samplenumber")
		}
		if (file_type=="loops_motifs") {
			rao_pref <- paste("Rao-loops-motifs", cell_type, sep="_")
			rao_names <- do.call(paste, list(rep(rao_pref,dim(rao_order)[1]), 1:dim(rao_order)[1], sep="_"))
			rao_bedpe <- data.frame(fixed_chrom1, rao_order[,2:3], fixed_chrom2, rao_order[,5:6], rao_names, empty_columns, empty_columns, empty_columns, empty_columns, rao_order[,21:24], rao_order[,26:29])
			colnames(rao_bedpe) <- c("chr1", "start1", "end1", "chr2", "start2", "end2", "name", "score", "strand1", "strand2", "samplenumber", "start_motif1", "end_motif1", "sequence_motif1", "orientation_motif1", "start_motif2", "end_motif2", "sequence_motif2", "orientation_motif2")
		}
		rownames(rao_bedpe) <- NULL

		# Fix naming issue with CH12 cells:
		if (cell_type=="CH12-LX") {
			rao_bedpe[,1] <- gsub("chrchr", "chr", rao_bedpe[,1], fixed=T)
			rao_bedpe[,4] <- gsub("chrchr", "chr", rao_bedpe[,4], fixed=T)
		}

		# Assign distances as score for loops, and size as score for domains:
		score <- NULL
		if (file_type=="domains") {
			score <- abs(rao_bedpe$end - rao_bedpe$start)
		} else {
			for (j in 1:dim(rao_bedpe)[1]) {
				if (rao_bedpe[j,1]==rao_bedpe[j,4]) {
					midp_1 <- (rao_bedpe[j,2]+rao_bedpe[j,3])/2
					midp_2 <- (rao_bedpe[j,5]+rao_bedpe[j,6])/2
					dist <- abs(midp_1 - midp_2)
				} else {
					dist <- 3036303846	# for inter-chromosomal, assign arbitrary high number, here the size of all hg19
				}
				score <- c(score, dist)
			}
		}
		rao_bedpe$score <- score

		# Perform liftover to translate from hg19 to hg38:
		#tmp <- bedpe_liftover(rao_bedpe, from="hg19", to="hg38")
		#rao_bedpe <- tmp[[1]]

	} else {
		rao_bedpe <- NULL
	}

	rao_bedpe
}
###############
# Note: Rao doesn't have strand info
###############
