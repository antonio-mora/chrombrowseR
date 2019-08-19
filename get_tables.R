get_tables <- function(file_location, txt_file, url, header=TRUE, comment.char="#", read=TRUE) {
	if (file.exists(txt_file) == TRUE) {
		cat("Reading available file...\n")
		geo_tab <- read.table(txt_file, header=header, comment.char=comment.char, sep='\t', quote="")
	} else {
		cat("Downloading file...\n")
		download.file(url, destfile=file_location)
		gunzip(file_location)
		cat("File has been saved as:\n")
		cat(paste(txt_file, "\n\n"))
		if (read==TRUE) {
			cat("Reading downloaded file...\n\n")
			geo_tab <- read.table(txt_file, header=header, comment.char=comment.char, sep='\t', quote="")
		} else {
			geo_tab <- NULL
		}
	}
	geo_tab
}

