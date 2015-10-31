full_range_comparison <- function(df_a, df_b) {
	cat("Intersection, differences and union of two dataframes of genomic ranges:")

	result <- list("a - b"=range_comparison(df_a, df_b, "difference"), "a&b"=range_comparison(df_a, df_b, "intersection"), "b&a"=range_comparison(df_b, df_a, "intersection"), "b - a"=range_comparison(df_b, df_a, "difference"), "a U b"=range_comparison(df_a, df_b, "union"))
}

