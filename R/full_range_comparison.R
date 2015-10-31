full_range_comparison <- function(df_a, df_b) {
	cat("Intersection, differences and union of two dataframes of genomic ranges:\n")

	result <- list("a - b"=range_comparison(df_a, df_b, "difference", info="off"),
			"a&b"=range_comparison(df_a, df_b, "intersection", info="off"), 
			"b&a"=range_comparison(df_b, df_a, "intersection", info="off"), 
			"b - a"=range_comparison(df_b, df_a, "difference", info="off"), 
			"a U b"=range_comparison(df_a, df_b, "union", info="off"))
	result
}

