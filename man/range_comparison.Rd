\name{range_comparison}
\alias{range_comparison}
\title{Intersection, Union and Difference of two data frames of genomic ranges}
\usage{
range_comparison(df_a, df_b, operation, ncols)
}
\arguments{
  \item{df_a}{A dataframe of a bed file.}

  \item{df_b}{A dataframe of a bed file.}

  \item{operation}{A string. Either "intersection", "union" or "difference". Default value: "intersection".}

  \item{ncols}{An integer. Number of columns to display in the output data frame. Currently, only accepting three columns: "chr", "start" and "end" (default value: 3).}
}
\value{
A dataframe with the intersection, union or difference between the input data frames \code{df_a} and \code{df_b}.
}
\description{
This function takes two data frames of genomic ranges and perform intersection, union or difference between them. Internally, this function uses the GenomicRanges package. The only goal of the function is to allow the user to directly work with data frames as input and output, avoiding the conversion to GRanges objects.
}
\examples{
df_1 <- data.frame( chr=rep("chr1",4), start=c(1,16,23,45), end=c(9,18,36,49), name=paste("region_",1:4,sep=""), score=rep(1,4), strand=rep("+",4) )
df_2 <- data.frame( chr=rep("chr1",4), start=c(4,19,28,49), end=c(9,21,39,55), name=paste("region_",1:4,sep=""), score=rep(1,4), strand=rep("+",4) )

range_comparison(df_1, df_2, operation="intersection")
range_comparison(df_1, df_2, operation="union")
range_comparison(df_1, df_2, operation="difference")
}
