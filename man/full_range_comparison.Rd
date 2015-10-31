\name{full_range_comparison}
\alias{full_range_comparison}
\title{Intersection, differences and union of two dataframes of genomic ranges}
\usage{
full_range_comparison(df_a, df_b)
}
\arguments{
  \item{df_a}{A dataframe of a bed file.}

  \item{df_b}{A dataframe of a bed file.}
}
\value{
A list of dataframes corresponding to: The difference \code{df_a} - \code{df_b}, the intersection \code{df_a} & \code{df_b}, the intersection \code{df_b} & \code{df_a}, the difference \code{df_b} - \code{df_a}, and the union \code{df_a} U \code{df_b}.
}
\description{
This function takes two dataframes of genomic ranges and perform all intersections, union and differences between them. The intersection is defined as those regions in the first dataframe that overlap at least 1 bp with regions in the second dataframe. Internally, this function uses the \code{GenomicRanges} package.
}
\examples{
df_1 <- data.frame( chr=rep("chr1",4), start=c(1,16,23,45), end=c(9,18,36,49), name=paste("region_",1:4,sep=""), score=rep(1,4), strand=rep("+",4) )
df_2 <- data.frame( chr=rep("chr1",4), start=c(4,19,28,49), end=c(9,21,39,55), name=paste("region_",1:4,sep=""), score=rep(1,4), strand=rep("+",4) )

full_range_comparison(df_1, df_2)
}
