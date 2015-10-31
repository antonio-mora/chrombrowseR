\name{range_annotator}
\alias{range_annotator}
\title{Annotation of a dataframe of genomic ranges}
\usage{
range_annotator(df_base, df_anno)
}
\arguments{
  \item{df_base}{A dataframe of a bed file of interest.}

  \item{df_anno}{A dataframe of anotated regions. For example, gene info.}
}
\value{
A dataframe with the bed file of interest plus an additional column corresponding to the annotation.
}
\description{
This function takes a dataframe of genomic ranges and a dataframe with annotated genomic regions, and assigns annotations to the original dataframe based on the overlap to annotated regions in the annotated dataframe.
}
\examples{
df_1 <- data.frame( chr=rep("chr1",4), start=c(1,16,23,45), end=c(9,18,36,49), name=paste("region_",1:4,sep=""), score=rep(1,4), strand=rep("+",4) )
df_2 <- data.frame( chr=rep("chr1",4), start=c(4,19,23,35), end=c(9,21,26,42), annotation=c("Gene1", "Gene2", "Gene3", "Gene4") )

range_annotator(df_1, df_2)
}
