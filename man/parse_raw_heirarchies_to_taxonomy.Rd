% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/old_taxa--taxonomy--internal.R
\name{parse_raw_heirarchies_to_taxonomy}
\alias{parse_raw_heirarchies_to_taxonomy}
\title{Infer edge list from hierarchies composed of character vectors}
\usage{
parse_raw_heirarchies_to_taxonomy(heirarchies, named_by_rank = FALSE)
}
\arguments{
\item{named_by_rank}{(`TRUE`/`FALSE`) If  `TRUE` and the input is a list of
vectors with each vector named by ranks, include that rank info in the
output object, so it can be accessed by `out$taxon_ranks()`. If `TRUE`,
taxa with different ranks, but the same name and location in the taxonomy,
will be considered different taxa.}
}
\value{
A list of character vectors.
}
\description{
Infer edge list and unique taxa from hierarchies.
}
\keyword{internal}
