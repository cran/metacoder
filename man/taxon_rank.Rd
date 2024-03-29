% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/old_taxa--taxon_rank.R
\name{taxon_rank}
\alias{taxon_rank}
\title{Taxon rank class}
\usage{
taxon_rank(name, database = NULL)
}
\arguments{
\item{name}{(character) rank name. required}

\item{database}{(character) database class object, optional}
}
\value{
An `R6Class` object of class `TaxonRank`
}
\description{
Stores the rank of a taxon. This is typically used to store where taxon
information came from in [taxon()] objects.
}
\examples{
taxon_rank("species")
taxon_rank("genus")
taxon_rank("kingdom")

(x <- taxon_rank(
  "species",
  database_list$ncbi
))
x$rank
x$database

# a null taxon_name object
taxon_name(NULL)
}
\seealso{
Other classes: 
\code{\link{hierarchies}()},
\code{\link{hierarchy}()},
\code{\link{taxa}()},
\code{\link{taxmap}()},
\code{\link{taxon}()},
\code{\link{taxon_database}()},
\code{\link{taxon_id}()},
\code{\link{taxon_name}()},
\code{\link{taxonomy}()}
}
\concept{classes}
