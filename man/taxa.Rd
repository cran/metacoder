% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/old_taxa--taxa.R
\name{taxa}
\alias{taxa}
\title{A class for multiple taxon objects}
\usage{
taxa(..., .list = NULL)
}
\arguments{
\item{...}{Any number of object of class [taxon()]}

\item{.list}{An alternate to the `...` input. Any number of object of class
[taxon()]. Cannot be used with `...`.}
}
\value{
An `R6Class` object of class `Taxon`
}
\description{
Stores one or more [taxon()] objects. This is just a thin wrapper for a list
of [taxon()] objects.
}
\details{
This is the documentation for the class called `taxa`. If you are looking for
the documentation for the package as a whole: [taxa-package].
}
\examples{
(a <- taxon(
  name = taxon_name("Poa annua"),
  rank = taxon_rank("species"),
  id = taxon_id(93036)
))
taxa(a, a, a)

# a null set
x <- taxon(NULL)
taxa(x, x, x)

# combo non-null and null
taxa(a, x, a)
}
\seealso{
Other classes: 
\code{\link{hierarchies}()},
\code{\link{hierarchy}()},
\code{\link{taxmap}()},
\code{\link{taxon}()},
\code{\link{taxon_database}()},
\code{\link{taxon_id}()},
\code{\link{taxon_name}()},
\code{\link{taxon_rank}()},
\code{\link{taxonomy}()}
}
\concept{classes}
