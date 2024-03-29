% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/old_taxa--taxmap--parsers.R
\name{extract_tax_data}
\alias{extract_tax_data}
\title{Extracts taxonomy info from vectors with regex}
\usage{
extract_tax_data(
  tax_data,
  key,
  regex,
  class_key = "taxon_name",
  class_regex = "(.*)",
  class_sep = NULL,
  sep_is_regex = FALSE,
  class_rev = FALSE,
  database = "ncbi",
  include_match = FALSE,
  include_tax_data = TRUE
)
}
\arguments{
\item{tax_data}{A vector from which to extract taxonomy information.}

\item{key}{(`character`) The identity of the capturing groups defined using
`regex`. The length of `key` must be equal to the number of capturing
groups specified in `regex`. Any names added to the terms will be used as
column names in the output. Only `"info"` can be used multiple times. Each
term must be one of those described below:
* `taxon_id`: A unique numeric id for a taxon for a particular `database`
(e.g. ncbi accession number). Requires an internet connection.
* `taxon_name`: The name of a taxon (e.g. "Mammalia" or "Homo sapiens").
Not necessarily unique, but interpretable by a particular `database`.
Requires an internet connection.
* `fuzzy_name`: The name of a taxon, but check for misspellings first.
Only use if you think there are misspellings. Using `"taxon_name"` is
faster.
* `class`: A list of taxon information that constitutes the full taxonomic
classification (e.g. "K_Mammalia;P_Carnivora;C_Felidae"). Individual
taxa are separated by the `class_sep` argument and the information is
parsed by the `class_regex` and `class_key` arguments.
* `seq_id`: Sequence ID for a particular database that is associated with a
taxonomic classification. Currently only works with the "ncbi" database.
* `info`: Arbitrary taxon info you want included in the output. Can be used
more than once.}

\item{regex}{(`character` of length 1) A regular expression with capturing
groups indicating the locations of relevant information. The identity of
the information must be specified using the `key` argument.}

\item{class_key}{(`character` of length 1) The identity of the capturing
groups defined using `class_regex`. The length of `class_key` must be equal
to the number of capturing groups specified in `class_regex`. Any names
added to the terms will be used as column names in the output. Only
`"info"` can be used multiple times. Each term must be one of those
described below:
* `taxon_name`: The name of a taxon. Not necessarily unique.
* `taxon_rank`: The rank of the taxon. This will be used to add rank info
into the output object that can be accessed by `out$taxon_ranks()`.
* `info`: Arbitrary taxon info you want included in the output. Can be used
more than once.}

\item{class_regex}{(`character` of length 1)
A regular expression with capturing groups indicating the locations of data
for each taxon in the `class` term in the `key` argument. The identity of
the information must be specified using the `class_key` argument. The
`class_sep` option can be used to split the classification into data for
each taxon before matching. If `class_sep` is `NULL`, each match of
`class_regex` defines a taxon in the classification.}

\item{class_sep}{(`character` of length 1)
Used with the `class` term in the `key` argument. The character(s) used to
separate individual taxa within a classification. After the string defined
by the `class` capture group in `regex` is split by `class_sep`, its
capture groups are extracted by `class_regex` and defined by `class_key`.
If `NULL`, every match of `class_regex` is used instead with first
splitting by `class_sep`.}

\item{sep_is_regex}{(`TRUE`/`FALSE`) Whether or not `class_sep` should be
used as a [regular expression](https://en.wikipedia.org/wiki/Regular_expression).}

\item{class_rev}{(`logical` of length 1)
Used with the `class` term in the `key` argument. If `TRUE`, the order of
taxon data in a classification is reversed to be specific to broad.}

\item{database}{(`character` of length 1) The name of the database that
patterns given in `parser` will apply to. Valid databases include "ncbi",
"itis", "eol", "col", "tropicos", "nbn", and "none". `"none"` will cause no
database to be queried; use this if you want to not use the internet. NOTE:
Only `"ncbi"` has been tested extensively so far.}

\item{include_match}{(`logical` of length 1) If `TRUE`, include the part of
the input matched by `regex` in the output object.}

\item{include_tax_data}{(`TRUE`/`FALSE`) Whether or not to include `tax_data`
as a dataset.}
}
\value{
Returns an object of type [taxmap()]
}
\description{
Convert taxonomic information in a character vector into a [taxmap()] object.
The location and identity of important information in the input is specified
using a [regular expression](https://en.wikipedia.org/wiki/Regular_expression)
with capture groups and a corresponding key. An object of type [taxmap()] is
returned containing the specified information. See the `key` option for
accepted sources of taxonomic information.
}
\section{Failed Downloads}{
 If you have invalid inputs or a download fails for
  another reason, then there will be a "unknown" taxon ID as a placeholder
  and failed inputs will be assigned to this ID. You can remove these using
  [filter_taxa()] like so: `filter_taxa(result, taxon_ids != "unknown")`. Add
  `drop_obs = FALSE` if you want the input data, but want to remove the
  taxon.
}

\examples{

\dontrun{

  # For demonstration purposes, the following example dataset has all the
  # types of data that can be used, but any one of them alone would work.
  raw_data <- c(
  ">id:AB548412-tid:9689-Panthera leo-tax:K_Mammalia;P_Carnivora;C_Felidae;G_Panthera;S_leo",
  ">id:FJ358423-tid:9694-Panthera tigris-tax:K_Mammalia;P_Carnivora;C_Felidae;G_Panthera;S_tigris",
  ">id:DQ334818-tid:9643-Ursus americanus-tax:K_Mammalia;P_Carnivora;C_Felidae;G_Ursus;S_americanus"
  )

  # Build a taxmap object from classifications
  extract_tax_data(raw_data,
                   key = c(my_seq = "info", my_tid = "info", org = "info", tax = "class"),
                   regex = "^>id:(.+)-tid:(.+)-(.+)-tax:(.+)$",
                   class_sep = ";", class_regex = "^(.+)_(.+)$",
                   class_key = c(my_rank = "info", tax_name = "taxon_name"))

  # Build a taxmap object from taxon ids
  # Note: this requires an internet connection
  extract_tax_data(raw_data,
                   key = c(my_seq = "info", my_tid = "taxon_id", org = "info", tax = "info"),
                   regex = "^>id:(.+)-tid:(.+)-(.+)-tax:(.+)$")

  # Build a taxmap object from ncbi sequence accession numbers
  # Note: this requires an internet connection
  extract_tax_data(raw_data,
                   key = c(my_seq = "seq_id", my_tid = "info", org = "info", tax = "info"),
                   regex = "^>id:(.+)-tid:(.+)-(.+)-tax:(.+)$")

  # Build a taxmap object from taxon names
  # Note: this requires an internet connection
  extract_tax_data(raw_data,
                   key = c(my_seq = "info", my_tid = "info", org = "taxon_name", tax = "info"),
                   regex = "^>id:(.+)-tid:(.+)-(.+)-tax:(.+)$")
}
}
\seealso{
Other parsers: 
\code{\link{lookup_tax_data}()},
\code{\link{parse_dada2}()},
\code{\link{parse_edge_list}()},
\code{\link{parse_greengenes}()},
\code{\link{parse_mothur_tax_summary}()},
\code{\link{parse_mothur_taxonomy}()},
\code{\link{parse_newick}()},
\code{\link{parse_phylo}()},
\code{\link{parse_phyloseq}()},
\code{\link{parse_qiime_biom}()},
\code{\link{parse_rdp}()},
\code{\link{parse_silva_fasta}()},
\code{\link{parse_tax_data}()},
\code{\link{parse_ubiome}()},
\code{\link{parse_unite_general}()}
}
\concept{parsers}
