% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parsers.R
\name{parse_silva_fasta}
\alias{parse_silva_fasta}
\title{Parse SILVA FASTA release}
\usage{
parse_silva_fasta(file = NULL, input = NULL, include_seqs = TRUE)
}
\arguments{
\item{file}{The path to a FASTA file containing sequences to use. Either
"input" or "file" must be supplied but not both.}

\item{input}{(\code{character}) One of the following: 
\describe{
  \item{A character vector of sequences}{See the example below for what this
  looks like. The parser \code{\link{read_fasta}} produces output like this.}
  \item{A list of character vectors}{Each vector should have one base per element.}
  \item{A "DNAbin" object}{This is the result of parsers like
  \code{\link[ape]{read.FASTA}}.}
  \item{A list of "SeqFastadna" objects}{This is the result of parsers like
  \code{\link[seqinr]{read.fasta}}.}
  Either "input" or "file" must be supplied but not both.
}}

\item{include_seqs}{(\code{logical} of length 1) If \code{TRUE}, include
sequences in the output object.}
}
\value{
\code{\link{taxmap}}
}
\description{
Parses an SILVA FASTA file that can be found at
\url{https://www.arb-silva.de/no_cache/download/archive/release_128/Exports/}.
}
\details{
The input file has a format like:

\preformatted{ >GCVF01000431.1.2369
Bacteria;Proteobacteria;Gammaproteobacteria;Oceanospiril...
CGUGCACGGUGGAUGCCUUGGCAGCCAGAGGCGAUGAAGGACGUUGUAGCCUGCGAUAAGCUCCGGUUAGGUGGCAAACA
ACCGUUUGACCCGGAGAUCUCCGAAUGGGGCAACCCACCCGUUGUAAGGCGGGUAUCACCGACUGAAUCCAUAGGUCGGU
... }
}
\seealso{
Other parsers: 
\code{\link{extract_tax_data}()},
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
\code{\link{parse_tax_data}()},
\code{\link{parse_ubiome}()},
\code{\link{parse_unite_general}()}
}
\concept{parsers}
