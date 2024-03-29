% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metacoder-package.R
\name{metacoder}
\alias{metacoder}
\title{Metacoder}
\description{
A package for planning and analysis of amplicon metagenomics research projects.
}
\details{
The goal of the \code{metacoder} package is to provide a set of tools for:

\itemize{
  \item Standardized parsing of taxonomic information from diverse resources.
  \item Visualization of statistics distributed over taxonomic classifications.
  \item Evaluating potential metabarcoding primers for taxonomic specificity.
  \item Providing flexible functions for analyzing taxonomic and abundance data.
}

To accomplish these goals, \code{metacoder} leverages resources from other R packages, interfaces with
external programs, and provides novel functions where needed to allow for entire analyses within R.
}
\section{Documentation}{


The full documentation can be found online at \url{https://grunwaldlab.github.io/metacoder_documentation/}.

There is also a short vignette included for offline use that can be accessed  by the following code:

\code{browseVignettes(package = "metacoder")}

\strong{Plotting:}

\itemize{
  \item \code{\link{heat_tree}}
  \item \code{\link{heat_tree_matrix}}
}

\strong{In silico PCR:}

\itemize{
  \item \code{\link{primersearch}}
}

\strong{Analysis:}

\itemize{
  \item \code{\link{calc_taxon_abund}}
  \item \code{\link{calc_obs_props}}
  \item \code{\link{rarefy_obs}}
  \item \code{\link{compare_groups}}
  \item \code{\link{zero_low_counts}}
  \item \code{\link{calc_n_samples}}
  \item \code{\link{filter_ambiguous_taxa}}
}

\strong{Parsers:}

\itemize{
  \item \code{\link{parse_greengenes}}
  \item \code{\link{parse_mothur_tax_summary}}
  \item \code{\link{parse_mothur_taxonomy}}
  \item \code{\link{parse_newick}}
  \item \code{\link{parse_phyloseq}}
  \item \code{\link{parse_phylo}}
  \item \code{\link{parse_qiime_biom}}
  \item \code{\link{parse_rdp}}
  \item \code{\link{parse_silva_fasta}}
  \item \code{\link{parse_unite_general}}
}

\strong{Writers:}

\itemize{
  \item \code{\link{write_greengenes}}
  \item \code{\link{write_mothur_taxonomy}}
  \item \code{\link{write_rdp}}
  \item \code{\link{write_silva_fasta}}
  \item \code{\link{write_unite_general}}
}

\strong{Database querying:}

\itemize{
  \item \code{\link{ncbi_taxon_sample}}
}
}

\section{Main classes}{


  These are the classes users would typically interact with:

  * [taxon]: A class used to define a single taxon. Many other classes in the
  `taxa`` package include one or more objects of this class.
  * : Stores one or more [taxon] objects. This is just a thin wrapper
  for a list of [taxon] objects.
  * [hierarchy]: A class containing an ordered list of [taxon] objects that
  represent a hierarchical classification.
  * [hierarchies]: A list of taxonomic classifications.  This is just a thin wrapper
  for a list of [hierarchy] objects.
  * [taxonomy]: A taxonomy composed of [taxon] objects organized in a tree
  structure. This differs from the [hierarchies] class in how the [taxon]
  objects are stored. Unlike a [hierarchies] object, each unique taxon is
  stored only once and the relationships between taxa are stored in an
  edgelist.
  * [taxmap]: A class designed to store a taxonomy and associated
  user-defined data. This class builds on the [taxonomy] class. User defined
  data can be stored in the list `obj$data`, where `obj` is a taxmap
  object. Any number of user-defined lists, vectors, or tables mapped
  to taxa can be manipulated in a cohesive way such that relationships
  between taxa and data are preserved.
}

\section{Minor classes}{


  These classes are mostly components for the larger classes above and would
  not typically be used on their own.

  * [taxon_database]: Used to store information about taxonomy databases.
  * [taxon_id]: Used to store taxon IDs, either arbitrary or from a
  particular taxonomy database.
  * [taxon_name]: Used to store taxon names, either arbitrary or from a
  particular taxonomy database.
  * [taxon_rank]: Used to store taxon ranks (e.g. species, family), either
  arbitrary or from a particular taxonomy database.
}

\section{Major manipulation functions}{


 These are some of the more important functions used to filter data in classes
 that store multiple taxa, like [hierarchies], [taxmap], and [taxonomy].

  * [filter_taxa]: Filter taxa in a [taxonomy] or [taxmap] object with a
  series of conditions. Relationships between remaining taxa and user-defined
  data are preserved (There are many options controlling this).
  * [filter_obs]: Filter user-defined data [taxmap] object with a series of
  conditions. Relationships between remaining taxa and user-defined data are
  preserved (There are many options controlling this);
  * [sample_n_taxa]: Randomly sample taxa. Has same abilities as
  [filter_taxa].
  * [sample_n_obs]: Randomly sample observations. Has same abilities as
  [filter_obs].
  * [mutate_obs]: Add datasets or columns to datasets in the `data` list of
  [taxmap] objects.
  * [pick]: Pick out specific taxa, while others are dropped in [hierarchy]
  and [hierarchies] objects.
  * [pop]: Pop out taxa (drop them) in [hierarchy] and [hierarchies] objects.
  * [span]: Select a range of taxa, either by two names, or relational
  operators in [hierarchy] and [hierarchies] objects.
}

\section{Mapping functions}{


  There are lots of functions for getting information for each taxon.

  * [subtaxa]: Return data for the subtaxa of each taxon in an [taxonomy] or
  [taxmap] object.
  * [supertaxa]: Return data for the supertaxa of each taxon in an [taxonomy]
  or [taxmap] object.
  * [roots]: Return data for the roots of each taxon in an [taxonomy] or
  [taxmap] object.
  * [leaves]: Return data for the leaves of each taxon in an [taxonomy] or
  [taxmap] object.
  * [obs]: Return user-specific data for each taxon and all of its subtaxa in
  an [taxonomy] or [taxmap] object.
}

\section{The kind of classes used}{


 Note, this is mostly of interest to developers and advanced users.

 The classes in the `taxa` package are mostly
 [R6](https://adv-r.hadley.nz/r6.html) classes ([R6Class]). A few of the
 simpler ones ( and [hierarchies]) are
 [S3](https://adv-r.hadley.nz/s3.html) instead. R6 classes are different than
 most R objects because they are
 [mutable](https://en.wikipedia.org/wiki/Immutable_object) (e.g. A function
 can change its input without returning it). In this, they are more similar
 to class systems in
 [object-oriented](https://en.wikipedia.org/wiki/Object-oriented_programming)
 languages like python. As in other object-oriented class systems, functions
 are thought to "belong" to classes (i.e. the data), rather than functions
 existing independently of the data. For example, the function `print` in R
 exists apart from what it is printing, although it will change how it prints
 based on what the class of the data is that is passed to it. In fact, a user
 can make a custom print method for their own class by defining a function
 called `print.myclassname`. In contrast, the functions that operate on R6
 functions are "packaged" with the data they operate on. For example, a print
 method of an object for an R6 class might be called like
 `my_data$print()` instead of `print(my_data)`.
}

\section{The two ways to call functions}{


 Note, you will need to read the previous section to fully understand this one.

 Since the R6 function syntax (e.g. `my_data$print()`) might be confusing to
 many R users, all functions in `taxa` also have S3 versions. For example,
 the [filter_taxa()] function can be called on a [taxmap] object called
 `my_obj` like `my_obj$filter_taxa(...)` (the R6 syntax) or
 `filter_taxa(my_obj, ...)` (the S3 syntax). For some functions, these two
 way of calling the function can have different effect. For functions that do
 not returned a modified version of the input (e.g. [subtaxa()]), the two ways have identical behavior.
 However, functions like [filter_taxa()], that modify their inputs, actually
 change the object passed to them as the first argument as well as returning that
 object. For example,

 `my_obj <- filter_taxa(my_obj, ...)`

 and

 `my_obj$filter_taxa(...)`

 and

 `new_obj <- my_obj$filter_taxa(...)`

 all replace `my_obj` with the filtered result, but

 `new_obj <- filter_taxa(my_obj, ...)`

 will not modify `my_obj`.
}

\section{Non-standard evaluation}{


  This is a rather advanced topic.

  Like packages such as `ggplot2` and [dplyr], the `taxa` package uses
  non-standard evaluation to allow code
  to be more readable and shorter. In effect, there are variables that only
  "exist" inside a function call and depend on what is passed to that function
  as the first parameter (usually a class object). For example, in the `dpylr`
  function [filter()], column names can be used as if they were independent
  variables. See `?dpylr::filter` for examples of this. The `taxa` package builds on this idea.

  For many functions that work on [taxonomy] or [taxmap] objects (e.g. [filter_taxa]),
  some functions that return per-taxon information (e.g. [taxon_names()]) can
  be referred to by just the name of the function. When one of these functions
  are referred to by name, the function is run on the relevant object and its
  value replaces the function name. For example,

  `new_obj <- filter_taxa(my_obj, taxon_names == "Bacteria")`

   is identical to:

  `new_obj <- filter_taxa(my_obj, taxon_names(my_obj) == "Bacteria")`

   which is identical to:

  `new_obj <- filter_taxa(my_obj, my_obj$taxon_names() == "Bacteria")`

   which is identical to:

  `my_names <- taxon_names(my_obj)`

  `new_obj <- filter_taxa(my_obj, my_names == "Bacteria")`

  For `taxmap` objects, you can also use names of user defined lists, vectors,
  and the names of columns in user-defined tables that are stored in the
  `obj$data` list. See [filter_taxa()] for examples. You can even add your own
  functions that are called by name by adding them to the `obj$funcs` list.
  For any object with functions that use non-standard evaluation, you can see
  what values can be used  with [all_names()] like `all_names(obj)`.
}

\section{Dependencies and inspiration}{


  Various elements of the `taxa` package were inspired by the [dplyr] and
  [taxize] packages. This package started as parts of the `metacoder` and
  `binomen` packages. There are also many dependencies that make `taxa`
  possible.
}

\section{Feedback and contributions}{


  Find a problem? Have a suggestion? Have a question? Please submit an issue
  at our [GitHub repository](https://github.com/ropensci/taxa):

  [https://github.com/ropensci/taxa/issues](https://github.com/ropensci/taxa/issues)

  A GitHub account is free and easy to set up. We welcome feedback! If you
  don't want to use GitHub for some reason, feel free to email us. We do
  prefer posting to github since it allows others that might have the same
  issue to see our conversation. It also helps us keep track of what problems
  we need to address.

  Want to contribute code or make a change to the code? Great, thank you!
  Please [fork](https://help.github.com/articles/fork-a-repo/) our GitHub
  repository and submit a [pull request](https://help.github.com/articles/about-pull-requests/).
}

\author{
Zachary Foster and Niklaus Grunwald
}
