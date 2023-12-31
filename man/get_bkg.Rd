% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_bkg.R
\name{get_bkg}
\alias{get_bkg}
\title{Calculate sequence background.}
\usage{
get_bkg(sequences, k = 1:3, as.prob = NULL, pseudocount = 0,
  alphabet = NULL, to.meme = NULL, RC = FALSE, list.out = NULL,
  nthreads = 1, merge.res = TRUE, window = FALSE, window.size = 0.1,
  window.overlap = 0)
}
\arguments{
\item{sequences}{\code{\link{XStringSet}} Input sequences. Note that if
multiple sequences are present, the results will be combined into one
(unless \code{merge.res = FALSE}).}

\item{k}{\code{integer} Size of k-let. Background can be calculated for any
k-let size.}

\item{as.prob}{Deprecated.}

\item{pseudocount}{\code{integer(1)} Add a count to each possible k-let. Prevents
any k-let from having 0 or 1 probabilities.}

\item{alphabet}{\code{character(1)} Provide a custom alphabet to calculate a
background for. If \code{NULL}, then standard letters will be assumed for
DNA, RNA and AA sequences, and all unique letters found will be used
for \code{BStringSet} type sequences. Note that letters which are not a part
of the standard DNA/RNA/AA alphabets or in the provided alphabet will
not be counted in the totals during probability calculations.}

\item{to.meme}{If not \code{NULL}, then \code{\link[=get_bkg]{get_bkg()}} will return the sequence
background in MEME Markov Background Model format. Input for this argument
will be used for \code{cat(..., file = to.meme)} within \code{\link[=get_bkg]{get_bkg()}}. See
\url{http://meme-suite.org/doc/bfile-format.html} for a description of
the format.}

\item{RC}{\code{logical(1)} Calculate the background of the reverse complement
of the input sequences as well. Only valid for DNA/RNA.}

\item{list.out}{Deprecated.}

\item{nthreads}{\code{numeric(1)} Run \code{\link[=get_bkg]{get_bkg()}} in parallel with \code{nthreads}
threads. \code{nthreads = 0} uses all available threads.
Note that no speed up will occur for jobs with only a single sequence.}

\item{merge.res}{\code{logical(1)} Whether to merge results from all sequences
or return background data for individual sequences.}

\item{window}{\code{logical(1)} Determine background in windows.}

\item{window.size}{\code{numeric} Window size. If a number between 0 and 1 is
provided, the value is calculated as the number multiplied by the sequence
length.}

\item{window.overlap}{\code{numeric} Overlap between windows. If a number
between 0 and 1 is provided, the value is calculated as the number
multiplied by the sequence length.}
}
\value{
If \code{to.meme = NULL}, a \code{DataFrame} with columns \code{klet}, \code{count},
and \code{probability}. If \code{merge.res = FALSE}, there will be an additional
\code{sequence} column. If \code{window = TRUE}, there will be an additional \code{start}
and \code{stop} columns.

If \code{to.meme} is not \code{NULL}, then \code{NULL} is returned, invisibly.
}
\description{
For a set of input sequences, calculate the overall sequence background for
any k-let size. For very large sequences DNA and RNA sequences (in the billions of bases),
please be aware of the much faster and more efficient
\code{\link[Biostrings:nucleotideFrequency]{Biostrings::oligonucleotideFrequency()}}.
\code{\link[=get_bkg]{get_bkg()}} can still be used in these cases, though it may take several seconds or
minutes to calculate the results (depending on requested k-let sizes).
}
\examples{
## Compare to Biostrings version
library(Biostrings)
seqs.DNA <- create_sequences()
bkg.DNA <- get_bkg(seqs.DNA, k = 3)
bkg.DNA2 <- oligonucleotideFrequency(seqs.DNA, 3, 1, as.prob = FALSE)
bkg.DNA2 <- colSums(bkg.DNA2)
all(bkg.DNA$count == bkg.DNA2)

## Create a MEME background file
get_bkg(seqs.DNA, k = 1:3, to.meme = stdout(), pseudocount = 1)

## Non-DNA/RNA/AA alphabets
seqs.QWERTY <- create_sequences("QWERTY")
bkg.QWERTY <- get_bkg(seqs.QWERTY, k = 1:2)

}
\references{
Bailey TL, Elkan C (1994). “Fitting a mixture model by expectation
maximization to discover motifs in biopolymers.” \emph{Proceedings of
the Second International Conference on Intelligent Systems for
Molecular Biology}, \strong{2}, 28-36.
}
\seealso{
\code{\link[=create_sequences]{create_sequences()}}, \code{\link[=scan_sequences]{scan_sequences()}}, \code{\link[=shuffle_sequences]{shuffle_sequences()}}
}
\author{
Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@uwaterloo.ca}
}
