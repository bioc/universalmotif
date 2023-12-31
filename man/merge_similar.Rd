% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/merge_similar.R
\name{merge_similar}
\alias{merge_similar}
\title{Identify and merge similar motifs within a collection of motifs (or simply
cluster motifs).}
\usage{
merge_similar(motifs, threshold = 0.95, threshold.type = "score.abs",
  method = "PCC", use.type = "PPM", min.overlap = 6, min.mean.ic = 0,
  tryRC = TRUE, relative_entropy = FALSE, normalise.scores = FALSE,
  min.position.ic = 0, score.strat.compare = "a.mean",
  score.strat.merge = "sum", nthreads = 1, return.clusters = FALSE)
}
\arguments{
\item{motifs}{See \code{\link[=convert_motifs]{convert_motifs()}} for acceptable motif formats.}

\item{threshold}{\code{numeric(1)} The minimum (for similarity metrics) or maximum (for
distance metrics) threshold score for merging.}

\item{threshold.type}{\code{character(1)} Type of score used for thresholding.
Currently unused.}

\item{method}{\code{character(1)} One of PCC, EUCL, SW, KL, BHAT, HELL,
SEUCL, MAN, WEUCL, WPCC. See \code{\link[=compare_motifs]{compare_motifs()}}. (The ALLR and ALLR_LL
methods cannot be used for distance matrix construction.)}

\item{use.type}{\code{character(1)} One of \code{'PPM'} and \code{'ICM'}.
The latter allows for taking into account the background
frequencies if \code{relative_entropy = TRUE}. Note that \code{'ICM'} is not
allowed when \code{method = c("ALLR", "ALLR_LL")}.}

\item{min.overlap}{\code{numeric(1)} Minimum overlap required when aligning the
motifs. Setting this to a number higher then the width of the motifs
will not allow any overhangs. Can also be a number between 0 and 1,
representing the minimum fraction that the motifs must overlap.}

\item{min.mean.ic}{\code{numeric(1)} Minimum mean information content between the
two motifs for an alignment to be scored. This helps prevent scoring
alignments between low information content regions of two motifs. Note that
this can result in some comparisons failing if no alignment passes the
mean IC threshold. Use \code{\link[=average_ic]{average_ic()}} to filter out low IC motifs to get around
this if you want to avoid getting \code{NA}s in your output.}

\item{tryRC}{\code{logical(1)} Try the reverse complement of the motifs as well,
report the best score.}

\item{relative_entropy}{\code{logical(1)} Change the ICM calculation affecting
\code{min.position.ic} and \code{min.mean.ic}. See \code{\link[=convert_type]{convert_type()}}.}

\item{normalise.scores}{\code{logical(1)} Favour alignments which leave fewer
unaligned positions, as well as alignments between motifs of similar length.
Similarity scores are multiplied by the ratio of
aligned positions to the total number of positions in the larger motif,
and the inverse for distance scores.}

\item{min.position.ic}{\code{numeric(1)} Minimum information content required between
individual alignment positions for it to be counted in the final alignment
score. It is recommended to use this together with \code{normalise.scores = TRUE},
as this will help punish scores resulting from only a fraction of an
alignment.}

\item{score.strat.compare}{\code{character(1)} The \code{score.strat} parameter used
by \code{\link[=compare_motifs]{compare_motifs()}}. For clustering purposes, the \code{"sum"} option cannot
be used.}

\item{score.strat.merge}{\code{character(1)} The \code{score.strat} parameter used
by \code{\link[=merge_motifs]{merge_motifs()}}. As discussed in \code{\link[=merge_motifs]{merge_motifs()}}, the \code{"sum"} option
is recommended over \code{"a.mean"} to maximize the overlap between motifs.}

\item{nthreads}{\code{numeric(1)} Run \code{\link[=compare_motifs]{compare_motifs()}} in parallel with \code{nthreads}
threads. \code{nthreads = 0} uses all available threads.}

\item{return.clusters}{\code{logical(1)} Return the clusters instead of
merging.}
}
\value{
See \code{\link[=convert_motifs]{convert_motifs()}} for available output formats.
}
\description{
Given a list of motifs, \code{\link[=merge_similar]{merge_similar()}} will identify similar motifs with
\code{\link[=compare_motifs]{compare_motifs()}}, and merge similar ones with \code{\link[=merge_motifs]{merge_motifs()}}.
}
\details{
See \code{\link[=compare_motifs]{compare_motifs()}} for more info on comparison parameters, and
\code{\link[=merge_motifs]{merge_motifs()}} for more info on motif merging.
}
\examples{
\dontrun{
library(MotifDb)
motifs <- filter_motifs(MotifDb, family = "bHLH")[1:50]
length(motifs)
motifs <- merge_similar(motifs)
length(motifs)
}

}
\seealso{
\code{\link[=compare_motifs]{compare_motifs()}}, \code{\link[=merge_motifs]{merge_motifs()}}
}
\author{
Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@uwaterloo.ca}
}
