% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_transfac.R
\name{write_transfac}
\alias{write_transfac}
\title{Export motifs in TRANSFAC format.}
\usage{
write_transfac(motifs, file, overwrite = FALSE, append = FALSE,
  name.tag = "ID", altname.tag = "NA")
}
\arguments{
\item{motifs}{See \code{\link[=convert_motifs]{convert_motifs()}} for acceptable formats.}

\item{file}{\code{character(1)} File name.}

\item{overwrite}{\code{logical(1)} Overwrite existing file.}

\item{append}{\code{logical(1)} Add to an existing file.}

\item{name.tag}{\code{character(1)} The tag to use when writing the motifs
name slot.}

\item{altname.tag}{\code{character(1)} The tag to use when writing the
motifs altname slot. Note that no tag will be written if the
slot is empty.}
}
\value{
\code{NULL}, invisibly.
}
\description{
Convert motifs to TRANSFAC format and write to file.
}
\details{
If the family slot of a motif is not empty, then its contents will
included using the HC tag. Similarly for the organism slot using the
tag OS. The default name and alternate name tags are ID and NA,
respectively, though these can be set manually.
}
\examples{
jaspar <- read_jaspar(system.file("extdata", "jaspar.txt",
                                    package = "universalmotif"))
write_transfac(jaspar, tempfile())

}
\references{
Wingender E, Dietze P, Karas H, Knuppel R (1996). “TRANSFAC: A
Database on Transcription Factors and Their DNA Binding Sites.”
\emph{Nucleic Acids Research}, \strong{24}, 238-241.
}
\seealso{
\code{\link[=read_transfac]{read_transfac()}}

Other write_motifs: 
\code{\link{write_homer}()},
\code{\link{write_jaspar}()},
\code{\link{write_matrix}()},
\code{\link{write_meme}()},
\code{\link{write_motifs}()}
}
\author{
Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@uwaterloo.ca}
}
\concept{write_motifs}
