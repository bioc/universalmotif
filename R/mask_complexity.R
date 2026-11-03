#' Mask low-complexity regions of DNA and RNA sequences.
#'
#' Homopolymers, microsatellites and other short repeats give motifs rich in
#' the same letters many overlapping hits, which can inflate scanning results
#' and enrichment. `mask_complexity()` finds these regions with the symmetric
#' DUST algorithm and replaces them with `letter`, so that [scan_sequences()],
#' [enrich_motifs()] and the like skip them.
#'
#' @param seqs `XStringSet` DNA or RNA sequences. Read a soft-masked FASTA
#'   file with [Biostrings::readBStringSet()] to keep its lowercase letters,
#'   since `DNAStringSet` objects have none.
#' @param window.size `integer(1)` The DUST window in bases.
#' @param threshold `integer(1)` The score above which a region is masked. A
#'   higher threshold masks less.
#' @param lowercase `logical(1)` Also mask lowercase letters, as soft-masked
#'   genomes mark their repeats.
#' @param letter `character(1)` The masking letter. It must be valid for the
#'   class of `seqs`, such as `"-"` or `"N"` for `DNAStringSet`.
#' @param return.ranges `logical(1)` Return the regions as a `GRanges`
#'   object, which [mask_ranges()] takes, instead of the masked sequences.
#'   Unnamed sequences are named by their place in `seqs`. Requires the
#'   \pkg{GenomicRanges} package.
#' @param nthreads `numeric(1)` Number of threads. Zero uses all available
#'   threads.
#'
#' @return The masked sequences, of the same class as `seqs`, or with
#'   `return.ranges = TRUE` a `GRanges` object of the regions.
#'
#' @details
#' Regions are found with symmetric DUST (Morgulis et al. 2006) as
#' implemented in sdust (Li 2016), and match those of NCBI's `dustmasker`
#' with the same `window.size` and `threshold`, 64 and 20 by default.
#' Letters other than A, C, G, T and U, in either case, split a sequence into
#' parts that are scored on their own, so that a gap of Ns or an earlier mask
#' does not join the sequence either side of it; `dustmasker` instead scores
#' across them, and caps the window at 64.
#'
#' When masking sequences before an enrichment test, mask the targets and the
#' background in the same way, and if the background is made by shuffling
#' the targets, mask after shuffling, or the masking letters will be
#' shuffled into the background.
#'
#' A `BStringSet` given with `lowercase = TRUE` comes back as a `BStringSet`,
#' which can be turned into a `DNAStringSet` for scanning with
#' [Biostrings::DNAStringSet()].
#'
#' @references
#' Li H (2016). "Minimap and miniasm: fast mapping and de novo assembly for
#' noisy long sequences." *Bioinformatics*, **32**, 2103-2110.
#'
#' Morgulis A, Gertz EM, Schaffer AA, Agarwala R (2006). "A fast and
#' symmetric DUST implementation to mask low-complexity DNA sequences."
#' *Journal of Computational Biology*, **13**, 1028-1040.
#'
#' @examples
#' seqs <- Biostrings::DNAStringSet(c(
#'   s1 = "GGTTGTACTGTACAGGAACTTATAGAATAATTCTGGATGTATATATATATATATATATATATATATATATTCCAAACG",
#'   s2 = "ACGTACGTAAAAAAAAAAAAAAAAAAAAAAAAACGTTGCAGTCGATCGATGCATGCAT"))
#' mask_complexity(seqs)
#'
#' ## Fewer, stronger regions with a higher threshold, as ranges:
#' if (requireNamespace("GenomicRanges", quietly = TRUE)) {
#'   mask_complexity(seqs, threshold = 30, return.ranges = TRUE)
#' }
#'
#' ## Soft-masked letters too:
#' soft <- Biostrings::BStringSet(c(s3 = "acgtacgtCAGTCGATCGATGCATGCATTTGACCAGT"))
#' mask_complexity(soft, lowercase = TRUE)
#'
#' @seealso [mask_seqs()], [mask_ranges()], [sequence_complexity()],
#'   [scan_sequences()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @export
mask_complexity <- function(seqs, window.size = 64, threshold = 20,
  lowercase = FALSE, letter = "-", return.ranges = FALSE, nthreads = 1) {

  if (!is(seqs, "XStringSet"))
    stop("`seqs` must be an `XStringSet` object", call. = FALSE)
  alph <- seqtype(seqs)
  if (alph == "AA")
    stop("`mask_complexity()` only works with DNA and RNA sequences",
      call. = FALSE)
  whole <- function(x, lo) is.numeric(x) && length(x) == 1 && !is.na(x) &&
    x == round(x) && x >= lo
  if (!whole(window.size, 4))
    stop("`window.size` must be a whole number above 3", call. = FALSE)
  if (!whole(threshold, 1))
    stop("`threshold` must be a whole number above 0", call. = FALSE)
  if (!is.logical(lowercase) || length(lowercase) != 1 || is.na(lowercase))
    stop("`lowercase` must be TRUE or FALSE", call. = FALSE)
  if (!is.character(letter) || length(letter) != 1 || nchar(letter) != 1)
    stop("`letter` must be a single character", call. = FALSE)
  if (!whole(nthreads, 0))
    stop("`nthreads` must be a whole number of 0 or more", call. = FALSE)

  strs <- as.character(seqs)

  if (return.ranges) {
    if (!requireNamespace("GenomicRanges", quietly = TRUE))
      stop("`return.ranges = TRUE` requires the GenomicRanges package",
        call. = FALSE)
    found <- dust_intervals_cpp(unname(strs), window.size, threshold,
      lowercase, nthreads)
    nms <- names(seqs)
    if (is.null(nms)) nms <- as.character(seq_along(seqs))
    n <- vapply(found, function(x) length(x$start), integer(1))
    return(GenomicRanges::GRanges(rep(nms, n), IRanges::IRanges(
      unlist(lapply(found, function(x) x$start)),
      unlist(lapply(found, function(x) x$end)))))
  }

  strs <- mask_complexity_cpp(unname(strs), window.size, threshold,
    lowercase, letter, nthreads)
  out <- switch(alph, DNA = DNAStringSet(strs), RNA = RNAStringSet(strs),
    BStringSet(strs))
  names(out) <- names(seqs)
  out

}
