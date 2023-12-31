#' Scan sequences for matches to input motifs.
#'
#' For sequences of any alphabet, scan them using the PWM matrices of
#' a set of input motifs.
#'
#' @param motifs See `convert_motifs()` for acceptable motif formats.
#' @param sequences \code{\link{XStringSet}} Sequences to scan. Alphabet
#'    should match motif.
#' @param threshold `numeric(1)` See details.
#' @param threshold.type `character(1)` One of `c('pvalue', 'qvalue',
#'    'logodds', 'logodds.abs')`. See details.
#' @param RC `logical(1)` If `TRUE`, check reverse complement of the input
#'    sequences. Only available for DNA/RNA.
#' @param use.freq `numeric(1)` The default, 1, uses the motif matrix (from
#'    the `motif['motif']` slot) to search for sequences. If a higher
#'    number is used, then the matching k-let matrix from the
#'    `motif['multifreq']` slot is used. See [add_multifreq()].
#' @param verbose `numeric(1)` Describe progress, from none (`0`) to
#'    verbose (`3`).
#' @param nthreads `numeric(1)` Run [scan_sequences()] in parallel with `nthreads`
#'    threads. `nthreads = 0` uses all available threads.
#'    Note that no speed up will occur for jobs with only a single motif and
#'    sequence.
#' @param motif_pvalue.k `numeric(1)` Control [motif_pvalue()] approximation.
#'    See [motif_pvalue()]. Only used if `motif_pvalue.method = "exhaustive"`.
#' @param use.gaps `logical(1)` Set this to `FALSE` to ignore motif gaps, if
#'    present.
#' @param allow.nonfinite `logical(1)` If `FALSE`, then apply a pseudocount if
#'    non-finite values are found in the PWM. Note that if the motif has a
#'    pseudocount greater than zero and the motif is not currently of type PWM,
#'    then this parameter has no effect as the pseudocount will be
#'    applied automatically when the motif is converted to a PWM internally. This
#'    value is set to `FALSE` by default in order to stay consistent with
#'    pre-version 1.8.0 behaviour. Also note that this parameter is not
#'    compatible with `motif_pvalue.method = "dynamic"`. A message will be printed
#'    if a pseudocount is applied. To disable this, set
#'    `options(pseudocount.warning=FALSE)`.
#' @param warn.NA `logical(1)` Whether to warn about the presence of non-standard
#'    letters in the input sequence, such as those in masked sequences.
#' @param calc.pvals `logical(1)` Calculate P-values for each hit. This is a
#'    convenience option which simply gives `motif_pvalue()` the input motifs
#'    and the scores of each hit. Be careful about setting this to `TRUE` if
#'    you anticipate getting thousands of hits and are using
#'    `motif_pvalue.method = "exhaustive"`: expect to wait a few seconds or
#'    minutes for the calculations to finish. Increasing the `nthreads` value
#'    can help greatly here. See Details for more information on P-value
#'    calculation. If `motif_pvalue.method = "dynamic"`, then this is usually
#'    not an issue.
#' @param return.granges `logical(1)` Return the results as a `GRanges` object.
#'    Requires the `GenomicRanges` package to be installed.
#' @param no.overlaps `logical(1)` Remove overlapping hits from the same motifs.
#'    Overlapping hits from different motifs are preserved. Please note that the
#'    current implementation of this feature can add significantly to the run
#'    time for large inputs.
#' @param no.overlaps.by.strand `logical(1)` Whether to discard overlapping hits
#'    from the opposite strand (`TRUE`), or to only discard overlapping hits on the
#'    same strand (`FALSE`).
#' @param no.overlaps.strat `character(1)` One of `c("score", "order")`.
#'    The former option keeps the highest scoring overlapping hit (and the first
#'    of these within ties), and the latter simply keeps the first overlapping hit.
#' @param respect.strand `logical(1)` If  motifs are DNA/RNA,
#'    then setting this option to `TRUE` will make `scan_sequences()` only
#'    scan the strands of the input sequences as indicated in the motif
#'    `strand` slot.
#' @param motif_pvalue.method `character(1)` One of `c("dynamic", "exhaustive")`.
#'    Algorithm used for calculating P-values. The `"exhaustive"` method
#'    involves finding all possible motif matches at or above the specified
#'    score using a branch-and-bound algorithm, which can be computationally
#'    intensive (Hartman et al., 2013). Additionally, the computation
#'    must be repeated for each hit. The `"dynamic"` method calculates the
#'    distribution of possible motif scores using a much faster dynamic
#'    programming algorithm, and can be recycled for multiple
#'    scores (Grant et al., 2011). The only
#'    disadvantage is the inability to use `allow.nonfinite = TRUE`.
#'    See [motif_pvalue()] for details.
#' @param calc.qvals `logical(1)` Whether to also calculate adjusted
#'    P-values. Only valid if `calc.pvals = TRUE`.
#' @param calc.qvals.method `character(1)` One of `c("fdr", "BH", "bonferroni")`.
#'    The method for calculating adjusted P-values. These are described in
#'    depth in the Sequence Searches vignette. Also see Noble (2009).
#'
#' @return `DataFrame`, `GRanges` with each row representing one hit. If the input
#'    sequences are \code{\link{DNAStringSet}} or \code{\link{RNAStringSet}},
#'    then an additional column with the strand is included. Function args are
#'    stored in the `metadata` slot. If `return.granges = TRUE`
#'    then a `GRanges` object is returned.
#'
#' @details
#'
#' ## Logodds scoring
#' Similar to [Biostrings::matchPWM()], the scanning method uses
#' logodds scoring. (To see the scoring matrix for any motif, simply
#' run `convert_type(motif, "PWM")`. For a `multifreq` scoring
#' matrix: `apply(motif["multifreq"][["2"]], 2, ppm_to_pwm)`). In order
#' to score a sequence, at each position within a sequence of length equal
#' to the length of the motif, the scores for each base are summed. If the
#' score sum is above the desired threshold, it is kept.
#'
#' ## Thresholds
#' If `threshold.type = 'logodds'`, then the `threshold` value is multiplied
#' by the maximum possible motif scores. To calculate the
#' maximum possible scores a motif (of type PWM) manually, run
#' `motif_score(motif, 1)`. If \code{threshold.type = 'pvalue'},
#' then threshold logodds scores are generated using [motif_pvalue()].
#' Finally, if \code{threshold.type = 'logodds.abs'}, then the exact values
#' provided will be used as thresholds. Finally, if `threshold.type = 'qvalue'`,
#' then the threshold is calculated as if `threshold.type = 'pvalue'` and the
#' final set of hits are filtered based on their calculated Q-value. (Note:
#' this means that the `thresh.score` column will be incorrect!) This is done
#' since most Q-values cannot be calculated prior to scanning. If you are
#' running a very large job, it may be wise to use a P-value threshold
#' followed by manually filtering by Q-value; this will avoid the scanning
#' have to parse the larger number of hits from the internally-lowered threshold.
#'
#' ## Non-standard letters
#' Non-standard letters (such as "N", "+", "-", ".", etc in \code{\link{DNAString}}
#' objects) will be safely ignored, resulting only in a warning and a very
#' minor performance cost. This can used to scan
#' masked sequences. See \code{\link[Biostrings:maskMotif]{Biostrings::mask()}}
#' for masking sequences
#' (generating \code{\link{MaskedXString}} objects), and [Biostrings::injectHardMask()]
#' to recover masked \code{\link{XStringSet}} objects for use with [scan_sequences()].
#' There is also a provided wrapper function which performs both steps: [mask_seqs()].
#'
#' @references
#'
#' Grant CE, Bailey TL, Noble WS (2011). "FIMO: scanning for occurrences
#' of a given motif." *Bioinformatics*, **27**, 1017-1018.
#'
#' Hartmann H, Guthohrlein EW, Siebert M, Soding SLJ (2013).
#' “P-value-based regulatory motif discovery using positional weight
#' matrices.” *Genome Research*, **23**, 181-194.
#'
#' Noble WS (2009). "How does multiple testing work?" *Nature Biotechnology*,
#' **27**, 1135-1137.
#'
#' @examples
#' ## any alphabet can be used
#' \dontrun{
#' set.seed(1)
#' alphabet <- paste(c(letters), collapse = "")
#' motif <- create_motif("hello", alphabet = alphabet)
#' sequences <- create_sequences(alphabet, seqnum = 1000, seqlen = 100000)
#' scan_sequences(motif, sequences)
#' }
#'
#' ## Sequence masking:
#' if (R.Version()$arch != "i386") {
#' library(Biostrings)
#' data(ArabidopsisMotif)
#' data(ArabidopsisPromoters)
#' seq <- mask_seqs(ArabidopsisPromoters, "AAAAA")
#' scan_sequences(ArabidopsisMotif, seq)
#' # A warning regarding the presence of non-standard letters will be given,
#' # but can be safely ignored in this case.
#' }
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @seealso [add_multifreq()], [Biostrings::matchPWM()],
#'    [enrich_motifs()], [motif_pvalue()]
#' @export
scan_sequences <- function(motifs, sequences, threshold = 0.0001,
  threshold.type = c("pvalue", "qvalue", "logodds", "logodds.abs"),
  RC = FALSE, use.freq = 1, verbose = 0,
  nthreads = 1, motif_pvalue.k = 8, use.gaps = TRUE, allow.nonfinite = FALSE,
  warn.NA = TRUE, calc.pvals = TRUE, return.granges = FALSE,
  no.overlaps = FALSE, no.overlaps.by.strand = FALSE,
  no.overlaps.strat = c("score", "order"),
  respect.strand = FALSE, motif_pvalue.method = c("dynamic", "exhaustive"),
  calc.qvals = calc.pvals, calc.qvals.method = c("fdr", "BH", "bonferroni")) {

  # TODO: add a flag to use the bkg probabilities from the actual input sequence
  # to be used in motif_pvalue() instead of using the bkgs from the motifs

  # param check --------------------------------------------
  args <- as.list(environment())
  all_checks <- character()
  num_check <- check_fun_params(list(threshold = args$threshold,
                                     use.freq = args$use.freq,
                                     verbose = args$verbose,
                                     nthreads = args$nthreads,
                                     motif_pvalue.k = args$motif_pvalue.k),
                                c(0, 1, 1, 1, 1), logical(), TYPE_NUM)
  logi_check <- check_fun_params(list(RC = args$RC, use.gaps = args$use.gaps,
                                      return.granges = args$return.granges,
                                      no.overlaps = args$no.overlaps,
                                      calc.qvals = args$calc.qvals,
                                      no.overlaps.by.strand = args$no.overlaps.by.strand),
                                 numeric(), logical(), TYPE_LOGI)
  s4_check <- check_fun_params(list(sequences = args$sequences), numeric(),
                               logical(), TYPE_S4)
  all_checks <- c(all_checks, num_check, logi_check, s4_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  motif_pvalue.method <- match.arg(motif_pvalue.method)
  calc.qvals.method <- match.arg(calc.qvals.method)
  threshold.type <- match.arg(threshold.type)
  no.overlaps.strat <- match.arg(no.overlaps.strat)

  if (motif_pvalue.method == "dynamic" && allow.nonfinite
      && (calc.pvals = TRUE || threshold.type %in% c("pvalue", "qvalue")))
    stop(wmsg("`motif_pvalue.method = \"dynamic\"` and `allow.nonfinite = TRUE` are ",
        "not compatible when `calc.pvals = TRUE` or ",
        "`threshold.type = c(\"pvalue\", \"qvalue\")`"), call. = FALSE)

  if (verbose > 2) {
    message(" * Input parameters")
    message("   * motifs:              ", deparse(substitute(motifs)))
    message("   * sequences:           ", deparse(substitute(sequences)))
    message("   * threshold:           ", ifelse(length(threshold) > 1, "...",
                                                 threshold))
    message("   * threshold.type:      ", threshold.type)
    message("   * RC:                  ", RC)
    message("   * respect.strand:      ", respect.strand)
    message("   * use.freq:            ", use.freq)
    message("   * use.gaps:            ", use.gaps)
    message("   * calc.pvals:          ", calc.pvals)
    message("   * calc.qvals:          ", calc.qvals)
    message("   * calc.qvals.method:   ", calc.qvals.method)
    message("   * no.overlaps:         ", no.overlaps)
    message("   * verbose:             ", verbose)
  }

  if (!no.overlaps.strat %in% c("score", "order"))
    stop("`no.overlaps.strat` must be \"score\" or \"order\"", call. = FALSE)

  if (missing(motifs) || missing(sequences)) {
    stop("need both motifs and sequences")
  }

  if (calc.qvals && !calc.pvals)
    message("`calc.qvals = TRUE` is ignored when `calc.pvals = FALSE`")

  if (RC && respect.strand)
    message(wmsg("Note: `RC=TRUE` is ignored when `respect.strand=TRUE`"))
  else if (respect.strand)
    RC <- TRUE

  if (verbose > 0) message(" * Processing motifs")

  if (verbose > 1) message(
    "   * Scanning ", length(motifs),
    ifelse(length(motifs) > 1, " motifs", " motif"), " in ", length(sequences),
    ifelse(length(sequences) > 1, " sequences", " sequence"),
    " of average size ", round(mean(width(sequences))))

  motifs <- convert_motifs(motifs)
  if (!is.list(motifs)) motifs <- list(motifs)
  motifs <- convert_type_internal(motifs, "PWM")
  needsfix <- vapply(motifs, function(x) any(is.infinite(x@motif)), logical(1))
  if (any(needsfix) && !allow.nonfinite) {
    warn_pseudo(paste0("Set `allow.nonfinite = TRUE` to prevent this behaviour ",
      "(when `motif_pvalue.method = \"exhaustive\", or `calc.pvals = FALSE`",
      " and `threshold.type = c(\"logodds\", \"logodds.abs\")`)."))
    for (i in which(needsfix)) {
      motifs[[i]] <- suppressMessages(normalize(motifs[[i]]))
    }
  }

  mot.names <- vapply(motifs, function(x) x@name, character(1))

  mot.gaps <- lapply(motifs, function(x) x@gapinfo)
  mot.hasgap <- vapply(mot.gaps, function(x) x@isgapped, logical(1))
  if (any(mot.hasgap) && use.gaps) {
    gapdat <- process_gapped_motifs(motifs, mot.hasgap)
  }

  mot.pwms <- lapply(motifs, function(x) x@motif)
  mot.alphs <- vapply(motifs, function(x) x@alphabet, character(1))
  if (length(unique(mot.alphs)) != 1) stop("can only scan using one alphabet")
  mot.alphs <- unique(mot.alphs)
  if (verbose > 1) message("   * Motif alphabet: ", mot.alphs)

  seq.names <- names(sequences)
  if (is.null(seq.names)) seq.names <- as.character(seq_len(length(sequences)))

  seq.alph <- seqtype(sequences)
  if (seq.alph != "B" && seq.alph != mot.alphs)
    stop("Motif and Sequence alphabets do not match")
  else if (seq.alph == "B")
    seq.alph <- mot.alphs
  if (respect.strand && !seq.alph %in% c("DNA", "RNA"))
    stop("`respect.strand = TRUE` is only valid for DNA/RNA motifs")
  if (RC && !seq.alph %in% c("DNA", "RNA")) {
    warning("`RC = TRUE` is only valid for DNA/RNA motifs, ignoring")
    RC <- FALSE
  }

  if (use.freq > 1) {
    if (any(mot.hasgap) && use.gaps)
      stop("use.freq > 1 cannot be used with gapped motifs")
    if (any(vapply(motifs, function(x) length(x@multifreq) == 0, logical(1))))
      stop("missing multifreq slots")
    check_multi <- vapply(motifs,
                          function(x) any(names(x@multifreq) %in%
                                          as.character(use.freq)),
                          logical(1))
    if (!any(check_multi)) stop("not all motifs have correct multifreqs")
  }

  if (use.freq == 1) {
    score.mats <- mot.pwms
  } else {
    score.mats <- lapply(motifs,
                         function(x) x@multifreq[[as.character(use.freq)]])
    for (i in seq_along(score.mats)) {
      score.mats[[i]] <- MATRIX_ppm_to_pwm(score.mats[[i]],
                                           nsites = motifs[[i]]@nsites,
                                           pseudocount = motifs[[i]]@pseudocount,
                                           bkg = motifs[[i]]@bkg[rownames(score.mats[[i]])])
    }
  }

  max.scores <- vapply(motifs, function(x)
    suppressMessages(motif_score(x, 1, use.freq, threshold.type = "fromzero",
        allow.nonfinite = allow.nonfinite)),
    numeric(1))
  if (!allow.nonfinite)
    min.scores <- vapply(motifs, function(x)
      suppressMessages(motif_score(x, 0, use.freq)), numeric(1))
  else
    min.scores <- vapply(motifs, function(x) motif_score_min(x, use.freq), numeric(1))

  if (threshold.type == "logodds") {

    thresholds <- max.scores * threshold

  } else if (threshold.type == "logodds.abs") {

    if (!length(threshold) %in% c(length(motifs), 1))
      stop(wmsg("for threshold.type = 'logodds.abs', a threshold must be provided for
                every single motif or one threshold recycled for all motifs"))

    if (length(threshold) == 1) threshold <- rep(threshold, length(motifs))
    thresholds <- threshold

  } else if (threshold.type %in% c("pvalue", "qvalue")) {

    if (threshold.type == "qvalue" && !calc.qvals) {
      stop("`calc.qvals` must be `TRUE` if `threshold.type = \"qvalue\"`")
    }

    if (verbose > 0)
      message(" * Converting P-values to logodds thresholds")
    thresholds <- motif_pvalue(motifs, pvalue = threshold, use.freq = use.freq,
                               method = motif_pvalue.method,
                               k = motif_pvalue.k, allow.nonfinite = allow.nonfinite)
    if (any(is.infinite(thresholds))) {
      stop(wmsg("Found -Inf values in threshold(s); try setting manual ",
          "thresholds with either `threshold.type=` \"logodds\" or ",
          "\"logodds.abs\" instead of \"pvalue\""),
        call. = FALSE)
    }
    for (i in seq_along(thresholds)) {
      if (thresholds[i] > max.scores[i]) thresholds[i] <- max.scores[i]
    }
    if (verbose > 3) {
      for (i in seq_along(thresholds)) {
        message("   * Motif ", mot.names[i], ": max.score = ", max.scores[i],
                ", threshold = ", round(thresholds[i], 3))
      }
    }
    thresholds <- unlist(thresholds)

  } else stop("unknown 'threshold.type'")

  for (i in seq_along(threshold)) {
    if (threshold[i] > max.scores[i])
      warning(wmsg("Threshold [", round(threshold[i], 3), "] for motif ", i,
          " is higher than the max possible threshold [", max.scores[i], "]"),
        immediate. = TRUE, call. = FALSE)
  }

  alph <- switch(seq.alph, "DNA" = "ACGT", "RNA" = "ACGU",
                 "AA" = collapse_cpp(AA_STANDARD2), seq.alph)
  sequences.original <- sequences
  sequences <- as.character(sequences)
  strands <- rep("+", length(score.mats))

  if (any(mot.hasgap) && use.gaps) {
    strands <- strands[gapdat$IDs]
    mot.names <- mot.names[gapdat$IDs]
    score.mats <- lapply(gapdat$motifs, function(x) x@motif)
    thresholds <- thresholds[gapdat$IDs]
    min.scores <- min.scores[gapdat$IDs]
    max.scores <- max.scores[gapdat$IDs]
  }

  score.mats.original <- score.mats

  if (RC || respect.strand) {
    if (respect.strand) {
      mot.strands <- vapply(motifs, function(x) x@strand, character(1))
      keep.pos <- rep(TRUE, length(motifs))
      keep.neg <- rep(TRUE, length(motifs))
      keep.pos[mot.strands == "-"] <- FALSE
      keep.neg[mot.strands == "+"] <- FALSE
    } else {
      keep.pos <- rep(TRUE, length(motifs))
      keep.neg <- rep(TRUE, length(motifs))
    }
    strands <- c(strands[keep.pos], rep("-", length(score.mats))[keep.neg])
    mot.names <- c(mot.names[keep.pos], mot.names[keep.neg])
    thresholds <- c(thresholds[keep.pos], thresholds[keep.neg])
    score.mats.rc <- lapply(score.mats,
                            function(x) matrix(rev(as.numeric(x)), nrow = nrow(x)))
    score.mats <- c(score.mats[keep.pos], score.mats.rc[keep.neg])
    min.scores <- c(min.scores[keep.pos], min.scores[keep.neg])
    max.scores <- c(max.scores[keep.pos], max.scores[keep.neg])
    mot.indices <- c(seq_along(motifs)[keep.pos], seq_along(motifs)[keep.neg])
    motifs <- c(motifs[keep.pos], motifs[keep.neg])
  }

  thresholds[thresholds == Inf] <- min_max_ints()$max / 1000
  thresholds[thresholds == -Inf] <- min_max_ints()$min / 1000

  if (allow.nonfinite) {
    for (i in seq_along(score.mats)) {
      if (any(is.infinite(score.mats[[i]]))) {
        min_val1 <- min_max_ints()$min / ncol(score.mats[[i]])
        min_val2 <- as.integer(log2(nrow(score.mats[[i]])) * ncol(score.mats[[i]])) * 1000
        min_val <- (min_val1 + min_val2) / 1000
        score.mats[[i]][is.infinite(score.mats[[i]])] <- min_val
      }
    }
  }

  if (verbose > 0) message(" * Scanning")

  res <- scan_sequences_cpp(score.mats, sequences, use.freq, alph, thresholds,
    nthreads, allow.nonfinite, warn.NA)

  if (verbose > 1) message("   * Number of matches: ", nrow(res))
  if (verbose > 0) message(" * Processing results")

  thresholds[thresholds <= min_max_ints()$min / 1000] <- -Inf
  thresholds[thresholds >= min_max_ints()$max / 1000] <- Inf

  res$thresh.score <- thresholds[res$motif]
  res$min.score <- min.scores[res$motif]
  res$max.score <- max.scores[res$motif]
  res$score.pct <- res$score / res$max.score * 100
  if (seq.alph %in% c("DNA", "RNA")) res$strand <- strands[res$motif]
  res$motif <- mot.names[res$motif]
  res$sequence <- seq.names[res$sequence]

  if (nrow(res) == 0) message("No hits found.")

  if (RC && nrow(res) > 0) res <- adjust_rc_hits(res, seq.alph)

  if (RC) res$motif.i <- mot.indices[res$motif.i]

  out <- as(res, "DataFrame")
  out@metadata <- list(
    args = args[-c(1:2)],
    seqlengths = structure(width(sequences), names = names(sequences))
  )

  if (nrow(out) && any(mot.hasgap) && use.gaps) {
    out$match <- add_gap_dots_cpp(out$match, gapdat$gaplocs[out$motif.i])
    out$motif.i <- gapdat$IDs[out$motif.i]
  }

  if (verbose > 1) message("   * Calculating P-values")
  if (nrow(out) && calc.pvals) {

    out$pvalue <- NA_real_

    if (motif_pvalue.method == "exhaustive") {
      out$pvalue <- motif_pvalue(motifs[out$motif.i], out$score, use.freq = use.freq,
        nthreads = nthreads, allow.nonfinite = allow.nonfinite, k = motif_pvalue.k,
        method = motif_pvalue.method)
    } else {
      # TODO: do multithreaded w/ c++
      for (i in unique(out$motif.i)) {
        which.rows <- which(out$motif.i == i)
        out$pvalue[which.rows] <- motif_pvalue(motifs[[i]], out$score[which.rows],
          use.freq = use.freq, nthreads = nthreads, allow.nonfinite = allow.nonfinite,
          k = motif_pvalue.k, method = motif_pvalue.method)
      }
    }

    if (verbose > 1) message("   * Calculating Q-values")
    if (calc.qvals) {

      max_hits <- function(m, seqs) {
        mLen <- ncol(m)
        mMax <- sum(width(seqs) - mLen + 1)
        if (RC) mMax * 2 else mMax
      }

      out$qvalue <- NA_real_

      # TODO: do multithreaded w/ c++
      if (calc.qvals.method == "fdr") {
        for (i in unique(out$motif.i)) {
          which.rows <- which(out$motif.i == i)
          out$qvalue[which.rows] <- calc_motif_fdr(
            max_hits(motifs[[i]], sequences), out$score[which.rows], out$pvalue[which.rows])
        }
      } else if (calc.qvals.method == "BH") {
        for (i in unique(out$motif.i)) {
          which.rows <- which(out$motif.i == i)
          out$qvalue[which.rows] <- calc_motif_bh(
            max_hits(motifs[[i]], sequences), out$pvalue[which.rows])
        }
      } else if (calc.qvals.method == "bonferroni") {
        for (i in unique(out$motif.i)) {
          which.rows <- which(out$motif.i == i)
          out$qvalue[which.rows] <- calc_motif_bonferroni(
            max_hits(motifs[[i]], sequences), out$pvalue[which.rows])
        }
      }

      if (threshold.type == "qvalue") {
        out <- out[out$qvalue <= threshold, ]
      }

    }

  } else if (!nrow(out) && calc.pvals) {
    out$pvalue <- numeric()
    if (calc.qvals) {
      out$qvalue <- numeric()
    }
  }

  if (nrow(out) && no.overlaps) {
    if (verbose > 1) message("   * Removing overlapping hits")
    # TODO: multithreaded c++
    if (RC && no.overlaps.by.strand) {
      row.indices.plus <- which(out$strand == "+")
      row.indices.minus <- which(out$strand == "-")
      row.indices.plus <- remove_masked_hits(out, row.indices.plus, no.overlaps.strat)
      row.indices.minus <- remove_masked_hits(switch_antisense_coords_cpp(out),
        row.indices.minus, no.overlaps.strat)
      row.indices <- c(row.indices.plus, row.indices.minus)
    } else if (RC) {
      row.indices <- remove_masked_hits(switch_antisense_coords_cpp(out),
        seq_len(nrow(out)), no.overlaps.strat)
    } else {
      row.indices <- seq_len(nrow(out))
      row.indices <- remove_masked_hits(out, seq_len(nrow(out)), no.overlaps.strat)
    }
    out <- out[row.indices, ]
  }

  if (verbose > 1) message(" * Final number of matches: ", nrow(out))

  if (return.granges) {
    if (verbose > 1) message("   * Processing results as GRanges")
    if (is.null(names(sequences))) {
      # warning(wmsg("Input sequences have no names, assigning names 1:",
      #     length(sequences)), call. = FALSE)
      names(sequences) <- 1:length(sequences)
    }
    colnames(out)[3] <- "seqname"
    if (RC) {
      out <- switch_antisense_coords_cpp(out)
    }
    out <- granges_fun(GenomicRanges::GRanges(out,
        seqlengths = structure(width(sequences), names = names(sequences))))
    sort(out)
  } else {
    out[order(out$motif.i, out$sequence, out$start), ]
  }

}

calc_motif_bonferroni <- function(mMax, pvals) {
  pmin(pvals * mMax, 1)
}

calc_motif_bh <- function(mMax, pvals) {
  pmin(pvals / ((rank(pvals) / mMax) * 100), 1)
  # pmin(pvals * (rank(pvals) / mMax), 1)
}

calc_motif_fdr <- function(mMax, scores, pvals) {
  scoreDF <- data.frame(OriginalOrder = seq_along(scores), Scores = scores,
    Pvals = pvals)
  scoreDF <- scoreDF[order(scoreDF$Scores), ]
  scoreDF$ObsHits <- rev(seq_len(nrow(scoreDF)))
  scoreDF$NulHits <- mMax * scoreDF$Pvals
  scoreDF$FDR <- scoreDF$NulHits / scoreDF$ObsHits
  scoreDF$FDR <- cummin(scoreDF$FDR)
  scoreDF$FDR <- pmin(scoreDF$FDR, 1)

  scoreDF$FDR[order(scoreDF$OriginalOrder)]
}

remove_masked_hits <- function(x, i = seq_len(nrow(x)), strat = "score") {
  if (!length(i)) return(i)
  y <- x[i, ]
  y$index.tokeep <- i
  switch(strat, score = remove_masked_hits_by_score(y),
    order = remove_masked_hits_by_order(y))
}

remove_masked_hits_by_order <- function(y) {
  sort(unlist(by(y, list(y$sequence, y$motif.i), function(z) {
    dedup_by_order(z, flatten_group_matrix(get_overlap_groups(z)))
  }, simplify = FALSE)))
}

remove_masked_hits_by_score <- function(y) {
  sort(unlist(by(y, list(y$sequence, y$motif.i), function(z) {
    dedup_by_score(z, flatten_group_matrix(get_overlap_groups(z)))
  }, simplify = FALSE)))
}

get_overlap_groups <- function(x) {
  y <- as.matrix(findOverlaps(IRanges(x$start, x$stop)))
  y <- xtabs(~queryHits + subjectHits, y)
  matrix(as.integer(y), nrow = nrow(x))
}

flatten_group_matrix <- function(x) {
  if (all(x == 1)) {
    # All overlapping
    rep(1, length(diag(x)))
  } else if (!sum(x[lower.tri(x)]) && !sum(x[upper.tri(x)])) {
    # None overlapping
    seq_along(diag(x))
  } else {
    # Cluster overlapping
    cutree(hclust(as.dist(1 - x)), h = 0.5)
  }
}

dedup_by_order <- function(x, i) {
  x$index.tokeep[!duplicated(i)]
}

dedup_by_score <- function(x, i) {
  unlist(by(x, i, function(y) {
    y$index.tokeep[which.max(y$score)]
  }, simplify = FALSE))
}

adjust_rc_hits <- function(res, alph) {
  rev.strand <- res$strand == "-"
  if (any(rev.strand)) {
    start <- res$stop[rev.strand]
    stop <- res$start[rev.strand]
    res$stop[rev.strand] <- stop
    res$start[rev.strand] <- start
    matches <- res$match[rev.strand]
    if (alph == "DNA")
      matches <- as.character(reverseComplement(DNAStringSet(matches)))
    else if (alph == "RNA")
      matches <- as.character(reverseComplement(RNAStringSet(matches)))
    res$match[rev.strand] <- matches
  }
  res
}

# Note: It's probably a lot faster to scan the individual submotifs and then
# process the gapped motifs afterwards, versus scanning all possible gapped
# motif combinations. Would need to think about how to score the submotifs
# though; so for now, go with the dumb and slow brute force option.

process_gapped_motifs <- function(motifs, hasgap) {
  motifs[hasgap] <- lapply(motifs[hasgap], ungap_single)
  motifs_gapped <- mapply(function(x, y) rep(x, length(y)), hasgap, motifs,
    SIMPLIFY = FALSE)
  IDs <- mapply(function(x, y) rep(x, length(y)), seq_along(motifs), motifs,
    SIMPLIFY = FALSE)
  out <- list(
    motifs = do.call(c, motifs),
    gapped = do.call(c, motifs_gapped),
    IDs = do.call(c, IDs)
  )
  out$gaplocs <- lapply(seq_along(out$motifs), function(x) integer())
  out$gaplocs[out$gapped] <- lapply(out$motifs[out$gapped], get_gaplocs)
  out
}

get_gaplocs <- function(x) {
  y <- strsplit(x@name, "/", fixed = TRUE)[[1]]
  npos <- seq_len(ncol(x@motif))
  lens <- vapply(y, function(x) strsplit(x, "_L", fixed = TRUE)[[1]][2], character(1))
  lens <- as.numeric(lens)
  lens <- lapply(lens, function(x) seq(1, x))
  lenslens <- cumsum(vapply(lens[-length(lens)], length, integer(1)))
  lens <- mapply(function(x, y) x + y, lens, c(0, lenslens), SIMPLIFY = FALSE)
  gapped <- grepl("BLANK", y)
  do.call(c, lens[gapped])
}

get_submotifs <- function(m) {
  n <- length(m@gapinfo@gaploc)
  mname <- m@name
  submotifs <- vector("list", n + 1)
  submotifs[[1]] <- subset(m, seq(1, m@gapinfo@gaploc[1]))
  submotifs[[length(submotifs)]] <- subset(
    m, seq(m@gapinfo@gaploc[n] + 1, ncol(m))
  )
  if (length(submotifs) > 2) {
    for (i in seq_along(submotifs)[-c(1, length(submotifs))]) {
      submotifs[[i]] <- subset(
        m, seq(m@gapinfo@gaploc[i - 1] + 1, m@gapinfo@gaploc[i])
      )
    }
  }
  for (i in seq_along(submotifs)) {
    submotifs[[i]]@name <- paste0("SUB_N", i, "_L", ncol(submotifs[[i]]@motif))
  }
  submotifs
}

make_blank_motif <- function(n, N, alph) {
  alphlen <- switch(alph, DNA = 4, RNA = 4, AA = 20, nchar(alph))
  mot <- matrix(0, nrow = alphlen, ncol = n)
  create_motif(mot, type = "PWM", alphabet = alph,
    name = paste0("BLANK_N", N, "_L", n))
}

ungap_single <- function(m) {
  gaplens <- mapply(
    seq, m@gapinfo@mingap, m@gapinfo@maxgap, SIMPLIFY = FALSE
  )
  gaplens <- expand.grid(gaplens)
  out <- vector("list", nrow(gaplens))
  submotifs <- get_submotifs(m)
  for (i in seq_along(out)) {
    tmp <- list(submotifs[[1]])
    for (j in seq_len(ncol(gaplens))) {
      if (gaplens[[j]][i] == 0) {
        tmp <- c(tmp, list(submotifs[[j + 1]]))
      } else {
        tmp <- c(tmp, list(make_blank_motif(gaplens[[j]][i], j, m@alphabet),
            submotifs[[j + 1]]))
      }
    }
    out[[i]] <- do.call(cbind, tmp)
  }
  out
}

motif_score_min <- function(x, use.freq) {
  if (any(is.infinite(x@motif)))
    -Inf
  else
    suppressMessages(motif_score(x, 0, use.freq))
}

granges_fun <- function(FUN, env = parent.frame()) {
  if (requireNamespace("GenomicRanges", quietly = TRUE)) {
    eval(substitute(FUN), envir = env)
  } else {
    stop(wmsg("The 'GenomicRanges' package must be installed for `return.granges=TRUE`. ",
        "[BiocManager::install(\"GenomicRanges\")]"), call. = FALSE)
  }
}
