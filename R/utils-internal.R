#' Utility functions.
#'
#' Utility functions have been split into two categories: those related to
#' motifs ?'utils-motif', and those related to sequences ?'utils-sequence'.
#'
#' @seealso [utils-motif], [utils-sequence]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @name utilities
NULL

# INTERNAL CONSTANTS -----------------------------------------------------------

DNA_DI <- c("AA", "AC", "AG", "AT",
            "CA", "CC", "CG", "CT",
            "GA", "GC", "GG", "GT",
            "TA", "TC", "TG", "TT")

# AA_STANDARD2 <- sort(AA_STANDARD)
AA_STANDARD2 <- safeExplode("ACDEFGHIKLMNPQRSTVWY")

# TYPE_NULL <- 0L
# TYPE_SYM  <- 1L
# TYPE_ENV  <- 4L
TYPE_LOGI <- 10L
# TYPE_INT  <- 13L
TYPE_NUM  <- 14L
# TYPE_COMP <- 15L
TYPE_CHAR <- 16L
# TYPE_DOT  <- 17L
# TYPE_ANY  <- 18L
TYPE_S4   <- 25L

UNIVERSALMOTIF_SLOTS <- c(

  "name",
  "altname",
  "family",
  "organism",
  "motif",
  "alphabet",
  "type",
  "icscore",
  "nsites",
  "pseudocount",
  "bkg",
  "bkgsites",
  "consensus",
  "strand",
  "pval",
  "qval",
  "eval",
  "multifreq",
  "extrainfo"

)

COMPARE_METRICS <- c("PCC", "EUCL", "SW", "KL", "WEUCL",
                     "ALLR", "BHAT", "HELL", "WPCC",
                     "SEUCL",  "MAN", "ALLR_LL")

# Credit to https://github.com/omarwagih/ggseqlogo/blob/master/R/col_schemes.r
# for the colours.
DNA_COLOURS <- c(A = "#109648", C = "#255C99", G = "#F7B32B", T = "#D62839")
RNA_COLOURS <- c(A = "#109648", C = "#255C99", G = "#F7B32B", U = "#D62839")
AA_COLOURS <- c(G = "#058644", S = "#058644", T = "#058644", Y = "#058644",
  C = "#058644", Q = "#720091", N = "#720091", K = "#0046C5", R = "#0046C5",
  H = "#0046C5", D = "#C5003E", E = "#C5003E", A = "#2E2E2E", V = "#2E2E2E",
  L = "#2E2E2E", I = "#2E2E2E", P = "#2E2E2E", W = "#2E2E2E", F = "#2E2E2E",
  M = "#2E2E2E")

# INTERNAL UTILITIES ----------------------------------------------------------- 

shrink_string <- function(name, maxLen = 5, suffix = "..") {
  if (nchar(name) > maxLen) {
    name <- paste0(substr(name, 1, maxLen), suffix)
  }
  name
}

warn_pseudo <- function(v = 1) {
  # Let's calm down on the warnings maybe...
  if (isTRUE(getOption("pseudocount.warning"))) {
    message(wmsg("Note: Added a pseudocount."))
    # if (v == 1) {
    #   message(wmsg("Note: found -Inf values in motif PWM, adding a pseudocount. ",
    #     "(To turn off this message: `options(pseudocount.warning=FALSE)`.) ",
    #     "Set `allow.nonfinite = TRUE` to prevent this behaviour."))
    # } else if (v == 2) {
    #   message(wmsg("Note: found -Inf values in motif PWM, adding a pseudocount. ",
    #     "(To turn off this message: `options(pseudocount.warning=FALSE)`.)"))
    # } else {
    #   message(wmsg("Note: found -Inf values in motif PWM, adding a pseudocount. ",
    #     "(To turn off this message: `options(pseudocount.warning=FALSE)`.) ", v))
    # }
  }
  invisible()
}

get_nsites <- function(motifs) {
  out <- numeric(length(motifs))
  for (i in seq_along(out)) {
    n <- motifs[[i]]@nsites
    out[i] <- ifelse(length(n) == 1 && n > 1, n, 100)
  }
  out
}

get_bkgs <- function(motifs, use.freq = 1) {

  if (use.freq == 1) {

    out <- lapply(motifs, function(x) x@bkg[seq_len(nrow(x@motif))])

  } else {

    out <- vector("list", length(motifs))
    for (i in seq_along(out)) {
      alph <- rownames(motifs[[i]]@motif)
      alph <- get_klets(alph, use.freq)
      bkg <- motifs[[i]]@bkg[alph]
      if (length(bkg) != nrow(motifs[[i]]@multifreq[[as.character(use.freq)]]))
        stop("Missing higher order background in motif: ", motifs[[i]]@name)
      out[[i]] <- bkg
    }

  }

  out

}

.internal_convert <- function(motifs, class = NULL) {

  if (is.null(class)) {

    CLASS <- class(motifs)
    CLASS_PKG <- attributes(CLASS)$package
    CLASS_IN <- collapse_cpp(c(CLASS_PKG, "-", CLASS))

    CLASS_IN

  } else {

    if (length(class) == 1 && class[1] != "universalmotif-universalmotif" &&
        class[1] != "MotifDb-MotifList") {

      tryCatch(motifs <- convert_motifs(motifs, class),
               error = function(e) message("motifs converted to class 'universalmotif'"))

    } else if (length(class) > 1 || class[1] == "MotifDb-MotifList")
      message("motifs converted to class 'universalmotif'")

    motifs

  }

}

# for a motif of length 4, the transition matrix is something like this:
#       bkg pos1 pos2 pos3 pos4
#  bkg    0    1    0    0    0
# pos1    0    0    1    0    0
# pos2    0    0    0    1    0
# pos3    0    0    0    0    1
# pos4    1    0    0    0    0

wmsg2 <- function(..., exdent = 0, indent = 0)
  paste0(strwrap(paste0(..., collapse = ""), exdent = exdent, indent = indent),
         collapse = "\n")

lapply_ <- function(X, FUN, ..., BP = FALSE, PB = FALSE) {

  FUN <- match.fun(FUN)

  if (!BP) {

    if (!PB) {

      out <- lapply(X, FUN, ...)

    } else {

      out <- vector("list", length(X))
      max <- length(X)
      print_pb(0)
      if (is.list(X)) {
        for (i in seq_along(X)) {
          out[[i]] <- do.call(FUN, list(X[[i]], ...))
          update_pb(i, max)
        }
      } else {
        for (i in seq_along(X)) {
          out[[i]] <- do.call(FUN, list(X[i], ...))
          update_pb(i, max)
        }
      }

    }

  } else {

    if (requireNamespace("BiocParallel", quietly = TRUE)) {
      out <- BiocParallel::bplapply(X, FUN, ...)
    } else {
      stop("'BiocParallel' is not installed")
    }
    # BPPARAM <- BiocParallel::bpparam()
    # if (PB) BPPARAM$progressbar <- TRUE
    # out <- BiocParallel::bplapply(X, FUN, ..., BPPARAM = BPPARAM)

  }

  out

}

mapply_ <- function(FUN, ..., MoreArgs = NULL, SIMPLIFY = TRUE,
                    USE.NAMES = TRUE, BP = FALSE, PB = FALSE) {

  FUN <- match.fun(FUN)

  if (!BP) {

    if (!PB) {

      out <- mapply(FUN, ..., MoreArgs = MoreArgs, SIMPLIFY = SIMPLIFY,
                    USE.NAMES = USE.NAMES)

    } else {

      # not sure how to implement USE.NAMES here, get error sometimes
      dots <- list(...)
      dots.len <- vapply(dots, length, numeric(1))
      dots.len.max <- max(dots.len)
      dots <- lapply(dots, rep, length.out = dots.len.max)
      out <- vector("list", dots.len.max)

      print_pb(0)
      for (i in seq_len(dots.len.max)) {
        dots.i <- mapply(function(dots, i) {
                           if (is.list(dots)) dots[[i]]
                           else dots[i]
                    }, dots, i, SIMPLIFY = FALSE)
        out[[i]] <- do.call(FUN, c(dots.i, MoreArgs))
        update_pb(i, dots.len.max)
      }

      if (SIMPLIFY && length(dots))
        out <- simplify2array(out, higher = (SIMPLIFY == "array"))

    }

  } else {

    if (requireNamespace("BiocParallel", quietly = TRUE)) {
      BPPARAM <- BiocParallel::bpparam()
      if (PB) BPPARAM$progressbar <- TRUE
      out <- BiocParallel::bpmapply(FUN, ..., MoreArgs = MoreArgs,
                                    SIMPLIFY = SIMPLIFY, USE.NAMES = USE.NAMES,
                                    BPPARAM = BPPARAM)
    } else {
      stop("'BiocParallel' is not installed")
    }

  }

  out

}
