#include <Rcpp.h>
#include <RcppThread.h>
#include <algorithm>
#include <cmath>
#include "types.h"

void deal_with_higher_k_NA(list_int_t &seq_ints, const int &k, const int &let_len) {

  int tmp = 0;
  for (std::size_t i = 0; i < seq_ints.size(); ++i) {
    for (std::size_t j = 0; j < seq_ints[i].size() - k + 1; ++j) {
      tmp = 0;
      for (int b = 0; b < k; ++b) {
        if (seq_ints[i][j + b] < 0) {
          tmp = -1;
          break;
        }
        tmp += seq_ints[i][j + b] * pow(let_len, k - b - 1);
      }
      seq_ints[i][j] = tmp;
    }
  }

}
void deal_with_higher_k(list_int_t &seq_ints, const int &k, const int &let_len) {

  int tmp = 0;
  for (std::size_t i = 0; i < seq_ints.size(); ++i) {
    for (std::size_t j = 0; j < seq_ints[i].size() - k + 1; ++j) {
      tmp = 0;
      for (int b = 0; b < k; ++b) {
        tmp += seq_ints[i][j + b] * pow(let_len, k - b - 1);
      }
      seq_ints[i][j] = tmp;
    }
  }

}


vec_int_t scan_single_seq_NA(const list_int_t &motif, const vec_int_t &sequence,
    const int &k) {

  vec_int_t result;
  result.reserve(sequence.size());

  int tmp = 0;
  for (std::size_t i = 0; i < sequence.size() - k + 1 - motif.size() + 1; ++i) {
    tmp = 0;
    for (std::size_t j = 0; j < motif.size(); ++j) {
      if (sequence[i + j] < 0)
        tmp += -999999;
      else
        tmp += motif[j][sequence[i + j]];
    }
    result.push_back(tmp);
  }

  return result;

}

vec_int_t scan_single_seq(const list_int_t &motif, const vec_int_t &sequence,
    const int &k) {

  vec_int_t result;
  result.reserve(sequence.size());

  int tmp = 0;
  for (std::size_t i = 0; i < sequence.size() - k + 1 - motif.size() + 1; ++i) {
    tmp = 0;
    for (std::size_t j = 0; j < motif.size(); ++j) {
      tmp += motif[j][sequence[i + j]];
    }
    result.push_back(tmp);
  }

  return result;

}

list_mat_t scan_sequences_cpp_internal(const list_mat_t &score_mats,
    const list_char_t &seq_vecs, const int &k, vec_char_t &alph,
    const int &nthreads, const bool &warnNA) {

  bool use_na_fun = false;
  list_int_t seq_ints(seq_vecs.size());

  vec_int_t na_index(seq_vecs.size(), 0);
  RcppThread::parallelFor(0, seq_vecs.size(),
      [&seq_ints, &alph, &na_index, &seq_vecs] (std::size_t i) {

        seq_ints[i].reserve(seq_vecs[i].size());
        for (std::size_t j = 0; j < seq_vecs[i].size(); ++j) {
          bool na_check = true;
          for (std::size_t a = 0; a < alph.size(); ++a) {
            if (seq_vecs[i][j] == alph[a]) {
              seq_ints[i].push_back(a);
              na_check = false;
              break;
            }
          }
          if (na_check) {
            seq_ints[i].push_back(-1);
            na_index[i] = 1;
          }
        }

      }, nthreads);

  if (std::accumulate(na_index.begin(), na_index.end(), 0) > 0) {
    use_na_fun = true;
    if (warnNA) {
      Rcpp::warning("Non-standard letters detected. These were ignored.");
    }
  }

  if (k > 1 && use_na_fun)
    deal_with_higher_k_NA(seq_ints, k, alph.size());
  else if (k > 1)
    deal_with_higher_k(seq_ints, k, alph.size());

  list_mat_t out(score_mats.size());

  if (use_na_fun) {

    RcppThread::parallelFor(0, out.size(),
        [&out, &score_mats, &seq_ints, &k] (std::size_t i) {
          out[i].reserve(seq_ints.size());
          for (std::size_t j = 0; j < seq_ints.size(); ++j) {
            out[i].push_back(scan_single_seq_NA(score_mats[i], seq_ints[j], k));
          }
        }, nthreads);

  } else {

    RcppThread::parallelFor(0, out.size(),
        [&out, &score_mats, &seq_ints, &k] (std::size_t i) {
          out[i].reserve(seq_ints.size());
          for (std::size_t j = 0; j < seq_ints.size(); ++j) {
            out[i].push_back(scan_single_seq(score_mats[i], seq_ints[j], k));
          }
        }, nthreads);

  }

  return out;

}

list_int_t format_results(const list_mat_t &out_pre, const vec_int_t &scores,
    const list_mat_t &motifs) {

  list_int_t res(5);

  /* not sure what to do about .reserve() here */

  for (std::size_t i = 0; i < out_pre.size(); ++i) {              // motif
    for (std::size_t j = 0; j < out_pre[i].size(); ++j) {         // sequence
      for (std::size_t b = 0; b < out_pre[i][j].size(); ++b) {    // position
        if (out_pre[i][j][b] >= scores[i]) {                      // score
          res[0].push_back(i + 1);                 // motif
          res[1].push_back(j + 1);                 // sequence
          res[2].push_back(b + 1);                 // start
          res[3].push_back(b + motifs[i].size());  // stop
          res[4].push_back(out_pre[i][j][b]);      // score
        }
      }
    }
  }

  return res;

}

vec_str_t get_matches(const list_int_t &res, const vec_str_t &seq_vecs,
    const list_mat_t &motifs) {

  vec_str_t out;
  out.reserve(res[0].size());

  for (std::size_t i = 0; i < res[0].size(); ++i) {
    out.push_back(seq_vecs[res[1][i] - 1].substr(res[2][i] - 1, motifs[res[0][i] - 1].size()));
  }

  return out;

}

void replace_gap_chars(str_t &seqstring, const vec_int_t &gaplocs) {
  for (std::size_t i = 0; i < gaplocs.size(); ++i) {
    seqstring.replace(gaplocs[i] - 1, 1, ".");
  }
}

/* C++ ENTRY ---------------------------------------------------------------- */

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector calc_hit_gc(const Rcpp::StringVector &hits, const bool ignoreN = false) {
  Rcpp::NumericVector res(hits.size());
  if (ignoreN) {
    int GC, AT;
    for (R_xlen_t i = 0; i < hits.size(); ++i) {
      GC = 0; AT = 0;
      for (R_xlen_t j = 0; j < hits[i].size(); ++j) {
        switch (hits[i][j]) {
          case 'S':
          case 'C':
          case 'G': GC++; break;
          case 'W':
          case 'A':
          case 'U':
          case 'T': AT++; break;
        }
      }
      res[i] = double(GC) / double(AT + GC);
    }
  } else {
    int GC;
    for (R_xlen_t i = 0; i < hits.size(); ++i) {
      GC = 0;
      for (R_xlen_t j = 0; j < hits[i].size(); ++j) {
        switch (hits[i][j]) {
          case 'S':
          case 'C':
          case 'G': GC++;
        }
      }
      res[i] = double(GC) / double(hits[i].size());
    }
  }
  return res;
}

// [[Rcpp::export(rng = false)]]
Rcpp::DataFrame switch_antisense_coords_cpp(const Rcpp::DataFrame &res) {
  Rcpp::DataFrame out = Rcpp::clone(res);
  Rcpp::CharacterVector col_strand = res["strand"];
  Rcpp::LogicalVector to_switch(col_strand.size());
  for (R_xlen_t i = 0; i < to_switch.size(); ++i) {
    to_switch[i] = col_strand[i] == "-";
  }
  Rcpp::IntegerVector col_start = res["start"];
  Rcpp::IntegerVector col_stop = res["stop"];
  Rcpp::IntegerVector col_start_new = Rcpp::clone(col_start);
  Rcpp::IntegerVector col_stop_new = Rcpp::clone(col_stop);
  for (R_xlen_t i = 0; i < to_switch.size(); ++i) {
    if (to_switch[i]) {
      col_start_new[i] = col_stop[i];
      col_stop_new[i] = col_start[i];
    }
  } 
  out["start"] = col_start_new;
  out["stop"] = col_stop_new;
  return out;
}

// [[Rcpp::export(rng = false)]]
std::vector<std::string> add_gap_dots_cpp(std::vector<std::string> seqs,
    const std::vector<std::vector<int>> &gaplocs) {
  for (std::size_t i = 0; i < seqs.size(); ++i) {
    if (gaplocs[i].size() > 0) {
      replace_gap_chars(seqs[i], gaplocs[i]);
    }
  }
  return seqs;
}

// [[Rcpp::export(rng = false)]]
Rcpp::DataFrame scan_sequences_cpp(const Rcpp::List &score_mats,
    const std::vector<std::string> &seq_vecs, const int &k, const std::string &alph,
    const std::vector<double> &min_scores, const int &nthreads,
    const bool &allow_nonfinite = false, const bool &warnNA = true) {

  vec_char_t alph2(alph.begin(), alph.end());

  vec_int_t min_scores2;
  min_scores2.reserve(min_scores.size());
  for (std::size_t i = 0; i < min_scores.size(); ++i) {
    min_scores2.push_back(min_scores[i] * 1000);
  }

  list_char_t seq2_vecs;
  seq2_vecs.reserve(seq_vecs.size());
  for (std::size_t i = 0; i < seq_vecs.size(); ++i) {
    seq2_vecs.push_back(vec_char_t(seq_vecs[i].begin(), seq_vecs[i].end()));
  }

  std::vector<int> motif_sizes(score_mats.size());
  std::vector<int> seq_sizes(seq_vecs.size());

  list_mat_t score2_mats(score_mats.size());
  for (R_xlen_t i = 0; i < score_mats.size(); ++i) {
    Rcpp::NumericMatrix tmp = score_mats[i];
    score2_mats[i].reserve(tmp.ncol());
    motif_sizes[i] = tmp.ncol();
    for (R_xlen_t j = 0; j < tmp.ncol(); ++j) {
      Rcpp::NumericVector tmp2 = tmp(Rcpp::_, j);
      tmp2 = tmp2 * 1000;
      score2_mats[i].push_back(vec_int_t(tmp2.begin(), tmp2.end()));
    }
  }

  for (std::size_t i = 0; i < seq_vecs.size(); ++i) {
    seq_sizes[i] = seq_vecs[i].size();
  }
  for (std::size_t i = 0; i < motif_sizes.size(); ++i) {
    for (std::size_t j = 0; j < seq_sizes.size(); ++j) {
      if (seq_sizes[j] < motif_sizes[i]) {
        Rcpp::stop("Found sequence(s) shorter than the width of the motif(s)");
      }
    }
  }

  list_mat_t out_pre = scan_sequences_cpp_internal(score2_mats, seq2_vecs, k,
      alph2, nthreads, warnNA);

  list_int_t res = format_results(out_pre, min_scores2, score2_mats);

  vec_num_t scores2 = vec_num_t(res[4].begin(), res[4].end());
  for (std::size_t i = 0; i < scores2.size(); ++i) {
    scores2[i] /= 1000;
  }

  vec_str_t matches = get_matches(res, seq_vecs, score2_mats);

  return Rcpp::DataFrame::create(
        Rcpp::_["motif"] = res[0],
        Rcpp::_["motif.i"] = res[0],
        Rcpp::_["sequence"] = res[1],
        Rcpp::_["start"] = res[2],
        Rcpp::_["stop"] = res[3],
        Rcpp::_["score"] = scores2,
        Rcpp::_["match"] = matches,
        Rcpp::_["stringsAsFactors"] = false
      );

}
