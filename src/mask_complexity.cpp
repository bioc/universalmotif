#include <Rcpp.h>
#include <RcppThread.h>
#include <deque>
#include <vector>
#include <string>
#include <algorithm>

// Symmetric DUST, ported from sdust.c of minimap2
// (https://github.com/lh3/minimap2), which is under the MIT License:
// Copyright (c) 2018- Dana-Farber Cancer Institute, 2017-2018 Broad
// Institute, Inc. It finds the same regions as dustmasker. Letters other than
// A, C, G, T and U split a sequence into runs scored on their own.

typedef std::pair<std::size_t, std::size_t> dust_interval_t;

static const std::size_t DUST_WORD = 3;

struct dust_perfect_t {
  std::size_t start, finish;
  long r, l;
};

class Dust {

  std::size_t window;
  long threshold;
  // The words of the window, oldest first, in a ring of `cap` slots.
  std::vector<int> ring;
  std::size_t cap, head = 0, nwords = 0;
  // Words of the window, and of its suffix of L words in which no word
  // appears too often, with their counts and scores.
  std::size_t L = 0;
  long rw = 0, rv = 0;
  long cw[64] = {0}, cv[64] = {0};
  std::deque<dust_perfect_t> perfect;

  inline int word(std::size_t k) const { return ring[(head + k) % cap]; }

  void shift_window(int t) {
    if (nwords > window - DUST_WORD) {
      int s = ring[head];
      head = (head + 1) % cap;
      --nwords;
      rw -= --cw[s];
      if (L > nwords) {
        --L;
        rv -= --cv[s];
      }
    }
    ring[(head + nwords) % cap] = t;
    ++nwords;
    ++L;
    rw += cw[t]++;
    rv += cv[t]++;
    if (cv[t] * 10 > threshold * 2) {
      int s;
      do {
        s = word(nwords - L);
        rv -= --cv[s];
        --L;
      } while (s != t);
    }
  }

  // Add the perfect intervals ending at the newest word, kept sorted by
  // descending start. sdust rescans the intervals from the first for every
  // word, but the best score so far is a running maximum, so the scan can
  // carry on from where it stopped.
  void find_perfect(std::size_t start) {
    long c[64];
    std::copy(cv, cv + 64, c);
    long r = rv, max_r = 0, max_l = 0;
    std::size_t n = nwords, j = 0;
    for (std::size_t k = n - L; k-- > 0; ) {
      int t = word(k);
      r += c[t]++;
      long new_r = r, new_l = (long) (n - k - 1);
      if (new_r * 10 > threshold * new_l) {
        while (j < perfect.size() && perfect[j].start >= k + start) {
          const dust_perfect_t &p = perfect[j];
          if (max_r == 0 || p.r * max_l > max_r * p.l) {
            max_r = p.r;
            max_l = p.l;
          }
          ++j;
        }
        if (max_r == 0 || new_r * max_l >= max_r * new_l) {
          max_r = new_r;
          max_l = new_l;
          dust_perfect_t p = {k + start, n + DUST_WORD - 1 + start, new_r, new_l};
          perfect.insert(perfect.begin() + j, p);
          ++j;
        }
      }
    }
  }

  // Move the perfect interval starting first into `out` once the window has
  // passed its start, and drop those it has left behind.
  void save(std::size_t start, std::vector<dust_interval_t> &out) {
    if (perfect.empty() || perfect.back().start >= start) return;
    const dust_perfect_t &p = perfect.back();
    if (!out.empty() && p.start <= out.back().second) {
      out.back().second = std::max(out.back().second, p.finish);
    } else {
      out.push_back(dust_interval_t(p.start, p.finish));
    }
    std::size_t keep = perfect.size();
    while (keep > 0 && perfect[keep - 1].start < start) --keep;
    perfect.resize(keep);
  }

public:

  Dust(std::size_t window, int threshold) : window(window), threshold(threshold),
    ring(window), cap(window) {}

  // Run over the bases of `run`, coded 0 to 3, which start at `offset`.
  void run(const std::vector<int> &codes, std::size_t offset,
      std::vector<dust_interval_t> &out) {
    int t = 0;
    for (std::size_t k = 0; k < codes.size(); ++k) {
      std::size_t l = k + 1;
      t = ((t << 2) | codes[k]) & 63;
      if (l >= DUST_WORD) {
        std::size_t start = (l > window ? l - window : 0) + offset;
        save(start, out);
        shift_window(t);
        if (rw * 10 > (long) L * threshold) find_perfect(start);
      }
    }
    // As in sdust, the flush starts one past where the window would.
    std::size_t n = codes.size() + 1;
    std::size_t start = (n > window ? n - window : 0) + offset + 1;
    while (!perfect.empty()) save(start++, out);
  }

};

static inline int dust_code(char c) {
  switch (c) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': case 'U': case 'u': return 3;
    default: return -1;
  }
}

std::vector<dust_interval_t> dust_one(const std::string &seq, std::size_t window,
    int threshold, bool lowercase) {

  std::vector<dust_interval_t> out;
  std::vector<int> codes;
  std::size_t i = 0, n = seq.size();
  while (i < n) {
    if (dust_code(seq[i]) < 0) {
      ++i;
      continue;
    }
    std::size_t begin = i;
    codes.clear();
    while (i < n && dust_code(seq[i]) >= 0) codes.push_back(dust_code(seq[i++]));
    Dust(window, threshold).run(codes, begin, out);
  }

  if (lowercase) {
    for (std::size_t k = 0; k < n; ) {
      if (seq[k] >= 'a' && seq[k] <= 'z') {
        std::size_t begin = k;
        while (k < n && seq[k] >= 'a' && seq[k] <= 'z') ++k;
        out.push_back(dust_interval_t(begin, k));
      } else {
        ++k;
      }
    }
    std::sort(out.begin(), out.end());
    std::vector<dust_interval_t> merged;
    for (std::size_t k = 0; k < out.size(); ++k) {
      if (!merged.empty() && out[k].first <= merged.back().second) {
        merged.back().second = std::max(merged.back().second, out[k].second);
      } else {
        merged.push_back(out[k]);
      }
    }
    out.swap(merged);
  }

  return out;

}

// [[Rcpp::export(rng = false)]]
Rcpp::List dust_intervals_cpp(const std::vector<std::string> &seqs,
    const int window, const int threshold, const bool lowercase,
    const int nthreads) {

  std::vector<std::vector<dust_interval_t>> found(seqs.size());
  RcppThread::parallelFor(0, seqs.size(),
      [&seqs, &found, window, threshold, lowercase] (std::size_t i) {
        found[i] = dust_one(seqs[i], window, threshold, lowercase);
      }, nthreads);

  Rcpp::List out(seqs.size());
  for (std::size_t i = 0; i < seqs.size(); ++i) {
    Rcpp::IntegerVector starts(found[i].size()), ends(found[i].size());
    for (std::size_t k = 0; k < found[i].size(); ++k) {
      // From 1, with the end included.
      starts[k] = found[i][k].first + 1;
      ends[k] = found[i][k].second;
    }
    out[i] = Rcpp::List::create(Rcpp::_["start"] = starts, Rcpp::_["end"] = ends);
  }
  return out;

}

// [[Rcpp::export(rng = false)]]
std::vector<std::string> mask_complexity_cpp(std::vector<std::string> seqs,
    const int window, const int threshold, const bool lowercase,
    const std::string &letter, const int nthreads) {

  const char l = letter[0];
  RcppThread::parallelFor(0, seqs.size(),
      [&seqs, window, threshold, lowercase, l] (std::size_t i) {
        std::vector<dust_interval_t> found = dust_one(seqs[i], window,
          threshold, lowercase);
        for (std::size_t k = 0; k < found.size(); ++k) {
          std::fill(seqs[i].begin() + found[k].first,
            seqs[i].begin() + found[k].second, l);
        }
      }, nthreads);
  return seqs;

}
