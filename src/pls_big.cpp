#include <Rcpp.h>

#if defined(__has_include)
#  if __has_include(<bigmemory/BigMatrix.h>) &&                                 \
      __has_include(<bigmemory/MatrixAccessor.hpp>)
#    define BIGPLS_HAS_BIGMEMORY 1
#  else
#    define BIGPLS_HAS_BIGMEMORY 0
#  endif
#else
#  define BIGPLS_HAS_BIGMEMORY 1
#endif

using namespace Rcpp;

#if BIGPLS_HAS_BIGMEMORY

#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>
#include <algorithm>
#include <cmath>
#include <numeric>

namespace {

inline void check_double_matrix(const bigmemory::BigMatrix &mat) {
  if (mat.matrix_type() != 8) {
    stop("Only double-precision big.matrix objects are supported");
  }
}

inline NumericVector compute_xt_vec(bigmemory::MatrixAccessor<double> &acc,
                                    std::size_t n, std::size_t p,
                                    const NumericVector &vec) {
  NumericVector result(p);
  for (std::size_t col = 0; col < p; ++col) {
    const double *col_ptr = acc[col];
    double sum = 0.0;
    for (std::size_t row = 0; row < n; ++row) {
      sum += col_ptr[row] * vec[row];
    }
    result[col] = sum;
  }
  return result;
}

inline NumericVector compute_x_vec(bigmemory::MatrixAccessor<double> &acc,
                                   std::size_t n, std::size_t p,
                                   const NumericVector &weights) {
  NumericVector result(n);
  for (std::size_t col = 0; col < p; ++col) {
    const double *col_ptr = acc[col];
    double coeff = weights[col];
    if (coeff == 0.0) {
      continue;
    }
    for (std::size_t row = 0; row < n; ++row) {
      result[row] += col_ptr[row] * coeff;
    }
  }
  return result;
}

inline NumericVector deflate_xt(const NumericMatrix &Tprev,
                                const NumericMatrix &Pprev,
                                const NumericVector &vec, int k) {
  if (k == 0) {
    return NumericVector(vec);
  }
  NumericVector result = clone(vec);
  NumericVector t_cross(k);
  std::size_t n = Tprev.nrow();
  for (int comp = 0; comp < k; ++comp) {
    double sum = 0.0;
    for (std::size_t row = 0; row < n; ++row) {
      sum += Tprev(row, comp) * vec[row];
    }
    t_cross[comp] = sum;
  }
  std::size_t p = Pprev.nrow();
  for (std::size_t col = 0; col < p; ++col) {
    double adjust = 0.0;
    for (int comp = 0; comp < k; ++comp) {
      adjust += Pprev(col, comp) * t_cross[comp];
    }
    result[col] -= adjust;
  }
  return result;
}

inline NumericVector deflate_xw(const NumericMatrix &Tprev,
                                const NumericMatrix &Pprev,
                                const NumericVector &vec, int k) {
  if (k == 0) {
    return NumericVector(vec);
  }
  NumericVector result = clone(vec);
  NumericVector p_cross(k);
  std::size_t p = Pprev.nrow();
  for (int comp = 0; comp < k; ++comp) {
    double sum = 0.0;
    for (std::size_t col = 0; col < p; ++col) {
      sum += Pprev(col, comp) * vec[col];
    }
    p_cross[comp] = sum;
  }
  std::size_t n = Tprev.nrow();
  for (std::size_t row = 0; row < n; ++row) {
    double adjust = 0.0;
    for (int comp = 0; comp < k; ++comp) {
      adjust += Tprev(row, comp) * p_cross[comp];
    }
    result[row] -= adjust;
  }
  return result;
}

inline NumericVector deflate_Yt(const NumericMatrix &Tprev,
                                const NumericMatrix &Qprev,
                                const NumericVector &vec, int k) {
  if (k == 0) {
    return NumericVector(vec);
  }
  NumericVector result = clone(vec);
  NumericVector t_cross(k);
  std::size_t n = Tprev.nrow();
  for (int comp = 0; comp < k; ++comp) {
    double sum = 0.0;
    for (std::size_t row = 0; row < n; ++row) {
      sum += Tprev(row, comp) * vec[row];
    }
    t_cross[comp] = sum;
  }
  std::size_t q = Qprev.nrow();
  for (std::size_t col = 0; col < q; ++col) {
    double adjust = 0.0;
    for (int comp = 0; comp < k; ++comp) {
      adjust += Qprev(col, comp) * t_cross[comp];
    }
    result[col] -= adjust;
  }
  return result;
}

inline NumericVector deflate_Yq(const NumericMatrix &Tprev,
                                const NumericMatrix &Qprev,
                                const NumericVector &vec, int k) {
  if (k == 0) {
    return NumericVector(vec);
  }
  NumericVector result = clone(vec);
  NumericVector q_cross(k);
  std::size_t q = Qprev.nrow();
  for (int comp = 0; comp < k; ++comp) {
    double sum = 0.0;
    for (std::size_t col = 0; col < q; ++col) {
      sum += Qprev(col, comp) * vec[col];
    }
    q_cross[comp] = sum;
  }
  std::size_t n = Tprev.nrow();
  for (std::size_t row = 0; row < n; ++row) {
    double adjust = 0.0;
    for (int comp = 0; comp < k; ++comp) {
      adjust += Tprev(row, comp) * q_cross[comp];
    }
    result[row] -= adjust;
  }
  return result;
}

inline NumericVector initial_u(const NumericMatrix &Y,
                               const NumericMatrix &Tprev,
                               const NumericMatrix &Qprev,
                               int k) {
  std::size_t n = Y.nrow();
  std::size_t q = Y.ncol();
  NumericVector candidate(n);
  for (std::size_t col = 0; col < q; ++col) {
    for (std::size_t row = 0; row < n; ++row) {
      candidate[row] = Y(row, col);
    }
    if (k > 0) {
      for (int comp = 0; comp < k; ++comp) {
        double qval = Qprev(col, comp);
        if (qval == 0.0) {
          continue;
        }
        for (std::size_t row = 0; row < n; ++row) {
          candidate[row] -= Tprev(row, comp) * qval;
        }
      }
    }
    double norm = std::sqrt(std::inner_product(candidate.begin(), candidate.end(), candidate.begin(), 0.0));
    if (norm > 0.0) {
      return candidate;
    }
  }
  stop("Unable to initialise response scores; Y appears to be constant");
}

} // anonymous namespace

// [[Rcpp::export]]
Rcpp::List pls_big_cpp(SEXP xpMat, const NumericMatrix &Y, int ncomp,
                       double tol = 1e-6, int max_iter = 500) {
  if (ncomp <= 0) {
    stop("ncomp must be positive");
  }
  if (tol <= 0) {
    stop("tol must be positive");
  }
  if (max_iter <= 0) {
    stop("max_iter must be positive");
  }
  XPtr<bigmemory::BigMatrix> xp(xpMat);
  bigmemory::BigMatrix &mat = *xp;
  check_double_matrix(mat);
  std::size_t n = mat.nrow();
  std::size_t p = mat.ncol();
  std::size_t q = Y.ncol();
  if (Y.nrow() != static_cast<int>(n)) {
    stop("Number of rows in X and Y must match");
  }
  if (n == 0 || p == 0 || q == 0) {
    stop("X and Y must have positive dimensions");
  }
  int comps = std::min<int>(ncomp, std::min<int>(n, p));
  bigmemory::MatrixAccessor<double> accessor(mat);

  NumericMatrix Tmat(n, comps);
  NumericMatrix Umat(n, comps);
  NumericMatrix Wmat(p, comps);
  NumericMatrix Pmat(p, comps);
  NumericMatrix Qmat(q, comps);
  NumericVector B(comps);

  for (int k = 0; k < comps; ++k) {
    NumericVector u = initial_u(Y, Tmat, Qmat, k);
    NumericVector tvec(n);
    NumericVector wvec(p);
    NumericVector qvec(q);

    NumericVector u_old(n);
    for (int iter = 0; iter < max_iter; ++iter) {
      std::copy(u.begin(), u.end(), u_old.begin());

      NumericVector xtu = compute_xt_vec(accessor, n, p, u);
      xtu = deflate_xt(Tmat, Pmat, xtu, k);
      double wnorm = std::sqrt(std::inner_product(xtu.begin(), xtu.end(), xtu.begin(), 0.0));
      if (wnorm == 0.0) {
        stop("Encountered zero variance direction when computing weights");
      }
      for (std::size_t idx = 0; idx < p; ++idx) {
        wvec[idx] = xtu[idx] / wnorm;
      }

      NumericVector Xw = compute_x_vec(accessor, n, p, wvec);
      Xw = deflate_xw(Tmat, Pmat, Xw, k);
      std::copy(Xw.begin(), Xw.end(), tvec.begin());
      double tnorm2 = std::inner_product(tvec.begin(), tvec.end(), tvec.begin(), 0.0);
      if (tnorm2 == 0.0) {
        stop("Score vector has zero norm; cannot continue");
      }

      NumericVector ytt(q);
      for (std::size_t col = 0; col < q; ++col) {
        double sum = 0.0;
        for (std::size_t row = 0; row < n; ++row) {
          sum += Y(row, col) * tvec[row];
        }
        ytt[col] = sum;
      }
      ytt = deflate_Yt(Tmat, Qmat, ytt, k);
      for (std::size_t idx = 0; idx < q; ++idx) {
        qvec[idx] = ytt[idx] / tnorm2;
      }

      NumericVector u_new(n);
      for (std::size_t col = 0; col < q; ++col) {
        double coeff = qvec[col];
        if (coeff == 0.0) {
          continue;
        }
        for (std::size_t row = 0; row < n; ++row) {
          u_new[row] += Y(row, col) * coeff;
        }
      }
      u_new = deflate_Yq(Tmat, Qmat, u_new, k);
      double qnorm2 = std::inner_product(qvec.begin(), qvec.end(), qvec.begin(), 0.0);
      if (qnorm2 > 0.0) {
        for (std::size_t row = 0; row < n; ++row) {
          u[row] = u_new[row] / qnorm2;
        }
      } else {
        std::copy(u_new.begin(), u_new.end(), u.begin());
      }

      double diff = 0.0;
      for (std::size_t row = 0; row < n; ++row) {
        double delta = u[row] - u_old[row];
        diff += delta * delta;
      }
      if (std::sqrt(diff) < tol) {
        break;
      }
      if (iter == max_iter - 1) {
        warning("NIPALS iterations did not converge within max_iter");
      }
    }

    NumericVector xtt = compute_xt_vec(accessor, n, p, tvec);
    xtt = deflate_xt(Tmat, Pmat, xtt, k);
    double tnorm2 = std::inner_product(tvec.begin(), tvec.end(), tvec.begin(), 0.0);
    for (std::size_t idx = 0; idx < p; ++idx) {
      Pmat(idx, k) = xtt[idx] / tnorm2;
      Wmat(idx, k) = wvec[idx];
    }

    double b = 0.0;
    for (std::size_t row = 0; row < n; ++row) {
      Tmat(row, k) = tvec[row];
      Umat(row, k) = u[row];
      b += tvec[row] * u[row];
    }
    for (std::size_t idx = 0; idx < q; ++idx) {
      Qmat(idx, k) = qvec[idx];
    }
    B[k] = b / tnorm2;
  }

  return List::create(_["scores"] = Tmat,
                      _["Yscores"] = Umat,
                      _["weights"] = Wmat,
                      _["loadings"] = Pmat,
                      _["Yloadings"] = Qmat,
                      _["coefficients"] = B);
}

#else  // BIGPLS_HAS_BIGMEMORY

// [[Rcpp::export]]
Rcpp::List pls_big_cpp(SEXP /*xpMat*/, const NumericMatrix & /*Y*/, int /*ncomp*/,
                       double /*tol*/ , int /*max_iter*/) {
  stop("bigmemory headers are not available. Install the bigmemory package "
       "to use pls_big(); the R fallback matrixpls_stream_bigmatrix() "
       "remains available.");
}

#endif  // BIGPLS_HAS_BIGMEMORY

