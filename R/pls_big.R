#' Partial least squares for bigmemory matrices
#'
#' @description
#' `pls_big()` fits a partial least squares (PLS) model using a NIPALS
#' algorithm implemented in C++ that operates directly on
#' [bigmemory::big.matrix] inputs. Both in-memory and file-backed matrices are
#' supported. `matrixpls_stream_bigmatrix()` is a pure R fallback that performs
#' the same computation by streaming over chunks of a file-backed
#' `big.matrix` without loading it fully into memory.
#'
#' @param X Either a `bigmemory::big.matrix`, a
#'   `bigmemory::big.matrix.descriptor`, or a character path to a delimited file
#'   that can be read with [bigmemory::read.big.matrix()].
#' @param Y Numeric response matrix with matching number of rows. Vectors are
#'   coerced to a one-column matrix.
#' @param ncomp Number of latent components to extract.
#' @param tol Convergence tolerance for the iterative updates.
#' @param max_iter Maximum number of iterations for the NIPALS inner loop.
#' @param stream Logical; when `TRUE`, force the chunk-wise R implementation
#'   via [matrixpls_stream_bigmatrix()].
#' @param num.rows.chunk Number of rows to load per chunk when
#'   `stream = TRUE`.
#' @param backingfile,backingpath,descriptorfile Optional arguments passed to
#'   [bigmemory::read.big.matrix()] when `X` is a file path.
#' @param type Storage mode to use when reading from a file path, defaulting to
#'   `"double"`.
#' @param ... Reserved for future extensions.
#'
#' @return A list containing
#'   * `scores`: X-score matrix (`T`)
#'   * `Yscores`: Y-score matrix (`U`)
#'   * `weights`: weight matrix (`W`)
#'   * `loadings`: X-loading matrix (`P`)
#'   * `Yloadings`: Y-loading matrix (`Q`)
#'   * `coefficients`: regression coefficients linking `T` and `U`
#'
#' @export
pls_big <- function(X, Y, ncomp = 2L, tol = 1e-6, max_iter = 500L,
                    stream = FALSE, num.rows.chunk = 1e6,
                    backingfile = NULL, backingpath = NULL,
                    descriptorfile = NULL, type = "double", ...) {
  Y <- as.matrix(Y)
  if (!is.numeric(Y)) {
    stop("Y must be numeric")
  }
  big_X <- resolve_big_matrix(X, backingfile = backingfile,
                              backingpath = backingpath,
                              descriptorfile = descriptorfile,
                              type = type)
  if (nrow(big_X) != nrow(Y)) {
    stop("X and Y must have matching numbers of rows")
  }
  if (stream) {
    return(matrixpls_stream_bigmatrix(big_X, Y, ncomp = ncomp, tol = tol,
                                      max_iter = max_iter,
                                      num.rows.chunk = num.rows.chunk))
  }
  addr <- bigmemory:::address(big_X)
  pls_big_cpp(addr, Y, ncomp = ncomp, tol = tol, max_iter = max_iter)
}

#' @rdname pls_big
#' @export
matrixpls_stream_bigmatrix <- function(X, Y, ncomp = 2L, tol = 1e-6,
                                       max_iter = 500L,
                                       num.rows.chunk = 1e6, ...) {
  Y <- as.matrix(Y)
  if (!is.numeric(Y)) {
    stop("Y must be numeric")
  }
  big_X <- resolve_big_matrix(X, ...)
  n <- nrow(big_X)
  p <- ncol(big_X)
  q <- ncol(Y)
  if (nrow(Y) != n) {
    stop("X and Y must have matching numbers of rows")
  }
  if (ncomp <= 0) {
    stop("ncomp must be positive")
  }
  ncomp <- min(ncomp, n, p)
  tol <- if (!is.numeric(tol) || length(tol) != 1) 1e-6 else tol
  max_iter <- if (!is.numeric(max_iter) || length(max_iter) != 1) 500L else max_iter
  if (tol <= 0) {
    stop("tol must be positive")
  }
  if (max_iter <= 0) {
    stop("max_iter must be positive")
  }
  if (num.rows.chunk <= 0) {
    stop("num.rows.chunk must be positive")
  }
  chunk_sizes <- c(0, rep(num.rows.chunk, floor(n/num.rows.chunk)))
  remainder <- n %% num.rows.chunk
  if (remainder > 0) {
    chunk_sizes <- c(chunk_sizes, remainder)
  }
  chunk_offsets <- cumsum(chunk_sizes)
  chunk_count <- length(chunk_sizes) - 1

  Tmat <- matrix(0, n, ncomp)
  Umat <- matrix(0, n, ncomp)
  Wmat <- matrix(0, p, ncomp)
  Pmat <- matrix(0, p, ncomp)
  Qmat <- matrix(0, q, ncomp)
  B <- numeric(ncomp)

  get_chunk_rows <- function(idx) {
    start <- chunk_offsets[idx] + 1
    end <- chunk_offsets[idx + 1]
    seq.int(start, end)
  }

  for (k in seq_len(ncomp)) {
    u <- initial_u_stream(Y, Tmat, Qmat, k)
    tvec <- numeric(n)
    wvec <- numeric(p)
    qvec <- numeric(q)
    u_old <- numeric(n)

    for (iter in seq_len(max_iter)) {
      u_old[] <- u
      xtu <- numeric(p)
      for (chunk in seq_len(chunk_count)) {
        rows <- get_chunk_rows(chunk)
        subX <- big_X[rows, , drop = FALSE]
        xtu <- xtu + crossprod(subX, u[rows])
      }
      if (k > 1) {
        xtu <- xtu - Pmat[, seq_len(k - 1), drop = FALSE] %*%
          crossprod(Tmat[, seq_len(k - 1), drop = FALSE], u)
      }
      wnorm <- sqrt(sum(xtu^2))
      if (wnorm == 0) {
        stop("Encountered zero variance direction when computing weights")
      }
      wvec <- as.numeric(xtu / wnorm)

      tvec <- numeric(n)
      for (chunk in seq_len(chunk_count)) {
        rows <- get_chunk_rows(chunk)
        subX <- big_X[rows, , drop = FALSE]
        tvec[rows] <- tvec[rows] + as.vector(subX %*% wvec)
      }
      if (k > 1) {
        tvec <- tvec - Tmat[, seq_len(k - 1), drop = FALSE] %*%
          (t(Pmat[, seq_len(k - 1), drop = FALSE]) %*% wvec)
      }
      tnorm2 <- sum(tvec^2)
      if (tnorm2 == 0) {
        stop("Score vector has zero norm; cannot continue")
      }

      ytt <- as.numeric(crossprod(Y, tvec))
      if (k > 1) {
        ytt <- ytt - Qmat[, seq_len(k - 1), drop = FALSE] %*%
          crossprod(Tmat[, seq_len(k - 1), drop = FALSE], tvec)
      }
      qvec <- ytt/tnorm2

      u_new <- as.numeric(Y %*% qvec)
      if (k > 1) {
        u_new <- u_new - Tmat[, seq_len(k - 1), drop = FALSE] %*%
          (t(Qmat[, seq_len(k - 1), drop = FALSE]) %*% qvec)
      }
      qnorm2 <- sum(qvec^2)
      if (qnorm2 > 0) {
        u <- u_new/qnorm2
      } else {
        u <- u_new
      }

      diff <- sqrt(sum((u - u_old)^2))
      if (diff < tol) {
        break
      }
      if (iter == max_iter) {
        warning("NIPALS iterations did not converge within max_iter")
      }
    }

    xtt <- numeric(p)
    for (chunk in seq_len(chunk_count)) {
      rows <- get_chunk_rows(chunk)
      subX <- big_X[rows, , drop = FALSE]
      xtt <- xtt + crossprod(subX, tvec[rows])
    }
    if (k > 1) {
      xtt <- xtt - Pmat[, seq_len(k - 1), drop = FALSE] %*%
        crossprod(Tmat[, seq_len(k - 1), drop = FALSE], tvec)
    }
    Pmat[, k] <- xtt/tnorm2
    Wmat[, k] <- wvec
    Tmat[, k] <- tvec
    Umat[, k] <- u
    Qmat[, k] <- qvec
    B[k] <- sum(tvec * u)/tnorm2
  }

  list(scores = Tmat, Yscores = Umat, weights = Wmat, loadings = Pmat,
       Yloadings = Qmat, coefficients = B)
}

resolve_big_matrix <- function(X, backingfile = NULL, backingpath = NULL,
                               descriptorfile = NULL, type = "double", ...) {
  if (bigmemory::is.big.matrix(X)) {
    return(X)
  }
  if (inherits(X, "big.matrix.descriptor")) {
    return(bigmemory::attach.big.matrix(X))
  }
  if (is.character(X) && length(X) == 1) {
    return(bigmemory::read.big.matrix(filename = X, sep = ",", skip = 0,
                                      header = TRUE, backingfile = backingfile,
                                      backingpath = backingpath,
                                      descriptorfile = descriptorfile,
                                      type = type, ...))
  }
  stop("Unsupported X specification for big.matrix input")
}

initial_u_stream <- function(Y, Tmat, Qmat, k) {
  n <- nrow(Y)
  q <- ncol(Y)
  if (k <= 1) {
    candidates <- seq_len(max(1, q))
  } else {
    candidates <- seq_len(q)
  }
  for (col in candidates) {
    u <- Y[, col]
    if (k > 1) {
      u <- u - Tmat[, seq_len(k - 1), drop = FALSE] %*% Qmat[col, seq_len(k - 1), drop = FALSE]
    }
    if (sqrt(sum(u^2)) > 0) {
      return(as.numeric(u))
    }
  }
  stop("Unable to initialise response scores; Y appears to be constant")
}

