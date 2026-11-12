# Internal optimisation routines for synthetic control data fusion.
#
# Notation (paper): reference-domain outcomes F (units x T), reference
# covariates Z (units x d_r), target-domain covariates X (units x d_t).
# The donor weights w live on the probability simplex.
#
#   NSE_F(w) = ||F_1 - F_0' w||^2 / T
#   NSE_Z(w) = ||Z_1 - Z_0' w||^2 / d_r
#   NSE_X(w) = ||X_1 - X_0' w||^2 / d_t
#
# Stage 1 (baselines): NSE_Z* and NSE_X* are the smallest attainable
#   covariate discrepancies on the simplex.
# Stage 2 (for each budget b = (b_F, b_Z, b_X)):
#   w(b) = argmin_w  b_F NSE_F(w) + b_Z NSE_Z(w) + b_X NSE_X(w)
#          s.t.      NSE_Z(w) <= (1 + eta_Z)(1 + NSE_Z*) - 1
#                    NSE_X(w) <= (1 + eta_X)(1 + NSE_X*) - 1
#                    w >= 0, sum(w) = 1
# Stage 3: the selected budget minimises NSE_F(w(b)) over the grid.
#
# The ECOS backend solves each problem as a second-order cone program (the
# objective is replaced by its square root, which has the same minimiser).
# The Gurobi backend reproduces the original research implementation
# (kernlab::ipop baselines + Gurobi QCQP) line by line.

# ---------------------------------------------------------------------------
# Problem container
# ---------------------------------------------------------------------------

# Bundle the treated/donor blocks for one fitting problem.
.fusion_problem <- function(F_treated, F_control, Z_treated, Z_control,
                            X_treated, X_control) {
  F_control <- as.matrix(F_control)
  Z_control <- as.matrix(Z_control)
  X_control <- as.matrix(X_control)
  storage.mode(F_control) <- "double"
  storage.mode(Z_control) <- "double"
  storage.mode(X_control) <- "double"
  list(
    F1 = as.numeric(F_treated), F0 = F_control,
    Z1 = as.numeric(Z_treated), Z0 = Z_control,
    X1 = as.numeric(X_treated), X0 = X_control,
    J = nrow(F_control), T = ncol(F_control),
    dr = ncol(Z_control), dt = ncol(X_control)
  )
}

# Normalised squared error of weights w for one block.
.nse <- function(w, treated, control) {
  sum((treated - as.numeric(crossprod(control, w)))^2) / length(treated)
}

# Zero-out tiny weights exactly as in the research code (no renormalisation).
.threshold_weights <- function(w, threshold) {
  w <- as.numeric(w)
  w[w < threshold] <- 0
  w
}

# ---------------------------------------------------------------------------
# ECOS helpers
# ---------------------------------------------------------------------------

.ecos_control <- function() {
  ECOSolveR::ecos.control(
    maxit = 200L, feastol = 1e-8, abstol = 1e-8, reltol = 1e-8
  )
}

# ECOSolveR equilibrates the problem data in place and only approximately
# restores it afterwards, so any R object handed to it can be altered (for
# example a literal `1` becomes 0.9999999999999999, which then corrupts every
# use of that constant). Always pass freshly allocated copies.
.ecos_solve <- function(c, G, h, dims, A, b, control) {
  G@x <- G@x + 0
  G@i <- G@i + 0L
  G@p <- G@p + 0L
  A@x <- A@x + 0
  A@i <- A@i + 0L
  A@p <- A@p + 0L
  ECOSolveR::ECOS_csolve(
    c = as.numeric(c) + 0, G = G, h = as.numeric(h) + 0,
    dims = lapply(dims, function(d) as.integer(d) + 0L),
    A = A, b = as.numeric(b) + 0, control = control
  )
}

.as_dgc <- function(M) {
  methods::as(methods::as(Matrix::Matrix(M, sparse = TRUE), "CsparseMatrix"),
              "generalMatrix")
}

# ECOS exit codes 0 (optimal) and 10 (optimal up to reduced accuracy).
.ecos_ok <- function(res) {
  !is.null(res$retcodes) && res$retcodes[["exitFlag"]] %in% c(0L, 10L)
}

# Simplex-constrained least squares:  min ||C' w - t||  s.t. w in simplex.
.ecos_simplex_ls <- function(treated, control) {
  J <- nrow(control)
  d <- length(treated)
  # variables x = (w, s); minimise s
  G <- rbind(
    cbind(-diag(J), 0),                   # w >= 0
    c(rep(0, J), -1),                     # SOC head: s
    cbind(-t(control), 0)                 # SOC tail: C' w - t
  )
  h <- c(rep(0, J), 0, -treated)
  res <- .ecos_solve(
    c = c(rep(0, J), 1),
    G = .as_dgc(G),
    h = h,
    dims = list(l = J, q = d + 1L, e = 0L),
    A = .as_dgc(matrix(c(rep(1, J), 0), nrow = 1L)),
    b = 1,
    control = .ecos_control()
  )
  if (!.ecos_ok(res)) {
    stop("ECOS failed to solve the baseline simplex least-squares problem.",
         call. = FALSE)
  }
  res$x[seq_len(J)]
}

# Build the (budget-independent part of the) SOCP for the constrained fit.
# Returns a closure that solves the problem for any budget vector.
.ecos_qcqp_factory <- function(prob, rad_Z, rad_X) {
  J <- prob$J; T <- prob$T; dr <- prob$dr; dt <- prob$dt
  # When there are more reference periods than donors, replace the F block by
  # its QR factor: ||F_0' w - F_1||^2 = ||R w - Q'F_1||^2 + const, so the
  # minimiser is unchanged while the cone gets smaller.
  F_rows <- t(prob$F0)
  F_rhs <- prob$F1
  if (T > J) {
    q <- qr(F_rows)
    R <- qr.R(q)[, order(q$pivot), drop = FALSE]
    F_rows <- R
    F_rhs <- qr.qty(q, prob$F1)[seq_len(nrow(R))]
  }
  nF <- nrow(F_rows)
  m_obj <- nF + dr + dt
  # Variables x = (w_1, ..., w_J, s).  Cones, in order:
  #   linear:  w >= 0                                     (J rows)
  #   SOC 1:   || A_b w - a_b || <= s                       (1 + m_obj rows)
  #   SOC 2:   || Z_0' w - Z_1 || <= rad_Z                  (1 + dr rows)
  #   SOC 3:   || X_0' w - X_1 || <= rad_X                  (1 + dt rows)
  G <- rbind(
    cbind(-diag(J), 0),
    c(rep(0, J), -1),
    cbind(-F_rows, 0),
    cbind(-t(prob$Z0), 0),
    cbind(-t(prob$X0), 0),
    0,
    cbind(-t(prob$Z0), 0),
    0,
    cbind(-t(prob$X0), 0)
  )
  h_base <- c(rep(0, J), 0, -F_rhs, -prob$Z1, -prob$X1,
              rad_Z, -prob$Z1, rad_X, -prob$X1)
  # group label of every row: 0 = unscaled, 1 = F, 2 = Z, 3 = X (objective)
  row_group <- c(rep(0L, J + 1L), rep(1L, nF), rep(2L, dr), rep(3L, dt),
                 rep(0L, dr + dt + 2L))
  Gs <- .as_dgc(G)
  base_x <- Gs@x
  nz_group <- row_group[Gs@i + 1L]
  dims <- list(l = J, q = c(1L + m_obj, 1L + dr, 1L + dt), e = 0L)
  A <- .as_dgc(matrix(c(rep(1, J), 0), nrow = 1L))
  cvec <- c(rep(0, J), 1)
  ctrl <- .ecos_control()

  function(b_F, b_Z, b_X) {
    sc <- c(1, sqrt(b_F / T), sqrt(b_Z / dr), sqrt(b_X / dt))
    Gs@x <- base_x * sc[nz_group + 1L]
    h <- h_base * sc[row_group + 1L]
    res <- .ecos_solve(c = cvec, G = Gs, h = h, dims = dims, A = A, b = 1,
                       control = ctrl)
    if (!.ecos_ok(res)) return(NULL)
    res$x[seq_len(J)]
  }
}

# ---------------------------------------------------------------------------
# Gurobi helpers (mirror the research implementation exactly)
# ---------------------------------------------------------------------------

.require_gurobi <- function() {
  if (!requireNamespace("gurobi", quietly = TRUE)) {
    stop(
      "solver = \"gurobi\" requires the 'gurobi' R package and a valid ",
      "Gurobi licence. Install it from your Gurobi distribution, or use ",
      "solver = \"ecos\" (the default).",
      call. = FALSE
    )
  }
  if (!requireNamespace("kernlab", quietly = TRUE)) {
    stop("solver = \"gurobi\" also requires the 'kernlab' package.",
         call. = FALSE)
  }
}

# kernlab::ipop baseline as in optimize_w_ipop() of the research code.
.ipop_simplex_ls <- function(treated, control, margin = 1e-4, sigf = 5,
                             bound = 10) {
  J <- nrow(control)
  d <- length(treated)
  Xt <- matrix(treated)
  Xc <- t(control)
  V <- diag(rep(1 / d, d), nrow = d)
  H <- t(Xc) %*% V %*% Xc
  cc <- -1 * t(Xt) %*% V %*% Xc
  res <- kernlab::ipop(
    c = cc, H = H, A = t(rep(1, J)), b = 1, l = rep(0, J), u = rep(1, J),
    r = 0, bound = bound, margin = margin, sigf = sigf
  )
  as.numeric(kernlab::primal(res))
}

.gurobi_qcqp_factory <- function(prob, NSE_Z_baseline, NSE_X_baseline,
                                 eta_Z, eta_X, time_limit = 100) {
  J <- prob$J; T <- prob$T; dr <- prob$dr; dt <- prob$dt
  F0 <- prob$F0; Z0 <- prob$Z0; X0 <- prob$X0
  F1 <- prob$F1; Z1 <- prob$Z1; X1 <- prob$X1
  FF <- F0 %*% t(F0); ZZ <- Z0 %*% t(Z0); XX <- X0 %*% t(X0)
  qF <- t(F1) %*% t(F0); qZ <- t(Z1) %*% t(Z0); qX <- t(X1) %*% t(X0)
  Q_Z <- ZZ
  L_Z <- -2 * t(matrix(Z1)) %*% t(Z0)
  rhs_Z <- dr * ((1 + eta_Z) * (1 + NSE_Z_baseline) - 1) - sum(Z1^2)
  Q_X <- XX
  L_X <- -2 * t(matrix(X1)) %*% t(X0)
  rhs_X <- dt * ((1 + eta_X) * (1 + NSE_X_baseline) - 1) - sum(X1^2)
  params <- list(OutputFlag = 0, TimeLimit = time_limit)

  function(b_F, b_Z, b_X) {
    model <- list(
      Q = b_F / T * FF + b_Z / dr * ZZ + b_X / dt * XX,
      obj = as.vector(-2 * (b_F / T * qF + b_Z / dr * qZ + b_X / dt * qX)),
      modelsense = "min",
      quadcon = list(
        list(Qc = Q_Z, q = as.vector(L_Z), rhs = rhs_Z),
        list(Qc = Q_X, q = as.vector(L_X), rhs = rhs_X)
      ),
      A = matrix(rep(1, J), nrow = 1L, ncol = J),
      rhs = c(1) + 0,
      sense = "=",
      lb = numeric(J),
      ub = rep(1, J)
    )
    res <- tryCatch(gurobi::gurobi(model, params = params),
                    error = function(e) NULL)
    if (is.null(res) || !is.list(res) || is.null(res$status) ||
        res$status != "OPTIMAL") {
      return(NULL)
    }
    x <- as.numeric(res$x)
    if (length(x) != J || anyNA(x)) return(NULL)
    x
  }
}

# ---------------------------------------------------------------------------
# Stage 1: baselines
# ---------------------------------------------------------------------------

.baseline_nse <- function(prob, solver, threshold) {
  solve_ls <- switch(solver,
                     ecos = .ecos_simplex_ls,
                     gurobi = .ipop_simplex_ls)
  wZ <- .threshold_weights(solve_ls(prob$Z1, prob$Z0), threshold)
  wX <- .threshold_weights(solve_ls(prob$X1, prob$X0), threshold)
  c(Z = .nse(wZ, prob$Z1, prob$Z0), X = .nse(wX, prob$X1, prob$X0))
}

# ---------------------------------------------------------------------------
# Stages 2-3: constrained fits over the budget grid and selection
# ---------------------------------------------------------------------------

# Fit the data-fusion weights for one treated unit / donor pool.
#
# Returns a list with the selected weights, budget, criterion values and
# baselines; `weights` is NULL when no budget produced a feasible solution.
.fit_fusion_weights <- function(prob, budgets, eta_Z, eta_X, solver,
                                threshold, keep_path = FALSE,
                                baseline = NULL) {
  if (is.null(baseline)) baseline <- .baseline_nse(prob, solver, threshold)
  rad2_Z <- prob$dr * ((1 + eta_Z) * (1 + baseline[["Z"]]) - 1)
  rad2_X <- prob$dt * ((1 + eta_X) * (1 + baseline[["X"]]) - 1)

  solve_b <- switch(
    solver,
    ecos = .ecos_qcqp_factory(prob, sqrt(max(rad2_Z, 0)),
                              sqrt(max(rad2_X, 0))),
    gurobi = .gurobi_qcqp_factory(prob, baseline[["Z"]], baseline[["X"]],
                                  eta_Z, eta_X)
  )

  nb <- nrow(budgets)
  path_loss <- rep(NA_real_, nb)
  best_loss <- Inf
  best_k <- NA_integer_
  best_w <- NULL
  for (k in seq_len(nb)) {
    w <- solve_b(budgets[k, 1L], budgets[k, 2L], budgets[k, 3L])
    if (is.null(w)) next
    w <- .threshold_weights(w, threshold)
    loss <- .nse(w, prob$F1, prob$F0)
    path_loss[k] <- loss
    if (!is.na(loss) && loss < best_loss) {
      best_loss <- loss
      best_k <- k
      best_w <- w
    }
  }

  out <- list(
    weights = best_w,
    budget = if (is.na(best_k)) NULL else budgets[best_k, ],
    loss = if (is.na(best_k)) NA_real_ else best_loss,
    baseline = baseline,
    n_feasible = sum(!is.na(path_loss)),
    n_budgets = nb
  )
  if (keep_path) out$path <- cbind(budgets, NSE_F = path_loss)
  out
}
