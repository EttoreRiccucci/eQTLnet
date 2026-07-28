## =============================================================================
## eqtlnet_regression.R
## -----------------------------------------------------------------------------
## Purpose : Infer gene regulatory networks (GRNs)
##           using six regularisation-based regression approaches:
##
##   1. Classical LASSO              			(glmnet, alpha = 1, uniform penalty)
##   2. Classical Ridge Regression   			(glmnet, alpha = 0, uniform penalty)
##   3. Weighted LASSO   aka eQTLnet-LASSO              (glmnet, alpha = 1, penalty.factor = p-values)
##   4. Weighted Ridge Regression aka eQTLnet-RR  	(glmnet, alpha = 0, penalty.factor = p-values)
##   5. Mixed OLS–LASSO  aka "eQTLnet-LASSO-OLS" (two separate fits per TG:   OLS via lm() on eQTL-supported TFs ,
##                         LASSO via cv.glmnet on unsupported TFs (weight == 1);
##   6. Mixed OLS–Ridge aka "eQTLnet-RR-OLS"   (same split as abpove, but with Ridge regression instead of LASSO for the
##                         unsupported TFs)
##
## Framework :
##   - For each target gene (TG), all TFs are predictors.
##     When the TG is itself a TF, it is removed from the predictor set because we do not address self-regulation.
##   - Optimal lambda is chosen by k-fold cross-validation (cv.glmnet).
##   - The final model is refit on ALL samples at lambda.min.
##   - Results are saved as coefficient matrices (TG × TF)
##
## Inputs  :
##   TGsexpr.csv        — (accessions × target genes) expression matrix (including transcription factors and non-TFs)
##   TFsexpr.csv        — (accessions × TFs) expression matrix for transcription factors only, same format
##   w_matrix_updated.csv — (TGs × TFs) matrix of eQTL p-values;
##                          row names = TG names, col names = TF names;
##                          1 = for the remaining TF as explain in the manuscript
##
## Outputs : Six RDS files saved to ./results/, one per method:
##   coef_LASSO.rds, coef_ridge.rds,
##   coef_weighted_LASSO.rds, coef_weighted_ridge.rds,
##   coef_mixed_LASSO.rds,    coef_mixed_ridge.rds
##
##   Each RDS contains a named list with:
##     $coef_matrix  : matrix (n_TG × n_TF), coefficient for each TF→TG edge
##     $lambda_min   : named vector, selected lambda per TG
##     $cv_errors    : list of cv.glmnet objects (one per TG), for diagnostics
##
## Dependencies : glmnet, doParallel, foreach
##               install.packages(c("glmnet", "doParallel", "foreach"))
##
## Hardware note (the computer used):
##   16 physical cores / 1 socket / 2 threads per core = 22 logical CPUs
##   N_CORES is set to 14 (physical cores − 2) — the optimal setting for
##   CPU-bound R work; hyperthreaded logical cores beyond the physical count
##   add little for glmnet and can cause memory-bandwidth contention.
##   Parallelism is across TARGET GENES (one TG per worker, all 6 models
##   fitted inside that worker) rather than across CV folds, which gives
##   far lower overhead.
## =============================================================================


# ── 0. Packages ───────────────────────────────────────────────────────────────
for (pkg in c("glmnet", "doParallel", "foreach")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Please install the '%s' package: install.packages('%s')", pkg, pkg))
  }
}
library(glmnet)
library(doParallel)
library(foreach)


# ═════════════════════════════════════════════════════════════════════════════
# ── 1. USER SETTINGS  (edit this block based on your requirements) 
# ═════════════════════════════════════════════════════════════════════════════

TG_EXPR_FILE   <- "TGsexpr.csv"         # TG expression matrix (accessions × TGs)
TF_EXPR_FILE   <- "TFsexpr.csv"         # TF expression matrix (accessions × TFs)
W_MATRIX_FILE  <- "w_matrix_updated.csv" # eQTL p-value weight matrix (TGs × TFs)

OUT_DIR        <- "results"             # where coefficient matrices are saved

N_FOLDS        <- 10                    # k for k-fold CV (10 for n=663 arabidopsis accessions)
LAMBDA_RULE    <- "lambda.min"          # "lambda.min" or "lambda.1se"
MIXED_CUTOFF   <- 1                     # p-value threshold for mixed OLS models:
#   weight == MIXED_CUTOFF → LASSO or Ridge
#   weight <  MIXED_CUTOFF → OLS (eQTL-supported)

N_CORES        <- 14L                   # physical cores to use for parallel
# execution across target genes.
# This machine: 16 physical cores,
# 22 logical (hyperthreaded).
# Using 14 = 16 physical − 2 keeps the
# system responsive and avoids memory-
# bandwidth contention from hyperthreading.
# Adjust downward on a shared or
# memory-constrained system.

MASTER_SEED    <- 2024L                 # reproducibility seed for CV fold assignment

# ═════════════════════════════════════════════════════════════════════════════


# ── 2. Helper  
start_cluster <- function(n_cores) {
  n_use <- min(n_cores, parallel::detectCores(logical = FALSE))
  message(sprintf(
    "Starting cluster: %d workers (physical cores available: %d, logical: %d)",
    n_use,
    parallel::detectCores(logical = FALSE),
    parallel::detectCores(logical = TRUE)
  ))
  cl <- parallel::makePSOCKcluster(n_use)
  doParallel::registerDoParallel(cl)
  parallel::clusterEvalQ(cl, library(glmnet))
  message(sprintf("Cluster ready: %d workers registered.", n_use))
  cl
}

# Fit one cv.glmnet and extract coefficients at lambda.min.
# Timing is measured with proc.time() inside the worker to reflects
# actual CPU work regardless of parallel scheduling
fit_one_tg <- function(y, X, alpha, penalty_factors = NULL,
                       n_folds = 10, lambda_rule = "lambda.min",
                       seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  cv_args <- list(
    x            = X,
    y            = y,
    alpha        = alpha,
    nfolds       = n_folds,
    standardize  = TRUE, 
    intercept    = TRUE
  )
  if (!is.null(penalty_factors)) {
    cv_args$penalty.factor <- penalty_factors
  }
  
  # Time the CV fitting 
  t0     <- proc.time()["elapsed"]
  cv_fit <- do.call(cv.glmnet, cv_args)
  elapsed <- as.numeric(proc.time()["elapsed"] - t0)
  
  lam_best <- cv_fit[[lambda_rule]]
  
  # Extract coefficients at best lambda (drop intercept row)
  coef_vec <- as.numeric(coef(cv_fit, s = lam_best))[-1]
  names(coef_vec) <- colnames(X)
  
  list(coefs = coef_vec, lambda_min = lam_best, elapsed_sec = elapsed)
}

# Fit the Mixed OLS–LASSO or Mixed OLS–Ridge model for one TG.
fit_mixed_tg <- function(y, X_LASSO, X_ols, tfs_LASSO, tfs_ols, all_tfs,
                         alpha, n_folds = 10, lambda_rule = "lambda.min",
                         seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  t0           <- proc.time()["elapsed"]
  beta         <- setNames(numeric(length(all_tfs)), all_tfs)
  lambda_best  <- NA_real_
  
  if (ncol(X_LASSO) > 0) {
    cv_fit      <- cv.glmnet(X_LASSO, y, alpha = alpha,
                             nfolds = n_folds, standardize = FALSE,
                             intercept = FALSE, family = "gaussian")
    lambda_best <- cv_fit[[lambda_rule]]
    # Refit at best lambda to get clean coefficients (same as your original code)
    LASSO_coeffs <- as.numeric(
      glmnet(X_LASSO, y, alpha = alpha, lambda = lambda_best,
             standardize = FALSE, intercept = FALSE,
             family = "gaussian")$beta
    )
    beta[tfs_LASSO] <- LASSO_coeffs
  }
  
  if (ncol(X_ols) > 0) {
    ols_fit         <- lm(y ~ 0 + X_ols)   # no intercept, matches your code
    beta[tfs_ols]   <- coef(ols_fit)
  }
  
  elapsed <- as.numeric(proc.time()["elapsed"] - t0)
  list(coefs = beta, lambda_min = lambda_best, elapsed_sec = elapsed)
}


# ── 3. Load data 
message("Loading expression and weight matrices ...")

# Expression matrices
TG_expr <- as.matrix(read.csv(TG_EXPR_FILE, row.names = 1, check.names = FALSE))
TF_expr <- as.matrix(read.csv(TF_EXPR_FILE, row.names = 1, check.names = FALSE))

# Weight matrix: rows = TGs, columns = TFs
W_mat   <- as.matrix(read.csv(W_MATRIX_FILE, row.names = 1, check.names = FALSE))

# Align accession order (such that rows match between TG and TF matrices)
common_acc  <- intersect(rownames(TG_expr), rownames(TF_expr))
stopifnot("No common accessions between TG and TF matrices" = length(common_acc) > 0)
TG_expr <- TG_expr[common_acc, , drop = FALSE]
TF_expr <- TF_expr[common_acc, , drop = FALSE]

# Working gene sets
all_TFs <- colnames(TF_expr)     # all TFs (predictors when not self)
all_TGs <- colnames(TG_expr)     # all TGs (response variables)

# Align weight matrix to working gene sets
# Rows = TGs present in W_mat ∩ all_TGs; cols = TFs present in W_mat ∩ all_TFs
TGs_with_weights <- intersect(all_TGs, rownames(W_mat))
TFs_with_weights <- intersect(all_TFs, colnames(W_mat))
W_mat <- W_mat[TGs_with_weights, TFs_with_weights, drop = FALSE]

message(sprintf("  Accessions (n)  : %d", length(common_acc)))
message(sprintf("  Target genes    : %d  (weight-matrix coverage: %d)",
                length(all_TGs), length(TGs_with_weights)))
message(sprintf("  Transcription factors : %d  (weight-matrix coverage: %d)",
                length(all_TFs), length(TFs_with_weights)))


# ── 4. Start parallel cluster 
cl <- start_cluster(N_CORES)
# Register a clean shutdown so the cluster is always stopped, even on error
on.exit(parallel::stopCluster(cl), add = TRUE)


# ── 6. Parallel loop over target genes
message(sprintf("\nFitting models for %d target genes on %d workers ...",
                length(all_TGs), N_CORES))
message(sprintf("  k-fold CV     : %d folds", N_FOLDS))
message(sprintf("  Lambda rule   : %s", LAMBDA_RULE))
message(sprintf("  Mixed cutoff  : weight < %g → OLS; weight == %g → LASSO/Ridge",
                MIXED_CUTOFF, MIXED_CUTOFF))

t_start <- proc.time()

# start_cluster 
results_list <- foreach(
  tg_idx = seq_along(all_TGs),
  .packages  = "glmnet",
  .export    = c("fit_one_tg", "fit_mixed_tg"),  
  .errorhandling = "pass"          
) %dopar% {

  tg <- all_TGs[tg_idx]

  # ── Build predictor matrix for this TG 
  # Remove this TG from predictors if it is itself a TF (no self-regulation)
  tfs_for_tg <- setdiff(all_TFs, tg)
  X          <- TF_expr[, tfs_for_tg, drop = FALSE]
  y          <- TG_expr[, tg]

  # ── Retrieve weight vector for this TG 
  if (tg %in% rownames(W_mat)) {
    tfs_in_w            <- intersect(tfs_for_tg, colnames(W_mat))
    pf_full             <- setNames(rep(1, length(tfs_for_tg)), tfs_for_tg)
    pf_full[tfs_in_w]   <- W_mat[tg, tfs_in_w]  # fill eQTL p-values
    # TFs absent from W_mat keep pf = 1 (maximum penalty)
  } else {
    pf_full <- setNames(rep(1, length(tfs_for_tg)), tfs_for_tg)
  }

  tfs_ols <- names(pf_full)[pf_full <  MIXED_CUTOFF]
  tfs_reg <- names(pf_full)[pf_full == MIXED_CUTOFF]

  X_ols   <- X[, tfs_ols, drop = FALSE]   # OLS predictor block
  X_reg   <- X[, tfs_reg, drop = FALSE]   # LASSO/Ridge predictor block

  # Per-TG seed: ensures reproducible CV fold splits across workers
  tg_seed <- MASTER_SEED + tg_idx

  # ── Fit the models 

  # 1. Classical LASSO  (uniform penalty)
  res_LASSO <- fit_one_tg(y, X, alpha = 1,
                          n_folds = N_FOLDS, lambda_rule = LAMBDA_RULE,
                          seed = tg_seed)

  # 2. Classical Ridge  (alpha = 0, uniform penalty)
  res_ridge <- fit_one_tg(y, X, alpha = 0,
                          n_folds = N_FOLDS, lambda_rule = LAMBDA_RULE,
                          seed = tg_seed)

  # 3. Weighted LASSO   (penalty.factor = eQTL p-values; smaller pval receive less shrinkage)
  res_wLASSO <- fit_one_tg(y, X, alpha = 1,
                           penalty_factors = pf_full,
                           n_folds = N_FOLDS, lambda_rule = LAMBDA_RULE,
                           seed = tg_seed)

  # 4. Weighted Ridge   
  res_wridge <- fit_one_tg(y, X, alpha = 0,
                           penalty_factors = pf_full,
                           n_folds = N_FOLDS, lambda_rule = LAMBDA_RULE,
                           seed = tg_seed)

  # 5. Mixed OLS–LASSO:
  res_mixed_LASSO <- fit_mixed_tg(y,
                                  X_LASSO   = X_reg,  X_ols = X_ols,
                                  tfs_LASSO = tfs_reg, tfs_ols = tfs_ols,
                                  all_tfs   = tfs_for_tg,
                                  alpha     = 1,
                                  n_folds   = N_FOLDS, lambda_rule = LAMBDA_RULE,
                                  seed      = tg_seed)

  # 6. Mixed OLS–Ridge:
  res_mixed_ridge <- fit_mixed_tg(y,
                                  X_LASSO   = X_reg,  X_ols = X_ols,
                                  tfs_LASSO = tfs_reg, tfs_ols = tfs_ols,
                                  all_tfs   = tfs_for_tg,
                                  alpha     = 0,
                                  n_folds   = N_FOLDS, lambda_rule = LAMBDA_RULE,
                                  seed      = tg_seed)

  # ── Return result for TG
  list(
    tg         = tg,
    tfs        = tfs_for_tg,
    coefs = list(
      LASSO       = res_LASSO$coefs,
      ridge       = res_ridge$coefs,
      wLASSO      = res_wLASSO$coefs,
      wridge      = res_wridge$coefs,
      mixed_LASSO = res_mixed_LASSO$coefs,
      mixed_ridge = res_mixed_ridge$coefs
    ),
    lambdas = list(
      LASSO       = res_LASSO$lambda_min,
      ridge       = res_ridge$lambda_min,
      wLASSO      = res_wLASSO$lambda_min,
      wridge      = res_wridge$lambda_min,
      mixed_LASSO = res_mixed_LASSO$lambda_min,
      mixed_ridge = res_mixed_ridge$lambda_min
    ),
    times = list(
      LASSO       = res_LASSO$elapsed_sec,
      ridge       = res_ridge$elapsed_sec,
      wLASSO      = res_wLASSO$elapsed_sec,
      wridge      = res_wridge$elapsed_sec,
      mixed_LASSO = res_mixed_LASSO$elapsed_sec,
      mixed_ridge = res_mixed_ridge$elapsed_sec
    )
  )
}

total_time <- round((proc.time() - t_start)["elapsed"], 1)
message(sprintf("\nParallel fitting complete in %.1f seconds (%.1f minutes).",
                total_time, total_time / 60))

# ── Check 
errors <- vapply(results_list, inherits, logical(1), "error")
if (any(errors)) {
  warning(sprintf(
    "%d target gene(s) failed during fitting. Check results_list[errors] for details.",
    sum(errors)
  ))
  message("Failed TGs: ", paste(all_TGs[errors], collapse = ", "))
}

# ── Assemble coefficient matrices from the list 
message("Assembling coefficient matrices ...")

# Initialise matrices: rows = TGs, cols = TFs
make_coef_matrix <- function(tgs, tfs) {
  matrix(0, nrow = length(tgs), ncol = length(tfs),
         dimnames = list(tgs, tfs))
}
coef_LASSO       <- make_coef_matrix(all_TGs, all_TFs)
coef_ridge       <- make_coef_matrix(all_TGs, all_TFs)
coef_wLASSO      <- make_coef_matrix(all_TGs, all_TFs)
coef_wridge      <- make_coef_matrix(all_TGs, all_TFs)
coef_mixed_LASSO <- make_coef_matrix(all_TGs, all_TFs)
coef_mixed_ridge <- make_coef_matrix(all_TGs, all_TFs)

lambda_LASSO       <- setNames(numeric(length(all_TGs)), all_TGs)
lambda_ridge       <- setNames(numeric(length(all_TGs)), all_TGs)
lambda_wLASSO      <- setNames(numeric(length(all_TGs)), all_TGs)
lambda_wridge      <- setNames(numeric(length(all_TGs)), all_TGs)
lambda_mixed_LASSO <- setNames(numeric(length(all_TGs)), all_TGs)
lambda_mixed_ridge <- setNames(numeric(length(all_TGs)), all_TGs)

# Per-TG timing vectors (seconds per cv.glmnet call, summed across all TGs
# later to give total CPU time per method)
time_LASSO       <- setNames(numeric(length(all_TGs)), all_TGs)
time_ridge       <- setNames(numeric(length(all_TGs)), all_TGs)
time_wLASSO      <- setNames(numeric(length(all_TGs)), all_TGs)
time_wridge      <- setNames(numeric(length(all_TGs)), all_TGs)
time_mixed_LASSO <- setNames(numeric(length(all_TGs)), all_TGs)
time_mixed_ridge <- setNames(numeric(length(all_TGs)), all_TGs)

for (res in results_list) {
  if (inherits(res, "error")) next  
  tg  <- res$tg
  tfs <- res$tfs                   
  coef_LASSO[tg,       tfs] <- res$coefs$LASSO
  coef_ridge[tg,       tfs] <- res$coefs$ridge
  coef_wLASSO[tg,      tfs] <- res$coefs$wLASSO
  coef_wridge[tg,      tfs] <- res$coefs$wridge
  coef_mixed_LASSO[tg, tfs] <- res$coefs$mixed_LASSO
  coef_mixed_ridge[tg, tfs] <- res$coefs$mixed_ridge
  lambda_LASSO[tg]       <- res$lambdas$LASSO
  lambda_ridge[tg]       <- res$lambdas$ridge
  lambda_wLASSO[tg]      <- res$lambdas$wLASSO
  lambda_wridge[tg]      <- res$lambdas$wridge
  lambda_mixed_LASSO[tg] <- res$lambdas$mixed_LASSO
  lambda_mixed_ridge[tg] <- res$lambdas$mixed_ridge

  time_LASSO[tg]       <- res$times$LASSO
  time_ridge[tg]       <- res$times$ridge
  time_wLASSO[tg]      <- res$times$wLASSO
  time_wridge[tg]      <- res$times$wridge
  time_mixed_LASSO[tg] <- res$times$mixed_LASSO
  time_mixed_ridge[tg] <- res$times$mixed_ridge
}
message("Assembly complete.")


# ──  save results
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# Helper to build the result list and save
save_result <- function(coef_mat, lambda_vec, method_name) {
  result <- list(
    coef_matrix = coef_mat,
    lambda_min  = lambda_vec,
    method      = method_name,
    n_accessions = length(common_acc),
    n_TG        = length(all_TGs),
    n_TF        = length(all_TFs)
  )
  out_path <- file.path(OUT_DIR, paste0("coef_", method_name, ".rds"))
  saveRDS(result, file = out_path)
  message(sprintf("  Saved: %s", out_path))
  invisible(result)
}

message("\nSaving coefficient matrices ...")
save_result(coef_LASSO,       lambda_LASSO,       "LASSO")
save_result(coef_ridge,       lambda_ridge,       "ridge")
save_result(coef_wLASSO,      lambda_wLASSO,      "weighted_LASSO")
save_result(coef_wridge,      lambda_wridge,      "weighted_ridge")
save_result(coef_mixed_LASSO, lambda_mixed_LASSO, "mixed_LASSO")
save_result(coef_mixed_ridge, lambda_mixed_ridge, "mixed_ridge")


message("\n── Non-zero edges per method ──")
count_edges <- function(mat, label) {
  n <- sum(mat != 0)
  message(sprintf("  %-25s : %d non-zero coefficients", label, n))
}
count_edges(coef_LASSO,       "Classical LASSO")
count_edges(coef_ridge,       "Classical Ridge")
count_edges(coef_wLASSO,      "Weighted LASSO")
count_edges(coef_wridge,      "Weighted Ridge")
count_edges(coef_mixed_LASSO, "Mixed OLS-LASSO")
count_edges(coef_mixed_ridge, "Mixed OLS-Ridge")


# ── 9.  time summary
method_labels <- c(
  "Classical LASSO",
  "Classical Ridge",
  "Weighted LASSO",
  "Weighted Ridge",
  "Mixed OLS-LASSO",
  "Mixed OLS-Ridge"
)

time_vectors <- list(
  time_LASSO, time_ridge,
  time_wLASSO, time_wridge,
  time_mixed_LASSO, time_mixed_ridge
)

# Get summary data frame
timing_summary <- data.frame(
  Method            = method_labels,
  N_TG_fitted       = vapply(time_vectors, function(v) sum(v > 0), integer(1)),
  Total_CPU_sec     = vapply(time_vectors, sum,  numeric(1)),
  Mean_sec_per_TG   = vapply(time_vectors, mean, numeric(1)),
  Min_sec_per_TG    = vapply(time_vectors, min,  numeric(1)),
  Max_sec_per_TG    = vapply(time_vectors, max,  numeric(1)),
  stringsAsFactors  = FALSE
)
timing_summary$Total_CPU_min <- round(timing_summary$Total_CPU_sec / 60, 2)
timing_summary$Mean_sec_per_TG <- round(timing_summary$Mean_sec_per_TG, 4)
timing_summary$Min_sec_per_TG  <- round(timing_summary$Min_sec_per_TG,  4)
timing_summary$Max_sec_per_TG  <- round(timing_summary$Max_sec_per_TG,  4)
timing_summary$Total_CPU_sec   <- round(timing_summary$Total_CPU_sec,   2)

# Print table if needed
message("\n── Computational time per method ──")
message(sprintf(
  "  Wall-clock time for parallel fitting: %.1f seconds (%.2f minutes) on %d workers",
  total_time, total_time / 60, N_CORES
))
message(sprintf("  %-25s  %10s  %12s  %14s  %12s  %12s",
                "Method", "N_TG", "CPU_sec", "CPU_min", "Mean_s/TG", "Max_s/TG"))
message(strrep("-", 90))
for (i in seq_len(nrow(timing_summary))) {
  r <- timing_summary[i, ]
  message(sprintf("  %-25s  %10d  %12.2f  %14.2f  %12.4f  %12.4f",
                  r$Method, r$N_TG_fitted,
                  r$Total_CPU_sec, r$Total_CPU_min,
                  r$Mean_sec_per_TG, r$Max_sec_per_TG))
}

# Save timing summary 
timing_csv <- file.path(OUT_DIR, "computational_time_summary.csv")
write.csv(timing_summary, file = timing_csv, row.names = FALSE, quote = FALSE)
message(sprintf("\n  Timing summary saved to: %s", timing_csv))

# save the full per-TG timing matrix
per_tg_timing <- data.frame(
  TG              = all_TGs,
  LASSO_sec       = time_LASSO,
  ridge_sec       = time_ridge,
  wLASSO_sec      = time_wLASSO,
  wridge_sec      = time_wridge,
  mixed_LASSO_sec = time_mixed_LASSO,
  mixed_ridge_sec = time_mixed_ridge,
  row.names       = NULL
)
per_tg_csv <- file.path(OUT_DIR, "computational_time_per_TG.csv")
write.csv(per_tg_timing, file = per_tg_csv, row.names = FALSE, quote = FALSE)
message(sprintf("  Per-TG timing saved to  : %s", per_tg_csv))

message("\n=== Done. Results are in: ", normalizePath(OUT_DIR), " ===")

