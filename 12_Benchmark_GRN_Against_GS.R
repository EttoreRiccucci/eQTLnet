###############################################################
# 12_Benchmark_GRN_Against_GS
#
# Description
#   Evaluate inferred regulatory interactions against a
#   maize gold-standard network using:
#
#     1. Balanced positive/negative sets
#     2. AUROC
#     3. AUPR
#     4. R² thresholds
#     5. Adjusted R² thresholds
#
###############################################################

############################
# Libraries
############################

library(data.table)
library(dplyr)
library(precrec)

############################
# Utility function
############################

loadRData <- function(file_name) {
  
  e <- new.env()
  
  nm <- load(
    file_name,
    envir = e
  )
  
  e[[nm[1]]]
}

############################
# Load inferred network
############################

beta_df <- loadRData("beta_dfnorm_nopv.RData")

beta_mx <- loadRData("beta_mx_nopv.RData")


rsq <- loadRData("rsq.RData")

adj_rsq <- loadRData("adj_rsq.RData")

beta_df$int <- paste0(
  beta_df$from,
  beta_df$to
)

############################
# Gold-standard network
############################

positive_interactions <- fread("AM_Final_Maize_GS.txt",  data.table = FALSE)

colnames(positive_interactions) <- c(
  "from",
  "to",
  "weight"
)

positive_interactions$int <- paste0(
  positive_interactions$from,
  positive_interactions$to
)

############################
# Generate negatives
############################

all_gene_pairs <- expand.grid(
  from = unique(
    positive_interactions$from
  ),
  to = unique(
    positive_interactions$to
  )
)

negative_interactions <-
  all_gene_pairs[
    !paste(
      all_gene_pairs$from,
      all_gene_pairs$to
    ) %in%
      positive_interactions$int,
  ]

negative_interactions$weight <- 0

############################
# Evaluation function
############################

evaluate_auc <- function(
    prediction_df,
    gs_df
) {
  
  shared <- intersect(
    prediction_df$int,
    gs_df$int
  )
  
  pred <- prediction_df[
    prediction_df$int %in% shared,
  ]
  
  truth <- gs_df[
    gs_df$int %in% shared,
  ]
  
  truth <- truth |>
    distinct(int, .keep_all = TRUE) |>
    arrange(int)
  
  pred <- pred |>
    arrange(int)
  
  stopifnot(
    identical(
      truth$int,
      pred$int
    )
  )
  
  scores <- evalmod(
    scores = pred$weight,
    labels = truth$weight
  )
  
  auc_tbl <- auc(scores)
  
  data.frame(
    AUROC = auc_tbl[1, 4],
    AUPR = auc_tbl[2, 4],
    N = length(shared)
  )
}

############################
# Balanced benchmarking
############################

benchmark_results <- list()

for (seed in 1:20) {
  
  set.seed(seed)
  
  neg_sample <-
    negative_interactions[
      sample(
        nrow(negative_interactions),
        nrow(positive_interactions)
      ),
    ]
  
  gs_balanced <- rbind(
    positive_interactions,
    neg_sample
  )
  
  gs_balanced$int <- paste0(
    gs_balanced$from,
    gs_balanced$to
  )
  
  benchmark_results[[seed]] <-
    cbind(
      seed = seed,
      evaluate_auc(
        beta_df,
        gs_balanced
      )
    )
}

benchmark_results <-
  bind_rows(
    benchmark_results
  )

write.table(benchmark_results,"balanced_benchmark_metrics.txt",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

############################
# Threshold analysis
############################

evaluate_threshold_curve <- function(
    metric_vector,
    prediction_df,
    positive_interactions
) {
  
  thresholds <- seq(
    0,
    1,
    by = 0.05
  )
  
  out <- data.frame(
    threshold = thresholds,
    AUROC = NA,
    AUPR = NA
  )
  
  for (i in seq_along(thresholds)) {
    
    selected_targets <-
      names(
        metric_vector[
          metric_vector >
            thresholds[i]
        ]
      )
    
    pred_subset <-
      prediction_df[
        prediction_df$to %in%
          selected_targets,
      ]
    
    auc_res <- evaluate_auc(
      pred_subset,
      positive_interactions
    )
    
    out$AUROC[i] <- auc_res$AUROC
    out$AUPR[i] <- auc_res$AUPR
  }
  
  out
}

############################
# R² thresholds
############################

curves_rsq <-
  evaluate_threshold_curve(
    metric_vector = rsq,
    prediction_df = beta_df,
    positive_interactions =
      positive_interactions
  )

write.table(curves_rsq,"curves_rsq.txt"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

############################
# Adjusted R² thresholds
############################

curves_adj_rsq <-
  evaluate_threshold_curve(
    metric_vector = adj_rsq,
    prediction_df = beta_df,
    positive_interactions =
      positive_interactions
  )

write.table(curves_adj_rsq,"curves_adj_rsq.txt",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
