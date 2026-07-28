###############################################################
# 11_Build_Beta_Network_Matrices
#
# Description:
#   1. Load regulator and regulated expression objects.
#   2. Load regression coefficients (beta values).
#   3. Remove constant-response genes.
#   4. Build edge-list representation of beta coefficients.
#   5. Normalize beta coefficients using min-max scaling.
#   6. Create regulator-target matrix.
#
# Inputs (results/):
#   - regulators_completenopv.RData
#   - regulated_completenopv.RData
#   - beta_cv_nopv.RData
#
# Outputs (results/):
#   - beta_df_nopv.RData
#   - beta_dfnorm_nopv.RData
#   - beta_mx_nopv.RData
#
###############################################################

############################
# Libraries
############################

library(data.table)
library(vctrs)
library(rlist)
library(dplyr)
library(tidyr)


############################
# Utility functions
############################

loadRData <- function(file_name) {
  
  e <- new.env()
  
  nm <- load(file_name, envir = e)
  
  e[[nm[1]]]
}

min_max_normalize <- function(x) {
  
  if (length(unique(x)) <= 1) {
    return(rep(0, length(x)))
  }
  
  abs(
    (x - min(x)) /
      (max(x) - min(x))
  )
}

############################
# Load data
############################

x_g <- loadRData("regulators_completenopv.RData")


y_g <- loadRData("regulated_completenopv.RData")

beta_cv <- loadRData("beta_cv_nopv.RData")

############################
# Center response variables
############################

y_scaled <- lapply(y_g,function(y) { scale(y, center = TRUE, scale = FALSE)})

names(y_scaled) <- names(y_g)

############################
# Remove constant responses
############################

n_unique <- sapply(y_scaled, function(x) length(unique(x)))

constant_ids <- which(n_unique == 1)

if (length(constant_ids) > 0) {
  
  y_scaled <- y_scaled[-constant_ids]
  
  cat(
    "Removed constant-response genes:",
    length(constant_ids),
    "\n"
  )
}

############################
# Build beta edge list
############################

beta_df <- vector(
  "list",
  length(y_scaled)
)

for (i in seq_along(y_scaled)) {
  
  weights <- unname(
    beta_cv[[i]]
  )
  
  beta_df[[i]] <- data.frame(
    from = colnames(x_g),
    to = rep(
      names(y_scaled)[i],
      ncol(x_g)
    ),
    weight = weights
  )
}

beta_df <- bind_rows(
  beta_df
)

cat(
  "Beta edge list dimensions:",
  dim(beta_df),
  "\n"
)

save(beta_df, "beta_df_nopv.RData")


############################
# Normalize beta values
############################

beta_norm <- lapply(
  beta_cv,
  function(x) {
    
    x <- unname(x)
    
    if (length(x) == 1) {
      
      return(x)
      
    } else {
      
      return(
        min_max_normalize(x)
      )
    }
  }
)

names(beta_norm) <- names(beta_cv)

############################
# Build normalized edge list
############################

beta_dfnorm <- vector(
  "list",
  length(y_scaled)
)

for (i in seq_along(y_scaled)) {
  
  beta_dfnorm[[i]] <- data.frame(
    from = colnames(x_g),
    to = rep(
      names(y_scaled)[i],
      ncol(x_g)
    ),
    weight = unname(
      beta_norm[[i]]
    )
  )
}

beta_dfnorm <- bind_rows(
  beta_dfnorm
)

cat(
  "Normalized edge list dimensions:",
  dim(beta_dfnorm),
  "\n"
)

save(
  beta_dfnorm,
  file = file.path(
    results_dir,
    "beta_dfnorm_nopv.RData"
  )
)

############################
# Create beta matrix
############################

beta_tb <- pivot_wider(
  beta_dfnorm,
  names_from = from,
  values_from = weight
)

row_names <- beta_tb$to

beta_tb <- beta_tb[, -1]

beta_mx <- apply(
  as.matrix(beta_tb),
  2,
  as.numeric
)

rownames(beta_mx) <- row_names

############################
# Save matrix
############################

save( beta_mx,"beta_mx_nopv.RData")
