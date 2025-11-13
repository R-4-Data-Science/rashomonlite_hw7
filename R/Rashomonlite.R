#' Fit multiple logistic regression models of a specified dimension
#'
#' This function fits all (or a sampled subset) of generalized linear models (GLMs)
#' of a given dimension \code{k} from the predictors in \code{X}. It computes model AIC,
#' regression coefficients, and p-values for each model.
#'
#' @param X A data frame or matrix of predictors. Columns should be numeric or
#'   factors suitable for inclusion in a logistic regression.
#' @param y A numeric vector of binary responses (0 = failure, 1 = success),
#'   with length equal to the number of rows in \code{X}.
#' @param k Integer specifying the number of predictors to include in each model.
#' @param m Integer specifying the maximum number of models to fit for this dimension.
#'   If the total number of possible combinations exceeds \code{m}, a random sample of
#'   size \code{m} is drawn (with fixed seed for reproducibility).
#'
#' @return A data frame with one row per fitted model and the following columns:
#' \describe{
#'   \item{vars}{Character vector of predictor names included in the model.}
#'   \item{aic}{Akaike Information Criterion (AIC) value for the model.}
#'   \item{coef}{Named numeric vector of model coefficients.}
#'   \item{pvals}{Named numeric vector of p-values for each coefficient.}
#' }
#' @author Prince Mensah Ansah,
#'         Mark Castilo Philip,
#'         Mohammad Al Srayheen
#' @export
fit_models_dim <- function(X, y, k, m) {
  all_vars <- colnames(X)
  n_vars <- length(all_vars)

  # Generate all combinations of size k
  if (k == 1) {
    combos <- as.list(all_vars)
  } else {
    combos <- combn(all_vars, k, simplify = FALSE)
  }

  # Sample m combinations if more than m exist
  if (length(combos) > m) {
    set.seed(123)
    combos <- combos[sample(length(combos), m)]
  }

  # Fit models
  models_list <- lapply(combos, function(vars) {
    formula_str <- paste("y ~", paste(vars, collapse = " + "))
    formula_obj <- as.formula(formula_str)

    fit <- suppressWarnings(
      stats::glm(formula_obj,
                 data = as.data.frame(cbind(y = y, X[, vars, drop = FALSE])),
                 family = stats::binomial,
                 control = stats::glm.control(maxit = 100))
    )

    coef_summary <- summary(fit)$coefficients
    coefs <- stats::coef(fit)
    pvals <- coef_summary[, 4]

    list(vars = vars, aic = stats::AIC(fit), coef = coefs, pvals = pvals)
  })

  data.frame(
    vars = I(lapply(models_list, function(x) x$vars)),
    aic = sapply(models_list, function(x) x$aic),
    coef = I(lapply(models_list, function(x) x$coef)),
    pvals = I(lapply(models_list, function(x) x$pvals)),
    stringsAsFactors = FALSE
  )
}


#' Select the top-performing models based on AIC
#'
#' This function takes a data frame of models (as returned by `fit_models_dim()`)
#' and returns the top fraction `alpha` of models sorted by lowest AIC.
#'
#' @param models_df A data frame returned by `fit_models_dim()`, containing AIC values.
#' @param alpha Numeric between 0 and 1 specifying the fraction of models to retain.
#'
#' @return A subset of the input data frame containing only the top `alpha` proportion
#'   of models sorted by lowest AIC.
#' @author Prince Mensah Ansah,
#'         Mark Castilo Philip,
#'         Mohammad Al Srayheen
#' @export
select_best <- function(models_df, alpha = 0.5) {
  n_select <- ceiling(alpha * nrow(models_df))
  sorted_idx <- order(models_df$aic)
  best_idx <- sorted_idx[1:n_select]
  models_df[best_idx, , drop = FALSE]
}



#' Expand top models to the next dimension
#'
#' @param best_df A data frame of best models from the previous dimension (as returned by [select_best()]).
#' @param all_vars Character vector of all available predictor names.
#' @param k_next Integer, size of the next dimension of models to grow to.
#' @param m Integer specifying the maximum number of new models to generate.
#'
#' @return A list of unique character vectors, where each vector represents the set
#'   of predictors for a candidate model of dimension `k_next`.
#' @author Prince Mensah Ansah,
#'         Mark Castilo Philip,
#'         Mohammad Al Srayheen
#'
#' @export
grow_from <- function(best_df, all_vars, k_next, m) {
  combos_list <- list()

  for (i in 1:nrow(best_df)) {
    base_vars <- best_df$vars[[i]]
    remaining_vars <- setdiff(all_vars, base_vars)

    for (new_var in remaining_vars) {
      new_combo <- sort(c(base_vars, new_var))
      combos_list[[length(combos_list) + 1]] <- new_combo
    }
  }

  combos_unique <- unique(combos_list)
  if (length(combos_unique) > m) {
    set.seed(123 + k_next)
    combos_unique <- combos_unique[sample(length(combos_unique), m)]
  }
  combos_unique
}


#' Run Rashomon model selection
#'
#' This function performs iterative model fitting and selection across increasing
#' model dimensions (from 1 up to \code{pmax}). At each step, it keeps the top-performing
#' models (based on AIC) and expands them to explore higher-dimensional models.
#'
#' @param X A data frame or matrix of predictors. All columns should be suitable
#'   for logistic regression (numeric or factors).
#' @param y A numeric vector of binary responses (0/1), with length matching \code{nrow(X)}.
#' @param pmax Integer specifying the maximum model dimension to explore.
#' @param m Integer specifying the maximum number of models per dimension.
#' @param alpha Numeric between 0 and 1 specifying the proportion of top models to retain at each dimension.
#'
#' @return A list containing:
#' \describe{
#'   \item{M}{List of top model data frames at each dimension.}
#'   \item{all_models}{List of all fitted models per dimension.}
#'   \item{predictor_counts}{Matrix showing the frequency of predictor usage by dimension.}
#'   \item{aic_summary}{Data frame summarizing AIC values per dimension.}
#' }
#'
#' @examples
#' # Example using Titanic data
#' library(titanic)
#' data("titanic_train")
#'
#' # Select relevant columns and remove missing values
#' titanic_data <- titanic_train[, c("Survived", "Pclass", "Sex", "Age", "SibSp", "Parch", "Fare", "Embarked")]
#' titanic_data <- na.omit(titanic_data)
#'
#' # Binary response variable (0 = did not survive, 1 = survived)
#' y <- as.numeric(titanic_data$Survived)
#'
#' # Create design matrix of predictors (excluding intercept)
#' X <- model.matrix(~ Pclass + Sex + Age + SibSp + Parch + Fare + Embarked,
#'                   data = titanic_data)[, -1]
#'
#' # Run Rashomon model selection
#' result <- run_rashomon(X, y, pmax = 3, m = 20, alpha = 0.5)
#'
#' # Inspect AIC summary to verify package works
#' head(result$aic_summary)
#' @author Prince Mensah Ansah,
#'         Mark Castilo Philip,
#'         Mohammad Al Srayheen
#'
#' @export
run_rashomon <- function(X, y, pmax = 5, m = 90, alpha = 0.5) {
  all_vars <- colnames(X)
  M <- list()
  all_models <- list()

  cat("Fitting dimension 1...\n")
  models_1 <- fit_models_dim(X, y, k = 1, m = length(all_vars))
  best_1 <- select_best(models_1, alpha = alpha)
  M[[1]] <- best_1
  all_models[[1]] <- models_1

  for (k in 2:pmax) {
    cat("Fitting dimension", k, "...\n")
    combos <- grow_from(M[[k - 1]], all_vars, k, m)

    models_list <- lapply(combos, function(vars) {
      formula_str <- paste("y ~", paste(vars, collapse = " + "))
      formula_obj <- as.formula(formula_str)
      fit <- suppressWarnings(
        stats::glm(formula_obj,
                   data = as.data.frame(cbind(y = y, X[, vars, drop = FALSE])),
                   family = stats::binomial,
                   control = stats::glm.control(maxit = 100))
      )
      coef_summary <- summary(fit)$coefficients
      coefs <- stats::coef(fit)
      pvals <- coef_summary[, 4]
      list(vars = vars, aic = stats::AIC(fit), coef = coefs, pvals = pvals)
    })

    models_df <- data.frame(
      vars = I(lapply(models_list, function(x) x$vars)),
      aic = sapply(models_list, function(x) x$aic),
      coef = I(lapply(models_list, function(x) x$coef)),
      pvals = I(lapply(models_list, function(x) x$pvals)),
      stringsAsFactors = FALSE
    )

    best_k <- select_best(models_df, alpha = alpha)
    M[[k]] <- best_k
    all_models[[k]] <- models_df
  }

  predictor_counts <- matrix(0, nrow = pmax, ncol = length(all_vars))
  colnames(predictor_counts) <- all_vars
  rownames(predictor_counts) <- paste0("Dim_", 1:pmax)

  for (k in 1:pmax) {
    models_k <- all_models[[k]]
    for (i in 1:nrow(models_k)) {
      vars <- models_k$vars[[i]]
      for (var in vars) {
        predictor_counts[k, var] <- predictor_counts[k, var] + 1
      }
    }
  }

  aic_summary_df <- do.call(rbind, lapply(1:pmax, function(k) {
    data.frame(dimension = k, aic = all_models[[k]]$aic)
  }))

  list(
    M = M,
    all_models = all_models,
    predictor_counts = predictor_counts,
    aic_summary = aic_summary_df
  )
}

