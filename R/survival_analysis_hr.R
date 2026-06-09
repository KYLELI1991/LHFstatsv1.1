# Load required packages
library(survival)
library(survminer)
library(dplyr)
library(tidyr)

#' Univariate Cox regression analysis with detailed result table (complete version)
#'
#' @param data Data frame containing PFS time, status, and other variables
#' @param time_var Name of PFS time variable (character)
#' @param status_var Name of PFS status variable (character, typically 1=event, 0=censored)
#' @param covariates Character vector of variable names to analyze
#' @param cut_method Method for dichotomizing continuous variables: "median" or "best" (maximally selected rank statistics)
#' @param add_continuous Whether to additionally include continuous variables directly (ungrouped, HR per 1-unit increase)
#' @param return_cutpoint Whether to return the best cutpoint value (only effective when cut_method="best")
#'
#' @return A data frame
#' @export
#'
#' @examples
#' library(survival)
#' data(lung)
#' lung <- lung %>%
#'   mutate(PFS = time,
#'          PFS_status = ifelse(status == 1, 0, 1))  # convert to 0=censored, 1=event
#' result <- uni_cox_table(lung, "PFS", "PFS_status", 
#'                         covariates = c("age", "sex", "ph.ecog"),
#'                         cut_method = "median",
#'                         add_continuous = TRUE)
#' print(result)
uni_cox_table <- function(data, time_var, status_var, covariates, 
                          cut_method = c("median", "best"), 
                          add_continuous = TRUE,
                          return_cutpoint = FALSE) {
  
  cut_method <- match.arg(cut_method)
  
  # ========== 1. Data validation ==========
  if (!is.data.frame(data)) {
    stop("data must be a data frame")
  }
  
  if (!all(c(time_var, status_var, covariates) %in% colnames(data))) {
    missing_vars <- c(time_var, status_var, covariates)[!c(time_var, status_var, covariates) %in% colnames(data)]
    stop("Missing columns: ", paste(missing_vars, collapse = ", "))
  }
  
  # Extract survival time and status
  surv_time <- data[[time_var]]
  surv_status <- data[[status_var]]
  
  # Check if status variable is 0/1
  if (!all(surv_status %in% c(0, 1))) {
    warning("Status variable should be 0=censored, 1=event. Attempting to convert...")
    if (all(surv_status %in% c(1, 2))) {
      surv_status <- surv_status - 1
      message("Converted status: 1->0, 2->1")
    } else {
      stop("Status variable must contain only 0/1 or 1/2 values")
    }
  }
  
  # Remove missing values
  complete_idx <- complete.cases(data[, c(time_var, status_var, covariates)])
  if (sum(!complete_idx) > 0) {
    message("Removed ", sum(!complete_idx), " rows with missing values")
    data <- data[complete_idx, ]
    surv_time <- data[[time_var]]
    surv_status <- data[[status_var]]
  }
  
  # ========== 2. Overall median PFS ==========
  overall_surv <- tryCatch({
    survfit(Surv(surv_time, surv_status) ~ 1, data = data)
  }, error = function(e) {
    stop("Failed to calculate overall survival: ", e$message)
  })
  
  overall_summary <- summary(overall_surv)
  overall_median <- overall_summary$table["median"]
  overall_ci_lower <- overall_summary$table["0.95LCL"]
  overall_ci_upper <- overall_summary$table["0.95UCL"]
  
  # Handle NA cases
  if (is.na(overall_median)) {
    overall_median_ci <- "NE"
  } else {
    overall_median_ci <- paste0(round(overall_median, 2), " (", 
                                round(overall_ci_lower, 2), "-", 
                                round(overall_ci_upper, 2), ")")
  }
  
  # ========== 3. Analyze each variable ==========
  res_list <- list()
  cutpoint_info <- list()
  
  for (var in covariates) {
    cat("Processing:", var, "\n")
    x <- data[[var]]
    
    # Determine if variable is continuous (numeric and >5 unique values)
    is_continuous <- is.numeric(x) && length(unique(na.omit(x))) > 5
    
    # ===== 3.1 Continuous variable raw analysis =====
    if (is_continuous && add_continuous) {
      # Direct Cox regression with continuous variable
      temp_data_cont <- data.frame(time = surv_time, status = surv_status, x = x)
      temp_data_cont <- temp_data_cont[complete.cases(temp_data_cont), ]
      
      cox_cont <- tryCatch({
        coxph(Surv(time, status) ~ x, data = temp_data_cont)
      }, error = function(e) {
        warning("Failed to fit Cox model for continuous variable ", var, ": ", e$message)
        return(NULL)
      })
      
      if (!is.null(cox_cont)) {
        cox_cont_sum <- summary(cox_cont)
        hr_cont <- exp(coef(cox_cont))
        ci_cont <- exp(confint(cox_cont))
        p_cont <- coef(summary(cox_cont))[, "Pr(>|z|)"]
        
        hr_ci_cont <- paste0(round(hr_cont, 3), " (", 
                             round(ci_cont[1], 3), "-", 
                             round(ci_cont[2], 3), ")")
        
        res_list[[length(res_list)+1]] <- data.frame(
          Variable = var,
          Type = "Continuous (per 1-unit increase)",
          Level = " - ",
          N = nrow(temp_data_cont),
          Events = sum(temp_data_cont$status),
          HR_95CI = hr_ci_cont,
          P_value = signif(p_cont, 3),
          Median_PFS_total = overall_median_ci,
          Median_PFS_subgroup = overall_median_ci,
          Cut_method = "continuous_raw",
          stringsAsFactors = FALSE
        )
      }
    }
    
    # ===== 3.2 Categorical/binary version =====
    if (is_continuous) {
      # Continuous variable: convert to categorical
      if (cut_method == "median") {
        cut_val <- median(x, na.rm = TRUE)
        group <- ifelse(x <= cut_val, paste0("<= ", round(cut_val, 2)), 
                        paste0("> ", round(cut_val, 2)))
        method_used <- paste0("median(", round(cut_val, 2), ")")
        var_type <- "Categorized (median split)"
      } else if (cut_method == "best") {
        # Best cutpoint, requires sufficient sample size
        if (sum(!is.na(x)) < 20) {
          warning(var, " has insufficient sample for best cutpoint, using median instead")
          cut_val <- median(x, na.rm = TRUE)
          group <- ifelse(x <= cut_val, paste0("<= ", round(cut_val, 2)), 
                          paste0("> ", round(cut_val, 2)))
          method_used <- paste0("median(", round(cut_val, 2), ")")
          var_type <- "Categorized (median split)"
        } else {
          cut_obj <- tryCatch({
            surv_cutpoint(data, time = time_var, event = status_var, variables = var)
          }, error = function(e) {
            warning("surv_cutpoint failed for ", var, ": ", e$message, ". Using median.")
            return(NULL)
          })
          
          if (!is.null(cut_obj)) {
            cut_val <- cut_obj$cutpoint[1, "cutpoint"]
            group <- ifelse(x <= cut_val, paste0("<= ", round(cut_val, 2)), 
                            paste0("> ", round(cut_val, 2)))
            method_used <- paste0("best_cutpoint(", round(cut_val, 2), ")")
            var_type <- "Categorized (optimal cutpoint)"
            cutpoint_info[[var]] <- cut_val
          } else {
            cut_val <- median(x, na.rm = TRUE)
            group <- ifelse(x <= cut_val, paste0("<= ", round(cut_val, 2)), 
                            paste0("> ", round(cut_val, 2)))
            method_used <- paste0("median(", round(cut_val, 2), ")")
            var_type <- "Categorized (median split)"
          }
        }
      }
      
      temp_data <- data.frame(time = surv_time, status = surv_status, group = group)
      temp_data <- temp_data[complete.cases(temp_data), ]
      
    } else {
      # Categorical variable
      group <- as.factor(x)
      temp_data <- data.frame(time = surv_time, status = surv_status, group = group)
      temp_data <- temp_data[complete.cases(temp_data), ]
      method_used <- "categorical"
      var_type <- "Categorical"
      cut_val <- NA
    }
    
    # Check if there are at least two groups with events after categorization
    if (length(unique(temp_data$group)) < 2) {
      warning(var, " has less than 2 groups after categorization, skipping")
      next
    }
    
    # Univariate Cox model
    cox_fit <- tryCatch({
      coxph(Surv(time, status) ~ group, data = temp_data)
    }, error = function(e) {
      warning("Failed to fit Cox model for ", var, ": ", e$message)
      return(NULL)
    })
    
    if (is.null(cox_fit)) {
      next
    }
    
    cox_sum <- summary(cox_fit)
    
    # Extract HR and CI
    coef_mat <- coef(summary(cox_fit))
    hr <- exp(coef(cox_fit))
    ci <- exp(confint(cox_fit))
    
    # Median PFS by group
    surv_by_group <- tryCatch({
      survfit(Surv(time, status) ~ group, data = temp_data)
    }, error = function(e) {
      warning("Failed to calculate survival for ", var, ": ", e$message)
      return(NULL)
    })
    
    if (is.null(surv_by_group)) {
      next
    }
    
    surv_summary <- summary(surv_by_group)
    group_median <- surv_summary$table[, "median"]
    group_lcl <- surv_summary$table[, "0.95LCL"]
    group_ucl <- surv_summary$table[, "0.95UCL"]
    
    # Calculate number of patients and events per group
    group_levels <- levels(as.factor(temp_data$group))
    group_n <- numeric(length(group_levels))
    group_events <- numeric(length(group_levels))
    
    for (i in seq_along(group_levels)) {
      idx <- which(temp_data$group == group_levels[i])
      group_n[i] <- length(idx)
      group_events[i] <- sum(temp_data$status[idx])
    }
    
    # Build output rows
    ref_level <- group_levels[1]
    for (i in seq_along(group_levels)) {
      lev <- group_levels[i]
      
      if (lev == ref_level) {
        hr_ci_str <- "1 (ref)"
        p_val_use <- NA
      } else {
        # Find the corresponding coefficient (handle possible different contrast formats)
        coef_row <- which(rownames(coef_mat) == paste0("group", lev))
        if (length(coef_row) == 0) {
          coef_row <- which(rownames(coef_mat) == lev)
        }
        if (length(coef_row) > 0) {
          hr_val <- exp(coef_mat[coef_row, "coef"])
          ci_low <- exp(coef_mat[coef_row, "coef"] - 1.96 * coef_mat[coef_row, "se(coef)"])
          ci_up <- exp(coef_mat[coef_row, "coef"] + 1.96 * coef_mat[coef_row, "se(coef)"])
          p_val_use <- coef_mat[coef_row, "Pr(>|z|)"]
          hr_ci_str <- paste0(round(hr_val, 3), " (", 
                              round(ci_low, 3), "-", 
                              round(ci_up, 3), ")")
        } else {
          hr_ci_str <- "NA"
          p_val_use <- NA
        }
      }
      
      # Median PFS string
      med_pfs <- ifelse(is.na(group_median[i]), NA, group_median[i])
      med_lcl <- ifelse(is.na(group_lcl[i]), NA, group_lcl[i])
      med_ucl <- ifelse(is.na(group_ucl[i]), NA, group_ucl[i])
      
      if (is.na(med_pfs)) {
        med_pfs_str <- "NE"
      } else {
        med_pfs_str <- paste0(round(med_pfs, 2), " (", 
                              round(med_lcl, 2), "-", 
                              round(med_ucl, 2), ")")
      }
      
      res_list[[length(res_list)+1]] <- data.frame(
        Variable = var,
        Type = var_type,
        Level = as.character(lev),
        N = group_n[i],
        Events = group_events[i],
        HR_95CI = hr_ci_str,
        P_value = ifelse(is.na(p_val_use), NA, signif(p_val_use, 3)),
        Median_PFS_total = overall_median_ci,
        Median_PFS_subgroup = med_pfs_str,
        Cut_method = method_used,
        stringsAsFactors = FALSE
      )
    }
  }
  
  # ========== 4. Combine results ==========
  if (length(res_list) == 0) {
    stop("No valid results generated")
  }
  
  result_df <- bind_rows(res_list)
  
  # Adjust column order
  result_df <- result_df %>%
    select(Variable, Type, Level, N, Events, HR_95CI, P_value, 
           Median_PFS_total, Median_PFS_subgroup, everything())
  
  # Add attributes
  attr(result_df, "total_patients") <- nrow(data)
  attr(result_df, "total_events") <- sum(surv_status)
  attr(result_df, "overall_median_pfs") <- overall_median_ci
  
  if (return_cutpoint && cut_method == "best" && length(cutpoint_info) > 0) {
    attr(result_df, "best_cutpoints") <- cutpoint_info
  }
  
  return(result_df)
}

#' Get variable names with p-value below a specified threshold
#'
#' @param result_table The result data frame from uni_cox_table function
#' @param p_threshold Numeric threshold for p-value (default: 0.05)
#' @param include_ref Whether to include reference levels (default: FALSE, only returns non-reference levels)
#' @param unique_vars Whether to return unique variable names (default: TRUE, removes duplicate variable names)
#'
#' @return A character vector of variable names (or a data frame with details if return_df=TRUE)
#' @export
#'
#' @examples
#' result <- uni_cox_table(lung_data, "PFS", "status", covariates = c("age", "sex", "treatment"))
#' sig_vars <- get_significant_vars(result, p_threshold = 0.05)
#' print(sig_vars)
get_significant_vars <- function(result_table, p_threshold = 0.05, 
                                 include_ref = FALSE, unique_vars = TRUE) {
  
  # Input validation
  if (!is.data.frame(result_table)) {
    stop("result_table must be a data frame")
  }
  
  if (!"P_value" %in% colnames(result_table)) {
    stop("result_table must contain a 'P_value' column")
  }
  
  if (!is.numeric(p_threshold) || p_threshold <= 0 || p_threshold >= 1) {
    stop("p_threshold must be a numeric value between 0 and 1")
  }
  
  # Filter rows based on p-value threshold
  # Exclude reference levels if include_ref = FALSE
  if (include_ref) {
    sig_rows <- result_table[!is.na(result_table$P_value) & result_table$P_value < p_threshold, ]
  } else {
    sig_rows <- result_table[!is.na(result_table$P_value) & 
                               result_table$P_value < p_threshold & 
                               result_table$HR_95CI != "1 (ref)", ]
  }
  
  if (nrow(sig_rows) == 0) {
    message("No significant variables found with p < ", p_threshold)
    return(character(0))
  }
  
  # Return unique variable names or full details
  if (unique_vars) {
    significant_vars <- unique(sig_rows$Variable)
    return(significant_vars)
  } else {
    return(sig_rows[, c("Variable", "Type", "Level", "HR_95CI", "P_value")])
  }
}

#' Prepare data for multivariate analysis with categorized variables
#' @export
#'
#' @param data Original data frame
#' @param time_var Name of time variable
#' @param status_var Name of status variable
#' @param sig_vars Character vector of significant variable names from univariate analysis
#' @param result_table Result table from uni_cox_table (contains cutpoint information)
#'
#' @return A list containing prepared data and formula
prepare_multivariate_data <- function(data, time_var, status_var, sig_vars, result_table) {
  
  # Extract cutpoint information for continuous variables
  cutpoint_info <- list()
  
  for (var in sig_vars) {
    # Check if this variable was categorized in univariate analysis
    var_rows <- result_table[result_table$Variable == var, ]
    
    if (nrow(var_rows) > 0) {
      cut_method <- unique(var_rows$Cut_method)
      
      # If it's a categorized variable (not continuous_raw and not categorical)
      if (any(grepl("cutpoint|median", cut_method)) && !grepl("continuous_raw", cut_method)) {
        # Extract the cut value from the level names
        level_names <- var_rows$Level
        cut_value <- NA
        
        for (lev in level_names) {
          if (grepl("<= ", lev)) {
            cut_value <- as.numeric(gsub("<= ", "", lev))
            break
          } else if (grepl("> ", lev)) {
            cut_value <- as.numeric(gsub("> ", "", lev))
            break
          }
        }
        
        if (!is.na(cut_value)) {
          cutpoint_info[[var]] <- cut_value
        }
      }
    }
  }
  
  # Create a new data frame with categorized variables
  multivariate_data <- data[, c(time_var, status_var)]
  
  for (var in sig_vars) {
    if (var %in% names(cutpoint_info)) {
      # Categorize continuous variable using the cutpoint from univariate analysis
      cut_val <- cutpoint_info[[var]]
      multivariate_data[[var]] <- ifelse(data[[var]] <= cut_val, 
                                         paste0("<= ", round(cut_val, 2)),
                                         paste0("> ", round(cut_val, 2)))
      multivariate_data[[var]] <- factor(multivariate_data[[var]], 
                                         levels = c(paste0("<= ", round(cut_val, 2)), 
                                                    paste0("> ", round(cut_val, 2))))
    } else {
      # Use original variable (categorical or continuous not categorized)
      multivariate_data[[var]] <- data[[var]]
    }
  }
  
  # Remove rows with missing values
  multivariate_data <- multivariate_data[complete.cases(multivariate_data), ]
  
  return(list(data = multivariate_data, cutpoints = cutpoint_info))
}

#' Perform multivariate Cox regression analysis
#'
#' @param data Original data frame
#' @param time_var Name of time variable
#' @param status_var Name of status variable
#' @param sig_vars Character vector of significant variable names from univariate analysis
#' @param result_table Result table from uni_cox_table (optional, used for maintaining consistent categorization)
#' @param method Method for variable selection: "enter" (enter all), "forward", "backward", or "stepwise"
#' @param p_enter P-value threshold for entry in stepwise selection (default: 0.05)
#' @param p_remove P-value threshold for removal in stepwise selection (default: 0.10)
#'
#' @return A list containing multivariate results
#' @export
#'
#' @examples
#' # First perform univariate analysis
#' uni_result <- uni_cox_table(data, "PFS", "status", covariates = c("age", "sex", "treatment"))
#' # Get significant variables
#' sig_vars <- get_significant_vars(uni_result, p_threshold = 0.05)
#' # Perform multivariate analysis
#' multi_result <- multivariate_cox(data, "PFS", "status", sig_vars, uni_result)
multivariate_cox <- function(data, time_var, status_var, sig_vars, 
                             result_table = NULL, 
                             method = c("enter", "forward", "backward", "stepwise"),
                             p_enter = 0.05, p_remove = 0.10) {
  
  method <- match.arg(method)
  
  if (length(sig_vars) == 0) {
    stop("No significant variables provided for multivariate analysis")
  }
  
  # Prepare data
  if (!is.null(result_table)) {
    prep_result <- prepare_multivariate_data(data, time_var, status_var, sig_vars, result_table)
    multi_data <- prep_result$data
  } else {
    # Use original variables without categorization
    multi_data <- data[, c(time_var, status_var, sig_vars)]
    multi_data <- multi_data[complete.cases(multi_data), ]
  }
  
  # Create survival object
  surv_obj <- Surv(multi_data[[time_var]], multi_data[[status_var]])
  
  # Create formula for multivariate model
  formula_str <- paste("surv_obj ~", paste(sig_vars, collapse = " + "))
  formula_obj <- as.formula(formula_str)
  
  # Fit multivariate Cox model based on selected method
  if (method == "enter") {
    # Enter all variables
    cox_model <- tryCatch({
      coxph(formula_obj, data = multi_data)
    }, error = function(e) {
      stop("Failed to fit multivariate Cox model: ", e$message)
    })
    
  } else {
    # Stepwise selection methods
    # Start with null model
    null_model <- coxph(surv_obj ~ 1, data = multi_data)
    full_model <- coxph(formula_obj, data = multi_data)
    
    if (method == "forward") {
      cox_model <- step(null_model, 
                        scope = formula_obj,
                        direction = "forward",
                        trace = 0,
                        k = qchisq(1 - p_enter, 1))
    } else if (method == "backward") {
      cox_model <- step(full_model, 
                        direction = "backward",
                        trace = 0,
                        k = qchisq(1 - p_remove, 1))
    } else { # stepwise
      cox_model <- step(null_model, 
                        scope = formula_obj,
                        direction = "both",
                        trace = 0,
                        k = qchisq(1 - p_enter, 1))
    }
  }
  
  # Extract results
  cox_summary <- summary(cox_model)
  coef_table <- as.data.frame(cox_summary$coefficients)
  conf_int <- exp(confint(cox_model))
  
  # Build result table
  result_df <- data.frame(
    Variable = rownames(coef_table),
    N = nrow(multi_data),
    Events = sum(multi_data[[status_var]]),
    HR = exp(coef_table$coef),
    HR_Lower = exp(conf_int[, 1]),
    HR_Upper = exp(conf_int[, 2]),
    HR_95CI = paste0(round(exp(coef_table$coef), 3), " (", 
                     round(exp(conf_int[, 1]), 3), "-", 
                     round(exp(conf_int[, 2]), 3), ")"),
    P_value = coef_table$`Pr(>|z|)`,
    stringsAsFactors = FALSE
  )
  
  # Add significance stars
  result_df$Significance <- ifelse(result_df$P_value < 0.001, "***",
                                   ifelse(result_df$P_value < 0.01, "**",
                                          ifelse(result_df$P_value < 0.05, "*", "")))
  
  # Add model diagnostics
  model_info <- list(
    model = cox_model,
    summary = cox_summary,
    result_table = result_df,
    n = nrow(multi_data),
    events = sum(multi_data[[status_var]]),
    method = method,
    formula = formula_str,
    concordance = cox_summary$concordance,
    loglik = cox_summary$loglik,
    AIC = AIC(cox_model)
  )
  
  # Test proportional hazards assumption
  ph_test <- tryCatch({
    cox.zph(cox_model)
  }, error = function(e) {
    warning("Failed to test proportional hazards assumption: ", e$message)
    return(NULL)
  })
  
  if (!is.null(ph_test)) {
    model_info$ph_test <- ph_test
  }
  
  attr(result_df, "model_info") <- model_info
  
  return(result_df)
}

#' Print comprehensive multivariate analysis results
#' @export
#'
#' @param multi_result Result from multivariate_cox function
#' @param show_diagnostics Whether to show model diagnostics (default: TRUE)
#'
#' @return No return value, prints results to console
print_multivariate_results <- function(multi_result, show_diagnostics = TRUE) {
  
  model_info <- attr(multi_result, "model_info")
  
  cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
  cat("MULTIVARIATE COX REGRESSION ANALYSIS\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n", sep = "")
  
  cat("Model Information:\n")
  cat("  Method:", model_info$method, "\n")
  cat("  Formula:", model_info$formula, "\n")
  cat("  Number of patients:", model_info$n, "\n")
  cat("  Number of events:", model_info$events, "\n")
  cat("  Concordance (C-index):", round(model_info$concordance[1], 4), 
      " (SE:", round(model_info$concordance[2], 4), ")\n")
  cat("  AIC:", round(model_info$AIC, 2), "\n")
  cat("  Log-likelihood:", round(model_info$loglik[2], 2), "\n\n")
  
  cat("Results:\n")
  print(multi_result[, c("Variable", "HR_95CI", "P_value", "Significance")], row.names = FALSE)
  
  if (show_diagnostics && !is.null(model_info$ph_test)) {
    cat("\n", paste(rep("-", 60), collapse = ""), "\n", sep = "")
    cat("Proportional Hazards Assumption Test (Schoenfeld residuals):\n")
    print(model_info$ph_test)
    
    if (min(model_info$ph_test$table[, "p"]) < 0.05) {
      cat("\n  WARNING: Proportional hazards assumption violated for some variables (p < 0.05)\n")
    } else {
      cat("\n  Proportional hazards assumption satisfied (all p > 0.05)\n")
    }
  }
  
  cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
}

#' Complete workflow: univariate -> multivariate analysis
#' @export
#'
#' @param data Data frame
#' @param time_var Time variable name
#' @param status_var Status variable name
#' @param covariates All covariates to test in univariate analysis
#' @param p_threshold P-value threshold for including variables in multivariate analysis
#' @param cut_method Method for dichotomizing continuous variables in univariate analysis
#' @param multivariate_method Method for multivariate analysis: "enter", "forward", "backward", "stepwise"
#' @param ... Additional arguments passed to uni_cox_table and multivariate_cox
#'
#' @return A list containing univariate results, significant variables, and multivariate results
complete_analysis <- function(data, time_var, status_var, covariates,
                              p_threshold = 0.05,
                              cut_method = "median",
                              multivariate_method = "enter",
                              ...) {
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
  cat("COMPLETE SURVIVAL ANALYSIS WORKFLOW\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n", sep = "")
  
  # Step 1: Univariate analysis
  cat("STEP 1: Univariate Cox Regression Analysis\n")
  cat(paste(rep("-", 70), collapse = ""), "\n", sep = "")
  
  uni_result <- uni_cox_table(data, time_var, status_var, covariates, 
                              cut_method = cut_method, ...)
  
  cat("\n✓ Univariate analysis completed\n")
  cat("  Total patients:", attr(uni_result, "total_patients"), "\n")
  cat("  Total events:", attr(uni_result, "total_events"), "\n")
  cat("  Overall median PFS:", attr(uni_result, "overall_median_pfs"), "\n\n")
  
  # Step 2: Get significant variables
  cat("STEP 2: Identifying Significant Variables\n")
  cat(paste(rep("-", 70), collapse = ""), "\n", sep = "")
  
  sig_vars <- get_significant_vars(uni_result, p_threshold = p_threshold)
  
  if (length(sig_vars) == 0) {
    cat("\n⚠ No significant variables found with p <", p_threshold, "\n")
    cat("Multivariate analysis will not be performed.\n")
    return(list(univariate = uni_result, 
                significant_vars = character(0),
                multivariate = NULL))
  }
  
  cat("\n✓ Found", length(sig_vars), "significant variable(s) with p <", p_threshold, ":\n")
  for (i in seq_along(sig_vars)) {
    cat("  ", i, ".", sig_vars[i], "\n")
  }
  cat("\n")
  
  # Step 3: Multivariate analysis
  cat("STEP 3: Multivariate Cox Regression Analysis\n")
  cat(paste(rep("-", 70), collapse = ""), "\n", sep = "")
  
  multi_result <- multivariate_cox(data, time_var, status_var, sig_vars,
                                   result_table = uni_result,
                                   method = multivariate_method)
  
  cat("\n✓ Multivariate analysis completed\n")
  
  # Display results
  print_multivariate_results(multi_result, show_diagnostics = TRUE)
  
  # Return all results
  return(list(
    univariate = uni_result,
    significant_vars = sig_vars,
    multivariate = multi_result,
    parameters = list(
      p_threshold = p_threshold,
      cut_method = cut_method,
      multivariate_method = multivariate_method
    )
  ))
}

# ====================== Complete usage example ======================
if(F){
  if (interactive()) {
    library(survival)
    library(dplyr)
    
    # Prepare data
    data(lung)
    test_data <- lung %>%
      mutate(PFS = time,
             PFS_status = ifelse(status == 1, 0, 1)) %>%  # convert to 0=censored, 1=event
      select(PFS, PFS_status, age, sex, ph.ecog, pat.karno, meal.cal)
    
    # Method 1: Step-by-step analysis
    
    # Step 1: Univariate analysis
    uni_result <- uni_cox_table(test_data, "PFS", "PFS_status",
                                covariates = c("age", "sex", "ph.ecog", "pat.karno", "meal.cal"),
                                cut_method = "median",
                                add_continuous = TRUE)
    
    # View univariate results
    print(uni_result)
    
    # Step 2: Get significant variables
    sig_vars <- get_significant_vars(uni_result, p_threshold = 0.1)
    cat("\nSignificant variables (p < 0.1):", paste(sig_vars, collapse = ", "), "\n")
    
    # Step 3: Multivariate analysis (enter all significant variables)
    multi_result <- multivariate_cox(test_data, "PFS", "PFS_status", 
                                     sig_vars, uni_result, 
                                     method = "enter")
    
    # Print multivariate results
    print_multivariate_results(multi_result)
    
    # Method 2: Stepwise selection
    multi_stepwise <- multivariate_cox(test_data, "PFS", "PFS_status",
                                       sig_vars, uni_result,
                                       method = "stepwise")
    
    print_multivariate_results(multi_stepwise)
    
    # Method 3: Complete automated workflow
    complete_results <- complete_analysis(test_data, "PFS", "PFS_status",
                                          covariates = c("age", "sex", "ph.ecog", "pat.karno", "meal.cal"),
                                          p_threshold = 0.1,
                                          cut_method = "median",
                                          multivariate_method = "stepwise")
    
    # Access individual components
    uni_table <- complete_results$univariate
    sig_vars_auto <- complete_results$significant_vars
    multi_table <- complete_results$multivariate
    
    # Get variable names with p < 0.05 (more stringent)
    sig_vars_05 <- get_significant_vars(uni_result, p_threshold = 0.05)
    cat("\nSignificant variables (p < 0.05):", paste(sig_vars_05, collapse = ", "), "\n")
  }
}
