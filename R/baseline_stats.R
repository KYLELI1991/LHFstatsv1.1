#' Baseline characteristics table (3-line table format)
#'
#' @param data A data frame
#' @param vars Character vector of variable names to summarize (optional, defaults to all columns)
#' @param var_labels Named vector for variable labels (optional, e.g., c("age"="Age (years)", "gender"="Gender"))
#' @param continuous_format Format for continuous variables: "median_iqr" (default), "mean_sd", or "median_sd"
#' @param decimal_places Number of decimal places for continuous variables (default = 2)
#' @param percentage_decimal_places Number of decimal places for percentages (default = 1)
#'
#' @return A data frame suitable for 3-line table output
#' @export
#'
#' @examples
#' df <- data.frame(
#'   age = c(25, 30, NA, 40, 35),
#'   gender = c("Male", "Female", "Male", NA, "Female"),
#'   sbp = c(120, 130, 140, NA, 125)
#' )
#' result <- baseline_table(df, continuous_format = "mean_sd", decimal_places = 1)
#' print_table3(result)
#'
baseline_table <- function(data, vars = NULL, var_labels = NULL, 
                           continuous_format = "median_iqr", 
                           decimal_places = 2,
                           percentage_decimal_places = 1) {
  # If variables not specified, use all columns
  if (is.null(vars)) {
    vars <- names(data)
  }
  
  # Validate continuous_format parameter
  if (!continuous_format %in% c("median_iqr", "mean_sd", "median_sd")) {
    stop("continuous_format must be 'median_iqr', 'mean_sd', or 'median_sd'")
  }
  
  # Store results
  results <- data.frame(
    Variable = character(),
    N = character(),
    Statistic = character(),
    stringsAsFactors = FALSE
  )
  
  for (var in vars) {
    x <- data[[var]]
    
    # Skip non-existent variables
    if (is.null(x)) {
      warning(paste("Variable", var, "does not exist in data, skipped"))
      next
    }
    
    # Get variable display name
    if (!is.null(var_labels) && var %in% names(var_labels)) {
      var_display <- var_labels[var]
    } else {
      var_display <- var
    }
    
    # Determine variable type
    if (is.factor(x) || is.character(x)) {
      # Categorical variable
      valid <- x[!is.na(x)]
      n_total <- length(valid)
      n_missing <- sum(is.na(x))
      
      if (length(valid) == 0) {
        stat_str <- "All missing"
        results <- rbind(results, data.frame(
          Variable = var_display,
          N = "",
          Statistic = stat_str,
          stringsAsFactors = FALSE
        ))
      } else {
        tbl <- table(valid)
        
        # Add variable name row (as section header)
        results <- rbind(results, data.frame(
          Variable = var_display,
          N = "",
          Statistic = "",
          stringsAsFactors = FALSE
        ))
        
        # One row per category
        for (lev in names(tbl)) {
          n <- tbl[lev]
          pct <- n / n_total * 100
          pct_formatted <- sprintf(paste0("%.", percentage_decimal_places, "f%%"), pct)
          stat_str <- pct_formatted
          results <- rbind(results, data.frame(
            Variable = paste0("  ", lev),
            N = as.character(n),
            Statistic = stat_str,
            stringsAsFactors = FALSE
          ))
        }
        
        # Add missing row if any missing values exist
        if (n_missing > 0) {
          missing_pct <- n_missing / (n_total + n_missing) * 100
          missing_pct_formatted <- sprintf(paste0("%.", percentage_decimal_places, "f%%"), missing_pct)
          results <- rbind(results, data.frame(
            Variable = "  Missing",
            N = as.character(n_missing),
            Statistic = missing_pct_formatted,
            stringsAsFactors = FALSE
          ))
        }
      }
      
    } else if (is.numeric(x) || is.integer(x)) {
      # Continuous variable
      valid <- x[!is.na(x)]
      n_total <- length(valid)
      n_missing <- sum(is.na(x))
      
      if (length(valid) == 0) {
        stat_str <- "All missing"
        results <- rbind(results, data.frame(
          Variable = var_display,
          N = "",
          Statistic = stat_str,
          stringsAsFactors = FALSE
        ))
      } else {
        mean_val <- mean(valid)
        sd_val <- sd(valid)
        median_val <- median(valid)
        q1 <- quantile(valid, 0.25)
        q3 <- quantile(valid, 0.75)
        
        # Format based on continuous_format parameter
        decimal_format <- paste0("%.", decimal_places, "f")
        
        if (continuous_format == "mean_sd") {
          # Format: Mean ± SD
          mean_formatted <- sprintf(decimal_format, mean_val)
          sd_formatted <- sprintf(decimal_format, sd_val)
          stat_str <- sprintf("%s ± %s", mean_formatted, sd_formatted)
        } else if (continuous_format == "median_sd") {
          # Format: Median ± SD
          median_formatted <- sprintf(decimal_format, median_val)
          sd_formatted <- sprintf(decimal_format, sd_val)
          stat_str <- sprintf("%s ± %s", median_formatted, sd_formatted)
        } else {
          # Format: Median (IQR) - default
          median_formatted <- sprintf(decimal_format, median_val)
          q1_formatted <- sprintf(decimal_format, q1)
          q3_formatted <- sprintf(decimal_format, q3)
          stat_str <- sprintf("%s (%s-%s)", median_formatted, q1_formatted, q3_formatted)
        }
        
        results <- rbind(results, data.frame(
          Variable = var_display,
          N = as.character(n_total),
          Statistic = stat_str,
          stringsAsFactors = FALSE
        ))
        
        # Add missing row if any missing values exist
        if (n_missing > 0) {
          missing_pct <- n_missing / (n_total + n_missing) * 100
          missing_pct_formatted <- sprintf(paste0("%.", percentage_decimal_places, "f%%"), missing_pct)
          results <- rbind(results, data.frame(
            Variable = "  Missing",
            N = as.character(n_missing),
            Statistic = missing_pct_formatted,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  # Add total sample size row at the beginning
  total_n <- nrow(data)
  results <- rbind(
    data.frame(Variable = "Total N", N = as.character(total_n), Statistic = "", stringsAsFactors = FALSE),
    results
  )
  
  return(results)
}

#' Print 3-line table
#'
#' @param tbl Data frame returned by baseline_table function
#' @param file Output file path (optional, e.g., "table1.doc")
#'
#' @return Prints 3-line table to console
#' @export
#'
print_table3 <- function(tbl, file = NULL) {
  # Calculate maximum width for each column
  var_width <- max(nchar(tbl$Variable), nchar("Variable"))
  n_width <- max(nchar(tbl$N), nchar("N"))
  stat_width <- max(nchar(tbl$Statistic), nchar("Statistic"))
  
  # Create separator line
  sep_line <- paste0("+", paste(rep("-", var_width + 2), collapse = ""),
                     "+", paste(rep("-", n_width + 2), collapse = ""),
                     "+", paste(rep("-", stat_width + 2), collapse = ""), "+")
  
  # Header
  header <- sprintf("| %-*s | %-*s | %-*s |", 
                    var_width, "Variable", 
                    n_width, "N",
                    stat_width, "Statistic")
  
  # Print to console
  cat(sep_line, "\n")
  cat(header, "\n")
  cat(sep_line, "\n")
  
  # Print each row
  for (i in 1:nrow(tbl)) {
    row <- sprintf("| %-*s | %-*s | %-*s |",
                   var_width, tbl$Variable[i],
                   n_width, tbl$N[i],
                   stat_width, tbl$Statistic[i])
    cat(row, "\n")
  }
  
  cat(sep_line, "\n")
  
  # Output to file if specified
  if (!is.null(file)) {
    sink(file)
    cat(sep_line, "\n")
    cat(header, "\n")
    cat(sep_line, "\n")
    for (i in 1:nrow(tbl)) {
      row <- sprintf("| %-*s | %-*s | %-*s |",
                     var_width, tbl$Variable[i],
                     n_width, tbl$N[i],
                     stat_width, tbl$Statistic[i])
      cat(row, "\n")
    }
    cat(sep_line, "\n")
    sink()
    message("Table saved to: ", file)
  }
}

#' Export 3-line table to Word or Excel
#'
#' @param tbl Data frame returned by baseline_table function
#' @param file Output file path (supports .csv, .xlsx, .txt)
#'
#' @export
#'
export_table3 <- function(tbl, file) {
  # Check file extension
  ext <- tolower(tools::file_ext(file))
  
  if (ext == "csv") {
    write.csv(tbl, file, row.names = FALSE, fileEncoding = "UTF-8")
    message("Table exported to: ", file)
  } else if (ext == "xlsx") {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("openxlsx package required: install.packages('openxlsx')")
    }
    openxlsx::write.xlsx(tbl, file, rowNames = FALSE)
    message("Table exported to: ", file)
  } else if (ext == "txt") {
    write.table(tbl, file, row.names = FALSE, sep = "\t", fileEncoding = "UTF-8")
    message("Table exported to: ", file)
  } else {
    stop("Unsupported file format. Please use .csv, .xlsx, or .txt")
  }
}

# ==================== Usage Examples ====================
if (F) {
  # Create example data
  set.seed(123)
  df <- data.frame(
    age = c(rnorm(80, 55, 10), rep(NA, 20)),
    gender = factor(c(rep("Male", 45), rep("Female", 40), rep(NA, 15))),
    bmi = c(rnorm(85, 24.5, 3), rep(NA, 15)),
    sbp = c(rnorm(90, 130, 15), rep(NA, 10))
  )
  
  # Set variable labels (optional)
  labels <- c(
    age = "Age (years)",
    gender = "Gender",
    bmi = "Body mass index (kg/m²)",
    sbp = "Systolic blood pressure (mmHg)"
  )
  
  # Example 1: Median (IQR) with 2 decimal places
  cat("\n=== Median (IQR), 2 decimal places ===\n")
  result1 <- baseline_table(df, var_labels = labels, 
                            continuous_format = "median_iqr", 
                            decimal_places = 2,
                            percentage_decimal_places = 1)
  print_table3(result1)
  
  # Example 2: Mean ± SD with 1 decimal place
  cat("\n=== Mean ± SD, 1 decimal place ===\n")
  result2 <- baseline_table(df, var_labels = labels, 
                            continuous_format = "mean_sd", 
                            decimal_places = 1,
                            percentage_decimal_places = 1)
  print_table3(result2)
  
  # Example 3: Median ± SD with 2 decimal places
  cat("\n=== Median ± SD, 2 decimal places ===\n")
  result3 <- baseline_table(df, var_labels = labels, 
                            continuous_format = "median_sd", 
                            decimal_places = 2,
                            percentage_decimal_places = 1)
  print_table3(result3)
  
  # Example 4: Custom percentage decimals (e.g., 2 decimal places for percentages)
  cat("\n=== Mean ± SD with 2 decimal percentages ===\n")
  result4 <- baseline_table(df, var_labels = labels, 
                            continuous_format = "mean_sd", 
                            decimal_places = 2,
                            percentage_decimal_places = 2)
  print_table3(result4)
  
  # Export to CSV
  export_table3(result1, "table1_median_iqr.csv")
}
