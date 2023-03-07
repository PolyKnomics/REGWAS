GwasWriter <- R6::R6Class(
  "GwasWriter",
  list(final_filename = "", tmpfile = "", setFilename = function(filename) {
    self$final_filename <- filename
  }, initialize = function() {
    self$tmpfile <- tempfile(pattern = ".regwass", tmpdir = ".")
  }, get_tmp_filename = function(chunk, prefix = "") {
    return(sprintf(paste0(self$tmpfile, "-", prefix, "%06d"), chunk))
  }, genotypes_list_to_Text = function(gen_list) {
    unlist(lapply(gen_list, function(l) {
      if (is.null(unlist(l))) {
        return("NA")
      } else {
        return(paste0(round(l, 1), sep = "", collapse = "|"))
      }
    }))
  }, write = function(results) {
    chunk <- attr(results, "chunk_number")


    if (chunk == 1) {
      outputheader <- TRUE
    } else {
      outputheader <- FALSE
    }

    # writting output and log is mostlty duplicated code to prevent
    # spaghetti code in the future get gz filename for output, for error
    # log a normal file
    outputfile <- gzfile(self$get_tmp_filename(chunk), "w")


    # format pvalue as scientific: all numeric (non int) are rounded
    # with 5 decimals
    scientificcolumns <- c("pvalue", "pen_chisq", "pen_df", "pen_pvalue", "LRT")
    available_sci_column <- scientificcolumns[scientificcolumns %in% names(results)]
    for (sci_column in available_sci_column) {
      results[, sci_column] <- format(results[, sci_column], scientific = TRUE)
    }

    numeric_columns <- sapply(results, class) == "numeric"
    results[, numeric_columns] <- round(results[, numeric_columns], 5)
    list_columns <- sapply(results, class) == "list"
    results[list_columns] <- apply(
      results[list_columns], 2,
      self$genotypes_list_to_Text
    )
    write.table(results, outputfile,
      quote = FALSE,
      row.names = FALSE, col.names = outputheader, sep = "\t"
    )
    close(outputfile)



    # getfilename for output, for error log a normal file

    outputfile <- file(self$get_tmp_filename(chunk, prefix = ".regwasslog"), "w")
    log_results <- attr(results, "log")
    numeric_columns <- sapply(log_results, class) == "numeric"
    log_results[, numeric_columns] <- round(log_results[, numeric_columns], 5)
    list_columns <- sapply(log_results, class) == "list"
    log_results[list_columns] <- apply(
      log_results[list_columns], 2,
      self$genotypes_list_to_Text
    )
    write.table(log_results, outputfile,
      quote = FALSE,
      row.names = FALSE, col.names = outputheader
    )
    close(outputfile)
  }, mergeresults = function(nchunks) {
    # merge result files
    temp_file_names <- self$get_tmp_filename(1:nchunks)
    if (sum(!file.exists(temp_file_names)) > 0) {
      stop(paste0(
        "Missing temporary file ",
        temp_file_names[!file.exists(temp_file_names)]
      ))
    }
    system2("cat",
      args = temp_file_names,
      stdout = self$final_filename, stderr = FALSE
    )
    system2("rm",
      args = temp_file_names,
      stdout = FALSE, stderr = FALSE
    )
    print(paste0("merging done for result file:", self$final_filename))

    temp_logfile_names <- self$get_tmp_filename(1:nchunks,
      prefix = ".regwasslog"
    )
    # merge log file
    log_file <- paste0(c(
      substr(
        self$final_filename, 1,
        nchar(self$final_filename) - 2
      ),
      "log"
    ),
    collapse = ""
    )

    if (sum(!file.exists(temp_logfile_names)) > 0) {
      stop(paste0(
        "Missing temporary file ",
        temp_logfile_names[!file.exists(temp_logfile_names)]
      ))
    }
    system2("cat", args = temp_logfile_names, stdout = log_file, stderr = FALSE)
    system2("rm", args = temp_logfile_names, stdout = FALSE, stderr = FALSE)
    print(paste0("merging done for log file:", log_file))
  })
)
