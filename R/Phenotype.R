Phenotype <- R6::R6Class("Phenotype", list(
  surv_object = "",
  samples = "",
  other_pheno = "",
  orignal_filename = "",
  pheno = "",
  loadPhenoType = function(phenotype) {
    self$orignal_filename <- phenotype
    if (substr(phenotype, (nchar(phenotype) + 1) - 4, nchar(phenotype)) == "bgen") {
      warning("You are trying to load a *.bgen file as phenotype file. Try to swap the position of the arguments")
    }
    rawpheno <- read.table(phenotype, header = TRUE)
    # remove rowa with NA from phenotype file
    self$pheno <- na.omit(rawpheno)

    # check if needed columns are availble (assuming the first column
    # is sampple the sample, but naming might be non standard: id,
    # sample, sampleid, s_id, etc)
    survcolumns <- c(names(self$pheno)[1], "tstop", "tstart", "status")
    for (columnname in survcolumns) {
      if (!(columnname %in% colnames(self$pheno))) {
        stop(paste0(
          "column named ",
          columnname,
          " is missing in phenotype file"
        ))
      }
    }
    other_pheno_index <- !(names(self$pheno) %in% survcolumns)
    self$other_pheno <- data.frame(self$pheno[, other_pheno_index])
    names(self$other_pheno) <- names(self$pheno)[other_pheno_index]
    self$samples <- as.character(self$pheno[, 1])
    # create survival object
    self$surv_object <- survival::Surv(
      self$pheno$tstart,
      self$pheno$tstop,
      self$pheno$status
    )
  },
  numberOfSamples = function() {
    return(length(unique(self$samples)))
  },
  samples_with_events = function() {
    return(unique(self$pheno[self$pheno$status == 1, 1]))
  },
  samples_without_events = function() {
    one <- unique(self$pheno[self$pheno$status == 1, 1])
    zero <- unique(self$pheno[self$pheno$status == 0, 1])
    return(setdiff(zero, one))
  }
))
