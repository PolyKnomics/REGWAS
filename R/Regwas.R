#' GWAS with Recurrent events
#'
#' Recurrent events
#'
#' @export
#'
#' @examples
#' # load files
#' pheno <- system.file("extdata", "big_phenotype.tsv", package = "REGWAS")
#' geno <- system.file("extdata", "example.8bits.bgen", package = "REGWAS")
#' regwas <- Regwas$new(geno, pheno)
#' # calculate each chunk
#' for (chunknr in c(1:regwas$g$getNumberChunks())) {
#'   regwas$gwaswriter$write(regwasworker(regwas, chunknr))
#' }
#' # merge all the results chunks of the regresion
#' regwas$gwaswriter$mergeresults(regwas$g$getNumberChunks())
Regwas <- R6::R6Class("Regwas",
  public = list(
    #' @field p phenotype object
    p = Phenotype$new(),
    #' @field g genotype object
    g = Genotype$new(),
    #' @field gwaswriter object to write results to disk
    gwaswriter = GwasWriter$new(),
    #' @field current_chunk number of current loaded chunk
    current_chunk = "",
    #' @field gen_data loaded data from a chunk
    gen_data = "",
    #' @field model object
    model = Model$new(),
    reg_data = "",
    NA_samples = "",
    variant = "",
    VariantMafFlip = "",

    #' @description
    #' Create a new Regwas object.
    #' @param genotype BGEN filename/path
    #' @param phenotype Phenotype filename/path
    #' @return A new `Regwas` object.
    initialize = function(genotype, phenotype, region = NA) {
      self$p <- Phenotype$new()
      self$p$loadPhenoType(phenotype)
      self$g <- Genotype$new()
      print("set genotype")
      self$g$setGenotype(genotype)
      self$checkSampleID()
      print("start gwas writer")
      self$gwaswriter <- GwasWriter$new()

      self$VariantMafFlip <- FALSE

      self$g$determ_chunks(self$p$numberOfSamples(), region = region)
      self$setFilename()

      self$model <- Model$new()
      self$model$set_covariants(colnames(self$p$other_pheno))
    },

    #' @description
    #' Check if sampleid in the phenotype file is also present in genotype file
    #'
    checkSampleID = function() {
      a <- DBI::dbConnect(RSQLite::SQLite(), self$g$indexfile)
      print(self$g$indexfile)
      first_variant <-
        DBI::dbGetQuery(
          a,
          "SELECT chromosome,position FROM Variant LIMIT 1;"
        )
      DBI::dbDisconnect(a)
      regions <-
        data.frame(
          chromosome = first_variant$chromosom,
          start = first_variant$position,
          end = first_variant$position + 1,
          stringsAsFactors = FALSE
        )
      first_variant_data <- rbgen::bgen.load(self$g$genotypefile, regions)

      missing <- self$p$samples %in% first_variant_data$samples
      if (sum(!missing) > 1) {
        missingsample <- paste0(unique(self$p$samples[!missing]),
          collapse = " "
        )
        stop(paste0(c(
          "could not find all sample ids from phenotype file in the genotype file. The following are missing",
          missingsample
        )))
      }
    },
    #' @description
    #' Load covariates for the regression.
    #' The covariates included coded genotype data (coding is based on
    #' set Model) and all other covariates given in phenotype file
    #' @param  variant variantnumber or name.
    SetCovariates = function(variant) {
      self$variant <- variant
      variantgenotype <- self$gen_data$data[variant, self$p$samples, ]
      if (self$model$genetic_model == "additive") {
        variant_coded <- 2 * variantgenotype[, 1] + variantgenotype[, 2]
      } else if (self$model$genetic_model == "dominant") {
        variant_coded <- variantgenotype[, 1] + variantgenotype[, 2]
      } else if (self$model$genetic_model == "recessive") {
        variant_coded <- variantgenotype[, 1]
      } else if (self$model$genetic_model == "overdominant") {
        variant_coded <- variantgenotype[, 2]
      } else if (self$model$genetic_model == "2df") {
        if (self$model$flipmaf == TRUE) {
          # get counts of genotypes
          n_allel1 <- sum(2 * variantgenotype[, 1] + variantgenotype[, 2])
          n_allel2 <- sum(2 * variantgenotype[, 3] + variantgenotype[, 2])
          if (n_allel1 < n_allel2) {
            self$VariantMafFlip <- TRUE
          }
          else {
            self$VariantMafFlip <- FALSE
          }
        }
        if (self$VariantMafFlip == FALSE) {
          variant_coded <- data.frame(
            snp1 = variantgenotype[, 2],
            snp2 = variantgenotype[, 1]
          )
        } else {
          variant_coded <- data.frame(
            snp1 = variantgenotype[, 2],
            snp2 = variantgenotype[, 3]
          )
        }
      } else {
        stop("No valid model selected")
      }
      # do not remove NA from this data frame since then this is not
      # anymore sync with the Surv object
      self$reg_data <- data.frame(variant_coded,
        self$p$other_pheno,
        samples = self$p$samples
      )
      self$NA_samples <- unique(self$regdata$samples[is.na(self$reg_data)])
    },

    #' @description
    #'
    setModel.additive = function() {
      self$model$set_additive()
    },
    setModel.2df = function() {
      self$model$set_2df()
    },

    #' @description
    #' Use frailty with a gaussian distrubution
    set_frailty.gaussian = function() {
      self$model$set_frailty.gaussian()
    },

    #' @description
    #' Use frailty with a gamma distrubution
    set_frailty.gamma = function() {
      self$model$set_frailty.gamma()
    },

    #' @description
    #' Use coxme
    set_coxme = function() {
      self$model$set_coxme()
    },

    #' @description
    #' Use non robust module
    set_default_model = function() {
      self$model$set_default_method()
    },
    set_2df_model = function() {
      self$model$set_2df()
    },
    getGenotypeCounts = function(variant, samples) {
      counts <- list(colSums(self$gen_data$data[variant, samples, ],
        na.rm = TRUE
      ))[1]
      if (self$VariantMafFlip) {
        fipcounts <- rev(counts)
        names(fipcounts) <- names(counts)
        counts <- fipcounts
      }
      return(counts)
    },
    getGenotypeCountsAll = function() {
      # variant,unique(regdata$samples)
      counts <- list(colSums(self$gen_data$data[
        self$variant,
        unique(self$reg_data$samples),
      ],
      na.rm = TRUE
      ))[1]
      if (self$VariantMafFlip) {
        fipcounts <- rev(counts)
        names(fipcounts) <- names(counts)
        counts <- fipcounts
      }
      return(counts)
    },

    getGenotypeCountsEvents = function() {
      counts <- list(colSums(self$gen_data$data[
        self$variant,
        self$p$samples_with_events(),
      ],
      na.rm = TRUE
      ))[1]
      if (self$VariantMafFlip) {
        fipcounts <- rev(counts)
        names(fipcounts) <- names(counts)
        counts <- fipcounts
      }
      return(counts)
    },

    getGenotypeCountsNonEvents = function() {
      counts <- list(colSums(self$gen_data$data[
        self$variant,
        self$p$samples_without_events(),
      ],
      na.rm = TRUE
      ))[1]
      if (self$VariantMafFlip) {
        fipcounts <- rev(counts)
        names(fipcounts) <- names(counts)
        counts <- fipcounts
      }
      return(counts)
    },
    #' @description
    #' get number of uniq samples used in the current regression
    getN = function() {
      return(length(unique(self$reg_data$samples)) - length(self$NA_samples))
    },
    #' @description Do the actual regression on current loaded variant
    regress = function() {
      surv_object <- self$p$surv_object
      cox <- stats::formula(self$model$formula)
      return(eval(self$model$regression))
    },
    regressnull = function() {
      surv_object <- self$p$surv_object
      cox <- stats::formula(self$model$formula)
      # remove variant coded and snp1 and snp2 if present
      cox <- update(cox, ~ . - variant_coded - snp1 - snp2)
      # find find NA in genotypes (only first column is also ok with 2df)
      non_na <- is.na(self$reg_data[, 1]) == FALSE
      # set a subset
      regnull <- eval(self$model$regressionnull)
      return(regnull$loglik[2])
    },
    #' @description
    #' return p-value of Hardy-Weinberg equilibrium based on chi square statistics
    #' @param gt genotype counts (must be 3 observations)
    HWEchisq = function(gt) {
      stopifnot(length(gt) == 3)
      n <- sum(gt)
      p <- (2 * gt[1] + gt[2]) / (2 * n)
      q <- 1 - p
      e <- c(p^2 * n, 2 * p * q * n, q^2 * n)
      return(pchisq(sum((gt - e)^2 / e), df = 1, lower.tail = FALSE))
    },
    #' @description
    #' return Minor allele frequency
    #' @param gt genotype counts (must be 3 observations)
    MAF = function(gt) {
      stopifnot(length(gt) == 3)
      maf <- (gt[1] + 0.5 * gt[2]) / sum(gt)
      return(min(maf, 1 - maf))
    },
    #' @description
    #' return the  MAF from samples and variant
    getMAF = function() {
      gt <- colSums(self$gen_data$data[
        self$variant,
        unique(self$reg_data$samples),
      ],
      na.rm = TRUE
      )
      return(self$MAF(gt))
    },
    #' @description
    #' return the  HWE from selected samples and variant
    getHWE = function() {
      gt <- colSums(self$gen_data$data[
        self$variant,
        unique(self$reg_data$samples),
      ],
      na.rm = TRUE
      )
      return(self$HWEchisq(gt))
    },

    #' @description
    #' get all information of the variants in current loaded chunk
    getInfoAllVariants = function() {
      return(self$gen_data$variant)
    },
    #' @description
    #' Load genotype data specified by chunknumber
    #' @param chunknr chunknumber
    loadChunk = function(chunknr) {
      self$gen_data <- rbgen::bgen.load(self$g$genotypefile,
        self$g$chunks[chunknr, ],
        samples = unique(self$p$samples)
      )
      self$current_chunk <- chunknr
    },
    #' @description
    #' Set resulting file name
    #' @param filename set filename, if not set, this will be a
    #' {phenotype filename}-{model}.tab.gz

    setFilename = function(filename = NA) {
      if (is.na(filename)) {
        self$gwaswriter$setFilename(
          paste0(
            tools::file_path_sans_ext(basename(self$p$orignal_filename)),
            "-",
            self$model$genetic_model,
            "_",
            self$model$current_method,
            ".tab.gz"
          )
        )
      } else {
        self$gwaswriter$setFilename(filename)
      }
    },
    #' @description
    #' Set flipmaf for 2df model Flip the allell coding if allel
    #' frequency of seccond allel is more than 0.5
    #' @param True of False

    setFlipMaf = function(flip) {
      if (is.logical(flip)) {
        self$model$flipmaf <- flip
      } else {
        stop("flip is not true or false")
      }
    },
    print = function(...) {
      cat("REGWAS current options: \n")
      cat("  Genetic Model:           ", self$model$genetic_model,
        "\n",
        sep = ""
      )
      cat("  Method:                  ", self$model$current_method,
        "\n",
        sep = ""
      )
      cat("  Model formula:           ", self$model$formula,
        "\n",
        sep = ""
      )

      cat("  FlipMAF (only 2df model):", self$model$flipmaf,
        "\n\n",
        sep = ""
      )

      cat(" Phenotype file:           ", self$p$orignal_filename,
        "\n",
        sep = ""
      )
      cat(" Genotype file:            ", self$g$genotypefile,
        "\n",
        sep = ""
      )

      cat(" Number of chunks:         ", self$g$getNumberChunks(),
        "\n",
        sep = ""
      )
      invisible(self)
    }
  )
)
