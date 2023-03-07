#' Function to perform regression
#'
#' @export
#'
#' @param gwas Regwas object
#' @param chunk_number chunk number to analyse
#' @param verbose print out stats of a variant before regression
#' @return data.frame with results of regression
#' @import survival
#' @importFrom stats pchisq
#'
regwasworker <- function(gwas, chunk_number, verbose = FALSE) {
  gwas$loadChunk(chunk_number)
  results <- gwas$getInfoAllVariants()
  attr(results, "chunk_number") <- chunk_number
  surv_object <- gwas$p$surv_object # nolint
  cox <- stats::formula(gwas$model$formula) # nolint


  log_colnames <- c(
    "chromosome", "position", "rsid", "number_of_alleles", "allele0",
    "allele1", "gt_all", "gt_events", "gt_non_events", "MAF", "HWE",
    "N", "type", "call", "message"
  )
  log_df <- data.frame(matrix(ncol = length(log_colnames), nrow = 0))
  colnames(log_df) <- log_colnames



  i <- 0
  for (variant in row.names(results)) {
    gwas$SetCovariates(variant)

    # get genotype count for different samples
    results[variant, "gt_all"] <- list(gwas$getGenotypeCountsAll())
    results[variant, "gt_events"] <- list(gwas$getGenotypeCountsEvents())
    results[variant, "gt_non_events"] <- list(gwas$getGenotypeCountsNonEvents())


    results[variant, "MAF"] <- gwas$getMAF()
    results[variant, "HWE"] <- gwas$getHWE()
    results[variant, "N"] <- gwas$getN()

    # swap the allel naming if there is a MAF flip
    if (gwas$VariantMafFlip) {
      temp_allel <- results[variant, "allele0"]
      results[variant, "allele0"] <- results[variant, "allele1"]
      results[variant, "allele1"] <- temp_allel
    }
    # do actual regression and catch error when regression fails
    r <- tryCatch(
      {
        r <- gwas$regress()
      },
      error = function(cond) {
        message(paste("run into error during regression"))
        log_df[variant, log_colnames] <<- c(
          results[
            variant,
            log_colnames[1:12]
          ],
          "error",
          toString(cond$call),
          cond$message
        )
        return(NA)
      },
      warning = function(cond) {
        message(paste("run into Warning during regression"))
        log_df[variant, log_colnames] <<- c(
          results[
            variant,
            log_colnames[1:12]
          ],
          "warning",
          toString(cond$call),
          cond$message
        )
        return(NA)
      },
      finally = {

      }
    )


    if (!is.na(r)[1]) {
      ndf <- 1
      if ("coxph" %in% class(r)) {
        frailty_model <- "frail" %in% names(r)
        if (gwas$model$genetic_model == "2df") {
          ndf <- 2
          covnumber <- 1
          beta_SNP_A1A1 <- r$coefficients[covnumber]
          results[variant, "beta_A1A1"] <- beta_SNP_A1A1


          if (frailty_model) {
            se <- sqrt(diag(r$var)[covnumber])
            results[variant, "se_A1A1"] <- se
            working_se <- se
          } else {
            se <- sqrt(diag(r$naive.var)[covnumber])
            results[variant, "se_A1A1"] <- se
            robust_se <- sqrt(diag(r$var)[covnumber])
            results[variant, "RobustSE_A1A1"] <- robust_se
            working_se <- robust_se
          }


          expbeta <- exp(beta_SNP_A1A1)
          results[variant, "HR_A1A1"] <- expbeta
          results[variant, "HRLower95_A1A1"] <-
            exp(beta_SNP_A1A1 - 1.96 * working_se)
          results[variant, "HRUpper95_A1A1"] <-
            exp(beta_SNP_A1A1 + 1.96 * working_se)


          covnumber <- 2
          beta_SNP_A1A2 <- r$coefficients[covnumber]
          results[variant, "beta_A1A2"] <- beta_SNP_A1A2

          if (frailty_model) {
            se <- sqrt(diag(r$var)[covnumber])
            results[variant, "se_A1A2"] <- se
            working_se <- se
          } else {
            se <- sqrt(diag(r$naive.var)[covnumber])
            results[variant, "se_A1A2"] <- se
            robust_se <- sqrt(diag(r$var)[covnumber])
            results[variant, "RobustSE_A1A2"] <- robust_se
            working_se <- robust_se
          }


          expbeta <- exp(beta_SNP_A1A2)
          results[variant, "HR_A1A2"] <- expbeta
          results[variant, "HRLower95_A1A2"] <-
            exp(beta_SNP_A1A2 - 1.96 * working_se)
          results[variant, "HRUpper95_A1A2"] <-
            exp(beta_SNP_A1A2 + 1.96 * working_se)
        } else {
          # non2df model
          covnumber <- 1
          beta <- r$coefficients[covnumber]
          results[variant, "beta"] <- beta
          if (frailty_model) {
            se <- sqrt(diag(r$var)[covnumber])
            results[variant, "se"] <- se
            working_se <- se
          } else {
            se <- sqrt(diag(r$naive.var)[covnumber])
            results[variant, "se"] <- se
            robust_se <- sqrt(diag(r$var)[covnumber])
            results[variant, "RobustSE"] <- robust_se
            working_se <- robust_se
          }
          z_value <- beta / working_se
          results[variant, "Z"] <- z_value
          pvalue <- stats::pchisq((beta / working_se)^2, 1, lower.tail = FALSE)
          results[variant, "pvalue"] <- pvalue
          expbeta <- exp(beta)
          results[variant, "HR"] <- expbeta
          results[variant, "HRLower95"] <- exp(beta - 1.96 * working_se)
          results[variant, "HRUpper95"] <- exp(beta + 1.96 * working_se)
        }
        logliknull <- gwas$regressnull()
        results[variant, "logliknull"] <- logliknull
        results[variant, "loglik"] <- r$loglik[2]

        results[variant, "LRT"] <- stats::pchisq(2 * (r$loglik[2] -
          logliknull),
        ndf,
        lower.tail = FALSE
        )
      }
      else {
        # for coxme (model is not coxph)
        if (gwas$model$genetic_model == "2df") {
          ndf <- 2
          covnumber <- 1
          beta_SNP_A1A1 <- r$coefficients[covnumber]
          results[variant, "beta_A1A1"] <- beta_SNP_A1A1
          n_frailty <- nrow(r$var) - length(r$coefficients)
          se_SNP_A1A1 <- sqrt(r$var[
            n_frailty + covnumber,
            n_frailty + covnumber
          ])
          results[variant, "se_A1A1"] <- se_SNP_A1A1
          z_value_SNP_A1A1 <- beta_SNP_A1A1 / se_SNP_A1A1
          results[variant, "z_value_A1_A1"] <- z_value_SNP_A1A1
          pvalue_SNP_A1A1 <- stats::pchisq((beta_SNP_A1A1 / se_SNP_A1A1)^2,
            1,
            lower.tail = FALSE
          )
          results[variant, "pvalue_SNP_A1A1"] <- pvalue_SNP_A1A1
          random_eff_stddev_SNP_A1A1 <- sqrt(r$var[covnumber, covnumber])
          results[variant, "random_eff_stddev_SNP_A1A1"] <-
            random_eff_stddev_SNP_A1A1

          covnumber <- 2
          beta_A1A2 <- r$coefficients[covnumber]
          results[variant, "beta_A1A2"] <- beta_A1A2

          se_A1A2 <- sqrt(r$var[n_frailty + covnumber, n_frailty + covnumber])
          results[variant, "se_A1A2"] <- se_A1A2
          z_value_A1A2 <- beta_A1A2 / se_A1A2
          results[variant, "z_value_A1A2"] <- z_value_A1A2
          pvalue_A1A2 <- stats::pchisq((beta_A1A2 / se_A1A2)^2,
            1,
            lower.tail = FALSE
          )
          results[variant, "pvalue_A1A2"] <- pvalue_A1A2
          random_eff_stddev_A1A2 <- sqrt(r$var[covnumber, covnumber])
          results[variant, "random_eff_stddev_A1A2"] <- random_eff_stddev_A1A2
        }
        else {
          ndf <- 1
          covnumber <- 1
          beta <- r$coefficients[covnumber]
          results[variant, "beta"] <- beta


          n_frailty <- nrow(r$var) - length(r$coefficients)
          se <- sqrt(r$var[n_frailty + covnumber, n_frailty + covnumber])
          results[variant, "se"] <- se
          z_value <- beta / se
          results[variant, "z_value"] <- z_value
          pvalue <- stats::pchisq((beta / se)^2, 1, lower.tail = FALSE)
          results[variant, "pvalue"] <- pvalue
          random_eff_stddev <- sqrt(r$vcoef[[1]][[1]])
          results[variant, "random_eff_stddev"] <- random_eff_stddev

          expbeta <- exp(beta)
          results[variant, "HR"] <- expbeta
          results[variant, "HRLower95"] <- exp(beta - 1.96 * se)
          results[variant, "HRUpper95"] <- exp(beta + 1.96 * se)
        }
        # stuff both needed for 2df and other models


        loglik <- r$loglik
        loglik[3] <- loglik[3] + r$penalty

        int_chisq <- abs(2 * (loglik[1] - loglik[2]))
        int_df <- r$df[1]
        int_pvalue <- pchisq(int_chisq, int_df, lower.tail = FALSE)

        pen_chisq <- abs(2 * (loglik[1] - loglik[3]))
        pen_df <- r$df[2]
        pen_pvalue <- pchisq(pen_chisq, pen_df, lower.tail = FALSE)

        results[variant, "int_chisq"] <- int_chisq
        results[variant, "int_df"] <- int_df
        results[variant, "int_pvalue"] <- int_pvalue

        results[variant, "pen_chisq"] <- pen_chisq
        results[variant, "pen_df"] <- pen_df
        results[variant, "pen_pvalue"] <- pen_pvalue

        logliknull <- gwas$regressnull()
        results[variant, "logliknull"] <- logliknull
        results[variant, "loglik"] <- r$loglik[2]
        results[variant, "LRT"] <- stats::pchisq(2 * (r$loglik[2] -
          logliknull),
        ndf,
        lower.tail = FALSE
        )
      }
    } # end R is NA

    # info is <0 when the fitting did not converge lrt<-2*(r$loglik[2]-logliknull)
    # N <- r$N



    # block is needed to give printed status updates during printing
    if (i %% 250 == 0) {
      print(paste0("regresions ", i, " done"))
      utils::flush.console()
    }
    i <- i + 1
  }
  attr(results, "log") <- log_df
  return(results)
}
