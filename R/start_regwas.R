#' Title
#'
#' @param bgen BGEN filename/path
#' @param pheno Phenotype filename/path
#' @param model genetic model possible options: "additive", "dominant", "recessive", "2df"
#' @param method  regression model to be used: "robustcoxph"(default), "frailty.gamma", "frailty.gaussian", "coxme"
#' @param flipmaf Flip MAF(only useful with 2df)
#' @param threads Amount of threads to use (default: 1)
#' @param output output file name
#' @param region region specify region to process. The notation is in samtools region format e.g. for chromosome 22 "22" or "chr22"  (depends on own annotation) and for a region of chromosome 22: "22:1000-3000"
#'
#' @return
#' @export
#' @import doParallel foreach
#' @importFrom parallel makeCluster
#' @examples
start_regwas <- function(bgen,
                         pheno,
                         model = "additive",
                         method = "AGmodel",
                         flipmaf = FALSE,
                         threads = 1,
                         output = NA,
                         region = NA) {
  methods <- c(
    "AGmodel", "frailty.gamma",
    "frailty.gaussian", "frailty.gaussian.coxme"
  )
  models <- c("additive", "dominant", "recessive", "2df")

  if (!method %in% methods) {
    stop(
      "Invalid 'method' value. Should be one of",
      paste(methods, sep = ",", collapse = " ")
    )
  }
  if (!model %in% models) {
    stop(
      "Invalid 'model' value. Should be one of :\n",
      paste(models, sep = ",", collapse = " ")
    )
  }



  regwas_object <- Regwas$new(bgen, pheno, region = region)

  switch(method,
    "AGmodel" = regwas_object$set_default_model(),
    "frailty.gamma" = regwas_object$set_frailty.gamma(),
    "frailty.gaussian" = regwas_object$set_frailty.gaussian(),
    "frailty.gaussian.coxme" = regwas_object$set_coxme()
  )

  switch(model,
    "additive" = regwas_object$setModel.additive(),
    "2df" = regwas_object$setModel.2df(),
    "dominant" = stop("dominant model not implemented"),
    "recessive" = stop("recessive model not implemented")
  )



  regwas_object$setFlipMaf(flipmaf)
  nchunks <- regwas_object$g$getNumberChunks()
  # do not spawn more threads then chunks are available: prevent overhead
  if (threads > nchunks) {
    threads <- nchunks
  }
  regwas_object$setFilename(output)
  print(regwas_object)

  # register the cores to use
  cl <- parallel::makeCluster(threads)
  doParallel::registerDoParallel(cl)
  chunknr <- 0
  foreach(chunknr = c(1:nchunks)) %dopar% {
    regwas_object$gwaswriter$write(REGWAS::regwasworker(
      regwas_object,
      chunknr
    ))
  }

  regwas_object$gwaswriter$mergeresults(nchunks)
}
