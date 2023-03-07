Model <- R6::R6Class("model", list(
  formula = "",
  genetic_model = "additive",
  current_method = "default",
  covariants = "",
  variant_name = "variant_coded",
  regression = "",
  regressionnull = "",
  flipmaf = FALSE,
  initialize = function() {
    self$setMethod(self$current_method)
  },
  set_additive = function() {
    self$genetic_model <- "additive"
    self$variant_name <- "variant_coded"
    self$setMethod(self$current_method)
  },
  set_2df = function() {
    self$genetic_model <- "2df"
    self$variant_name <- "snp1+snp2"
    self$setMethod(self$current_method)
  },
  setMethod = function(method) {
    validmethods <- c("default", "frailty.gamma", "frailty.gaussian", "coxme")
    if (method %in% validmethods) {
      if (method == "default") {
        self$set_default_method()
      }
      else if (method == "frailty.gamma") {
        self$set_frailty.gamma()
      }
      else if (method == "frailty.gaussian") {
        self$set_frailty.gaussian()
      }
      else if (method == "coxme") {
        self$set_coxme()
      }
    } else {
      warning(paste0(
        method,
        " is non vaild method: should be one of ",
        validmethods
      ))
    }
  },
  set_covariants = function(covars) {
    self$covariants <- paste0("+", (paste(covars, sep = "+", collapse = NULL)))
    # correct covars if it is NA
    if (length(covars) == 1) {
      if (is.na(covars)) {
        self$covariants <- ""
      }
    }
    # in case no column names a character0 object is given
    if (identical(covars, character(0))) {
      self$covariants <- ""
    }

    # update formula etc with new covariates
    self$setMethod(self$current_method)
  },

  set_frailty.gaussian = function() {
    self$formula <- paste(c(
      "surv_object~",
      self$variant_name,
      self$covariants,
      "+frailty.gaussian(samples)"
    ),
    collapse = ""
    )
    self$regression <- quote(survival::coxph(cox, data = self$reg_data))
    self$regressionnull <- quote(survival::coxph(cox,
      data = self$reg_data,
      subset = non_na
    ))
    self$current_method <- "frailty.gaussian"
  },
  set_frailty.gamma = function() {
    self$formula <- paste(c(
      "surv_object~",
      self$variant_name,
      self$covariants,
      "+frailty.gamma(samples)"
    ),
    collapse = ""
    )
    self$regression <- quote(survival::coxph(cox, data = self$reg_data))
    self$regressionnull <- quote(survival::coxph(cox,
      data = self$reg_data,
      subset = non_na
    ))
    self$current_method <- "frailty.gamma"
  },
  set_default_method = function() {
    self$formula <- paste(c(
      "surv_object~",
      self$variant_name,
      self$covariants
    ),
    collapse = ""
    )

    self$regression <- quote(survival::coxph(cox,
      data = self$reg_data,
      cluster = samples
    ))
    self$regressionnull <- quote(survival::coxph(cox,
      data = self$reg_data,
      cluster = samples,
      subset = non_na
    ))
    self$current_method <- "default"
  },
  set_coxme = function() {
    self$formula <- paste(c(
      "surv_object~",
      self$variant_name,
      self$covariants,
      "+(1|samples)"
    ),
    collapse = ""
    )
    self$regression <- quote(coxme::coxme(cox, data = self$reg_data))
    self$regressionnull <- quote(coxme::coxme(cox,
      data = self$reg_data,
      subset = non_na
    ))

    self$current_method <- "coxme"
  }
))
