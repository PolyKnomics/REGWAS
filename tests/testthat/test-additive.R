regwas <- Regwas$new("../snps_goudenstandaard_n85.bgen", "../pheno_goudenstandaard_n85_agesex.txt")

test_that("model_encoding", {
  regwas$loadChunk(1)
  regwas$SetCovariates("rs16984366")

  goudstandaard <- read.table("../bladder2_SNPs_covar.txt", header = T)
  variant_coded <- as.integer(goudstandaard$NBCS_SNP_genotyped == 2) * 2 + as.integer(goudstandaard$NBCS_SNP_genotyped == 1)

  goudstandaard <- cbind(goudstandaard, variant_coded)
  expect_equal(goudstandaard$variant_coded, regwas$reg_data$variant_coded)
})


test_that("overview normal", {
  regwas$set_default_model()
  re_results <- regwasworker(regwas, 1)
  expect_equal(class(re_results), "data.frame")
  goudstandaard <- read.table("../bladder2_SNPs_covar.txt", header = T)
  z <- summary(coxph(formula = Surv(tstart, tstop, status) ~ NBCS_SNP_genotyped +
    +age + sex, data = goudstandaard, cluster = id))

  expect_equal(z$coefficients["NBCS_SNP_genotyped", "Pr(>|z|)"], re_results["rs16984366", "pvalue"])
  expect_equal(z$coefficients["NBCS_SNP_genotyped", "Pr(>|z|)"], 0.0817, tolerance = .0001)


  expect_equal(z$conf.int["NBCS_SNP_genotyped", "lower .95"], re_results["rs16984366", "HRLower95"], tolerance = .00001)
  expect_equal(z$conf.int["NBCS_SNP_genotyped", "lower .95"], 0.9550, tolerance = .0001)

  expect_equal(z$conf.int["NBCS_SNP_genotyped", "upper .95"], re_results["rs16984366", "HRUpper95"], tolerance = .00001)
  expect_equal(z$conf.int["NBCS_SNP_genotyped", "upper .95"], 2.178, tolerance = .001)

  expect_equal(z$conf.int["NBCS_SNP_genotyped", "exp(coef)"], re_results["rs16984366", "HR"], tolerance = .00001)
  expect_equal(z$conf.int["NBCS_SNP_genotyped", "exp(coef)"], 1.4423, tolerance = .0001)
  # Value is not present in output but can be calculated with beta
  expect_equal(z$conf.int["NBCS_SNP_genotyped", "exp(-coef)"], exp(-1 * (re_results["rs16984366", "beta"])), tolerance = .00001)
  expect_equal(z$conf.int["NBCS_SNP_genotyped", "exp(-coef)"], 0.6933, tolerance = .0001)

  m <- coxph(formula = Surv(tstart, tstop, status) ~ NBCS_SNP_genotyped + age + sex, data = goudstandaard, cluster = id)
  n <- coxph(formula = Surv(tstart, tstop, status) ~ age + sex, data = goudstandaard, cluster = id)

  #
  expect_equal(-457.50245, re_results["rs16984366", "logliknull"], tolerance = .00001)

  expect_equal(-454.05985, re_results["rs16984366", "loglik"], tolerance = .00001)

  expect_equal(m$loglik[2], re_results["rs16984366", "loglik"], tolerance = .00001)
  expect_equal(anova(m, n)$`P(>|Chi|)`[2], re_results["rs16984366", "LRT"], tolerance = .00001)

  expect_equal(anova(m, n)$`P(>|Chi|)`[2], 0.008691, tolerance = .00001)
})



test_that("overview normal_no_pheno", {
  # phenotype file without aditional covariates
  regwas_no_pheno <- Regwas$new("../sig_snp_n85.bgen", "../pheno_goudenstandaard_n85_no_pheno.txt")
  re_results <- regwasworker(regwas_no_pheno, 1)
  expect_equal(class(re_results), "data.frame")
})




test_that("overview normal_sig", {
  regwas_sig <- Regwas$new("../sig_snp_n85.bgen", "../pheno_goudenstandaard_n85_agesex.txt")
  regwas_sig$set_default_model()
  re_results <- regwasworker(regwas_sig, 1)
  expect_equal(class(re_results), "data.frame")
  goudstandaard <- read.table("../bladder2_SNPs_covar.txt", header = T)
  z <- summary(coxph(formula = Surv(tstart, tstop, status) ~ NBCS_SNP_genotyped +
    +age + sex, data = goudstandaard, cluster = id))


  expect_equal(re_results["rs142088393", "pvalue"], 0.000000000000272583, tolerance = .0000000001)


  expect_equal(3.32737, re_results["rs142088393", "HRLower95"], tolerance = .00001)

  expect_equal(8.03239, re_results["rs142088393", "HRUpper95"], tolerance = .00001)

  expect_equal(5.16979, re_results["rs142088393", "HR"], tolerance = .00001)
  # Value is not present in output but can be calculated with beta
  expect_equal(exp(-1 * 1.64283), exp(-1 * (re_results["rs142088393", "beta"])), tolerance = .00001)

  expect_equal(0.21124, re_results["rs142088393", "se"], tolerance = .00001)
  expect_equal(0.22482, re_results["rs142088393", "RobustSE"], tolerance = .00001)

  expect_equal(-457.50245, re_results["rs142088393", "logliknull"], tolerance = .00001)

  expect_equal(-417.37407, re_results["rs142088393", "loglik"], tolerance = .00001)
})

test_that("overview coxme", {
  regwas$set_coxme()
  re_results <- regwasworker(regwas, 1)
  expect_equal(class(re_results), "data.frame")

  goudstandaard <- read.table("../bladder2_SNPs_covar.txt", header = T)
  z <- (coxme::coxme(Surv(tstart, tstop, status) ~ NBCS_SNP_genotyped + age + sex + (1 | id), data = goudstandaard))
  beta <- z$coefficients
  nvar <- length(beta)
  nfrail <- nrow(z$var) - nvar
  se <- sqrt(diag(as.matrix(z$var))[nfrail + 1:nvar])
  names(se) <- names(beta)
  p <- 1 - pchisq((beta / se)^2, 1)

  # p-value
  expect_equivalent(p["NBCS_SNP_genotyped"], re_results["rs16984366", "pvalue"])
  expect_equivalent(p["NBCS_SNP_genotyped"], 0.064, tolerance = 0.0005)

  # c("coef", "exp(coef)","se(coef)", "z", "p"
  # beta, exp(beta), se, round(beta/se,2),1 - pchisq((beta/ se)^2, 1)


  # se(coef)
  expect_equivalent(se["NBCS_SNP_genotyped"], re_results["rs16984366", "se"], tolerance = .00001)
  expect_equivalent(se["NBCS_SNP_genotyped"], 0.23945113, tolerance = .00001)

  expect_equivalent(beta["NBCS_SNP_genotyped"], re_results["rs16984366", "beta"], tolerance = .00001)
  expect_equivalent(beta["NBCS_SNP_genotyped"], 0.442973116, tolerance = .000001)


  m <- coxme::coxme(Surv(tstart, tstop, status) ~ NBCS_SNP_genotyped + age + sex + (1 | id), data = goudstandaard)
  n <- coxme::coxme(Surv(tstart, tstop, status) ~ age + sex + (1 | id), data = goudstandaard)

  expect_equal(anova(m, n)$`P(>|Chi|)`[2], re_results["rs16984366", "LRT"], tolerance = .00001)
  expect_equivalent(anova(m, n)$`P(>|Chi|)`[2], 0.04641, tolerance = .00001)

  expect_equal(sqrt(m$vcoef[[1]][[1]]), re_results["rs16984366", "random_eff_stddev"])
})


test_that("overview frailty gaussion", {
  regwas <- Regwas$new("../sig_snp_n85.bgen", "../pheno_goudenstandaard_n85_agesex.txt")

  regwas$set_frailty.gaussian()
  re_results <- regwasworker(regwas, 1)
  expect_equal(class(re_results), "data.frame")
  print(re_results["rs16984366", "RobustSE"])
  expect_equivalent(NULL, re_results["rs142088393", "RobustSE"])
  expect_equivalent(0.24814, re_results["rs142088393", "se"], tolerance = .00001)
})

test_that("overview frailty gamma", {
  regwas$set_frailty.gamma()
  expect_equal(class(regwasworker(regwas, 1)), "data.frame")
})
