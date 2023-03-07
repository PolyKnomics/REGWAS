regwas <- Regwas$new("../snps_goudenstandaard_n85.bgen", "../pheno_goudenstandaard_n85_agesex.txt")


test_that("model_encoding", {
  regwas$set_2df_model()
  regwas$loadChunk(1)
  regwas$SetCovariates("rs16984366")

  goudstandaard <- read.table("../bladder2_SNPs_covar.txt", header = T)
  snp1 <- as.integer(goudstandaard$NBCS_SNP_genotyped == 1)
  snp2 <- as.integer(goudstandaard$NBCS_SNP_genotyped == 2)
  goudstandaard <- cbind(goudstandaard, snp1, snp2)
  expect_equal(goudstandaard$snp1, regwas$reg_data$snp1)
  expect_equal(goudstandaard$snp2, regwas$reg_data$snp2)

  regwas$setFlipMaf(TRUE)
  regwas$SetCovariates("rs16984366")
  # MAF is to big to flip, so coding stays the same
  # goudstandaard <- read.table("../bladder2_SNPs_covar.txt", header = T)
  # snp1 <- as.integer(goudstandaard$NBCS_SNP_genotyped == 0) * 2 + as.integer(goudstandaard$NBCS_SNP_genotyped == 1)
  # snp2 <- as.integer(goudstandaard$NBCS_SNP_genotyped == 2) * 2

  goudstandaard <- cbind(goudstandaard, snp1, snp2)

  expect_equal(goudstandaard$snp1, regwas$reg_data$snp1)
  expect_equal(goudstandaard$snp2, regwas$reg_data$snp2)
})

test_that("overview normal_2df", {
  regwas$set_2df_model()
  # regwas$setFlipMaf(TRUE)
  re_results <- regwasworker(regwas, 1)
  expect_equal(class(re_results), "data.frame")
  goudstandaard <- read.table("../bladder2_SNPs_covar.txt", header = T)

  snp1 <- as.integer(goudstandaard$NBCS_SNP_genotyped == 2)

  snp2 <- as.integer(goudstandaard$NBCS_SNP_genotyped == 1)

  # snp3 <- as.integer(goudstandaard$NBCS_SNP_genotyped == 0)
  goudstandaard <- cbind(goudstandaard, snp1, snp2)
  z <- summary(coxph(formula = Surv(tstart, tstop, status) ~ snp1 + snp2 +
    +age + sex, data = goudstandaard, cluster = id))

  beta <- z$coefficients
  expect_equivalent(beta["snp2", "coef"], re_results["rs16984366", "beta_A1A1"], tolerance = .00001)
  # expect_equivalent(beta["snp1","coef"], 0.442973116, tolerance = .000001)


  expect_equivalent(beta["snp1", "coef"], re_results["rs16984366", "beta_A1A2"], tolerance = .00001)
  # expect_equivalent(beta["snp2","coef"], 0.442973116, tolerance = .000001)


  m <- coxph(formula = Surv(tstart, tstop, status) ~ snp1 + snp2 + age + sex, data = goudstandaard, cluster = id)
  n <- coxph(formula = Surv(tstart, tstop, status) ~ age + sex, data = goudstandaard, cluster = id)
  #
  expect_equal(-457.50245, re_results["rs16984366", "logliknull"], tolerance = .00001)
  expect_equal(n$loglik[2], re_results["rs16984366", "logliknull"], tolerance = .00001)

  expect_equal(-453.3615, re_results["rs16984366", "loglik"], tolerance = .00001)
  expect_equal(m$loglik[2], re_results["rs16984366", "loglik"], tolerance = .00001)



  expect_equal(anova(m, n)$`P(>|Chi|)`[2], re_results["rs16984366", "LRT"], tolerance = .00001)
  expect_equal(anova(m, n)$`P(>|Chi|)`[2], 0.01590744, tolerance = .00001)
})

test_that("2df coxme", {
  regwas$set_coxme()
  regwas$set_2df_model()
  re_results <- regwasworker(regwas, 1)
  expect_equal(class(re_results), "data.frame")

  goudstandaard <- read.table("../bladder2_SNPs_covar.txt", header = T)

  snp1 <- as.integer(goudstandaard$NBCS_SNP_genotyped == 1)
  snp2 <- as.integer(goudstandaard$NBCS_SNP_genotyped == 2)
  goudstandaard <- cbind(goudstandaard, snp1, snp2)

  z <- (coxme::coxme(Surv(tstart, tstop, status) ~ snp1 + snp2 + age + sex + (1 | id), data = goudstandaard))
  beta <- z$coefficients
  nvar <- length(beta)
  nfrail <- nrow(z$var) - nvar
  se <- sqrt(diag(as.matrix(z$var))[nfrail + 1:nvar])
  names(se) <- names(beta)
  p <- 1 - pchisq((beta / se)^2, 1)

  # p-value
  # expect_equivalent(p["NBCS_SNP_genotyped"], re_results["rs16984366", "pvalue"])
  # expect_equivalent(p["NBCS_SNP_genotyped"], 0.064, tolerance = 0.0005)

  # c("coef", "exp(coef)","se(coef)", "z", "p"
  # beta, exp(beta), se, round(beta/se,2),1 - pchisq((beta/ se)^2, 1)


  # se(coef)
  # expect_equivalent(se["NBCS_SNP_genotyped"], re_results["rs16984366", "se"], tolerance = .00001)
  # expect_equivalent(se["NBCS_SNP_genotyped"], 0.23945113, tolerance = .00001)

  # expect_equivalent(beta["NBCS_SNP_genotyped"], re_results["rs16984366", "beta"], tolerance = .00001)
  # expect_equivalent(beta["NBCS_SNP_genotyped"], 0.442973116, tolerance = .000001)


  m <- coxme::coxme(Surv(tstart, tstop, status) ~ snp1 + snp2 + age + sex + (1 | id), data = goudstandaard)
  n <- coxme::coxme(Surv(tstart, tstop, status) ~ age + sex + (1 | id), data = goudstandaard)

  # expect_equal(anova(m, n)$`P(>|Chi|)`[2], re_results["rs16984366", "LRT"], tolerance = .00001)
  # expect_equivalent(anova(m, n)$`P(>|Chi|)`[2], 0.04641, tolerance = .00001)
})


test_that("2df frailty gaussian", {
  regwas$set_frailty.gaussian()
  regwas$set_2df_model()
  expect_equal(class(regwasworker(regwas, 1)), "data.frame")
})

test_that("2df frailty gamma", {
  regwas$set_frailty.gamma()
  regwas$set_2df_model()
  expect_equal(class(regwasworker(regwas, 1)), "data.frame")
})
