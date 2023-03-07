

regwas <- Regwas$new("../snps_goudenstandaard_n85.bgen", "../pheno_goudenstandaard_n85_agesex.txt")

test_that("HWE", {
  expect_equal(regwas$HWEchisq(c(1469, 138, 5)), 1 - 0.637727, 0.0003)
  expect_equal(regwas$HWEchisq(c(40, 34, 11)), 0.3830583, tolerance = .0000001)
})


test_that("MAF", {
  expect_equal(regwas$MAF(c(0, 0, 100)), 0)
  expect_equal(regwas$MAF(c(100, 0, 0)), 0)
  expect_equal(regwas$MAF(c(100, 0, 1)), 0.0099, tolerance = .000001)
  expect_equal(regwas$MAF(c(1, 0, 100)), 0.0099, tolerance = .000001)
  expect_equal(regwas$MAF(c(100, 0, 100)), 0.5)
  expect_equal(regwas$MAF(c(0, 100, 0)), 0.5)
  expect_equal(regwas$MAF(c(12.6, 40.7, 31.6)), 0.3881037, tolerance = .0000001)
})
