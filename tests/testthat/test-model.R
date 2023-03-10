model <- Model$new()
test_that("default", {
  expect_equal(model$genetic_model, "additive")
  expect_equal(model$current_method, "default")
  expect_equal(model$formula, "surv_object~variant_coded")
  model$set_2df()
  expect_equal(model$formula, "surv_object~snp1+snp2")
  model$set_covariants(c("sex", "age"))
  expect_equal(model$formula, "surv_object~snp1+snp2+sex+age")
  model$set_covariants(NA)
})

test_that("coxme", {
  model$set_coxme()
  model$set_additive()
  expect_equal(model$current_method, "coxme")
  expect_equal(model$formula, "surv_object~variant_coded+(1|samples)")
  model$set_2df()
  expect_equal(model$genetic_model, "2df")
  expect_equal(model$formula, "surv_object~snp1+snp2+(1|samples)")
  model$set_additive()
  expect_equal(model$genetic_model, "additive")
  expect_equal(model$formula, "surv_object~variant_coded+(1|samples)")
})
