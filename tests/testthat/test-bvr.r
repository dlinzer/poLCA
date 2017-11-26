library(poLCA)
context("BVRs for latent class models")


test_that("BVRs for example model are as expected (no covariates)", {
  set.seed(201708)
  
  data(values)
  f <- cbind(A,B,C,D)~1
  M0 <- poLCA(f,values,nclass=1, verbose = FALSE) # log-likelihood: -543.6498
  M1 <- poLCA(f,values,nclass=2, verbose = FALSE) # log-likelihood: -504.4677
  
  bvr_0 <- bvr(M0)
  bvr_1 <- bvr(M1)
  
  expect_that(all(abs(bvr_0 - c(12.378947,  5.704402,  8.308752, 15.586100, 23.562456, 23.779559)) < 1e-6), is_true() )
  expect_that(all(abs(bvr_1 - c(1.172161691, 0.003303647, 0.010097938, 0.037913535, 0.016503607, 0.055858509)) < 1e-6), is_true() )
  
})

test_that("BVRs for example model are as expected (with covariates)", {
  set.seed(201708)
  
  data(election)
  
  f2a <- cbind(MORALG, CARESG, DISHONG, INTELG) ~ PARTY
  
  fit_nes <-poLCA(f2a, election, nclass = 3, nrep = 5, verbose = T)  
  
  bvr_nes <- bvr(fit_nes)
  
  expect_that(all(abs(bvr_nes - c(5.773705,  6.799946,  7.369901,  4.684608,  5.205810, 11.022726)) < 1e-6), is_true() )
})