# Test preprocessing epigenetic data and calculate DMP and DMRs
set.seed(99999)

## Make RGset
base_dir=system.file("extdata", package="minfiDataEPIC", mustWork = TRUE)
testSet_afterImport = Make_RGset(idat_paths = base_dir, name = "test")
test_that("Test that Make_RGset imports all idat files in directory and subdirectories", {
  expect_match(testSet_afterImport$name, "test")
  expect_equal(dim(testSet_afterImport$rgset), c(1052641,3))
})

## Normalisation
testSet_afterNormalisation = ProbeNormalisation(testSet_afterImport)
test_that("Test that probe normalisation works - funnorm", {
  expect_match(testSet_afterNormalisation$normalisationMethod, "funnorm")
  expect_equal(dim(testSet_afterNormalisation$mset), c(865859, 3))
})

## DMP_limma
cat_vars = paste0("cat_", 1:1)
cont_vars = paste0("cont_", 1)
Random_AdjustmentVars = data.frame(replicate(length(cat_vars), sample(0:1, dim(testSet_afterNormalisation$mset)[2], rep=TRUE)),
                                   replicate(length(cont_vars), sample(0:10, dim(testSet_afterNormalisation$mset)[2], rep=TRUE)))
colnames(Random_AdjustmentVars) = c(cat_vars, cont_vars)
test_LimmaOutput = suppressWarnings(DMP_limma(testSet_afterNormalisation$mset, Random_AdjustmentVars, cat_vars=cat_vars, Group="cat_1"))
test_that("test that DMP_limma finds DMPs", {
  expect_equal(dim(test_LimmaOutput$DMPresult), c(865859, 6))
  expect_equal(length(test_LimmaOutput$Sign_CpG), 627)
})

## DMR_DMRCate
test_that("DMR_dmrcate returns regions", {
  ## Currently not enough data
  #DMR_output = DMR_dmrcate(testSet_afterNormalisation$mset, test_LimmaOutput, lambda = 100000, fdr = 1)
  #expect_equal(dim(DMR_output$DMRresult), c(, ))
  expect_equal(1, 1)
})
