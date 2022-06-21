



test_that("Gene set handling", {

  genesets <- getGeneSets(c("C2|CP:KEGG"))

  sum(grepl("KEGG", names(genesets))) > 100

  df <- convertGeneSets(genesets, from = "list", to = "dataframe")
  genesets2 <- convertGeneSets(df, from = "dataframe", to = "list")

  expect_equal(genesets, genesets2)
})


test_that("ID mapping", {

  g1 <- convertGeneIDs(c("3098", "3099"), from = "ENTREZID", to = "SYMBOL")
  g2 <- convertGeneIDs(c("MTOR", "MYC"), from = "SYMBOL", to = "ENTREZID")

  expect_true(all(g1 == c("HK1", "HK2")))
  expect_true(all(g2 == c("2475", "4609")))

})


test_that("DE analysis", {

  tmp <- getTestData(n = 12, nrow = 100)

  # DESeq2
  ds2 <- runDESeq2(data = tmp$data, design = tmp$design, formula = ~ group, contrasts = tmp$contrasts, vst = FALSE)
  up <- rownames(subset(ds2$results[[1]], log2FC > 3 & padj <= 0.01))
  down <- rownames(subset(ds2$results[[1]], log2FC < -3 & padj <= 0.01))
  expect_type(ds2, "list")
  expect_equal(class(ds2$results[[1]]), "data.frame")
  expect_true(mean(up %in% tmp$up) > 0.9)
  expect_true(mean(down %in% tmp$down) > 0.9)


  # Limma
  lm <- runLIMMA(data = log2(ds2$normcounts + 1), design = tmp$design, formula = ~ group, contrasts = tmp$contrasts)
  up.lm <- rownames(subset(lm[[1]], log2FC > 3 & padj <= 0.01))
  down.lm <- rownames(subset(lm[[1]], log2FC < -3 & padj <= 0.01))
  expect_true(mean(up.lm %in% tmp$up) > 0.9)
  expect_true(mean(down.lm %in% tmp$down) > 0.9)


  # GSVA
  GS <- getTestGenesets(rownames(tmp$data))
  GS$up <- tmp$up
  GS$down <- tmp$down
  gsva <- runGSVA(log2(ds2$normcounts+1), genesets = GS)
  expect_equal(nrow(gsva), length(GS))
  expect_true(mean(as.numeric(gsva["up",tmp$design$group == 2])) > 0)
  expect_true(mean(as.numeric(gsva["down",tmp$design$group == 2])) < 0)


  # fGSEA
  fgsea <- runfGSEA(ds2$results[[1]], genesets = GS)
  expect_true(subset(fgsea, padj <= 0.01 & term == "up")[,"NES"] > 0.9)
  expect_true(subset(fgsea, padj <= 0.01 & term == "down")[,"NES"] < -0.9)

  # GSEA
  gsea <- runGSEA(ds2$results[[1]], genesets = GS, as.df = TRUE)
  expect_true(subset(fgsea, padj <= 0.01 & term == "up")[,"NES"] > 0.9)
  expect_true(subset(fgsea, padj <= 0.01 & term == "down")[,"NES"] < -0.9)

})


test_that("getRanks", {

  df <- data.frame(stat = runif(100))
  expect_equal(unname(getRanks(df, type = "+")), sort(df$stat, decreasing = TRUE))

  df <- data.frame(stat = runif(100))
  expect_equal(unname(getRanks(df, type = "p")), -log10(sort(df$stat, decreasing = FALSE)))

  df <- data.frame(stat = runif(100)-0.5)
  expect_equal(unname(getRanks(df, type = "")), sort(df$stat, decreasing = TRUE))

})



