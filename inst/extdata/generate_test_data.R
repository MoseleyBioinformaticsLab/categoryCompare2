library(jsonlite)
library(categoryCompare2)
library(org.Hs.eg.db)
library(GO.db)
library(magrittr)
library(dplyr)


library("affy")
library("hgu95av2.db")
library("estrogen")
library("limma")
library("genefilter")

data_dir <- system.file("extdata", package = "estrogen")
pheno_data <- read.AnnotatedDataFrame(file.path(data_dir, "estrogen.txt"),
                                      header = TRUE, sep = "", row.names = 1)
estrogen_raw <- ReadAffy(filenames = file.path(data_dir, rownames(pData(pheno_data))),
                         phenoData = pheno_data, verbose = TRUE)

estrogen_rma <- rma(estrogen_raw)

affy_ids = featureNames(estrogen_rma)
affy_other = AnnotationDbi::select(hgu95av2.db, keys = affy_ids, columns = c("ENTREZID", "SYMBOL", "UNIPROT")) |>
  dplyr::mutate(affy = PROBEID)

e10 <- estrogen_rma[, estrogen_rma$time.h == 10]
e10 <- nsFilter(e10, remove.dupEntrez=TRUE, var.filter=FALSE, 
                feature.exclude="^AFFX")$eset

e10$estrogen <- factor(e10$estrogen)
d10 <- model.matrix(~0 + e10$estrogen)
colnames(d10) <- unique(e10$estrogen)
fit10 <- lmFit(e10, d10)
c10 <- makeContrasts(present - absent, levels=d10)
fit10_2 <- contrasts.fit(fit10, c10)
eB10 <- eBayes(fit10_2)
table10 <- topTable(eB10, number=nrow(e10), p.value=1, adjust.method="BH")
table10$affy <- rownames(table10)
table10 = dplyr::left_join(table10, affy_other, by = "affy")

e48 <- estrogen_rma[, estrogen_rma$time.h == 48]
e48 <- nsFilter(e48, remove.dupEntrez=TRUE, var.filter=FALSE, 
                feature.exclude="^AFFX")$eset

e48$estrogen <- factor(e48$estrogen)
d48 <- model.matrix(~0 + e48$estrogen)
colnames(d48) <- unique(e48$estrogen)
fit48 <- lmFit(e48, d48)
c48 <- makeContrasts(present - absent, levels=d48)
fit48_2 <- contrasts.fit(fit48, c48)
eB48 <- eBayes(fit48_2)
table48 <- topTable(eB48, number=nrow(e48), p.value=1, adjust.method="BH")
table48$affy <- rownames(table48)

table48 = dplyr::left_join(table48, affy_other, by = "affy")
adj_cut <- 0.05

t48_entrez <- dplyr::filter(table48, adj.P.Val <= adj_cut) %>% dplyr::pull("ENTREZID") %>% unique()
t10_entrez <- dplyr::filter(table10, adj.P.Val <= adj_cut) %>% dplyr::pull("ENTREZID") %>% unique()
universe_entrez <- unique(table10$ENTREZID)

cat(t48_entrez, sep = "\n", file = "executable_related/test_data/48_entrez.txt")
cat(t10_entrez, sep = "\n", file = "executable_related/test_data/10_entrez.txt")
cat(universe_entrez, sep = "\n", file = "executable_related/test_data/universe_entrez.txt")

t48_symbol <- dplyr::filter(table48, adj.P.Val <= adj_cut) %>% extract2("SYMBOL") %>% unique()
t10_symbol <- dplyr::filter(table10, adj.P.Val <= adj_cut) %>% extract2("SYMBOL") %>% unique()
universe_symbol <- unique(table10$SYMBOL)

cat(t48_symbol, sep = "\n", file = "executable_related/test_data/48_symbol.txt")
cat(t10_symbol, sep = "\n", file = "executable_related/test_data/10_symbol.txt")
cat(universe_symbol, sep = "\n", file = "executable_related/test_data/universe_symbol.txt")

t48_uniprot <- dplyr::filter(table48, adj.P.Val <= adj_cut) %>% extract2("UNIPROT") %>% unique()
t10_uniprot <- dplyr::filter(table10, adj.P.Val <= adj_cut) %>% extract2("UNIPROT") %>% unique()
universe_uniprot <- unique(table10$UNIPROT)

cat(t48_uniprot, sep = "\n", file = "executable_related/test_data/48_uniprot.txt")
cat(t10_uniprot, sep = "\n", file = "executable_related/test_data/10_uniprot.txt")
cat(universe_uniprot, sep = "\n", file = "executable_related/test_data/universe_uniprot.txt")

