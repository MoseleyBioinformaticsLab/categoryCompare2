# generating much smaller test data

library(GO.db)
library(GOstats)
library(org.Hs.eg.db)
library(magrittr)

set.seed("1234")

# to choose genes, we will actually choose a set of GO terms first
all_hs_GO <- as.list(org.Hs.egGO2ALLEGS)
all_hs_GO <- lapply(all_hs_GO, unique)

len_GO <- sapply(all_hs_GO, length)

keep_GO <- (len_GO >= 2) & (len_GO <= 200)

all_hs_GO <- all_hs_GO[keep_GO]

init_GO <- sample(names(all_hs_GO), 200) %>% all_hs_GO[.]
init_GO_len <- sapply(init_GO, length)

diff_GO <- which((init_GO_len >= 30) & (init_GO_len <= 50))[1]
keep_genes <- init_GO[[diff_GO]]
n_keep <- length(keep_genes)

all_genes <- unique(unlist(init_GO, use.names = FALSE))

sample_genes <- all_genes[!(all_genes %in% keep_genes)]

universe_genes <- c(keep_genes, sample(sample_genes, 200 - n_keep))

c1_diff <- sample(keep_genes, 20)
c2_diff <- sample(keep_genes, 20)

c1_hyper <- new("GOHyperGParams", geneIds = c1_diff,
                 universeGeneIds = universe_genes,
                 annotation = "org.Hs.eg.db",
                 ontology = "BP",
                 pvalueCutoff = 1.0)
c2_hyper <- new("GOHyperGParams", geneIds = c2_diff,
                 universeGeneIds = universe_genes,
                 annotation = "org.Hs.eg.db",
                 ontology = "BP",
                 pvalueCutoff = 1.0)

c1_enrich <- hyperGTest(c1_hyper)
c2_enrich <- hyperGTest(c2_hyper)


all_hs_go <- keys(org.Hs.eg.db, keytype = "GOALL")
all_hs_gotype <- select(GO.db, keys = all_hs_go, columns = "ONTOLOGY")
all_hs_bp <- all_hs_gotype$GOID[all_hs_gotype$ONTOLOGY == "BP"]

all_hs_bp2gene <- select(org.Hs.eg.db, keys = all_hs_bp, keytype = "GOALL", columns = c("ENTREZID"))

bp_annotation <- split(all_hs_bp2gene$ENTREZID, all_hs_bp2gene$GOALL)
bp_annotation <- lapply(bp_annotation, unique)

c1_enrichment <- list(features = c1_diff,
                               universe = universe_genes,
                               annotation = bp_annotation,
                               enrich = list(pvalues = pvalues(c1_enrich),
                                             odds_ratio = oddsRatios(c1_enrich),
                                             counts = geneCounts(c1_enrich),
                                             stat_names = names(pvalues(c1_enrich))))

c2_enrichment <- list(features = c2_diff,
                               universe = universe_genes,
                               annotation = bp_annotation,
                               enrich = list(pvalues = pvalues(c2_enrich),
                                             odds_ratio = oddsRatios(c2_enrich),
                                             counts = geneCounts(c2_enrich),
                                             stat_names = names(pvalues(c2_enrich))))

save(c1_enrichment, c2_enrichment, bp_annotation, file = "test_data2.RData")
