# orgdb annotations

Generate an annotation object for genes based on an "org.\*.db" object,
and pulling information from it.

## Usage

``` r
get_db_annotation(
  orgdb = "org.Hs.eg.db",
  features = NULL,
  feature_type = "ENTREZID",
  annotation_type = "GO"
)
```

## Arguments

- orgdb:

  the name of the org.\*.db object

- features:

  which features to get annotations for

- feature_type:

  which type of IDs to map (see details)

- annotation_type:

  the type of annotation to grab (see details)

## Value

annotation object

## Details

This function generates a `categoryCompare2` annotation object from a
Bioconductor "org.\*.db" object. Even though different gene identifiers
can be used, almost all of the mappings are via ENTREZID.

The set of feature or gene keys that can be used to create the
annotations include:

- ENTREZID: ENTREZ gene ids

- ACCNUM: genbank accession numbers

- SYMBOL: gene symbols, eg ABCA1

- GENENAME: gene names, eg "ATP binding cassette subfamily A member 1"

- ENSEMBL: the ensembl gene ids (all start with ENSG...)

- ENSEMBLPROT: ensembl protein ids (ENSP...)

- ENSEMBLTRANS: ensemlb transcript ids (ENST...)

- REFSEQ: reference sequence IDs, NM, NP, NR, XP, etc

- UNIGENE: gene ids from UNIPROT eg Hs.88556

- UNIPROT: protein ids from UNIPROT eg P80404

The set of annotations that can be mapped to features include:

- GO: annotations from gene ontology

- PATH: KEGG Pathway identifiers (not updated since 2011!)

- CHRLOC: location on the chromosome

- OMIM: mendelian inheritance in man identifiers

- PMID: pubmed identifiers

- PROSITE

- PFAM: protein family identifiers

- IPI: protein-protein interactions

For GO annotations, it is also possible to pass `GO` to use all 3
sub-ontologies simultaneously, or any combination of `BP`, `MF`, and
`CC`.
