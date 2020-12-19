library(xml2)
library(rWikiPathways)
library(BridgeDbR)

wikipathways_release_date <- "20201210"

wikipathways_species <- c("Anopheles gambiae", "Arabidopsis thaliana", "Bacillus subtilis",
    "Bos taurus", "Caenorhabditis elegans", "Canis familiaris", "Danio rerio",
    "Drosophila melanogaster", "Equus caballus", "Escherichia coli",
    "Gallus gallus", "Gibberella zeae", "Homo sapiens", "Hordeum vulgare",
    "Mus musculus", "Mycobacterium tuberculosis", "Oryza sativa",
    "Pan troglodytes", "Plasmodium falciparum", "Populus trichocarpa",
    "Rattus norvegicus", "Saccharomyces cerevisiae", "Solanum lycopersicum",
    "Sus scrofa", "Zea mays")

for (species in wikipathways_species) {
    downloadPathwayArchive(date=wikipathways_release_date, organism=species, format="gpml")
}
