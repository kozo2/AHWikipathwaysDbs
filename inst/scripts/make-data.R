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

# metabolite_bridge_url <- "https://ndownloader.figshare.com/files/25388453"
# download.file(metabolite_bridge_url, "metabolites_20201104.bridge")

file <- "metabolites_20201104.bridge"
download.file(
  "https://ndownloader.figshare.com/files/25388453",
  location
)
location = normalizePath(file)
mapper <- loadDatabase(location)

# see https://bridgedb.github.io/pages/system-codes.html
# map(mapper, "HMDB0001058", source="Ch", target="Ck")
# map(mapper, "HMDB0001058", source="Ch", target="Ce")
# map(mapper, "HMDB0001058", source="Ch", target="Ca")
# map(mapper, "HMDB0001058", source="Ch", target="Dr")
# map(mapper, "HMDB0001058", source="Ch", target="Ik")
