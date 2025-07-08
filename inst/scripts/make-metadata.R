library(httr2)
library(jsonlite)
library(taxize)

url <- "https://www.wikipathways.org/json/listOrganisms.json"
resp <- request(url) |> req_perform()
json_data <- resp_body_string(resp)
parsed_data <- fromJSON(json_data)
wikipathways_species <- parsed_data$organisms

# Get NCBI taxonomy IDs for the species using taxize
wikipathways_taxonomyids <- as.character(get_uid(wikipathways_species))

titles <- paste("wikipathways", gsub(" ", "_", wikipathways_species),
                "metabolites.rda", sep = "_")
descriptions <- paste('Metabolite names linked to Wikipathways',
    wikipathways_species,
    'pathways (includes HMDB, KEGG, ChEBI, Drugbank, PubChem compound,',
    'ChemSpider, KNApSAcK, and Wikidata IDs + CAS + InChI Key)')

meta <- data.frame(
    Title = titles,
    Description = descriptions,
    BiocVersion = "3.21",
    Genome = NA,
    SourceType = "XML",
    SourceUrl = "https://data.wikipathways.org/20250610/gpml/",
    SourceVersion = "Jun 10 2025",
    Species = wikipathways_species,
    TaxonomyId = wikipathways_taxonomyids,
    Coordinate_1_based = NA,
    DataProvider = "WikiPathways",
    Maintainer = "Kozo Nishida <kozo.nishida@gmail.com>",
    RDataClass = "Tibble",
    DispatchClass = "Rda",
    RDataPath = paste0('AHWikipathwaysDbs/', titles),
    Tags = "Pathways:metabolites:HMDB:KEGG:ChEBI:Drugbank:CAS:ChemSpider"
)

write.csv(meta, file="metadata.csv", row.names=FALSE)
