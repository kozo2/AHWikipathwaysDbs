wikipathways_species <- c("Anopheles gambiae", "Arabidopsis thaliana",
    "Bacillus subtilis", "Bos taurus", "Caenorhabditis elegans",
    "Canis familiaris", "Danio rerio", "Drosophila melanogaster",
    "Equus caballus", "Escherichia coli", "Gallus gallus", "Gibberella zeae",
    "Homo sapiens", "Hordeum vulgare", "Mus musculus",
    "Mycobacterium tuberculosis", "Oryza sativa", "Pan troglodytes",
    "Plasmodium falciparum", "Populus trichocarpa", "Rattus norvegicus",
    "Saccharomyces cerevisiae", "Solanum lycopersicum", "Sus scrofa",
    "Zea mays")

wikipathways_taxonomyids <-c("7165", "3702",
    "1423", "9913", "6239",
    "9615", "7955", "7227",
    "9796", "562", "9031", "5518",
    "9606", "4513", "10090",
    "1773", "4530", "9598",
    "5833", "3694", "10116",
    "4932", "4081", "9823",
    "4577")

titles <- paste("wikipathways", gsub(" ", "_", wikipathways_species),
                "metabolites.rda", sep = "_")
descriptions <- paste('Metabolite names linked to Wikipathways',
    wikipathways_species,
    'pathways (includes HMDB, KEGG, ChEBI, Drugbank, PubChem compound,',
    'ChemSpider, KNApSAcK, and Wikidata IDs + CAS + InChI Key)')

meta <- data.frame(
    Title = titles,
    Description = descriptions,
    BiocVersion = "3.18",
    Genome = NA,
    SourceType = "XML",
    SourceUrl = "https://www.wikipathways.org/index.php/Download_Pathways",
    SourceVersion = "Jun 10 2023",
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
