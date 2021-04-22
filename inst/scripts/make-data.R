library(xml2)
library(stringr)
library(tibble)
library(BridgeDbR)
library(rWikiPathways)
library(dplyr)

wikipathways_release_date <- "20210410"

wikipathways_species <- c("Anopheles gambiae", "Arabidopsis thaliana",
    "Bacillus subtilis", "Bos taurus", "Caenorhabditis elegans",
    "Canis familiaris", "Danio rerio", "Drosophila melanogaster",
    "Equus caballus", "Escherichia coli", "Gallus gallus", "Gibberella zeae",
    "Homo sapiens", "Hordeum vulgare", "Mus musculus",
    "Mycobacterium tuberculosis", "Oryza sativa", "Pan troglodytes",
    "Plasmodium falciparum", "Populus trichocarpa", "Rattus norvegicus",
    "Saccharomyces cerevisiae", "Solanum lycopersicum", "Sus scrofa",
    "Zea mays")

# metabolite_bridge_url <- "https://ndownloader.figshare.com/files/25388453"
# download.file(metabolite_bridge_url, "metabolites_20201104.bridge")

file <- "metabolites_20210109.bridge"
# download.file(
#     "https://ndownloader.figshare.com/files/26001794",
#     location
# )
# location = normalizePath(file)
# mapper <- loadDatabase(location)
mapper <- loadDatabase(file)

for (species in wikipathways_species) {
    filename = downloadPathwayArchive(date=wikipathways_release_date,
                                    organism=species, format="gpml")
    gpmlpaths = utils::unzip(filename)
    df <- tibble(wpid = character(), wpid_version = character(),
                pathway_name = character(),
                metabolite_name = character(),
                HMDB_ID = character(), KEGG_ID = character(),
                ChEBI_ID = character(), DrugBank_ID = character(),
                PubChem_CID = character(), ChemSpider = character(),
                KNaPSAcK_ID = character(), Wikidata_ID = character(),
                CAS = character(), InChI_Key =character())
    for (gpmlpath in gpmlpaths) {
        print(paste0('reading', gpmlpath))
        df = readGpml(gpmlpath, df)
    }
    filename = paste("wikipathways", gsub(" ", "_", species),
                     "metabolites.rda", sep = "_")
    df = distinct(df)
    save(df, file = filename, compress = "xz")
}

readGpml <- function(gpmlpath, df) {
    x <- read_xml(gpmlpath)
    filename <- tools::file_path_sans_ext(gpmlpath)
    wpidinfo <- tail(str_split(filename, "_", simplify = TRUE)[1,], n=2)
    wpid <- wpidinfo[1]
    wpid_version <- wpidinfo[2]
    pathway_name <- xml_attr(x, "Name")
    pathway_description = ""
    pathway_category = ""
    for (xml_node in xml_children(x)) {
        if (xml_name(xml_node) == "DataNode") {
            if (!is.na(xml_attr(xml_node, "Type") == "Metabolite")) {
                metabolite_name = xml_attr(xml_node, "TextLabel")
                childs = xml_children(xml_node)
                xref_node = childs[length(childs)]
                database = xml_attr(xref_node, "Database")
                mol_id = xml_attr(xref_node, "ID")
                if (database == "HMDB") {
                    df = addrow_from_hmdb(df, wpid, wpid_version,
                                      pathway_name, metabolite_name, mol_id)
                } else if (database == "ChEBI") {
                    df = addrow_from_chebi(df, wpid, wpid_version,
                                      pathway_name, metabolite_name, mol_id)
                } else if (database == "KEGG Compound") {
                    df = addrow_from_kegg(df, wpid, wpid_version,
                                      pathway_name, metabolite_name, mol_id)
                } else if (database == "PubChem-compound") {
                    df = addrow_from_pubchem(df, wpid, wpid_version,
                                      pathway_name, metabolite_name, mol_id)
                } else if (database == "CAS") {
                    df = addrow_from_cas(df, wpid, wpid_version,
                                      pathway_name, metabolite_name, mol_id)
                } else if (database == "Chemspider") {
                  df = addrow_from_chemspider(df, wpid, wpid_version,
                                      pathway_name, metabolite_name, mol_id)
                } else if (database == "Wikidata") {
                  df = addrow_from_wikidata(df, wpid, wpid_version,
                                      pathway_name, metabolite_name, mol_id)
                } else if (database == "KNApSAcK") {
                  df = addrow_from_knapsack(df, wpid, wpid_version,
                                      pathway_name, metabolite_name, mol_id)
                } else if (is.null(database) || database == "") {
                    #print(paste0(metabolite_name, " does not have Xref"))
                } else {
                    #print(paste0(database, " is needed to addrow_from"))
                }
            }
        }
    }
    return(df)
}

addrow_from_hmdb <- function(df, wpid, wpid_version, pathway_name,
                        metabolite_name, mol_id){
    system_code = "Ch"
    Ck = get_mapped_id(mapper, mol_id, source = system_code, target = "Ck")
    Ce = get_mapped_id(mapper, mol_id, source = system_code, target = "Ce")
    Cpc = get_mapped_id(mapper, mol_id, source = system_code, target = "Cpc")
    Dr = get_mapped_id(mapper, mol_id, source = system_code, target = "Dr")
    Ca = get_mapped_id(mapper, mol_id, source = system_code, target = "Ca")
    Ik = get_mapped_id(mapper, mol_id, source = system_code, target = "Ik")
    Cs = get_mapped_id(mapper, mol_id, source = system_code, target = "Cs")
    Wd = get_mapped_id(mapper, mol_id, source = system_code, target = "Wd")
    Cks = get_mapped_id(mapper, mol_id, source = system_code, target = "Cks")
    df = add_row(df, wpid = wpid, wpid_version = wpid_version,
            pathway_name = pathway_name,
            metabolite_name = metabolite_name,
            HMDB_ID = mol_id, KEGG_ID = Ck, ChEBI_ID = Ce,
            PubChem_CID = Cpc, DrugBank_ID = Dr, CAS = Ca, InChI_Key = Ik,
            ChemSpider = Cs, Wikidata_ID = Wd, KNaPSAcK_ID = Cks)
    return(df)
}

addrow_from_kegg <- function(df, wpid, wpid_version, pathway_name,
                        metabolite_name, mol_id){
  system_code = "Ck"
  Ch = get_mapped_id(mapper, mol_id, source = system_code, target = "Ch")
  Ce = get_mapped_id(mapper, mol_id, source = system_code, target = "Ce")
  Cpc = get_mapped_id(mapper, mol_id, source = system_code, target = "Cpc")
  Dr = get_mapped_id(mapper, mol_id, source = system_code, target = "Dr")
  Ca = get_mapped_id(mapper, mol_id, source = system_code, target = "Ca")
  Ik = get_mapped_id(mapper, mol_id, source = system_code, target = "Ik")
  Cs = get_mapped_id(mapper, mol_id, source = system_code, target = "Cs")
  Wd = get_mapped_id(mapper, mol_id, source = system_code, target = "Wd")
  Cks = get_mapped_id(mapper, mol_id, source = system_code, target = "Cks")
  df = add_row(df, wpid = wpid, wpid_version = wpid_version,
            pathway_name = pathway_name,
            metabolite_name = metabolite_name,
            HMDB_ID = Ch, KEGG_ID = mol_id, ChEBI_ID = Ce,
            PubChem_CID = Cpc, DrugBank_ID = Dr, CAS = Ca, InChI_Key = Ik,
            ChemSpider = Cs, Wikidata_ID = Wd, KNaPSAcK_ID = Cks)
  return(df)
}

addrow_from_chebi <- function(df, wpid, wpid_version, pathway_name,
                        metabolite_name, mol_id){
    system_code = "Ce"
    Ch = get_mapped_id(mapper, mol_id, source = system_code, target = "Ch")
    Ck = get_mapped_id(mapper, mol_id, source = system_code, target = "Ck")
    Cpc = get_mapped_id(mapper, mol_id, source = system_code, target = "Cpc")
    Dr = get_mapped_id(mapper, mol_id, source = system_code, target = "Dr")
    Ca = get_mapped_id(mapper, mol_id, source = system_code, target = "Ca")
    Ik = get_mapped_id(mapper, mol_id, source = system_code, target = "Ik")
    Cs = get_mapped_id(mapper, mol_id, source = system_code, target = "Cs")
    Wd = get_mapped_id(mapper, mol_id, source = system_code, target = "Wd")
    Cks = get_mapped_id(mapper, mol_id, source = system_code, target = "Cks")
    df = add_row(df, wpid = wpid, wpid_version = wpid_version,
            pathway_name = pathway_name,
            metabolite_name = metabolite_name,
            HMDB_ID = Ch, KEGG_ID = Ck, ChEBI_ID = mol_id,
            PubChem_CID = Cpc, DrugBank_ID = Dr, CAS = Ca, InChI_Key = Ik,
            ChemSpider = Cs, Wikidata_ID = Wd, KNaPSAcK_ID = Cks)
    return(df)
}

addrow_from_pubchem <- function(df, wpid, wpid_version, pathway_name,
                        metabolite_name, mol_id){
    system_code = "Cpc"
    Ch = get_mapped_id(mapper, mol_id, source = system_code, target = "Ch")
    Ck = get_mapped_id(mapper, mol_id, source = system_code, target = "Ck")
    Ce = get_mapped_id(mapper, mol_id, source = system_code, target = "Ce")
    Dr = get_mapped_id(mapper, mol_id, source = system_code, target = "Dr")
    Ca = get_mapped_id(mapper, mol_id, source = system_code, target = "Ca")
    Ik = get_mapped_id(mapper, mol_id, source = system_code, target = "Ik")
    Cs = get_mapped_id(mapper, mol_id, source = system_code, target = "Cs")
    Wd = get_mapped_id(mapper, mol_id, source = system_code, target = "Wd")
    Cks = get_mapped_id(mapper, mol_id, source = system_code, target = "Cks")
    df = add_row(df, wpid = wpid, wpid_version = wpid_version,
            pathway_name = pathway_name,
            metabolite_name = metabolite_name,
            HMDB_ID = Ch, KEGG_ID = Ck, ChEBI_ID = Ce,
            PubChem_CID = mol_id, DrugBank_ID = Dr, CAS = Ca, InChI_Key = Ik,
            ChemSpider = Cs, Wikidata_ID = Wd, KNaPSAcK_ID = Cks)
    return(df)
}

addrow_from_cas <- function(df, wpid, wpid_version, pathway_name,
                        metabolite_name, mol_id = mol_id){
    system_code = "Ca"
    Ch = get_mapped_id(mapper, mol_id, source = system_code, target = "Ch")
    Ck = get_mapped_id(mapper, mol_id, source = system_code, target = "Ck")
    Ce = get_mapped_id(mapper, mol_id, source = system_code, target = "Ce")
    Cpc = get_mapped_id(mapper, mol_id, source = system_code, target = "Cpc")
    Dr = get_mapped_id(mapper, mol_id, source = system_code, target = "Dr")
    Ik = get_mapped_id(mapper, mol_id, source = system_code, target = "Ik")
    Cs = get_mapped_id(mapper, mol_id, source = system_code, target = "Cs")
    Wd = get_mapped_id(mapper, mol_id, source = system_code, target = "Wd")
    Cks = get_mapped_id(mapper, mol_id, source = system_code, target = "Cks")
    df = add_row(df, wpid = wpid, wpid_version = wpid_version,
            pathway_name = pathway_name,
            metabolite_name = metabolite_name,
            HMDB_ID = Ch, KEGG_ID = Ck, ChEBI_ID = Ce,
            PubChem_CID = Cpc, DrugBank_ID = Dr, CAS = mol_id, InChI_Key = Ik,
            ChemSpider = Cs, Wikidata_ID = Wd, KNaPSAcK_ID = Cks)
    return(df)
}

addrow_from_chemspider <- function(df, wpid, wpid_version, pathway_name,
                        metabolite_name, mol_id = mol_id){
    system_code = "Cs"
    Ch = get_mapped_id(mapper, mol_id, source = system_code, target = "Ch")
    Ck = get_mapped_id(mapper, mol_id, source = system_code, target = "Ck")
    Ce = get_mapped_id(mapper, mol_id, source = system_code, target = "Ce")
    Cpc = get_mapped_id(mapper, mol_id, source = system_code, target = "Cpc")
    Dr = get_mapped_id(mapper, mol_id, source = system_code, target = "Dr")
    Ca = get_mapped_id(mapper, mol_id, source = system_code, target = "Ca")
    Ik = get_mapped_id(mapper, mol_id, source = system_code, target = "Ik")
    Wd = get_mapped_id(mapper, mol_id, source = system_code, target = "Wd")
    Cks = get_mapped_id(mapper, mol_id, source = system_code, target = "Cks")
    df = add_row(df, wpid = wpid, wpid_version = wpid_version,
            pathway_name = pathway_name,
            metabolite_name = metabolite_name,
            HMDB_ID = Ch, KEGG_ID = Ck, ChEBI_ID = Ce,
            PubChem_CID = Cpc, DrugBank_ID = Dr, CAS = Ca, InChI_Key = Ik,
            ChemSpider = mol_id, Wikidata_ID = Wd, KNaPSAcK_ID = Cks)
    return(df)
}

addrow_from_wikidata <- function(df, wpid, wpid_version, pathway_name,
                        metabolite_name, mol_id = mol_id){
    system_code = "Wd"
    Ch = get_mapped_id(mapper, mol_id, source = system_code, target = "Ch")
    Ck = get_mapped_id(mapper, mol_id, source = system_code, target = "Ck")
    Ce = get_mapped_id(mapper, mol_id, source = system_code, target = "Ce")
    Cpc = get_mapped_id(mapper, mol_id, source = system_code, target = "Cpc")
    Dr = get_mapped_id(mapper, mol_id, source = system_code, target = "Dr")
    Ca = get_mapped_id(mapper, mol_id, source = system_code, target = "Ca")
    Ik = get_mapped_id(mapper, mol_id, source = system_code, target = "Ik")
    Cs = get_mapped_id(mapper, mol_id, source = system_code, target = "Cs")
    Cks = get_mapped_id(mapper, mol_id, source = system_code, target = "Cks")
    df = add_row(df, wpid = wpid, wpid_version = wpid_version,
            pathway_name = pathway_name,
            metabolite_name = metabolite_name,
            HMDB_ID = Ch, KEGG_ID = Ck, ChEBI_ID = Ce,
            PubChem_CID = Cpc, DrugBank_ID = Dr, CAS = Ca, InChI_Key = Ik,
            ChemSpider = Cs, Wikidata_ID = mol_id, KNaPSAcK_ID = Cks)
    return(df)
}

addrow_from_knapsack <- function(df, wpid, wpid_version, pathway_name,
                        metabolite_name, mol_id = mol_id){
    system_code = "Cks"
    Ch = get_mapped_id(mapper, mol_id, source = system_code, target = "Ch")
    Ck = get_mapped_id(mapper, mol_id, source = system_code, target = "Ck")
    Ce = get_mapped_id(mapper, mol_id, source = system_code, target = "Ce")
    Cpc = get_mapped_id(mapper, mol_id, source = system_code, target = "Cpc")
    Dr = get_mapped_id(mapper, mol_id, source = system_code, target = "Dr")
    Ca = get_mapped_id(mapper, mol_id, source = system_code, target = "Ca")
    Ik = get_mapped_id(mapper, mol_id, source = system_code, target = "Ik")
    Cs = get_mapped_id(mapper, mol_id, source = system_code, target = "Cs")
    Wd = get_mapped_id(mapper, mol_id, source = system_code, target = "Wd")
    df = add_row(df, wpid = wpid, wpid_version = wpid_version,
                 pathway_name = pathway_name,
                 metabolite_name = metabolite_name,
                 HMDB_ID = Ch, KEGG_ID = Ck, ChEBI_ID = Ce,
                 PubChem_CID = Cpc, DrugBank_ID = Dr, CAS = Ca, InChI_Key = Ik,
                 ChemSpider = Cs, Wikidata_ID = Wd, KNaPSAcK_ID = mol_id)
    return(df)
}

get_mapped_id <- function(mapper, mol_id, source, target){
    mapped_id = ""
    m = map(mapper, mol_id, source = source, target = target)
    if (dim(m)[1] > 0) {
      mapped_id = m[1,][["mapping"]]
    }
    return(mapped_id)
}
