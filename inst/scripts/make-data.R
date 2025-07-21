library(xml2)
library(stringr)
library(tibble)
library(dplyr)
library(BridgeDbR)
library(rWikiPathways)

wikipathways_release_date <- "20250710"

# Vertebrates
chicken_url <- paste0("https://zenodo.org/records/15905886/files/wikipathways-", wikipathways_release_date, "-gpml-Gallus_gallus.zip")
chimp_url <- paste0("https://zenodo.org/records/15905897/files/wikipathways-", wikipathways_release_date, "-gpml-Pan_troglodytes.zip")
cow_url <- paste0("https://zenodo.org/records/15905900/files/wikipathways-", wikipathways_release_date, "-gpml-Bos_taurus.zip")
dog_url <- paste0("https://zenodo.org/records/15905903/files/wikipathways-", wikipathways_release_date, "-gpml-Canis_lupus_familiaris.zip")
horse_url <- paste0("https://zenodo.org/records/15905906/files/wikipathways-", wikipathways_release_date, "-gpml-Equus_caballus.zip")
human_url <- paste0("https://zenodo.org/records/15905909/files/wikipathways-", wikipathways_release_date, "-gpml-Homo_sapiens.zip")
mouse_url <- paste0("https://zenodo.org/records/15905912/files/wikipathways-", wikipathways_release_date, "-gpml-Mus_musculus.zip")
pig_url <- paste0("https://zenodo.org/records/15905915/files/wikipathways-", wikipathways_release_date, "-gpml-Sus_scrofa.zip")
rat_url <- paste0("https://zenodo.org/records/15905918/files/wikipathways-", wikipathways_release_date, "-gpml-Rattus_norvegicus.zip")
zebrafish_url <- paste0("https://zenodo.org/records/15905921/files/wikipathways-", wikipathways_release_date, "-gpml-Danio_rerio.zip")

# Plants
arabidopsis_url <- paste0("https://zenodo.org/records/15905924/files/wikipathways-", wikipathways_release_date, "-gpml-Arabidopsis_thaliana.zip")
barley_url <- paste0("https://zenodo.org/records/15905925/files/wikipathways-", wikipathways_release_date, "-gpml-Hordeum_vulgare.zip")
common_wheat_url <- paste0("https://zenodo.org/records/15905926/files/wikipathways-", wikipathways_release_date, "-gpml-Triticum_aestivum.zip")
japanese_rice_url <- paste0("https://zenodo.org/records/15905928/files/wikipathways-", wikipathways_release_date, "-gpml-Oryza_sativa_Japonica.zip")
maize_url <- paste0("https://zenodo.org/records/15905927/files/wikipathways-", wikipathways_release_date, "-gpml-Zea_mays.zip")
poplar_url <- paste0("https://zenodo.org/records/15905929/files/wikipathways-", wikipathways_release_date, "-gpml-Populus_trichocarpa.zip")
tomato_url <- paste0("https://zenodo.org/records/15905936/files/wikipathways-", wikipathways_release_date, "-gpml-Solanum_lycopersicum.zip")
wine_grape_url <- paste0("https://zenodo.org/records/15905937/files/wikipathways-", wikipathways_release_date, "-gpml-Vitis_vinifera.zip")

# Parasites
mararia_parasite_url <- paste0("https://zenodo.org/records/15905930/files/wikipathways-", wikipathways_release_date, "-gpml-Plasmodium_falciparum.zip")

# Invertebrates
c_elegans_url <- paste0("https://zenodo.org/records/15905931/files/wikipathways-", wikipathways_release_date, "-gpml-Caenorhabditis_elegans.zip")
fruitfly_url <- paste0("https://zenodo.org/records/15905932/files/wikipathways-", wikipathways_release_date, "-gpml-Drosophila_melanogaster.zip")
mosquito_url <- paste0("https://zenodo.org/records/15905933/files/wikipathways-", wikipathways_release_date, "-gpml-Aedes_aegypti.zip")

# Fungi
f_graminearum_url <- paste0("https://zenodo.org/records/15905934/files/wikipathways-", wikipathways_release_date, "-gpml-Fusarium_graminearum.zip")
yeast_url <- paste0("https://zenodo.org/records/15905935/files/wikipathways-", wikipathways_release_date, "-gpml-Saccharomyces_cerevisiae.zip")

# Bacteria
a_woodii_url <- paste0("https://zenodo.org/records/15905938/files/wikipathways-", wikipathways_release_date, "-gpml-Aquifex_aeolicus.zip")
b_subtilis_url <- paste0("https://zenodo.org/records/15905939/files/wikipathways-", wikipathways_release_date, "-gpml-Bacillus_subtilis.zip")
c_bibrioides_url <- paste0("https://zenodo.org/records/15905940/files/wikipathways-", wikipathways_release_date, "-gpml-Campylobacter_jejuni.zip")
e_coli_url <- paste0("https://zenodo.org/records/15905941/files/wikipathways-", wikipathways_release_date, "-gpml-Escherichia_coli.zip")
tuberculosis_url <- paste0("https://zenodo.org/records/15905942/files/wikipathways-", wikipathways_release_date, "-gpml-Mycobacterium_tuberculosis.zip")

# bridgedb_url <- "https://figshare.com/ndownloader/files/51201419"
destfile <- "metabolites_20241215.bridge"
# bridgedb_resp <- request(bridgedb_url) |> req_perform()
# writeBin(resp_body_raw(bridgedb_resp), destfile)

mapper <- loadDatabase(destfile)

createWikipathwaysMetabolitesDb <- function() {
    for (species in wikipathways_species) {
        filename = downloadPathwayArchive(date=wikipathways_release_date,
                                          organism=species, format="gpml")
        gpmlpaths = utils::unzip(filename)
        df <- tibble(wpid = character(), wpid_version = character(),
                     pathway_name = character(),
                     metabolite_name = character(),
                     HMDB_ID = character(), KEGG_ID = character(),
                     ChEBI_ID = character(), DrugBank_ID = character(),
                     PubChem_CID = character(), ChemSpider_ID = character(),
                     KNApSAcK_ID = character(), Wikidata_ID = character(),
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
            ChemSpider_ID = Cs, Wikidata_ID = Wd, KNApSAcK_ID = Cks)
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
              ChemSpider_ID = Cs, Wikidata_ID = Wd, KNApSAcK_ID = Cks)
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
            ChemSpider_ID = Cs, Wikidata_ID = Wd, KNApSAcK_ID = Cks)
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
            ChemSpider_ID = Cs, Wikidata_ID = Wd, KNApSAcK_ID = Cks)
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
            ChemSpider_ID = Cs, Wikidata_ID = Wd, KNApSAcK_ID = Cks)
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
            ChemSpider_ID = mol_id, Wikidata_ID = Wd, KNApSAcK_ID = Cks)
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
            ChemSpider_ID = Cs, Wikidata_ID = mol_id, KNApSAcK_ID = Cks)
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
                 ChemSpider_ID = Cs, Wikidata_ID = Wd, KNApSAcK_ID = mol_id)
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
