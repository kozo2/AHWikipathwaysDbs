library(xml2)
library(stringr)
library(tibble)
library(BridgeDbR)

wikipathways_release_date <- "20210310"

wikipathways_species <- c("Anopheles gambiae", "Arabidopsis thaliana",
    "Bacillus subtilis", "Bos taurus", "Caenorhabditis elegans",
    "Canis familiaris", "Danio rerio", "Drosophila melanogaster",
    "Equus caballus", "Escherichia coli", "Gallus gallus", "Gibberella zeae",
    "Homo sapiens", "Hordeum vulgare", "Mus musculus",
    "Mycobacterium tuberculosis", "Oryza sativa", "Pan troglodytes",
    "Plasmodium falciparum", "Populus trichocarpa", "Rattus norvegicus",
    "Saccharomyces cerevisiae", "Solanum lycopersicum", "Sus scrofa",
    "Zea mays")

for (species in wikipathways_species) {
    filename = downloadPathwayArchive(date=wikipathways_release_date,
                                    organism=species, format="gpml")
    gpmlpaths = utils::unzip(filename)
    for (gpmlpath in gpmlpaths) {
        readGpml(gpmlpath)
    }
}

readGpml <- function(gpmlpath) {
    x <- read_xml(gpmlpath)
    filename <- tools::file_path_sans_ext(gpmlpath)
    wpidinfo <- tail(str_split(filename, "_", simplify = TRUE)[1,], n=2)
    wpid <- wpidinfo[1]
    wpid_version <- wpidinfo[2]
    pathway_name <- xml_attr(x, "Name")
    df <- tibble(wpid = character(), wpid_version = character(),
                 pathway_name = character(), pathway_description = character(),
                 pathway_category = character(), metabolite_name = character(),
                 HMDB_ID = character(), KEGG_ID = character(),
                 ChEBI_ID = character(), DrugBank_ID = character(),
                 PubChem_CID = character(), CAS = character(),
                 InChI_Key =character())
                 # database = character(), mol_id = character())
    pathway_description = ""
    pathway_category = ""
    for (xml_node in xml_children(x)) {
        if (xml_name(xml_node) == "Comment") {
            if (xml_attr(xml_node, "Source") == "WikiPathways-description") {
              pathway_description = xml_text(xml_node)
            } else if (xml_attr(xml_node, "Source") == "WikiPathways-category"){
              pathway_category = xml_text(xml_node)
            }
        } else if (xml_name(xml_node) == "DataNode") {
            if (xml_attr(xml_node, "Type") == "Metabolite") {
                metabolite_name = xml_attr(xml_node, "TextLabel")
                id_node = xml_children(xml_node)[2]
                database = xml_attr(id_node, "Database")
                mol_id = xml_attr(id_node, "ID")
                if (database == "HMDB") {
                    df = addrow_from_hmdb(df, wpid, wpid_version,
                                          pathway_name, pathway_description,
                                          pathway_category, metabolite_name,
                                          mol_id)
                } else if (database == "ChEBI") {
                    df = addrow_from_chebi(df, wpid, wpid_version,
                                          pathway_name, pathway_description,
                                          pathway_category, metabolite_name,
                                          mol_id)
                } else if (database == "KEGG Compound") {
                    df = addrow_from_kegg(df, wpid, wpid_version,
                                         pathway_name, pathway_description,
                                         pathway_category, metabolite_name,
                                         mol_id)
                } else if (database == "PubChem-compound") {
                    df = addrow_from_pubchem(df, wpid, wpid_version,
                                         pathway_name, pathway_description,
                                         pathway_category, metabolite_name,
                                         mol_id)
                } else {
                    print(paste(database, "is needed to addrow_from"))
                }
            # } else if (xml_attr(xml_node, "Type") == "Protein") {
            #   protein_name = xml_attr(xml_node, "TextLabel")
            #   print(protein_name)
            #   id_node = xml_children(xml_node)[2]
            #   database = xml_attr(id_node, "Database")
            #   print(database)
            #   mol_id = xml_attr(id_node, "ID")
            #   print(mol_id)
            # } else if (xml_attr(xml_node, "Type") == "GeneProduct") {
            #   gene_name = xml_attr(xml_node, "TextLabel")
            #   print(gene_name)
            #   id_node = xml_children(xml_node)[2]
            #   database = xml_attr(id_node, "Database")
            #   print(database)
            #   mol_id = xml_attr(id_node, "ID")
            #   print(mol_id)
            }
        }
    }
    print(df)
}

addrow_from_hmdb <- function(df, wpid, wpid_version, pathway_name,
                             pathway_description, pathway_category,
                             metabolite_name, mol_id = mol_id){
    system_code = "Ch"
    kegg_id = get_mapped_id(mapper, mol_id, source = system_code, target = "Ck")
    chebi_id = get_mapped_id(mapper, mol_id, source = system_code, target = "Ce")
    pubchem_id = get_mapped_id(mapper, mol_id, source = system_code, target = "Cpc")
    drugbank_id = get_mapped_id(mapper, mol_id, source = system_code, target = "Dr")
    cas = get_mapped_id(mapper, mol_id, source = system_code, target = "Ca")
    inchi_key = get_mapped_id(mapper, mol_id, source = system_code, target = "Ik")
    df = add_row(df, wpid = wpid, wpid_version = wpid_version,
                 pathway_name = pathway_name,
                 pathway_description = pathway_description,
                 pathway_category = pathway_category,
                 metabolite_name = metabolite_name,
                 HMDB_ID = mol_id,
                 KEGG_ID = kegg_id,
                 ChEBI_ID = chebi_id,
                 PubChem_CID = pubchem_id,
                 DrugBank_ID = drugbank_id,
                 CAS = cas,
                 InChI_Key = inchi_key)
    return(df)
}

addrow_from_chebi <- function(df, wpid, wpid_version, pathway_name,
                              pathway_description, pathway_category,
                              metabolite_name, mol_id = mol_id){
    system_code = "Ce"
    kegg_id = get_mapped_id(mapper, mol_id, source = system_code, target = "Ck")
    hmdb_id = get_mapped_id(mapper, mol_id, source = system_code, target = "Ch")
    pubchem_id = get_mapped_id(mapper, mol_id, source = system_code, target = "Cpc")
    drugbank_id = get_mapped_id(mapper, mol_id, source = system_code, target = "Dr")
    cas = get_mapped_id(mapper, mol_id, source = system_code, target = "Ca")
    inchi_key = get_mapped_id(mapper, mol_id, source = system_code, target = "Ik")
    df = add_row(df, wpid = wpid, wpid_version = wpid_version,
                 pathway_name = pathway_name,
                 pathway_description = pathway_description,
                 pathway_category = pathway_category,
                 metabolite_name = metabolite_name,
                 HMDB_ID = hmdb_id,
                 KEGG_ID = kegg_id,
                 ChEBI_ID = mol_id,
                 PubChem_CID = pubchem_id,
                 DrugBank_ID = drugbank_id,
                 CAS = cas,
                 InChI_Key = inchi_key)
    return(df)
}

addrow_from_kegg <- function(df, wpid, wpid_version, pathway_name,
                             pathway_description, pathway_category,
                             metabolite_name, mol_id = mol_id){
    system_code = "Ck"
    chebi_id = get_mapped_id(mapper, mol_id, source = system_code, target = "Ce")
    hmdb_id = get_mapped_id(mapper, mol_id, source = system_code, target = "Ch")
    pubchem_id = get_mapped_id(mapper, mol_id, source = system_code, target = "Cpc")
    drugbank_id = get_mapped_id(mapper, mol_id, source = system_code, target = "Dr")
    cas = get_mapped_id(mapper, mol_id, source = system_code, target = "Ca")
    inchi_key = get_mapped_id(mapper, mol_id, source = system_code, target = "Ik")
    df = add_row(df, wpid = wpid, wpid_version = wpid_version,
                 pathway_name = pathway_name,
                 pathway_description = pathway_description,
                 pathway_category = pathway_category,
                 metabolite_name = metabolite_name,
                 HMDB_ID = hmdb_id,
                 KEGG_ID = mol_id,
                 ChEBI_ID = chebi_id,
                 PubChem_CID = pubchem_id,
                 DrugBank_ID = drugbank_id,
                 CAS = cas,
                 InChI_Key = inchi_key)
    return(df)
}

addrow_from_pubchem <- function(df, wpid, wpid_version, pathway_name,
                                pathway_description, pathway_category,
                                metabolite_name, mol_id = mol_id){
    system_code = "Cpc"
    chebi_id = get_mapped_id(mapper, mol_id, source = system_code, target = "Ce")
    hmdb_id = get_mapped_id(mapper, mol_id, source = system_code, target = "Ch")
    kegg_id = get_mapped_id(mapper, mol_id, source = system_code, target = "Ck")
    drugbank_id = get_mapped_id(mapper, mol_id, source = system_code, target = "Dr")
    cas = get_mapped_id(mapper, mol_id, source = system_code, target = "Ca")
    inchi_key = get_mapped_id(mapper, mol_id, source = system_code, target = "Ik")
    df = add_row(df, wpid = wpid, wpid_version = wpid_version,
                 pathway_name = pathway_name,
                 pathway_description = pathway_description,
                 pathway_category = pathway_category,
                 metabolite_name = metabolite_name,
                 HMDB_ID = hmdb_id,
                 KEGG_ID = kegg_id,
                 ChEBI_ID = chebi_id,
                 PubChem_CID = mol_id,
                 DrugBank_ID = drugbank_id,
                 CAS = cas,
                 InChI_Key = inchi_key)
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

# metabolite_bridge_url <- "https://ndownloader.figshare.com/files/25388453"
# download.file(metabolite_bridge_url, "metabolites_20201104.bridge")

file <- "metabolites_20210109.bridge"
# download.file(
#     "https://ndownloader.figshare.com/files/26001794",
#     location
# )
# location = normalizePath(file)
mapper <- loadDatabase(location)

# see https://bridgedb.github.io/pages/system-codes.html
# map(mapper, "HMDB0001058", source="Ch", target="Ck")
# map(mapper, "HMDB0001058", source="Ch", target="Ce")
# map(mapper, "HMDB0001058", source="Ch", target="Ca")
# map(mapper, "HMDB0001058", source="Ch", target="Dr")
# map(mapper, "HMDB0001058", source="Ch", target="Ik")
