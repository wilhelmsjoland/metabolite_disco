library(CompoundDb)

## --- ChEBI ---
chebi_sdf <- "annotation_databases/chebi_lite_3_stars.sdf.gz"   # your downloaded file
chebi_tbl <- compound_tbl_sdf(chebi_sdf, collapse = "|")

chebi_meta <- make_metadata(
    source = "ChEBI",
    url = "https://www.ebi.ac.uk/chebi",
    source_version = "248",
    source_date = "2026-01-13",
    organism = NA_character_
)

chebi_db_file <- createCompDb(chebi_tbl, metadata = chebi_meta, path = ".")
chebi_db <- CompDb(chebi_db_file)

## --- PubChem (subset SDF you downloaded) ---
pubchem_sdf <- "pubchem_subset.sdf.gz"
pubchem_tbl <- compound_tbl_sdf(pubchem_sdf)

pubchem_meta <- make_metadata(
    source = "PubChem",
    url = "https://pubchem.ncbi.nlm.nih.gov",
    source_version = "YOUR_PUBCHEM_RELEASE_OR_EXPORT_TAG",
    source_date = "YYYY-MM",
    organism = NA_character_
)

pubchem_db_file <- createCompDb(pubchem_tbl, metadata = pubchem_meta, path = ".")
pubchem_db <- CompDb(pubchem_db_file)
