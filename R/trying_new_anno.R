library(CompoundDb)


# Use FULL SDF (has formula + monoisotopic mass). 3-star is smaller/faster.
chebi_sdf <- file.path("annotation_databases", "chebi_3_stars.sdf.gz")
chebi_sdf_readme <- file.path("annotation_databases", "chebi_README.txt")

if (!file.exists(chebi_sdf)) {
    curl::curl_download(
        url = "https://ftp.ebi.ac.uk/pub/databases/chebi/SDF/chebi_3_stars.sdf.gz",
        destfile = chebi_sdf
    )
}

if (!file.exists(chebi_sdf_readme)) {
    curl::curl_download(
        url = "https://ftp.ebi.ac.uk/pub/databases/chebi/SDF/README",
        destfile = chebi_sdf_readme
    )
}

cmps <- CompoundDb::compound_tbl_sdf(chebi_sdf)

# From the ChEBI SDF README (release + date)
metad <- CompoundDb::make_metadata(
    source = "ChEBI",
    url = "https://www.ebi.ac.uk/chebi/",
    source_version = "247",
    source_date = "2025-12-03",
    organism = NA
)

db_file <- CompoundDb::createCompDb(cmps, metadata = metad, path = "annotation_databases")
cdb_chebi <- CompoundDb::CompDb(db_file)
cdb_chebi
