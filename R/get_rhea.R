library(httr2)
library(readr)

rhea_id <- 25218

tsv <- request("https://www.rhea-db.org/rhea/") |>
    req_url_query(
        query   = paste0("RHEA:", rhea_id),
        columns = "rhea-id,equation,chebi-id,chebi",
        format  = "tsv",
        limit   = 1
    ) |>
    req_user_agent("your-script/1.0 (your.email@domain)") |>
    req_perform() |>
    resp_body_string()

df <- readr::read_tsv(I(tsv), show_col_types = FALSE)
df
<z