# %%
import os
from pathlib import Path
import polars as pl
import altair as alt
import seaborn as sns
import matplotlib.pyplot as plt

alt.data_transformers.enable("vegafusion")
# %%
data = pl.read_parquet(
    source=Path(
        f"/Volumes/bluecub/aglycone_release_100um_24h/combined_results/"
        f"aglycone_release_100um_24h_hits.parquet"
    )
)

# %%
counts = (
    data.filter(
        pl.col("sim") > 0.8,
        pl.col("rtime") > 30
        # pl.col("source") == "bio"
    )
    .select(["experiment", "inchikey", "sim"])
    .group_by("experiment", "inchikey")
    .agg(
        pl.len().alias("n")
    )
)

# %%
counts

# %%
data
# %%
alt.Chart(counts).mark_rect().encode(
    x=alt.X("experiment:N"),
    y=alt.Y("inchikey:N"),
    color=alt.Color("n:Q", scale=alt.Scale(scheme="viridis")),
    tooltip=["inchikey", "experiment", "n"],
)

# %%