# %%
import os
import re
import argparse
import polars as pl
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from pathlib import Path

dir_path = Path(
    "/Users/wilhelm/Proton/01_juniper/01_arbete/01_projekt/03_psm/retrorules"
)
retros = os.listdir(dir_path)
retro_paths = [f for f in retros if re.search(r".csv$", f)]
metanetx = (
    pl.read_csv(source=dir_path.joinpath("metanetx.csv"))
    .with_columns(pl.col("SCORE").cast(pl.Float64))
)
rhea = (
    pl.read_csv(source=dir_path.joinpath("rhea.csv"))
    .with_columns(pl.col("SCORE").cast(pl.Float64))
)
uspto = (
    pl.read_csv(source=dir_path.joinpath("uspto.csv"))
    .with_columns(pl.col("SCORE").cast(pl.Float64))
)
comb = pl.concat([metanetx, rhea, uspto], how="diagonal")

# Filter to relevant EC classes, valid templates, and a reasonable radius
rules = comb.filter(
    # pl.col("ECS").str.contains(r"^(3\.2\.1|1\.13\.11)") # &
    pl.col("VALID")
    & (pl.col("RADIUS_MAX") <= 4)
)

rules.write_parquet(file=dir_path.joinpath("rxn_rules.parquet"))

# %%
