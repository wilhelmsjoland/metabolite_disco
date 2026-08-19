# %%
import os
import re
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
    .with_columns(pl.lit("metanetx").alias("SOURCE"))
)
rhea = (
    pl.read_csv(source=dir_path.joinpath("rhea.csv"))
    .with_columns(pl.col("SCORE").cast(pl.Float64))
    .with_columns(pl.lit("rhea").alias("SOURCE"))
)
uspto = (
    pl.read_csv(source=dir_path.joinpath("uspto.csv"))
    .with_columns(pl.col("SCORE").cast(pl.Float64))
    .with_columns(pl.lit("uspto").alias("SOURCE"))
)
comb = pl.concat([metanetx, rhea, uspto], how="diagonal")
# %%
# Filter to relevant EC classes, valid templates, and a reasonable radius
rules = comb.filter(
    # pl.col("ECS").str.contains(r"^(3\.2\.1|1\.13\.11)") # &
    pl.col("VALID")
    & (pl.col("RADIUS_MAX") <= 4)
)

# %% 
def mine_reactions(df, smiles):
    substrate = Chem.MolFromSmiles(smiles)
    substrate_mass = Descriptors.ExactMolWt(substrate)
    results = []
    for row in df.iter_rows(named=True):
        try:
            rxn = AllChem.ReactionFromSmarts(row["TEMPLATE"])
            for product_set in rxn.RunReactants((substrate,)):
                for p in product_set:
                    try:
                        Chem.SanitizeMol(p)
                        smi = Chem.MolToSmiles(p)
                        mass = Descriptors.ExactMolWt(p)
                        results.append({
                            "rule_id": row["TEMPLATE_ID"],
                            "ec": row["ECS"],
                            "score": row["SCORE"],
                            "product_smiles": smi,
                            "mass_shift": mass - substrate_mass,
                            "source": row["SOURCE"]
                        })
                    except Exception:
                        continue
        except Exception:
            continue

    products_df = pl.DataFrame(results).unique(subset="product_smiles")
    return(products_df)


# %%
all_rxn = mine_reactions(rules, "C1=CC(=CC=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)O)O")

# %%
all_rxn
# %%
