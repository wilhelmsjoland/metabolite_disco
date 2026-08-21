# %%
import polars as pl
from rdkit import Chem
from microberx import MetabolitePredictor

smiles_input = "kaempferol,C1=CC(=CC=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O;afzelin,C[C@H]1[C@@H]([C@H]([C@H]([C@@H](O1)OC2=C(OC3=CC(=CC(=C3C2=O)O)O)C4=CC=C(C=C4)O)O)O)O"

# %%
smiles_list = smiles_input.split(";")
ref_smiles = dict(x.split(",") for x in smiles_list)

# %%
human_frames = []
for name, smi in ref_smiles.items():
    predictor = MetabolitePredictor(
        Chem.MolFromSmiles(smi),
        query_name=name,
        cut_off=0.6,
        biosystem="human",
    )
    predictor.run_prediction()
    human_frames.append(
        pl.from_pandas(predictor.predicted_metabolites)
        .with_columns(
            pl.lit(name).alias("compound"),
            pl.lit("microberx_human").alias("source"),
        )
    )

# %%
gutmicrobes_frames = []
for name, smi in ref_smiles.items():
    predictor = MetabolitePredictor(
        Chem.MolFromSmiles(smi),
        query_name=name,
        cut_off=0.6,
        biosystem="gutmicrobes",
    )
    predictor.run_prediction()
    gutmicrobes_frames.append(
        pl.from_pandas(predictor.predicted_metabolites)
        .with_columns(
            pl.lit(name).alias("compound"),
            pl.lit("microberx_gutbacteria").alias("source"),
        )
    )

# %%
results = pl.concat(human_frames + gutmicrobes_frames)

# %%
results

# %%
(
    results
    .filter(pl.col("similarity_products") == 1)
)

# %%
# The output of main_product vs secondary_products_smiles sucks it seems
(
    results
    .filter(pl.col("similarity_products") == 1)
)["secondary_products_smiles"][0]
# %%
