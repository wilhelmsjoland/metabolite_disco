# %%
from pathlib import Path
import pandas as pd
import os
import rdkit
from rdkit.Chem import rdFingerprintGenerator
import argparse

# %%
parser = argparse.ArgumentParser(
    description="Calculate similarity vs one more more compounds"
)
parser.add_argument(
    "-i", "--input",
    required=True,
    help="Path to metabolite_disco output directory"
)
parser.add_argument(
    "-s", "--smiles",
    required=True,
    help="Names and SMILES of one or more reference compounds, "
         "e.g. 'name1,smiles1;name2,smiles2'"
)
parser.add_argument(
    "-r", "--radius",
    type=int,
    default=3,
    help= "radius for morganfingerprint"
)
parser.add_argument(
    "-f", "--fpsize",
    type=int,
    default=2048,
    help= "fpsisze for morganfingerprint"
)
parser.add_argument(
    "-o", "--output",
    default=None,
    help="Path to output folder [default: <input>/report/]"
)

args = parser.parse_args()

# %%
# This should be an argparse
# exp_dir = Path(
#     f"/Volumes/bluecub/aglycone_release_100um_24h/output/"
#     f"afzelin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon"
# )
# ref_smiles = {
#     "kaempferol": "C1=CC(=CC=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O",
#     "afzelin": "C[C@H]1[C@@H]([C@H]([C@H]([C@@H](O1)OC2=C(OC3=CC(=CC(=C3C2=O)O)O)C4=CC=C(C=C4)O)O)O)O"
# }
# output_dir = Path(exp_dir).joinpath("report")

exp_dir = Path(args.input)
smiles_str = args.smiles
output_dir = Path(args.output) if args.output else Path(exp_dir).joinpath("report")
smiles_list = smiles_str.split(";")
ref_smiles = dict(x.split(",") for x in smiles_list)
radius_arg = args.radius
fpsize_arg = args.fpsize

# %%
def compute_similarity(df, reference_smiles, smiles_col):
    fps_df = (
        df
        .dropna(subset=[smiles_col])
        .reset_index(drop=True)
        .assign(
            fps=lambda df: df[smiles_col].apply(
                lambda x: fpgen.GetFingerprint(rdkit.Chem.MolFromSmiles(x))
            )
        )
    )

    sim_df = pd.concat(
        [
            fps_df.assign(
                sim_mol=name,
                sim=rdkit.DataStructs.BulkTanimotoSimilarity(
                    fpgen.GetFingerprint(rdkit.Chem.MolFromSmiles(smiles)),
                    fps_df["fps"]
                ),
            )
            for name, smiles in reference_smiles.items()
        ],
        ignore_index=True,
    ).drop(columns=["fps"])

    return sim_df


# %% 
# radius=1~ECFP2, radius=2~ECFP4, radius=3~ECFP6
fpgen = rdkit.Chem.rdFingerprintGenerator.GetMorganGenerator(
    radius=radius_arg,
    fpSize=fpsize_arg
)

biot_path = Path(exp_dir).joinpath("tables", "biotransformer_predictions.csv")
anno_path = Path(exp_dir).joinpath("tables", "annotation_predictions.csv")

biot = (
    pd.read_csv(filepath_or_buffer=biot_path)
    .drop(columns=["input"])
)
anno = pd.read_csv(filepath_or_buffer=anno_path)

anno_sims = compute_similarity(
    df=anno,
    reference_smiles=ref_smiles,
    smiles_col="smiles"
)
biot_sims = compute_similarity(
    df=biot,
    reference_smiles=ref_smiles,
    smiles_col="SMILES"
)

# %%
os.makedirs(Path(output_dir), exist_ok=True)
biot_sims.to_parquet(
    path=(
        output_dir.joinpath("biotransformer_similarities.parquet")
    ),
    index=False
)
anno_sims.to_parquet(
    path=(
        output_dir.joinpath("annotation_similarities.parquet")
    ),
    index=False
)
