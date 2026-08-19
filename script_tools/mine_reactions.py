# %%
import argparse
import polars as pl
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from pathlib import Path
# %%
parser = argparse.ArgumentParser(
    description="Filter and report on candidate hits vs one or more reference compounds"
)
parser.add_argument(
    "-i", "--input",
    required=True,
    help="Reaction rules parquet"
)
parser.add_argument(
    "-s", "--smiles",
    required=True,
    help=(
        "Names and SMILES of one or more reference compounds, "
        "e.g. 'name1,smiles1;name2,smiles2'"
    )
)
parser.add_argument(
    "-o", "--output",
    default=None,
    help="Path to output folder [default: <input>/report/]"
)

# %%
args = parser.parse_args()
smiles_list = args.smiles.split(";")

# smiles_arg = "aspirin,CC(=O)OC1=CC=CC=C1C(=O)O;caffeine,CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
# smiles_list = smiles_arg.split(";")

ref_smiles = dict(x.split(",") for x in smiles_list)
rxn_path = Path(args.input)
output = (
    Path(args.output).joinpath("rxns.parquet")
    if args.output
    else rxn_path.parent.joinpath("rxns.parquet")
)
rxn_rules = pl.read_parquet(source=rxn_path)

# %%

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
                            "source": row["DATASETS"]
                        })
                    except Exception:
                        continue
        except Exception:
            continue

    products_df = pl.DataFrame(results).unique(subset="product_smiles")
    return(products_df)

# %%
all_frames = []
for i, k in ref_smiles.items():
    rxns_it = (
        mine_reactions(rxn_rules, k)
        .with_columns(pl.lit(i).alias("compound"))
        .with_columns(pl.col("mass_shift").abs().alias("abs_mass_shift"))
        .sort("abs_mass_shift", descending=False)
        .filter(pl.col("abs_mass_shift") > 0)
    )
    all_frames.append(rxns_it)

all_rxn = pl.concat(all_frames)

output.parent.mkdir(parents=True, exist_ok=True)
all_rxn.write_parquet(output)
