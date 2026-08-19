# %%
import os
import re
from pathlib import Path
import polars as pl
import altair as alt
import seaborn as sns
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw
import base64
from io import BytesIO
from scipy.cluster.hierarchy import linkage, leaves_list

alt.data_transformers.enable("vegafusion")

# %%
# def smiles_to_data_uri(smiles: str, size: int = 150) -> str:
#     mol = Chem.MolFromSmiles(smiles)
#     options = Draw.MolDrawOptions()
#     options.padding = 0.0
#     img = Draw.MolToImage(mol, size=(size, size), options=options)
#     buf = BytesIO()
#     img.save(buf, format="PNG")
#     b64 = base64.b64encode(buf.getvalue()).decode()
#     return f"data:image/png;base64,{b64}"

def smiles_to_data_uri(smiles: str, size: int = 150) -> str:
    mol = Chem.MolFromSmiles(smiles)
    drawer = Draw.rdMolDraw2D.MolDraw2DSVG(size, size)
    drawer.drawOptions().padding = 0.01
    drawer.drawOptions().clearBackground = False
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    b64 = base64.b64encode(svg.encode()).decode()
    return f"data:image/svg+xml;base64,{b64}"

def get_smiles_from_skeleton(skeleton):
    return Chem.MolToSmiles(
        Chem.MolFromSmiles(
            data.filter(pl.col("skeleton_key") == skeleton)
            .get_column("smiles")
            .first()
        ),
        isomericSmiles=False,
    )

def cluster_order(
    df: pl.DataFrame,
    x: str,
    y: str,
    values: str) -> tuple[list, list]:
    matrix = df.pivot(index=y, on=x, values=values).fill_null(0)
    y_labels = matrix[y].to_list()
    mat_values = matrix.drop(y).to_numpy()
    y_order = [y_labels[i] for i in leaves_list(linkage(mat_values, method="average"))]
    x_labels = matrix.drop(y).columns
    x_order = [x_labels[i] for i in leaves_list(linkage(mat_values.T, method="average"))]
    return x_order, y_order

# %%
data = pl.read_parquet(
    source=Path(
        f"/Volumes/bluecub/aglycone_release_100um_24h/combined_results/"
        f"aglycone_release_100um_24h_hits.parquet"
    )
)
cols = data.columns
cols.remove("experiment")
cols.insert(cols.index("feature") + 1, "experiment")
data = data.select(cols)

# %%
exp_info_path = Path(
    "/Users/wilhelm/Proton/01_juniper/01_arbete/01_projekt/03_psm",
    "molekyler/information/aglycone_release_100um_24h",
    "experiment_matrix.xlsx"
)
exp_info = (
    pl.read_excel(exp_info_path)
    .drop(
        "sample",
        "group",
        "experiment_id",
        "unclean_strain",
        "unclean_condition")
    .filter(pl.col("condition") != "ycfa_glucose")
    .group_by("experiment")
    .agg(pl.all().exclude("experiment").unique())
    .explode("condition")
    .with_columns(
        pl.col("strain").list.sort().list.join("_and_")
    )
)

# %%
# This is not entirely what I want because on experiment that
# gets a bunch of the same scaffold biases this somewhat
################################################################################
# counts -----------------------------------------------------------------------
################################################################################
counts = (
    data
    .filter(pl.col("sim") > 0.3)
    .group_by("experiment", "scaffold_smiles")
    .agg(pl.len().alias("n"))
    .with_columns(
        pl.len().over("scaffold_smiles").alias("n_experiments"),
        pl.col("n").sum().over("scaffold_smiles").alias("n_total"),
    )
    .sort(["n_experiments", "n"], descending=True)
)

alt.Chart(counts).mark_rect().encode(
    x=alt.X("experiment:N"),
    y=alt.Y("scaffold_smiles:N"),
    color=alt.Color("n:Q", scale=alt.Scale(scheme="viridis")),
    tooltip=["scaffold_smiles", "experiment", "n"],
)
# .properties(width=1200, height=1200)

# %%
################################################################################
# presence/absence SMILES - ----------------------------------------------------
################################################################################
# pull out the interesting scaffolds from this
scaffold_prevalence = (
    data
    .filter(pl.col("source") == "biotransformer")
    .filter(pl.col("sim") > 0.2)
    .with_columns(
        pl.col("scaffold_smiles").map_elements(
            lambda smi: Chem.MolToSmiles(
                mol,
                isomericSmiles=False
            ) if (mol := Chem.MolFromSmiles(smi)) is not None else None,
            return_dtype=pl.String,
        )
    )
    .group_by("scaffold_smiles")
    .agg(
        pl.col("experiment").n_unique().alias("n_experiments"),
    )
    .sort("n_experiments", descending=True)
    # .filter(pl.col("n_experiments") > 5)
)

scaffold_prevalence_img = scaffold_prevalence.with_columns(
    pl.col("scaffold_smiles")
    .map_elements(
        smiles_to_data_uri,
        return_dtype=pl.String
    )
    .alias("structure")
)

base = alt.Chart(scaffold_prevalence_img).encode(
    x=alt.X(
        "scaffold_smiles:N",
        sort=alt.EncodingSortField(field="n_experiments", order="descending"),
        axis=alt.Axis(labelAngle=-45),
    ),
).properties(width=alt.Step(50))

points = base.mark_point().encode(
    y=alt.Y("n_experiments:Q"),
    tooltip=["scaffold_smiles", "n_experiments"],
)
structures = base.mark_image(width=50, height=50).encode(
    url="structure:N",
    x=alt.X(
        "scaffold_smiles:N",
        sort=alt.EncodingSortField(
            field="n_experiments",
            order="descending"),
            axis=None
        ),
    y=alt.value(20),
).properties(height=60)

(
    alt.vconcat(structures, points)
    .resolve_scale(x="shared")
    .configure_view(strokeWidth=0)
)



# %%
################################################################################
# presence/absence skeleton across all experiments -----------------------------
################################################################################
skeleton_prevalence = (
    data
    .filter(
        # pl.col("source") == "biotransformer",
        # pl.col("rtime") > 30,
        # pl.col("sim") > 0.1,
        # pl.col("sim") != 1,

        pl.col("source") == "biotransformer",
        pl.col("sim") > 0.1,
        pl.col("sim") != 1,
        # Filter for confidence - if others exists that's fine
        # but without them having [M-H]- it's almost impossible to know
        pl.col("adduct") == "[M-H]-",
    )
    .group_by("skeleton_key")
    .agg(
        pl.col("experiment").n_unique().alias("n_experiments"),
        pl.col("smiles").first().alias("smiles"),
    )
    .sort("n_experiments", descending=True)
    .filter(pl.col("n_experiments") > 5)
)

skeleton_prevalence_img = skeleton_prevalence.with_columns(
    pl.col("smiles")
    .map_elements(smiles_to_data_uri, return_dtype=pl.String)
    .alias("structure")
)

# %%
base = alt.Chart(skeleton_prevalence_img).encode(
    x=alt.X(
        "skeleton_key:N",
        sort=alt.EncodingSortField(field="n_experiments", order="descending"),
        axis=alt.Axis(labelAngle=-45),
    ),
).properties(width=alt.Step(50))

points = base.mark_point().encode(
    y=alt.Y("n_experiments:Q"),
    tooltip=["skeleton_key", "n_experiments", "smiles"],
)
structures = base.mark_image(width=50, height=50).encode(
    url="structure:N",
    x=alt.X(
        "skeleton_key:N",
        sort=alt.EncodingSortField(field="n_experiments", order="descending"),
        axis=None,
    ),
    y=alt.value(20),
).properties(height=60)

(
    alt.vconcat(structures, points)
    .resolve_scale(x="shared")
    .configure_view(strokeWidth=0)
)

# ################################################################################
# # presence/absence skeleton across strains -------------------------------------
# ################################################################################
# # %%
# skeleton_prevalence = (
#     data
#     .join(
#         other=exp_info,
#         on="experiment",
#         how="left"
#     )
#     .filter(
#         # pl.col("source") == "biotransformer",
#         # pl.col("rtime") > 30,
#         # pl.col("sim") > 0.1,
#         # pl.col("sim") != 1,

#         pl.col("source") == "biotransformer",
#         pl.col("sim") > 0.1,
#         # pl.col("sim") != 1,
#         # Filter for confidence - if others exists that's fine
#         # but without them having [M-H]- it's almost impossible to know
#         pl.col("adduct") == "[M-H]-",
#     )
#     .group_by("skeleton_key", "condition")
#     .agg(
#         pl.col("experiment").n_unique().alias("n_experiments"),
#         pl.col("smiles").first().alias("smiles"),
#     )
#     .sort("n_experiments", descending=True)
#     # .filter(pl.col("n_experiments") > 2)
# )

# x_order, y_order = cluster_order(
#     skeleton_prevalence,
#     x="skeleton_key",
#     y="condition",
#     values="n_experiments"
# )

# # %%
# skeleton_prevalence_img = skeleton_prevalence.with_columns(
#     pl.col("smiles")
#     .map_elements(smiles_to_data_uri, return_dtype=pl.String)
#     .alias("structure")
# )

# base = alt.Chart(skeleton_prevalence_img).encode(
#     x=alt.X("skeleton_key:N", sort=x_order, axis=alt.Axis(labelAngle=-45)),
# ).properties(width=alt.Step(50))

# heatmap = base.mark_rect().encode(
#     y=alt.Y("condition:N", sort=y_order),
#     color=alt.Color(
#         "n_experiments:O",
#         scale=alt.Scale(
#             domain=[1, 2, 3, 4],
#             range=["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"],
#         ),
#         legend=alt.Legend(title="Putative production across strains"),
#     ),
#     tooltip=["skeleton_key", "condition", "n_experiments", "smiles"],
# ).properties(height=alt.Step(30))

# structures = base.mark_image(width=50, height=50).encode(
#     url="structure:N",
#     x=alt.X("skeleton_key:N", sort=x_order, axis=None),
#     y=alt.value(20),
# ).properties(height=60)

# (
#     alt.vconcat(structures, heatmap)
#     .resolve_scale(x="shared")
#     .configure_view(strokeWidth=0)
# )

# ################################################################################
# # presence/absence skeleton across glycoside -----------------------------------
# ################################################################################
# # %%
# skeleton_prevalence = (
#     data
#     .join(
#         other=exp_info,
#         on="experiment",
#         how="left"
#     )
#     .filter(
#         # pl.col("source") == "biotransformer",
#         # pl.col("rtime") > 30,
#         # pl.col("sim") > 0.1,
#         # pl.col("sim") != 1,

#         pl.col("source") == "biotransformer",
#         pl.col("sim") > 0.1,
#         # pl.col("sim") != 1,
#         # Filter for confidence - if others exists that's fine
#         # but without them having [M-H]- it's almost impossible to know
#         pl.col("adduct") == "[M-H]-",
#     )
#     .group_by("skeleton_key", "strain")
#     .agg(
#         pl.col("experiment").n_unique().alias("n_experiments"),
#         pl.col("smiles").first().alias("smiles"),
#     )
#     .sort("n_experiments", descending=True)
#     # .filter(pl.col("n_experiments") > 2)
# )

# skeleton_prevalence_img = skeleton_prevalence.with_columns(
#     pl.col("smiles")
#     .map_elements(smiles_to_data_uri, return_dtype=pl.String)
#     .alias("structure")
# )

# x_order, y_order = cluster_order(
#     skeleton_prevalence,
#     x="skeleton_key",
#     y="strain",
#     values="n_experiments"
# )

# # %%
# base = alt.Chart(skeleton_prevalence_img).encode(
#     x=alt.X("skeleton_key:N", sort=x_order, axis=alt.Axis(labelAngle=-45)),
# ).properties(width=alt.Step(50))

# heatmap = base.mark_rect().encode(
#     y=alt.Y("strain:N", sort=y_order),
#     color=alt.Color(
#         "n_experiments:Q",
#         legend=alt.Legend(title="Putative production across glycosides"),
#     ),
#     tooltip=["skeleton_key", "strain", "n_experiments", "smiles"],
# ).properties(height=alt.Step(30))

# structures = base.mark_image(width=50, height=50).encode(
#     url="structure:N",
#     x=alt.X("skeleton_key:N", sort=x_order, axis=None),
#     y=alt.value(20),
# ).properties(height=60)

# (
#     alt.vconcat(structures, heatmap)
#     .resolve_scale(x="shared")
#     .configure_view(strokeWidth=0)
# )

################################################################################
# Interesting candidate hits ---------------------------------------------------
################################################################################

# %%
# sulfated glycoside
get_smiles_from_skeleton("UBZZAJXRNULYDE")
# This is dopamine so obviously interesting
get_smiles_from_skeleton("VYFYYTLLBUKUHU")
# DL-tyrosine
get_smiles_from_skeleton("OUYCCCASQSFEME")
# Vanillylmandelic Acid
get_smiles_from_skeleton("CGQCWMIAEPEHNQ")
# Sulfated molecules
get_smiles_from_skeleton("OGSRMTMBMOFFQC")
# Sulfated molecule with [M-H]-
get_smiles_from_skeleton("MGUNHDBGRGWJMR")

# %%
get_smiles_from_skeleton("MKUXAQIIEYXACX")

# %%
# Check also specifically for N and S and other elements in the string
(
    skeleton_prevalence
    .filter(pl.col("smiles").str.contains("S|N"))
)

# %%
pl.Config.set_tbl_rows(-1)
(
    data.filter(pl.col("skeleton_key") == "MKUXAQIIEYXACX")
    .sort(["adduct", "experiment"])
)
# %%
pl.Config.restore_defaults()
# Check the peaks for the hit for dopamine across experiments
# Need a script for this....

# %%

# %%
