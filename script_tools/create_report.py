
# %%
# TODO Needs cleaning!
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import DataStructs, inchi, rdFingerprintGenerator, rdFMCS, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import PandasTools
import re
import argparse
import os
import sys
# from IPython.display import HTML

# this RDKit build has no Cairo backend (MolDraw2DCairo unavailable), so use
# SVG rendering instead of the PNG-via-Cairo default, for any dataframe with
# a molecule column rendered via PandasTools from here on
PandasTools.molRepresentation = "svg"
PandasTools.RenderImagesInAllDataFrames(images=True)
# %%

################################################################################
# Argparse arguments -----------------------------------------------------------
################################################################################
# %%
parser = argparse.ArgumentParser(
    description="Filter and report on candidate hits vs one or more reference compounds"
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
    "-c", "--chromatogram",
    help="Path to chromatogram SVGs parquet [default: <input>/report/features.parquet]"
)
parser.add_argument(
    "-t", "--similarity_cutoff",
    type=float,
    default=0.1,
    help="Minimum Tanimoto similarity to keep [default: %(default)s]"
)
parser.add_argument(
    "-m", "--mcs_cutoff",
    type=float,
    default=0,
    help="Minimum MCS fraction of reference to keep [default: %(default)s]"
)
parser.add_argument(
    "-o", "--output",
    default=None,
    help="Path to output folder [default: <input>/report/]"
)

args = parser.parse_args()

smiles_list = args.smiles.split(";")
ref_smiles = dict(x.split(",") for x in smiles_list)
references = list(ref_smiles.items())
# %%

################################################################################
# Setting variables ------------------------------------------------------------
################################################################################
# %%
output_dir = args.input
report_dir = f"{output_dir}/report"
tables_dir = f"{output_dir}/tables"
chromatogram_svgs_path = (
    args.chromatogram if args.chromatogram
    else f"{output_dir}/report/features.parquet"
)
report_dir = args.output if args.output else f"{output_dir}/report"
os.makedirs(report_dir, exist_ok=True)

similarity_cutoff = args.similarity_cutoff
mcs_cutoff = args.mcs_cutoff

chromatogram_svgs = pd.read_parquet(chromatogram_svgs_path)

# an experiment with zero features retained upstream (filter_intensity.R)
# produces a genuinely columnless 0x0 parquet here (not just 0 rows) - treat
# that as "nothing to report" rather than crashing on a missing "feature" col
if chromatogram_svgs.empty or "feature" not in chromatogram_svgs.columns:
    print(f"No retained features for {output_dir} - nothing to report, skipping.")
    sys.exit(0)

feature_levels = chromatogram_svgs["feature"]
annotation_similarities = pd.read_csv(
    f"{report_dir}/annotation_similarities.csv",
    on_bad_lines = "warn"
)
biotransformer_similarities = pd.read_csv(
    f"{report_dir}/biotransformer_similarities.csv"
)
# %%

################################################################################
# Old interactive --------------------------------------------------------------
################################################################################
# %%
# html_export_path = "/Users/wilhelm/Desktop/similar_hits.html"
# output_dir = (
#     "/Volumes/bluecub/aglycone_release_100um_24h/output/"
#     "afzelin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon"
# )
# tables_dir = f"{output_dir}/tables"
# report_dir = f"{output_dir}/report"
# annotation_similarities = pd.read_csv(
#     f"{report_dir}/annotation_similarities.csv",
#     on_bad_lines = "warn"
# )
# biotransformer_similarities = pd.read_csv(
#     f"{report_dir}/biotransformer_similarities.csv"
# )
# # kaempferol, the afzelin aglycone
# reference_smiles = "C1=CC(=CC=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O"
# # Created a parquet with svg information
# chromatogram_svgs = pd.read_parquet("/Users/wilhelm/Desktop/features.parquet")
# similarity_cutoff = 0.2
# mcs_cutoff = 0.5
# %%

# %%
hits = pd.concat(
    [
        annotation_similarities[
            [
                "feature",
                "peak_id",
                "target_smiles",
                "sim",
                "sim_mol",
                "ppm_error",
                "mz",
                "rtime",
                "adduct"
            ]
        ].rename(
            columns={
                "feature": "feature_pred",
                "peak_id": "feature",
                "target_smiles": "smiles"
            }
        ).assign(source="mz_annotation"),
        biotransformer_similarities[
            [
                "feature",
                "SMILES",
                "sim",
                "sim_mol",
                "mz",
                "mzmed",
                "rtmed",
                "adduct"
            ]
        ].rename(
            columns={"SMILES": "smiles", "rtmed": "rtime"}
        ).assign(
            feature_pred=lambda df: df["feature"],
            # no direct mass-error column for biotransformer hits, so compute it
            # the same way (observed feature m/z vs. candidate's predicted m/z)
            ppm_error=lambda df: (df["mzmed"] - df["mz"]) / df["mz"] * 1e6,
            source="biotransformer",
        )
        # "mz" here was the candidate's own theoretical m/z (only needed for
        # the ppm_error calc above) - "mzmed" is the observed feature m/z,
        # same meaning as annotation_similarities' "mz" column, so rename to match
        .drop(columns=["mz"])
        .rename(columns={"mzmed": "mz"}),
    ],
    ignore_index=True,
).dropna(subset=["smiles"])

hits["mol"] = hits["smiles"].apply(Chem.MolFromSmiles)
hits = hits[hits["mol"].notna()].reset_index(drop=True)

# Different stereoisomers of the same skeleton (e.g. regiochemistry-preserving
# biotransformations applied at different positions) can independently match
# the same feature and show up as visually near-identical rows - the default
# Morgan fingerprint used for `sim` ignores stereochemistry too, so their
# similarity scores are literally identical. Collapse to one representative
# per feature/source/skeleton, using the InChIKey connectivity layer (first
# 14 chars, which by design excludes stereo/isotope/charge) as the "same
# skeleton" key.
hits["skeleton_key"] = hits["mol"].apply(lambda mol: inchi.MolToInchiKey(mol)[:14])
hits = hits.drop_duplicates(
    subset=["feature", "source", "sim_mol", "skeleton_key"]
).reset_index(drop=True)

# Restrict to features already flagged as interesting in
# sirius_add_to_pipeline.R (feature_levels there) - higher peak area in the
# glycoside-cultivated samples vs. glucose control, among other criteria.
# Written out from that script as a plain one-column CSV.
hits = hits[hits["feature"].isin(feature_levels)].reset_index(drop=True)

# pandas' to_html() runs a separate whitespace-preservation pass on string
# cell values (independent of escape=False, which only covers HTML-special-
# character escaping): real newlines become a literal "\n", and any run of
# 2+ consecutive spaces gets each replaced with a literal "&nbsp;". Both R
# SVGs embedded per reference below rely on this being fixed up first - see
# process_r_svg() for the full explanation.
chromatogram_svgs["chromatogram"] = [
    re.sub(r" {2,}", " ", svg.replace("\n", ""))
    for svg in chromatogram_svgs["chromatogram"]
]
# %%


################################################################################
# Per-reference analysis and report ---------------------------------------------
################################################################################
html_export_path = f"{report_dir}/similar_hits.html"


def process_reference(reference_smiles, ref_name):
    pcoa_export_path = f"{report_dir}/pcoa_{ref_name}.svg"

    # aglycone
    reference_mol = Chem.MolFromSmiles(reference_smiles)
    reference_num_atoms = reference_mol.GetNumAtoms()

    # hits holds one row per (candidate, reference) pair (from calc_similarity.py
    # scoring every candidate against every reference) - without the sim_mol
    # filter here, every call to process_reference() would see every reference's
    # rows pooled together, and would label rows scored against some *other*
    # reference with this call's reference_structure/reference_smiles instead.
    similar_hits = hits[
        (hits["sim"] >= similarity_cutoff) & (hits["sim_mol"] == ref_name)
    ].reset_index(drop=True)

    def mcs_result_vs_reference(mol):
        return rdFMCS.FindMCS([mol, reference_mol], timeout=10)

    def leftover_scaffold_smiles(mol, mcs_result):
        # Remove the atoms shared with the reference (the matched core) and
        # return the SMILES of the largest remaining piece - i.e. "what's
        # actually attached", instead of describing the whole molecule (which
        # is dominated by the shared core and hides what's actually different).
        if mcs_result.numAtoms == 0:
            return None
        core_pattern = Chem.MolFromSmarts(mcs_result.smartsString)
        if not mol.HasSubstructMatch(core_pattern):
            return None
        leftover = Chem.ReplaceCore(mol, core_pattern, labelByIndex=True)
        if leftover is None or leftover.GetNumAtoms() == 0:
            return None
        frags = Chem.GetMolFrags(leftover, asMols=True, sanitizeFrags=False)
        if not frags:
            return None
        largest_frag = max(frags, key=lambda f: f.GetNumAtoms())
        return Chem.MolToSmiles(largest_frag)

    def murcko_scaffold_smiles(mol):
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(scaffold)

    def murcko_generic_scaffold_smiles(mol):
        # Same Bemis-Murcko scaffold, but reduced to a bare carbon-skeleton
        # topology (all atoms -> C, all bonds -> single). This groups molecules
        # that share the same ring/linker connectivity even if they differ in
        # heteroatoms, aromaticity, or bond order - the exact scaffold above is
        # too strict for that and fragments genuinely-related structures apart.
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        generic_scaffold = MurckoScaffold.MakeScaffoldGeneric(scaffold)
        return Chem.MolToSmiles(generic_scaffold)

    reference_scaffold_smiles = murcko_scaffold_smiles(reference_mol)
    reference_generic_scaffold_smiles = murcko_generic_scaffold_smiles(
        reference_mol
    )

    # Mol objects (not just SMILES) of the reference's own scaffold, needed for
    # substructure containment checks below - a glycosylated/conjugated
    # candidate's *whole* scaffold will never equal kaempferol's bare scaffold
    # (it has extra rings from whatever's attached), but it can still *contain*
    # it. Equality was the bug; containment is the actual question.
    reference_scaffold_mol = MurckoScaffold.GetScaffoldForMol(reference_mol)
    reference_generic_scaffold_mol = MurckoScaffold.MakeScaffoldGeneric(
        reference_scaffold_mol
    )

    unique_smiles = similar_hits["smiles"].unique()
    unique_mols = {smi: Chem.MolFromSmiles(smi) for smi in unique_smiles}

    mcs_result_by_smiles = {smi: mcs_result_vs_reference(mol) for smi, mol in unique_mols.items()}
    mcs_by_smiles = {smi: result.numAtoms for smi, result in mcs_result_by_smiles.items()}
    scaffold_by_smiles = {smi: murcko_scaffold_smiles(mol) for smi, mol in unique_mols.items()}
    generic_scaffold_by_smiles = {
        smi: murcko_generic_scaffold_smiles(mol) for smi, mol in unique_mols.items()
    }
    leftover_by_smiles = {
        smi: leftover_scaffold_smiles(mol, mcs_result_by_smiles[smi])
        for smi, mol in unique_mols.items()
    }

    similar_hits["mcs_atoms_vs_reference"] = similar_hits["smiles"].map(mcs_by_smiles)
    similar_hits["mcs_fraction_of_reference"] = (
        similar_hits["mcs_atoms_vs_reference"] / reference_num_atoms
    )

    similar_hits["leftover_smiles"] = similar_hits["smiles"].map(leftover_by_smiles)

    similar_hits["scaffold_smiles"] = similar_hits["smiles"].map(scaffold_by_smiles)
    similar_hits["same_scaffold_as_reference"] = similar_hits["smiles"].map(
        lambda smi: unique_mols[smi].HasSubstructMatch(reference_scaffold_mol)
    )

    # Tanimoto similarity between each candidate's own scaffold and the
    # reference's scaffold specifically (not the whole molecule, that's `sim`) -
    # a graded measure rather than the boolean same_scaffold_as_reference.
    scaffold_morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    reference_scaffold_fp = scaffold_morgan_gen.GetFingerprint(reference_scaffold_mol)
    scaffold_similarity_by_smiles = {
        smi: DataStructs.TanimotoSimilarity(
            scaffold_morgan_gen.GetFingerprint(Chem.MolFromSmiles(scaffold_smi)),
            reference_scaffold_fp,
        )
        for smi, scaffold_smi in scaffold_by_smiles.items()
    }
    similar_hits["scaffold_similarity_to_reference"] = similar_hits["smiles"].map(
        scaffold_similarity_by_smiles
    )

    similar_hits["generic_scaffold_smiles"] = similar_hits["smiles"].map(generic_scaffold_by_smiles)
    similar_hits["same_generic_scaffold_as_reference"] = similar_hits["smiles"].map(
        lambda smi: unique_mols[smi].HasSubstructMatch(reference_generic_scaffold_mol)
    )

    # Compact group id for which scaffold family a hit belongs to (1 = most
    # common generic scaffold among these hits, 2 = next, etc.) - easier to
    # scan/group by than comparing the raw generic_scaffold_smiles strings.
    scaffold_group_by_smiles = (
        similar_hits["generic_scaffold_smiles"]
        .value_counts()
        .rank(method="first", ascending=False)
        .astype(int)
    )
    similar_hits["scaffold_group"] = similar_hits["generic_scaffold_smiles"].map(
        scaffold_group_by_smiles
    )

    # Conjugates: reference core substantially retained, but the candidate is
    # meaningfully bigger than what the MCS matched - i.e. something else is
    # attached, rather than just a small in-scaffold decoration (hydroxylation
    # etc). That's the opposite tail from same_scaffold_as_reference.
    similar_hits["candidate_num_atoms"] = similar_hits["smiles"].map(
        lambda smi: unique_mols[smi].GetNumAtoms()
    )
    similar_hits["extra_atoms_vs_reference"] = (
        similar_hits["candidate_num_atoms"] - similar_hits["mcs_atoms_vs_reference"]
    )

    def leftover_profile(leftover_smi):
        # Coarse composition of the attached piece: ring count + which
        # heteroatoms it contains. Two attachments with different exact
        # connectivity (e.g. one sugar ring vs. two linked sugar rings) still
        # get treated as "the same kind of thing" if this profile matches -
        # a genuinely more general grouping than exact scaffold topology, which
        # splits on ring *count* even when the underlying chemistry (a sugar-like
        # conjugate) is the same idea.
        if leftover_smi is None:
            return (0, False, False, False)
        frag = Chem.MolFromSmiles(leftover_smi)
        if frag is None:
            return (0, False, False, False)
        num_rings = rdMolDescriptors.CalcNumRings(frag)
        symbols = [atom.GetSymbol() for atom in frag.GetAtoms()]
        return (num_rings, "O" in symbols, "N" in symbols, "S" in symbols)

    leftover_profile_by_smiles = {
        smi: leftover_profile(leftover_by_smiles[smi]) for smi in unique_smiles
    }

    similar_hits["leftover_num_rings"] = similar_hits["smiles"].map(
        lambda smi: leftover_profile_by_smiles[smi][0]
    )
    similar_hits["leftover_has_o"] = similar_hits["smiles"].map(
        lambda smi: leftover_profile_by_smiles[smi][1]
    )
    similar_hits["leftover_has_n"] = similar_hits["smiles"].map(
        lambda smi: leftover_profile_by_smiles[smi][2]
    )
    similar_hits["leftover_has_s"] = similar_hits["smiles"].map(
        lambda smi: leftover_profile_by_smiles[smi][3]
    )

    def attachment_type(row):
        if row["leftover_num_rings"] == 0 and row["extra_atoms_vs_reference"] == 0:
            return "none"
        heteroatoms = "".join(
            symbol
            for symbol, present in (
                ("N", row["leftover_has_n"]),
                ("S", row["leftover_has_s"]),
                ("O", row["leftover_has_o"]),
            )
            if present
        ) or "none"
        return f"rings={row['leftover_num_rings']}_hetero={heteroatoms}"

    similar_hits["attachment_type"] = similar_hits.apply(attachment_type, axis=1)
    similar_hits["abs_ppm_error"] = similar_hits["ppm_error"].abs()

    # Filter with cutoff for MCS
    mcs_survivors = similar_hits[
        similar_hits["mcs_fraction_of_reference"] >= mcs_cutoff
    ].reset_index(drop=True)

    # Info on counts (grouping by the *generic* scaffold - see note above on why
    # the exact scaffold over-fragments closely-related structures)
    scaffold_counts = mcs_survivors["generic_scaffold_smiles"].value_counts()
    print(
        f"[{ref_name}] {len(scaffold_counts)} distinct scaffolds among",
        f"{len(mcs_survivors)} MCS survivors"
    )

    # Group by what's actually *attached* (leftover_smiles) rather than by the
    # whole-molecule scaffold, which is dominated by the shared core. Recurring
    # leftover fragments across different features/candidates are the strongest
    # "same kind of conjugate showing up more than once" signal so far.
    leftover_counts = (
        mcs_survivors.dropna(subset=["leftover_smiles"])["leftover_smiles"].value_counts()
    )
    print(
        f"[{ref_name}] {len(leftover_counts)} distinct attached fragments among",
        f"{leftover_counts.sum()} candidates"
    )

    # PCoA of every individual candidate assignment, vs. the reference itself -
    # shows how the actual proposed assignments cluster in chemical space
    # relative to the reference, using every hit rather than a reduced/grouped set.
    morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

    def pcoa(distance_matrix, n_components=2):
        n = distance_matrix.shape[0]
        d2 = distance_matrix**2
        centering = np.eye(n) - np.ones((n, n)) / n
        b = -0.5 * centering @ d2 @ centering

        eigenvalues, eigenvectors = np.linalg.eigh(b)
        order = np.argsort(eigenvalues)[::-1]
        eigenvalues = eigenvalues[order]
        eigenvectors = eigenvectors[:, order]

        coords = eigenvectors[:, :n_components] * np.sqrt(
            np.clip(eigenvalues[:n_components], 0, None)
        )
        explained_variance_ratio = (
            eigenvalues[:n_components] / eigenvalues[eigenvalues > 0].sum()
        )
        return coords, explained_variance_ratio

    hit_fps = [morgan_gen.GetFingerprint(mol) for mol in similar_hits["mol"]]
    reference_fp = morgan_gen.GetFingerprint(reference_mol)
    all_fps = hit_fps + [reference_fp]

    all_tanimoto_matrix = np.array(
        [DataStructs.BulkTanimotoSimilarity(fp, all_fps) for fp in all_fps]
    )
    all_distance_matrix = 1 - all_tanimoto_matrix

    all_pcoa_coords, all_pcoa_explained = pcoa(all_distance_matrix)
    hit_coords = all_pcoa_coords[:-1]
    reference_coords = all_pcoa_coords[-1]

    source_colors = {"mz_annotation": "#0072B2", "biotransformer": "#D55E00"}

    fig, ax = plt.subplots(figsize=(7, 6))
    for source, idx in similar_hits.groupby("source").indices.items():
        ax.scatter(
            hit_coords[idx, 0],
            hit_coords[idx, 1],
            s=15,
            color=source_colors[source],
            label=source,
            edgecolor="none",
        )

    ax.scatter(
        [reference_coords[0]],
        [reference_coords[1]],
        s=250,
        color="black",
        marker="*",
        label="reference compound",
        zorder=5,
    )

    ax.set_xlabel(f"PCo1 ({all_pcoa_explained[0]:.1%})")
    ax.set_ylabel(f"PCo2 ({all_pcoa_explained[1]:.1%})")
    ax.set_title(f"PCoA of individual candidate assignments vs {ref_name}")
    ax.legend(title="Source")
    fig.tight_layout()
    fig.savefig(pcoa_export_path, format="svg")
    plt.close(fig)

    # Resorting columns
    similar_hits["reference_name"] = ref_name
    similar_hits["reference_smiles"] = reference_smiles
    PandasTools.AddMoleculeColumnToFrame(
        similar_hits, smilesCol="reference_smiles", molCol="reference_structure"
    )
    PandasTools.AddMoleculeColumnToFrame(
        similar_hits, smilesCol="scaffold_smiles", molCol="scaffold_structure"
    )

    cols = list(similar_hits.columns)
    for col in ["reference_structure", "scaffold_structure", "reference_name", "reference_smiles", "scaffold_smiles", "mz"]:
        cols.remove(col)

    mol_idx = cols.index("mol")
    cols[mol_idx:mol_idx] = ["reference_structure", "scaffold_structure"]

    feature_idx = cols.index("feature") + 1
    cols[feature_idx:feature_idx] = ["reference_name", "reference_smiles", "scaffold_smiles"]

    rtime_idx = cols.index("rtime")
    cols[rtime_idx:rtime_idx] = ["mz"]

    similar_hits = similar_hits[cols]

    # Chromatogram SVGs for every feature (not just one), exported from
    # check_hit_xcms.R into chromatogram_svgs (columns: feature, chromatogram -
    # the raw SVG text, already whitespace-fixed up above). Merge onto
    # similar_hits by feature.
    _, mol_height = PandasTools.molSize

    def process_r_svg(svg_content, suffix):
        vb_width, vb_height = (
            float(x)
            for x in re.search(
                r"viewBox='0 0 ([\d.]+) ([\d.]+)'", svg_content
            ).groups()
        )
        target_width = round(vb_width * mol_height / vb_height)

        # only the root <svg> element's own width/height (the first
        # occurrence) - inner elements (clip rects etc.) must stay untouched
        svg_content = re.sub(
            r"width='[^']*'", f"width='{target_width}px'", svg_content, count=1
        )
        svg_content = re.sub(
            r"height='[^']*'", f"height='{mol_height}px'", svg_content, count=1
        )

        # explicit white background, drawn first/behind everything else -
        # these SVGs have no background rect of their own (RDKit's own SVGs
        # always draw one), so "transparent" areas show through to whatever's
        # behind them otherwise
        background_rect = (
            f"<rect width='{vb_width}' height='{vb_height}' style='fill:#FFFFFF;'/>"
        )
        svg_content = re.sub(
            r"(<svg[^>]*>)", rf"\1{background_rect}", svg_content, count=1
        )

        # clipPath ids are deterministic per plot layout, so different rows
        # embedding different features' plots can still collide if two plots
        # share the same panel geometry - suffix every id (and its references,
        # since they share the same string) with something unique per row.
        # Prefixed with ref_name too, since the same feature's chromatogram is
        # now embedded once per reference report.
        ids = set(re.findall(r"id='([^']+)'", svg_content))
        for id_value in ids:
            svg_content = svg_content.replace(id_value, f"{id_value}_{ref_name}_{suffix}")

        return svg_content

    ref_chromatogram_svgs = chromatogram_svgs.assign(
        chromatogram=[
            process_r_svg(svg, suffix=i)
            for i, svg in enumerate(chromatogram_svgs["chromatogram"])
        ]
    )

    similar_hits = similar_hits.merge(
        ref_chromatogram_svgs[["feature", "chromatogram"]], on="feature", how="left"
    )

    cols = list(similar_hits.columns)
    cols.remove("chromatogram")
    cols.insert(cols.index("mol") + 1, "chromatogram")
    similar_hits = similar_hits[cols]

    return similar_hits


similar_hits = pd.concat(
    [process_reference(reference_smiles, ref_name) for ref_name, reference_smiles in references],
    ignore_index=True,
)

# Sort: 1. overall similarity (sim, highest first), 2. feature, 3. scaffold
# similarity to the reference (highest first). Rows for different reference
# compounds are interleaved by this same ordering - reference_name is not a
# sort key, just a column to filter/scan by.
similar_hits = similar_hits.sort_values(
    ["sim", "feature", "scaffold_similarity_to_reference"],
    ascending=[False, True, False],
).reset_index(drop=True)

# Standalone HTML export - everything here is self-contained (RDKit SVGs +
# base64-embedded PNG), so it opens in any browser outside VS Code, no
# Python/RDKit needed to view it.
# Plain (non f-) string so all the JS's own {}'s don't need escaping.
toggle_script = """
document.addEventListener("DOMContentLoaded", function () {
  var table = document.querySelector("table");
  var headerCells = table.querySelectorAll("thead tr th");
  var allColumns = [];
  headerCells.forEach(function (th, idx) {
    // pandas adds a blank leading <th> for the row index - skip it, nothing
    // useful to toggle there
    if (th.textContent.trim() !== "") {
      allColumns.push({ idx: idx, name: th.textContent });
    }
  });

  var toolbar = document.getElementById("column-toggles");
  allColumns.forEach(function (col) {
    var label = document.createElement("label");
    label.style.marginRight = "12px";
    var checkbox = document.createElement("input");
    checkbox.type = "checkbox";
    checkbox.checked = true;
    checkbox.addEventListener("change", function () {
      var display = checkbox.checked ? "" : "none";
      table.querySelectorAll("tr").forEach(function (row) {
        var cell = row.children[col.idx];
        if (cell) cell.style.display = display;
      });
    });
    label.appendChild(checkbox);
    label.appendChild(document.createTextNode(" " + col.name));
    toolbar.appendChild(label);
  });
});
"""

with open(html_export_path, "w") as f:
    f.write(
        f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>similar_hits</title>
<style>
body {{ font-family: sans-serif; }}
table {{ border-collapse: collapse; }}
th, td {{ border: 1px solid #ccc; padding: 4px; text-align: center; vertical-align: middle; }}
th {{ position: sticky; top: 0; background: #f5f5f5; }}
#column-toggles {{ margin-bottom: 10px; padding: 8px; background: #f0f0f0; border-radius: 4px; }}
</style>
</head>
<body>
<div id="column-toggles"><strong>Toggle columns:</strong></div>
{similar_hits.to_html(escape=False)}
<script>
{toggle_script}
</script>
</body>
</html>
"""
    )

print(f"Saved to {html_export_path}")
