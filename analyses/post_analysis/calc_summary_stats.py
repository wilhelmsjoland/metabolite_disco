# %%
chart_path = Path(
    "/Volumes/bluecub/aglycone_release_100um_24h/combined_results/hit_plots"
)

data_with_info = (
    data
    .join(
        other=exp_info,
        on="experiment",
        how="left"
    )
    .filter(
        pl.col("source") == "biotransformer",
        # pl.col("adduct") == "[M-H]-",
    )
)

aglycone_exps = [
    "avicularin",
    "quercetin_3_o_b_d_glucose_7_o_b_d_gentiobioside",
    "prunin",
    "naringenin_7_o_b_d_glucuronide",
    "luteolin_7_o_rutinoside",
    "luteolin_7_o_glucuronide",
    "clitorin",
    "kaempferol_7_neohesperidoside",
    "kaempferol_3_glucuronide",
    "kaempferol_7_o_glucoside",
    "isorhamnetin_3_o_neohesperidoside",
    # don't have hesperetin-7-o-glucoside
    "sophoricoside",
    # nothing from apiin
    "isorhoifolin"
]

aglycone_exps = [
    "avicularin_b_uniformis_atcc_8492_and_bu_gsh_d_ggh_c_gsh_g",
    "quercetin_3_o_b_d_glucose_7_o_b_d_gentiobioside_b_uniformis_atcc_8492_and_bu_gsh_d_ggh_c_gsh_g",
    "prunin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon",
    "prunin_b_thetaiotaomicron_vpi_5482",
    "prunin_b_uniformis_atcc_8492_and_bu_gsh_d_ggh_c_gsh_g",
    "prunin_p_copri_i_ak263",
    "naringenin_7_o_b_d_glucuronide_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon",
    "naringenin_7_o_b_d_glucuronide_b_thetaiotaomicron_vpi_5482",
    "naringenin_7_o_b_d_glucuronide_b_uniformis_atcc_8492_and_bu_gsh_d_ggh_c_gsh_g",
    "luteolin_7_o_rutinoside_b_uniformis_atcc_8492_and_bu_gsh_d_ggh_c_gsh_g",
    "luteolin_7_o_glucuronide_b_uniformis_atcc_8492_and_bu_gsh_d_ggh_c_gsh_g",
    "luteolin_7_o_glucuronide_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon",
    "luteolin_7_o_glucuronide_b_thetaiotaomicron_vpi_5482",
    "clitorin_b_uniformis_atcc_8492_and_bu_gsh_d_ggh_c_gsh_g",
    "clitorin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon",
    "clitorin_b_thetaiotaomicron_vpi_5482",
    "kaempferol_7_neohesperidoside_b_uniformis_atcc_8492_and_bu_gsh_d_ggh_c_gsh_g",
    "kaempferol_7_neohesperidoside_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon",
    "kaempferol_7_neohesperidoside_b_thetaiotaomicron_vpi_5482",
    "kaempferol_3_glucuronide_b_uniformis_atcc_8492_and_bu_gsh_d_ggh_c_gsh_g",
    "kaempferol_3_glucuronide_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon",
    "kaempferol_3_glucuronide_b_thetaiotaomicron_vpi_5482",
    "kaempferol_7_o_glucoside_b_uniformis_atcc_8492_and_bu_gsh_d_ggh_c_gsh_g",
    "kaempferol_7_o_glucoside_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon",
    "kaempferol_7_o_glucoside_b_thetaiotaomicron_vpi_5482",
    "kaempferol_7_o_glucoside_p_copri_i_ak263",
    "isorhamnetin_3_o_neohesperidoside_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon",
    "hesperetin_7_o_glucoside_b_uniformis_atcc_8492_and_bu_gsh_d_ggh_c_gsh_g",
    "hesperetin_7_o_glucoside_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon",
    "hesperetin_7_o_glucoside_b_thetaiotaomicron_vpi_5482",
    "hesperetin_7_o_glucoside_p_copri_i_ak263",
    "sophoricoside_b_uniformis_atcc_8492_and_bu_gsh_d_ggh_c_gsh_g",
    "sophoricoside_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon",
    "sophoricoside_b_thetaiotaomicron_vpi_5482",
    "sophoricoside_p_copri_i_ak263",
    "apiin_b_uniformis_atcc_8492_and_bu_gsh_d_ggh_c_gsh_g",
    "apigenin_7_o_2g_rhamnosyl_gentiobioside_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon",
    "apigenin_7_o_2g_rhamnosyl_gentiobioside_b_thetaiotaomicron_vpi_5482",
    "isorhoifolin_b_uniformis_atcc_8492_and_bu_gsh_d_ggh_c_gsh_g",
    "isorhoifolin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon",
    "isorhoifolin_b_thetaiotaomicron_vpi_5482",
    "isorhoifolin_p_copri_i_ak263",
]

# %%
unique_feats_per_exp = (
    data_with_info
    .group_by("experiment")
    .agg(
        pl.col("feature").n_unique().alias("feats"),
        pl.col("condition").first(),
        pl.col("strain").first()
    )
    .filter(pl.col("experiment").is_in(aglycone_exps))
)

# Features per experiment
all_hits_points = alt.Chart(unique_feats_per_exp).mark_circle(size = 100).encode(
    x=alt.X(
        "condition:N",
        title="Supplemented glycoside",
        axis=alt.Axis(
            labelLimit=300,
            # titlePadding=50,
            labelAngle=45
        ),
        sort=alt.EncodingSortField(
            field="feats",
            op="max",
            order="descending",
        )
    ),
    y=alt.Y(
        "feats:Q",
        title="Features per experiment"
    ),
    fill=alt.Fill(
        "strain:N",
        legend=alt.Legend(
            labelLimit=400,
            title="Strains"
        )
    ),
).properties(padding={"bottom": 150})

all_hits_points.save(fp=chart_path.joinpath("all_hits.pdf"), format="pdf")

# %%
# CONDITION
all_exps = unique_feats_per_exp["experiment"].unique().sort().to_list()

for i in all_exps[:]:
        tink = (
            data_with_info
            .filter(
                pl.col("experiment") == i,
                # pl.col("sim") > 0.1
            )
            .with_columns(
                pl.col("skeleton_key").n_unique().over("feature").alias("total_hits")
            )
            .unique(subset=["feature", "skeleton_key"])
            .with_columns(
                pl.col("smiles")
                .map_elements(smiles_to_data_uri, return_dtype=pl.String)
                .alias("structure")
            )
            .with_columns(
                (
                    pl.col("feature")
                    + "\n(n="
                    + pl.col("total_hits").cast(pl.String)
                    + ", mz="
                    + pl.col("mz").round(2).cast(pl.String)
                    + ", rt="
                    + pl.col("rtime").round(1).cast(pl.String)
                    + ")"
                ).alias("feature_label")
            )
            .with_columns(
                (pl.col("sim").rank(method="ordinal", descending=True).over("feature") == 1).alias("is_top_hit")
            )
        )

        circle = alt.Chart(tink).mark_circle(size = 50).encode(
            x=alt.X(
                "sim:Q",
                title = "Tanimoto similarity",
                sort=alt.EncodingSortField(field="sim", order="descending"),
                scale=alt.Scale(domain = [0, 1], padding=0)
            ),
            y=alt.Y(
                "feature_label:N",
                title=None,
                scale=alt.Scale(padding=0),
                axis=alt.Axis(
                    labelExpr="split(datum.label, '\\n')",
                    labelAlign="center",
                    labelPadding=70,
                ),
            ),
            color=alt.condition(
                alt.datum.is_top_hit,
                alt.value("crimson"),
                alt.value("lightgray"),
            ),
            tooltip=["skeleton_key", "smiles"],
        ).properties(width = 80, height=alt.Step(40))

        structures = circle.mark_image(width=60, height=60).transform_filter(
            alt.datum.is_top_hit
        ).encode(
            url="structure:N",
            y=alt.Y(
                "feature_label:N",
                sort=alt.EncodingSortField(field="sim", order="descending"),
                axis=None,
                scale=alt.Scale(padding=0)
            ),
            x=alt.value(20),
        ).properties(width=60, height=alt.Step(40))

        full_plot = (
            alt.hconcat(circle, structures, spacing=5)
            .resolve_scale(y="shared")
            .configure_view(strokeWidth=0)
        )

        full_plot.save(fp=chart_path.joinpath(f"{i}.pdf"), format="pdf")

# %%

# %%
top_skeletons_base = (
    data_with_info
    .filter(
        pl.col("experiment").is_in(aglycone_exps)
        # pl.col("sim") > 0.1
    )
    .group_by("skeleton_key", "strain")
    .agg(
        pl.col("experiment").n_unique().alias("n_experiments"),
        pl.col("smiles").first().alias("smiles"),
    )
)

top_skeleton_keys = (
    top_skeletons_base
    .group_by("skeleton_key")
    .agg(pl.col("n_experiments").max().alias("best_strain_n"))
    .sort("best_strain_n", descending=True)
    .head(12)
    .get_column("skeleton_key")
)

top_skeletons = (
    top_skeletons_base
    .filter(pl.col("skeleton_key").is_in(top_skeleton_keys))
    .with_columns(
        pl.col("smiles")
        .map_elements(smiles_to_data_uri, return_dtype=pl.String)
        .alias("structure")
    )
)

comb_circle = alt.Chart(top_skeletons).mark_circle(size = 100).encode(
    x=alt.X(
        "skeleton_key:N",
        sort=alt.EncodingSortField(field="n_experiments", order="descending"),
        axis=alt.Axis(title="Unique metabolite")
        ),
    y=alt.Y(
        "n_experiments:Q",
        title="Hits across experiments"
    ),
    color=alt.Color(
        "strain:N"
    ),
    tooltip=["skeleton_key", "smiles", "n_experiments"],
).properties(height = 150, width=alt.Step(60))

comb_struct = comb_circle.mark_image(width=60, height=60).encode(
    url="structure:N",
    x=alt.X(
        "skeleton_key:N",
        sort=alt.EncodingSortField(field="n_experiments", order="descending"),
        axis=None,
    ),
    y=alt.value(20),
).properties(height=60, width=alt.Step(60))

full_n_plot = (
    alt.vconcat(comb_struct, comb_circle, spacing=6)
    .resolve_scale(x="shared")
    .configure_view(strokeWidth=0)
)

full_n_plot.save(fp=chart_path.joinpath(f"combined_hits_strain.pdf"), format="pdf")
# %%

# Combined_all
# %%
combined_hits = (
    data_with_info
    .filter(
        pl.col("experiment").is_in(aglycone_exps)
        # pl.col("sim") > 0.1
    )
    .group_by("skeleton_key")
    .agg(
        pl.col("experiment").n_unique().alias("n_experiments"),
        pl.col("smiles").first().alias("smiles"),
    )
    .sort("n_experiments", descending=True)
    .head(12)
    .with_columns(
        pl.col("smiles")
        .map_elements(smiles_to_data_uri, return_dtype=pl.String)
        .alias("structure")
    )
)

comb_circle = alt.Chart(combined_hits).mark_circle(size = 100).encode(
    x=alt.X(
        "skeleton_key:N",
        sort=alt.EncodingSortField(field="n_experiments", order="descending"),
        axis=alt.Axis(title="Unique metabolite")
        ),
    y=alt.Y(
        "n_experiments:Q",
        title="Hits across experiments"
    ),
    tooltip=["skeleton_key", "smiles", "n_experiments"],
).properties(height = 150, width=alt.Step(60))

comb_struct = comb_circle.mark_image(width=60, height=60).encode(
    url="structure:N",
    x=alt.X(
        "skeleton_key:N",
        sort=alt.EncodingSortField(field="n_experiments", order="descending"),
        axis=None,
    ),
    y=alt.value(20),
).properties(height=60, width=alt.Step(60))

full_n_plot = (
    alt.vconcat(comb_struct, comb_circle, spacing=6)
    .resolve_scale(x="shared")
    .configure_view(strokeWidth=0)
)

full_n_plot.save(fp=chart_path.joinpath(f"combined_hits.pdf"), format="pdf")


# %%
n_feats = []
for i in all_exps:
    report_path = Path(
        f"/Volumes/bluecub/aglycone_release_100um_24h/output/"
        f"{i}/report/features.parquet"
    )
    n_feat = pl.read_parquet(report_path)["feature"].unique().len()
    n_feats.append({"experiment": i, "n_feat": n_feat})

n_good_features = pl.DataFrame(
    n_feats,
    schema={
        "experiment": pl.Utf8,
        "n_feat": pl.Int64
    }
)

ratio_annotated_features = (
    unique_feats_per_exp
    .join(
        other=n_good_features,
        on="experiment",
        how="left"
    )
    .with_columns(
        (pl.col("feats") / pl.col("n_feat")*100).alias("ratio_annot")
    )
)

strains = ratio_annotated_features["strain"].unique().to_list()
strain_pattern = "|".join(re.escape(f"_{s}") for s in strains)

annotated_plot = alt.Chart(ratio_annotated_features).mark_bar().encode(
    x=alt.X(
        "experiment:N",
        sort=alt.EncodingSortField(
            field="ratio_annot",
            order="descending",
        ),
        axis=alt.Axis(
            labelLimit=300,
            title="Experiment",
            labelAngle=45,
            labelExpr=f"replace(datum.value, regexp('{strain_pattern}', 'g'), '')"
        )
    ),
    y=alt.Y(
        "ratio_annot:Q",
        axis=alt.Axis(title="Proportion annotated features (%)")
    ),
    color=alt.Color(
        "strain:N",
        legend=alt.Legend(
            labelLimit=300,
            title="Strain"
        )
    )
)

annotated_plot.save(fp=chart_path.joinpath(f"annotated.pdf"), format="pdf")

# %%
svg_exp="hesperetin_7_o_glucoside_b_uniformis_atcc_8492_and_bu_gsh_d_ggh_c_gsh_g"
svg_path=Path(
    f"/Volumes/bluecub/aglycone_release_100um_24h/output/"
    f"{svg_exp}/report/features.parquet"
)
svgs = pl.read_parquet(svg_path)
svg_string = svgs.filter(pl.col("feature") == "FT07629")["chromatogram"].item()
pdf_path = chart_path.joinpath("chromatogram.pdf")

subprocess.run(
    ["rsvg-convert", "-f", "pdf", "-o", str(pdf_path)],
    input=svg_string.encode(),
    check=True
)

# %%
svgs.filter(pl.col("feature") == "FT05181")["chromatogram"].item()
# %%
