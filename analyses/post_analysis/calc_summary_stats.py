# %%
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
    "luteolin_7_o_rutinoside_b_uniformis_atcc_8492_and_bu_gsh_d_ggh_c_gsh_g"
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
    "hesperetin_7_o_glucoside_b_uniformis_atcc_8492_and_bu_gsh_d_ggh_c_gsh_g",
    "hesperetin_7_o_glucoside_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon",
    "hesperetin_7_o_glucoside_b_thetaiotaomicron_vpi_5482",
    "hesperetin_7_o_glucoside_p_copri_i_ak263",
    "sophoricoside_b_uniformis_atcc_8492_and_bu_gsh_d_ggh_c_gsh_g",
    "sophoricoside_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon",
    "sophoricoside_b_thetaiotaomicron_vpi_5482",
    "sophoricoside_p_copri_i_ak263",
    "apiin_b_uniformis_atcc_8492_and_bu_gsh_d_ggh_c_gsh_g",
    "apigenin_7_o_2g_rhamnosyl_gentiobioside_b_thetaiotaomicron_vpi_5482",
    "isorhoifolin_b_uniformis_atcc_8492_and_bu_gsh_d_ggh_c_gsh_g",
    "isorhoifolin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon",
    "isorhoifolin_b_thetaiotaomicron_vpi_5482",
    "isorhoifolin_p_copri_i_ak263"
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
        # sort=alt.EncodingSortField(
        #     field="feats",
        #     order="descending"
        # )
    ),
    y=alt.Y("feats:Q"),
    fill=alt.Fill("strain:N"),
)



# %%
# CONDITION
all_exps = unique_feats_per_exp["experiment"].unique().sort().to_list()
chart_path = Path(
    "/Volumes/bluecub/aglycone_release_100um_24h/combined_results/hit_plots"
)
all_hits_points.save(fp=chart_path.joinpath("all_hits.pdf"), format="pdf")

# FIX THIS SO IT DOESNT PLOT EVERY STRAIN FOR EVERY CONDITION
# FIX SO IT ONLY PLOTS THE ONES IN THE aglycone_exps list!!!!!!!

for i in all_exps[:]:
    strains_in_condition = (
        data_with_info
        .filter(pl.col("condition") == i)
        .get_column("strain")
        .unique()
        .sort()
        .to_list()
    )

    for s in strains_in_condition[:]:
        tink = (
            data_with_info
            .filter(
                pl.col("condition") == i,
                pl.col("strain") == s,
                # pl.col("sim") > 0.1
            )
            .with_columns(
                pl.col("skeleton_key").n_unique().over("feature").alias("total_hits")
            )
            # .filter(
            #     pl.col("sim").rank(method="ordinal", descending=True).over("feature") <= 5
            # )
            .unique(subset=["feature", "skeleton_key"])
        )

        tink_img = tink.with_columns(
            pl.col("smiles")
            .map_elements(smiles_to_data_uri, return_dtype=pl.String)
            .alias("structure")
        )

        tink_img = tink_img.with_columns(
            (
                pl.col("feature")
                + " (n="
                + pl.col("total_hits").cast(pl.String)
                + ", mz="
                + pl.col("mz").round(2).cast(pl.String)
                + ", rt="
                + pl.col("rtime").round(1).cast(pl.String)
                + ")"
            ).alias("feature_label")
        )

        tink_img = tink_img.with_columns(
            (pl.col("sim").rank(method="ordinal", descending=True).over("feature") == 1).alias("is_top_hit")
        )

        circle = alt.Chart(tink_img).mark_circle(size = 100).encode(
            x=alt.X(
                "sim:Q",
                sort=alt.EncodingSortField(field="sim", order="descending"),
                scale=alt.Scale(domain = [0, 1])
            ),
            y=alt.Y("feature_label:N", title=None),
            color=alt.condition(
                alt.datum.is_top_hit,
                alt.value("crimson"),
                alt.value("lightgray"),
            ),
            tooltip=["skeleton_key", "smiles"],
        ).properties(height=alt.Step(50))

        structures = circle.mark_image(width=50, height=50).transform_filter(
            alt.datum.is_top_hit
        ).encode(
            url="structure:N",
            y=alt.Y(
                "feature_label:N",
                sort=alt.EncodingSortField(field="sim", order="descending"),
                axis=None,
            ),
            x=alt.value(20),
        ).properties(width=50, height=alt.Step(50))


        full_plot = (
            alt.hconcat(circle, structures, spacing=0)
            .resolve_scale(y="shared")
            .configure_view(strokeWidth=0)
        )

        full_plot.save(fp=chart_path.joinpath(f"{i}_{s}.pdf"), format="pdf")

# %%
full_plot

# %%
tink = (
    data_with_info
    .filter(
        pl.col("experiment").is_in(aglycone_exps)
        # pl.col("sim") > 0.1
    )
    .with_columns(
        pl.col("skeleton_key").n_unique().over("feature").alias("total_hits"),
        pl.col("experiment").n_unique().over("skeleton_key").alias("n_experiments"),
    )
    .filter(
        pl.col("sim").rank(method="ordinal", descending=True).over("feature") <= 5
    )
    .unique(subset=["feature", "skeleton_key"])
    .filter(pl.col("n_experiments") > 5)
)

tink_img = tink.with_columns(
    pl.col("smiles")
    .map_elements(smiles_to_data_uri, return_dtype=pl.String)
    .alias("structure")
)

tink_img = tink_img.with_columns(
    (pl.col("feature") + " (" + pl.col("total_hits").cast(pl.String) + ")").alias("feature_label")
)

circle = alt.Chart(tink_img).mark_circle(size = 100).encode(
    x=alt.X(
        "skeleton_key:N",
        sort=alt.EncodingSortField(field="n_experiments", order="descending"),
        axis=None
        ),
    y=alt.Y("n_experiments:Q", title=None),
    tooltip=["skeleton_key", "smiles", "n_experiments"],
)

structures = circle.mark_image(width=50, height=50).encode(
    url="structure:N",
    x=alt.X(
        "skeleton_key:N",
        sort=alt.EncodingSortField(field="n_experiments", order="descending"),
        axis=None,
    ),
    y=alt.value(20),
).properties(height=50, width=alt.Step(50))

full_n_plot = (
    alt.vconcat(circle, structures, spacing=0)
    .resolve_scale(x="shared")
    .configure_view(strokeWidth=0)
)

full_n_plot.save(fp=chart_path.joinpath(f"combined_hits.pdf"), format="pdf")

# %%
(
    tink
    .filter(pl.col("skeleton_key") == "VYFYYTLLBUKUHU")
)

# %%
