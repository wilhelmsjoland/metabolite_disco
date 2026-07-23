cd /Users/wilhelm/Proton/01_juniper/01_arbete/01_projekt/03_psm

source "$(conda info --base)/etc/profile.d/conda.sh"

# Rscript script_tools/filter_intensity.R \
#	-i /Volumes/bluecub/aglycone_release_100um_24h/output/afzelin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon \
#	-s 0.2 \
#	-f 5
#
#Rscript script_tools/chromatogram_plots.R \
#	-i /Volumes/bluecub/aglycone_release_100um_24h/output/afzelin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon \
#	-f /Volumes/bluecub/aglycone_release_100um_24h/output/afzelin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon/report/retained_features.csv \
#
conda run -n psm_chem python script_tools/create_report.py \
#	-i /Volumes/bluecub/aglycone_release_100um_24h/output/afzelin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon \
#	-s "C1=CC(=CC=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O" \
#	-c /Volumes/bluecub/aglycone_release_100um_24h/output/afzelin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon/report/features.parquet \
#	-t 0.2 \
#	-m 0
	
