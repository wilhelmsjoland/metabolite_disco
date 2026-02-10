cd /mnt/c/Users/wilhelm/Documents/MEGA/01_juniper/01_arbete/01_projekt/03_psm/predict_metabolism/biotransformer3.0jar

java -jar BioTransformer3.0_20230525.jar \
	-k pred \
	-b superbio \
	-isdf ../apigenin.sdf \
	-ocsv ../apigenin.csv \
	-a
