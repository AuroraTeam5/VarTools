#create temporary file 
TMPFILE=`mktemp -p /tmp`

#use tabix tool to determine specific chromosome. | use bcftools norm to split multiallelic variants in separate rows.
tabix -h $1 $2 | sed 's/##INFO=<ID=AC_Het,Number=A/##INFO=<ID=AC_Het,Number=./g'| bcftools norm -m - -o $TMPFILE

#making fifo to apply piping 
mkfifo pipe pipe1 pipe2 pipe4 pipe5 pipe6
chmod 666 pipe pipe1 pipe2 pipe4 pipe5 pipe6


mysql -u root -p201666 -e "SET FOREIGN_KEY_CHECKS = 0;"

#determine columns of variant table and load this data into mysql database.
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\n' $TMPFILE > pipe1 & mysql -u root -p201666 -e "LOAD DATA LOCAL INFILE 'pipe1' INTO TABLE DataSet4_db.variants_variant FIELDS TERMINATED BY '\t' (Chromosome_Number, Position, Var_ID, Reference, Alternate, Quality, Filter) SET id = NULL;"

#copy data of population table and annotation table in 2 pipes.
bcftools query -f '%INFO/AC\t%INFO/AC_AFR\t%INFO/AC_AMR\t%INFO/AC_Adj\t%INFO/AC_EAS\t%INFO/AC_FIN\t%INFO/AC_Hemi\t%INFO/AC_Het\t%INFO/AC_Hom\t%INFO/AC_NFE\t%INFO/AC_OTH\t%INFO/AC_SAS\t%INFO/AF\t%INFO/AN\t%INFO/AN_AFR\t%INFO/AN_AMR\t%INFO/AN_Adj\t%INFO/AN_EAS\t%INFO/AN_FIN\t%INFO/AN_NFE\t%INFO/AN_OTH\t%INFO/AN_SAS\t%INFO/DP\t%INFO/GQ_MEAN\t%INFO/GQ_STDDEV\t%INFO/Hemi_AFR\t%INFO/Hemi_AMR\t%INFO/Hemi_EAS\t%INFO/Hemi_FIN\t%INFO/Hemi_NFE\t%INFO/Hemi_OTH\t%INFO/Hemi_SAS\t%INFO/Het_AFR\t%INFO/Het_AMR\t%INFO/Het_EAS\t%INFO/Het_FIN\t%INFO/Het_NFE\t%INFO/Het_OTH\t%INFO/Het_SAS\t%INFO/Hom_AFR\t%INFO/Hom_AMR\t%INFO/Hom_EAS\t%INFO/Hom_FIN\t%INFO/Hom_NFE\t%INFO/Hom_OTH\t%INFO/Hom_SAS\t%INFO/InbreedingCoeff\t%INFO/DP_HIST\t%INFO/GQ_HIST\t%INFO/CSQ\n' $TMPFILE > pipe &


#determine columns of population table and use awk command to add id column.
cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49 -d$'\t' pipe |awk -F'\t' '{$1=++i FS $1;}1' OFS='\t' > pipe4 &



#load population data into its table in mysql database.
mysql -u root -p201666 -e "LOAD DATA LOCAL INFILE 'pipe4' INTO TABLE DataSet4_db.variants_population FIELDS TERMINATED BY '\t' (Variant_id,Allele_Count, AC_African_American, AC_American, AC_Adjusted, AC_East_Asian, AC_Finnish, AC_Hemizygous, AC_Heterozygous, AC_Homozygous, AC_Non_Finnish_European, AC_Other, AC_South_Asian, Allele_Frequency, Allele_Number, AN_African_American, AN_American, AN_Adjusted, AN_East_Asian, AN_Finnish, AN_Non_Finnish, AN_Other, AN_South_Asian, Read_Depth, Genotype_Quality_MEAN, GQ_Standard_Deviation, Hemi_African_American, Hemi_American, Hemi_East_Asian, Hemi_Finnish, Hemi_Non_Finnish_European, Hemi_Other, Hemi_South_Asian, Het_African_American, Het_American, Het_East_Asian, Het_Finnish, Het_Non_Finnish_European, Het_Other, Het_South_Asian, Hom_African_American, Hom_American, Hom_East_Asian, Hom_Finnish, Hom_Non_Finnish_European_Homozygous, Hom_Other, Hom_South_Asian, InbreedingCoef, DP_HIST, GQ_HIST) SET id = NULL;" 

#copy data of population table and annotation table in 2 pipes.
bcftools query -f '%INFO/AC\t%INFO/AC_AFR\t%INFO/AC_AMR\t%INFO/AC_Adj\t%INFO/AC_EAS\t%INFO/AC_FIN\t%INFO/AC_Hemi\t%INFO/AC_Het\t%INFO/AC_Hom\t%INFO/AC_NFE\t%INFO/AC_OTH\t%INFO/AC_SAS\t%INFO/AF\t%INFO/AN\t%INFO/AN_AFR\t%INFO/AN_AMR\t%INFO/AN_Adj\t%INFO/AN_EAS\t%INFO/AN_FIN\t%INFO/AN_NFE\t%INFO/AN_OTH\t%INFO/AN_SAS\t%INFO/DP\t%INFO/GQ_MEAN\t%INFO/GQ_STDDEV\t%INFO/Hemi_AFR\t%INFO/Hemi_AMR\t%INFO/Hemi_EAS\t%INFO/Hemi_FIN\t%INFO/Hemi_NFE\t%INFO/Hemi_OTH\t%INFO/Hemi_SAS\t%INFO/Het_AFR\t%INFO/Het_AMR\t%INFO/Het_EAS\t%INFO/Het_FIN\t%INFO/Het_NFE\t%INFO/Het_OTH\t%INFO/Het_SAS\t%INFO/Hom_AFR\t%INFO/Hom_AMR\t%INFO/Hom_EAS\t%INFO/Hom_FIN\t%INFO/Hom_NFE\t%INFO/Hom_OTH\t%INFO/Hom_SAS\t%INFO/InbreedingCoeff\t%INFO/DP_HIST\t%INFO/GQ_HIST\t%INFO/CSQ\n' $TMPFILE  > pipe1 &

#cut annotations column | use awk to add id column | use awk to split annotations into separate rows.
cut -f50 -d$'\t' pipe1 |awk -F'\t' '{$1=++i FS $1;}1' OFS=,|awk '{n=split($2,s,",");for (i=1;i<=n;i++) {$2=s[i];print}}' |tr ' ' '|'|awk -F'\t' '{$1=i+1 FS $1;}1' OFS=, |tr '	' '|' > pipe6 &

#load annotation data into its table in mysql database.
mysql -u root -p201666 -e "LOAD DATA LOCAL INFILE 'pipe6' INTO TABLE DataSet4_db.variants_annotation FIELDS TERMINATED BY '|' (Population_id,Variant_id,Allele, Consequence, Impact, Gene_Symbol, Gene, Feature, Feature_type, Biotype, Exon, Intron, HGVS_DNA, HGVS_Protein, cDNA_position, CDS_position, Protein_position, Amino_acids, Codons, Existing_variation, ALLELE_NUM, DISTANCE, STRAND, VARIANT_CLASS, MINIMISED, SYMBOL_SOURCE, HGNC_ID, CANONICAL, TSL, CCDS, ENSP, SWISSPROT, TREMBL, UNIPARC, SIFT, PolyPhen, DOMAINS, HGVS_OFFSET, GMAF, AFR_MAF, AMR_MAF, ASN_MAF, EAS_MAF, EUR_MAF, SAS_MAF, AA_MAF, EA_MAF, CLIN_SIG, SOMATIC, PHENO, PUBMED, MOTIF_NAME, MOTIF_POS, HIGH_INF_POS, MOTIF_SCORE_CHANGE, LoF_info, LoF_flags, LoF_filter, LoF, context, ancestral) SET id = NULL  ;" 


mysql -u root -p201666 -e "update DataSet4_db.variants_annotation set Population_id = Variant_id where Population_id=1;"

mysql -u root -p201666 -e "SET FOREIGN_KEY_CHECKS = 1;"

#remove piping files and temp file
rm pipe pipe1 pipe2 pipe4 pipe5 pipe6
rm $TMPFILE


