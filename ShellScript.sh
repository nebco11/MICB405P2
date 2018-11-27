#!/bin/bash

#
#
#  Created by Benjamin Willox on 2018-11-25.
#


DATA_DIR=/projects/micb405/resources/project_2/2018/SaanichInlet_10m
RESOURCE_DIR=/projects/micb405/resources/project_2/2018/
PROKKA_HELPER_PATH=/projects/micb405/resources/project_2/2018/SaanichInlet_10m/MetaBAT2_SaanichInlet_10m/gtdbtk_output
WORK_DIR=/home/bwillox_mb18/p2
MAGS=/projects/micb405/resources/project_2/2018/SaanichInlet\_10m/MetaBAT2\_SaanichInlet_10m/MedQPlus_MAGs
META_T_DIR=/projects/micb405/resources/project_2/2018/Metatranscriptomes/
ASSIGNED_DEPTH=10m
RPKM=/projects/micb405/resources/project_2/2018/rpkm
THREADS=8


echo

echo "------Starting the script------"

echo 

if [[ -n "$WORK_DIR" && -n "$DATA_DIR" ]]; then
	echo "working directory is $WORK_DIR"
	echo "data directory is $DATA_DIR"
else
	echo "ERROR: directories don't exist"
	exit 1
fi

cd $WORK_DIR



echo

echo "------PROKKA------"

echo

mkdir -p $WORK_DIR/prokka_output

for f in $MAGS/*.fa
do
bin=${f#*.};
bin=${bin%.*};

taxi=$(grep -w $bin $PROKKA_HELPER_PATH/gtdbtk.*.classification_pplacer.tsv | awk '{ print $2 }' | awk -F";" '{ print $1 }'| sed 's/d__//g');
prokka --kingdom $taxi --outdir $WORK_DIR/prokka_output/${bin}/ --force --prefix ${bin}_SaanichInlet_MAG_ORFs $f
echo
echo
echo $taxi
echo
echo
done


echo "Mapping prokka IDs to MAGs..."
for f in ~/p2/prokka_output/*/*faa
do
prokka_id=$( head -1 $f | awk -F_ '{ print $1 }' | sed 's/^>//g' )
name=$(basename $f)
mag_id=$( echo $name | sed 's/.faa//g' | sed 's/_SaanichInlet_MAG_ORFs//g')
echo $prokka_id,"SaanichInlet_10m."$mag_id
done > Prokka_MAG_map.csv




echo
echo "concatenating all ORFs (.faa files)..."


    cat $WORK_DIR/prokka_output/*/*.faa >> "$WORK_DIR/SaanichInlet_10m_all_MAGs_ORFs.faa"
echo "concatenating all .ffn files..."
    cat $WORK_DIR/prokka_output/*/*.ffn >> "$WORK_DIR/SaanichInlet_10m_all_ref.ffn"
echo "done"

echo "You should now upload this file to the KAAS servers and then download and clean"

echo

echo "------METATRANSCRIPTOME Time------"

echo

 mkdir -p $WORK_DIR/align/index $WORK_DIR/align/sam $WORK_DIR/align/logs $WORK_DIR/align/bam
 mv $WORK_DIR/SaanichInlet_10m_all_ref.ffn $WORK_DIR/align/index

echo "creating index from reference..."
$BWA index -p $WORK_DIR/align/index/SaanichInlet_10m_all_ref $WORK_DIR/align/index/SaanichInlet_10m_all_ref.ffn

for f in $META_T_DIR/*"$ASSIGNED_DEPTH"*
do
 	name=$(basename $f)
 	name=${name//.gz}
 	name=${name//.fastq}
    echo "aligning file: $name"
 	$BWA mem -t $THREADS -p $WORK_DIR/align/index/SaanichInlet_10m_all_ref $f \
 		1> $WORK_DIR/align/sam/"$name".sam 2> $WORK_DIR/align/logs/"$name"_log.txt
 done

echo "creating csvs from rpkm"
mkdir -p $WORK_DIR/rpkm
for f in $WORK_DIR/align/sam/*; do
	name=$(basename $f)
	name=$(echo $name | cut -d. -f1)
    echo "rpkm: $name"
	$RPKM -c $WORK_DIR/align/index/SaanichInlet_10m_all_ref.ffn -a $f -o $WORK_DIR/rpkm/$name.rpkm.csv
done

echo

echo "------DONE------"
