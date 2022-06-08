#!/bin/sh
#SBATCH --time=24:00:00
#SBATCH -N 2
#SBATCH -n 12
#SBATCH --mem=64000


###############################################################################
# Experiment details
###############################################################################


###############################################################################
# Environment variables
###############################################################################
WORKDIR=$HOME/CAQ/
MAP=$WORKDIR/map/mycobiome.txt
MYCOBIOME=$WORKDIR/mycobiome_single
MYCOUT=$MYCOBIOME/output/rdp
MYCONFIG=$WORKDIR/config/qiime_parameters.txt
MYSAMPLES=$WORKDIR/samples/single
GG=$HOME/.bin/qiime_dbs/greengenes/gg_13_8_otus/rep_set/97_otus.fasta
UNITY=$HOME/.bin/qiime_dbs/unite/2/sh_refs_qiime_ver7_97_01.12.2017.fasta


source activate qiime

echo "###############################"
echo "###############################"
echo "STARTING PRE-PROCESSING..."
echo "###############################"
echo "###############################"


###############################################################################
# Trimming
###############################################################################

echo " "
echo "WARNING: STARTING THE LIBRARY TRIMMING PROCESS..."

mkdir -p $MYCOUT/trimming

cd $MYSAMPLES

for i in $(ls | rev | cut -c 11- | rev | uniq)

do

seqtk trimfq -b 0 -e 0 ${i}_001.fastq > $MYCOUT/trimming/${i}_001.fastq

done

echo "FINISHED TRIMMING PROCESS."
echo " "

###############################################################################
# Renaming folders
###############################################################################

cd $MYCOUT/trimming/

for sample_single in $(ls)

do

sample_rename=`echo ${sample_single} | cut -d "_" -f1`
mkdir ${sample_rename}
mv ${sample_single} ${sample_rename}/fastqjoin.join.fastq

done

###############################################################################

###############################################################################
# Demultiplex .fastq sequence data. 
# "Turning off" filter parameters, and storing the demultiplexed .fastq file. 
###############################################################################

echo "WARNING: STARTING SPLITTING PROCESS ..." 

multiple_split_libraries_fastq.py \
                -i $MYCOUT/trimming/ \
                -o $MYCOUT/splitting/ \
                -p $MYCONFIG \
                --remove_filepath_in_name \
                --include_input_dir_path \
                --demultiplexing_method sampleid_by_file


echo "SPLITTING PROCESS ENDED."
echo " "

###############################################################################

###############################################################################
# Split FASTA file
###############################################################################

echo "WARNING: DIVIDING FILE FASTA IN 3 PARTS TO REDUCE THE SIZE AND ALLOW THE EXECUTION BY USEARCH61..."  


mkdir -p $MYCOUT/split_fasta
cd $MYCOUT/split_fasta

#pyfasta split -n3 $MYCOUT/splitting/seqs.fna
#cp $MYCOUT/fastx_filter/reads.*.fa .

fasta-splitter --n-parts 3 $MYCOUT/splitting/seqs.fna --out-dir $MYCOUT/split_fasta

echo "FINISHED DIVISION."

echo " "

###############################################################################

###############################################################################
# Identify and filter chimeric seqs
###############################################################################

echo "WARNING: IDENTIFYING AND REMOVING CHIMERIC SEQUENCES ..."

for sample_paired in $(ls)

do

mkdir -p $MYCOUT/no_chimera/$sample_paired

identify_chimeric_seqs.py \
            -m usearch61 \
            -i $MYCOUT/split_fasta/$sample_paired \
            -r $UNITY \
            -o $MYCOUT/no_chimera/$sample_paired/

filter_fasta.py \
            -f $MYCOUT/split_fasta/$sample_paired \
            -o $MYCOUT/no_chimera/$sample_paired/seqs_chimeras_filtered.fna \
            -s $MYCOUT/no_chimera/$sample_paired/chimeras.txt -n

echo "CHIMERIC SEQUENCES REMOVED."
echo " "

done

###############################################################################

###############################################################################
# Joined FASTA file
###############################################################################

mkdir -p $MYCOUT/joined_fasta
cd $MYCOUT/joined_fasta

cat $MYCOUT/no_chimera/*/seqs_chimeras_filtered.fna > $MYCOUT/joined_fasta/seqs_chimeras_filtered.fna

###############################################################################

echo "###############################"
echo "###############################"
echo "PRE-PROCESSING COMPLETED..."
echo "###############################"
echo "###############################"

echo " "

echo "###############################"
echo "###############################"
echo "STARTING PROCESSING..."
echo "###############################"
echo "###############################"

###############################################################################
# PICKING OTU
###############################################################################

echo " "
echo "WARNING: STARTING PICK OTU PROCESS USING SORTMERNA..."

mkdir -p $MYCOUT/otu_picking/

pick_open_reference_otus.py \
                -m sortmerna_sumaclust \
                -i $MYCOUT/joined_fasta/seqs_chimeras_filtered.fna \
                -o $MYCOUT/otu_picking/rdp \
                -r $UNITY \
                -p $MYCONFIG \
                --suppress_align_and_tree \
                --force

echo "PICK OTU ENDED."
echo " "

echo "WARNING: ADJUSTING OTU TABLE..."

cp $MYCOUT/otu_picking/rdp/otu_table_mc2_w_tax.biom \
         $MYCOUT/otu_picking/rdp/otu_table_sortmerna.biom


filter_taxa_from_otu_table.py \
                -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna.biom \
                -o $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered.biom \
                -p k__Fungi

filter_otus_from_otu_table.py \
                -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered.biom \
                -o $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton.biom \
                -n 2

filter_otus_from_otu_table.py \
                -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered.biom \
                -o $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001.biom \
                --min_count_fraction 0.0001

filter_fasta.py \
                -f $MYCOUT/otu_picking/rdp/rep_set.fna \
                -b $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001.biom \
                -o $MYCOUT/otu_picking/rdp/rep_set_0001.fna

assign_taxonomy.py \
                -o $MYCOUT/otu_picking/rdp/rdp_assigned_taxonomy_0001 \
                -i $MYCOUT/otu_picking/rdp/rep_set_0001.fna \
                -p $MYCONFIG

filter_otus_from_otu_table.py \
                -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton.biom \
                -o $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001.biom \
                --min_count_fraction 0.00001

filter_fasta.py \
                -f $MYCOUT/otu_picking/rdp/rep_set.fna \
                -b $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001.biom \
                -o $MYCOUT/otu_picking/rdp/rep_set_00001.fna

assign_taxonomy.py \
                -o $MYCOUT/otu_picking/rdp/rdp_assigned_taxonomy_00001 \
                -i $MYCOUT/otu_picking/rdp/rep_set_00001.fna \
                -p $MYCONFIG


echo "OTU TABLE ADJUSTMENTS FINISHED."
echo " "

echo "WARNING: CONVERTING OTU TABLE TO TSV..."

biom convert \
                -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001.biom \
                -o $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001.txt  \
                -m $MAP  \
                --header-key=taxonomy --output-metadata-id="Consensus Lineage" \
                --process-obs-metadata=taxonomy \
                --table-type="OTU table" \
                --to-tsv

biom convert \
                -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001.biom \
                -o $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001.txt  \
                -m $MAP  \
                --header-key=taxonomy --output-metadata-id="Consensus Lineage" \
                --process-obs-metadata=taxonomy \
                --table-type="OTU table" \
                --to-tsv

biom convert \
                -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton.biom \
                -o $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton.txt  \
                -m $MAP  \
                --header-key=taxonomy --output-metadata-id="Consensus Lineage" \
                --process-obs-metadata=taxonomy \
                --table-type="OTU table" \
                --to-tsv

biom convert \
                -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered.biom \
                -o $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered.txt  \
                -m $MAP  \
                --header-key=taxonomy --output-metadata-id="Consensus Lineage" \
                --process-obs-metadata=taxonomy \
                --table-type="OTU table" \
                --to-tsv

biom convert \
                -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna.biom \
                -o $MYCOUT/otu_picking/rdp/otu_table_sortmerna.txt  \
                -m $MAP  \
                --header-key=taxonomy --output-metadata-id="Consensus Lineage" \
                --process-obs-metadata=taxonomy \
                --table-type="OTU table" \
                --to-tsv


echo "CONVERSION OF OTU TABLE TO TSV COMPLETED."
echo " "

echo "WARNING: GENERATING OTU TABLE REPORT..."

biom summarize-table \
                -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001.biom \
                -o $MYCOUT/otu_picking/rdp/biom_summary_0001.txt

biom summarize-table \
                -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001.biom \
                -o $MYCOUT/otu_picking/rdp/biom_summary_00001.txt


biom summarize-table \
                -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton.biom \
                -o $MYCOUT/otu_picking/rdp/biom_summary_singleton.txt

biom summarize-table \
                -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered.biom \
                -o $MYCOUT/otu_picking/rdp/biom_summary_filtered.txt


biom summarize-table \
                -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna.biom \
                -o $MYCOUT/otu_picking/rdp/biom_summary.txt




echo "FINISHED REPORT."
echo " "

echo "###############################"
echo "###############################"
echo "PROCESSING COMPLETED..."
echo "###############################"
echo "###############################"

echo " "

###############################################################################
# Microbiome Analyst format
###############################################################################

echo "#################################################"
echo "#################################################"
echo "FORMATTING OTU TABLE FOR THE MICROBIOME ANALYST..."
echo "#################################################"
echo "#################################################"


echo " "
echo " "

mkdir -p $MYCOUT/otu_picking/rdp/microbiomeanalyst

sed "s/k__Fungi/k__Bacteria/g" \
        $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001.txt \
        > $MYCOUT/otu_picking/rdp/microbiomeanalyst/otu_table_sortmerna_filtered_singleton_abundance_0001_manalyst.txt

sed "s/k__Fungi/k__Bacteria/g" \
        $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001.txt \
        > $MYCOUT/otu_picking/rdp/microbiomeanalyst/otu_table_sortmerna_filtered_singleton_abundance_00001_manalyst.txt

sed "s/k__Fungi/k__Bacteria/g" \
        $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton.txt \
        > $MYCOUT/otu_picking/rdp/microbiomeanalyst/otu_table_sortmerna_filtered_singleton_manalyst.txt

sed "s/k__Fungi/k__Bacteria/g" \
        $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered.txt \
        > $MYCOUT/otu_picking/rdp/microbiomeanalyst/otu_table_sortmerna_filtered_manalyst.txt

sed "s/k__Fungi/k__Bacteria/g" \
	$MYCOUT/otu_picking/rdp/otu_table_sortmerna.txt \
	> $MYCOUT/otu_picking/rdp/microbiomeanalyst/otu_table_sortmerna_manalyst.txt

biom convert \
	-i $MYCOUT/otu_picking/rdp/microbiomeanalyst/otu_table_sortmerna_filtered_singleton_abundance_0001_manalyst.txt  \
	-o $MYCOUT/otu_picking/rdp/microbiomeanalyst/otu_table_sortmerna_filtered_singleton_abundance_0001_manalyst.biom \
	--to-hdf5 \
	--table-type="OTU table" \
	--process-obs-metadata taxonomy

biom convert \
        -i $MYCOUT/otu_picking/rdp/microbiomeanalyst/otu_table_sortmerna_filtered_singleton_abundance_00001_manalyst.txt  \
        -o $MYCOUT/otu_picking/rdp/microbiomeanalyst/otu_table_sortmerna_filtered_singleton_abundance_00001_manalyst.biom \
        --to-hdf5 \
        --table-type="OTU table" \
        --process-obs-metadata taxonomy

biom convert \
        -i $MYCOUT/otu_picking/rdp/microbiomeanalyst/otu_table_sortmerna_filtered_singleton_manalyst.txt  \
        -o $MYCOUT/otu_picking/rdp/microbiomeanalyst/otu_table_sortmerna_filtered_singleton_manalyst.biom \
        --to-hdf5 \
        --table-type="OTU table" \
        --process-obs-metadata taxonomy

biom convert \
        -i $MYCOUT/otu_picking/rdp/microbiomeanalyst/otu_table_sortmerna_filtered_manalyst.txt  \
        -o $MYCOUT/otu_picking/rdp/microbiomeanalyst/otu_table_sortmerna_filtered_manalyst.biom \
        --to-hdf5 \
        --table-type="OTU table" \
        --process-obs-metadata taxonomy

biom convert \
        -i $MYCOUT/otu_picking/rdp/microbiomeanalyst/otu_table_sortmerna_manalyst.txt  \
        -o $MYCOUT/otu_picking/rdp/microbiomeanalyst/otu_table_sortmerna_manalyst.biom \
        --to-hdf5 \
        --table-type="OTU table" \
        --process-obs-metadata taxonomy

###############################################################################
# Core Microbiome
###############################################################################


echo "###############################"
echo "###############################"
echo "CREATE CORE MICROBIOME..."
echo "###############################"
echo "###############################"


echo " "
echo " "

mkdir -p $MYCOUT/otu_picking/rdp/core

### 0001 abundance

filter_samples_from_otu_table.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001.biom \
        -o $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001_CAQ-B.biom \
        -m $MAP \
        -s 'ALL:CAQ-B'

filter_samples_from_otu_table.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001.biom \
        -o $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001_SUP-B.biom \
        -m $MAP \
        -s 'ALL:SUP-B'

filter_samples_from_otu_table.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001.biom \
        -o $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001_TB-B.biom \
        -m $MAP \
        -s 'ALL:TB-B'

filter_samples_from_otu_table.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001.biom \
        -o $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001_CAQ-A.biom \
        -m $MAP \
        -s 'ALL:CAQ-A'

filter_samples_from_otu_table.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001.biom \
        -o $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001_SUP-A.biom \
        -m $MAP \
        -s 'ALL:SUP-A'

filter_samples_from_otu_table.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001.biom \
        -o $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001_TB-A.biom \
        -m $MAP \
        -s 'ALL:TB-A'

compute_core_microbiome.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001_CAQ-B.biom \
        -o $MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_CAQ-B

compute_core_microbiome.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001_SUP-B.biom \
        -o $MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_SUP-B

compute_core_microbiome.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001_TB-B.biom \
        -o $MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_TB-B

compute_core_microbiome.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001_CAQ-A.biom \
        -o $MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_CAQ-A

compute_core_microbiome.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001_SUP-A.biom \
        -o $MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_SUP-A

compute_core_microbiome.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001_TB-A.biom \
        -o $MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_TB-A

merge_otu_tables.py \
        -i $MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_CAQ-A/core_table_100.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_CAQ-B/core_table_100.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_SUP-A/core_table_100.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_SUP-B/core_table_100.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_TB-A/core_table_100.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_TB-B/core_table_100.biom \
        -o $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_100.biom

merge_otu_tables.py \
        -i $MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_CAQ-A/core_table_90.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_CAQ-B/core_table_90.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_SUP-A/core_table_90.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_SUP-B/core_table_90.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_TB-A/core_table_90.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_TB-B/core_table_90.biom \
        -o $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_90.biom

merge_otu_tables.py \
        -i $MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_CAQ-A/core_table_80.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_CAQ-B/core_table_80.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_SUP-A/core_table_80.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_SUP-B/core_table_80.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_TB-A/core_table_80.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_TB-B/core_table_80.biom \
        -o $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_80.biom

merge_otu_tables.py \
        -i $MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_CAQ-A/core_table_70.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_CAQ-B/core_table_70.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_SUP-A/core_table_70.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_SUP-B/core_table_70.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_TB-A/core_table_70.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_TB-B/core_table_70.biom \
        -o $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_70.biom


biom convert \
	-i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_100.biom \
	-o $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_100.txt \
	-m $MAP \
	--header-key=taxonomy \
	--output-metadata-id="Consensus Lineage" \
	--process-obs-metadata=taxonomy \
	--table-type="OTU table" \
	--to-tsv

biom convert \
        -i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_90.biom \
        -o $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_90.txt \
        -m $MAP \
        --header-key=taxonomy \
        --output-metadata-id="Consensus Lineage" \
        --process-obs-metadata=taxonomy \
        --table-type="OTU table" \
        --to-tsv
biom convert \
        -i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_80.biom \
        -o $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_80.txt \
        -m $MAP \
        --header-key=taxonomy \
        --output-metadata-id="Consensus Lineage" \
        --process-obs-metadata=taxonomy \
        --table-type="OTU table" \
        --to-tsv

biom convert \
        -i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_70.biom \
        -o $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_70.txt \
        -m $MAP \
        --header-key=taxonomy \
        --output-metadata-id="Consensus Lineage" \
        --process-obs-metadata=taxonomy \
        --table-type="OTU table" \
        --to-tsv


### 00001 abundance

filter_samples_from_otu_table.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001.biom \
        -o $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001_CAQ-B.biom \
        -m $MAP \
        -s 'ALL:CAQ-B'

filter_samples_from_otu_table.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001.biom \
        -o $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001_SUP-B.biom \
        -m $MAP \
        -s 'ALL:SUP-B'

filter_samples_from_otu_table.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001.biom \
        -o $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001_TB-B.biom \
        -m $MAP \
        -s 'ALL:TB-B'


filter_samples_from_otu_table.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001.biom \
        -o $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001_CAQ-A.biom \
        -m $MAP \
        -s 'ALL:CAQ-A'

filter_samples_from_otu_table.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001.biom \
        -o $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001_SUP-A.biom \
        -m $MAP \
        -s 'ALL:SUP-A'

filter_samples_from_otu_table.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001.biom \
        -o $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001_TB-A.biom \
        -m $MAP \
        -s 'ALL:TB-A'


compute_core_microbiome.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001_CAQ-B.biom \
        -o $MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_CAQ-B

compute_core_microbiome.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001_SUP-B.biom \
        -o $MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_SUP-B


compute_core_microbiome.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001_TB-B.biom \
        -o $MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_TB-B

compute_core_microbiome.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001_CAQ-A.biom \
        -o $MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_CAQ-A

compute_core_microbiome.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001_SUP-A.biom \
        -o $MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_SUP-A

compute_core_microbiome.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001_TB-A.biom \
        -o $MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_TB-A

merge_otu_tables.py \
        -i $MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_CAQ-A/core_table_100.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_CAQ-B/core_table_100.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_SUP-A/core_table_100.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_SUP-B/core_table_100.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_TB-A/core_table_100.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_TB-B/core_table_100.biom \
        -o $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_100.biom

merge_otu_tables.py \
        -i $MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_CAQ-A/core_table_90.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_CAQ-B/core_table_90.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_SUP-A/core_table_90.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_SUP-B/core_table_90.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_TB-A/core_table_90.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_TB-B/core_table_90.biom \
        -o $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_90.biom

merge_otu_tables.py \
        -i $MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_CAQ-A/core_table_80.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_CAQ-B/core_table_80.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_SUP-A/core_table_80.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_SUP-B/core_table_80.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_TB-A/core_table_80.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_TB-B/core_table_80.biom \
        -o $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_80.biom

merge_otu_tables.py \
        -i $MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_CAQ-A/core_table_70.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_CAQ-B/core_table_70.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_SUP-A/core_table_70.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_SUP-B/core_table_70.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_TB-A/core_table_70.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_TB-B/core_table_70.biom \
        -o $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_70.biom


biom convert \
        -i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_100.biom \
        -o $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_100.txt \
        -m $MAP \
        --header-key=taxonomy \
        --output-metadata-id="Consensus Lineage" \
        --process-obs-metadata=taxonomy \
        --table-type="OTU table" \
        --to-tsv

biom convert \
        -i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_90.biom \
        -o $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_90.txt \
        -m $MAP \
        --header-key=taxonomy \
        --output-metadata-id="Consensus Lineage" \
        --process-obs-metadata=taxonomy \
        --table-type="OTU table" \
        --to-tsv
biom convert \
        -i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_80.biom \
        -o $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_80.txt \
        -m $MAP \
        --header-key=taxonomy \
        --output-metadata-id="Consensus Lineage" \
        --process-obs-metadata=taxonomy \
        --table-type="OTU table" \
        --to-tsv

biom convert \
        -i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_70.biom \
        -o $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_70.txt \
        -m $MAP \
        --header-key=taxonomy \
        --output-metadata-id="Consensus Lineage" \
        --process-obs-metadata=taxonomy \
        --table-type="OTU table" \
        --to-tsv

### Singletons

### singleton abundance

filter_samples_from_otu_table.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton.biom \
        -o $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_CAQ-B.biom \
        -m $MAP \
        -s 'ALL:CAQ-B'

filter_samples_from_otu_table.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton.biom \
        -o $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_SUP-B.biom \
        -m $MAP \
        -s 'ALL:SUP-B'

filter_samples_from_otu_table.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton.biom \
        -o $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_TB-B.biom \
        -m $MAP \
        -s 'ALL:TB-B'

filter_samples_from_otu_table.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton.biom \
        -o $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_CAQ-A.biom \
        -m $MAP \
        -s 'ALL:CAQ-A'

filter_samples_from_otu_table.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton.biom \
        -o $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_SUP-A.biom \
        -m $MAP \
        -s 'ALL:SUP-A'

filter_samples_from_otu_table.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton.biom \
        -o $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_TB-A.biom \
        -m $MAP \
        -s 'ALL:TB-A'

compute_core_microbiome.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_CAQ-B.biom \
        -o $MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_CAQ-B

compute_core_microbiome.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_SUP-B.biom \
        -o $MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_SUP-B

compute_core_microbiome.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_TB-B.biom \
        -o $MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_TB-B

compute_core_microbiome.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_CAQ-A.biom \
        -o $MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_CAQ-A

compute_core_microbiome.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_SUP-A.biom \
        -o $MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_SUP-A

compute_core_microbiome.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_TB-A.biom \
        -o $MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_TB-A


merge_otu_tables.py \
        -i $MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_CAQ-A/core_table_100.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_CAQ-B/core_table_100.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_SUP-A/core_table_100.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_SUP-B/core_table_100.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_TB-A/core_table_100.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_TB-B/core_table_100.biom \
        -o $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_100.biom

merge_otu_tables.py \
        -i $MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_CAQ-A/core_table_90.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_CAQ-B/core_table_90.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_SUP-A/core_table_90.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_SUP-B/core_table_90.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_TB-A/core_table_90.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_TB-B/core_table_90.biom \
        -o $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_90.biom

merge_otu_tables.py \
        -i $MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_CAQ-A/core_table_80.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_CAQ-B/core_table_80.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_SUP-A/core_table_80.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_SUP-B/core_table_80.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_TB-A/core_table_80.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_TB-B/core_table_80.biom \
        -o $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_80.biom

merge_otu_tables.py \
        -i $MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_CAQ-A/core_table_70.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_CAQ-B/core_table_70.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_SUP-A/core_table_70.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_SUP-B/core_table_70.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_TB-A/core_table_70.biom,$MYCOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_TB-B/core_table_70.biom \
        -o $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_70.biom


biom convert \
        -i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_100.biom \
        -o $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_100.txt \
        -m $MAP \
        --header-key=taxonomy \
        --output-metadata-id="Consensus Lineage" \
        --process-obs-metadata=taxonomy \
        --table-type="OTU table" \
        --to-tsv

biom convert \
        -i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_90.biom \
        -o $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_90.txt \
        -m $MAP \
        --header-key=taxonomy \
        --output-metadata-id="Consensus Lineage" \
        --process-obs-metadata=taxonomy \
        --table-type="OTU table" \
        --to-tsv
biom convert \
        -i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_80.biom \
        -o $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_80.txt \
        -m $MAP \
        --header-key=taxonomy \
        --output-metadata-id="Consensus Lineage" \
        --process-obs-metadata=taxonomy \
        --table-type="OTU table" \
        --to-tsv

biom convert \
        -i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_70.biom \
        -o $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_70.txt \
        -m $MAP \
        --header-key=taxonomy \
        --output-metadata-id="Consensus Lineage" \
        --process-obs-metadata=taxonomy \
        --table-type="OTU table" \
        --to-tsv


###############################################################################
# Microbiome Analyst format
###############################################################################


echo "#################################################"
echo "#################################################"
echo "FORMATING CORE FOR THE ANALYST MICROBIOME..."
echo "#################################################"
echo "#################################################"


echo " "
echo " "

mkdir -p $MYCOUT/otu_picking/rdp/core/microbiomeanalyst

### 0001 abundance

sed "s/k__Fungi/k__Bacteria/g" \
        $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_100.txt \
        > $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_abundance_0001_100_manalyst.txt

sed "s/k__Fungi/k__Bacteria/g" \
        $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_90.txt \
        > $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_abundance_0001_90_manalyst.txt

sed "s/k__Fungi/k__Bacteria/g" \
        $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_80.txt \
        > $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_abundance_0001_80_manalyst.txt

sed "s/k__Fungi/k__Bacteria/g" \
        $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_70.txt \
        > $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_abundance_0001_70_manalyst.txt

biom convert \
        -i $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_abundance_0001_100_manalyst.txt \
        -o $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_abundance_0001_100_manalyst.biom \
        --to-hdf5 \
        --table-type="OTU table" \
        --process-obs-metadata taxonomy

biom convert \
        -i $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_abundance_0001_90_manalyst.txt \
        -o $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_abundance_0001_90_manalyst.biom \
        --to-hdf5 \
        --table-type="OTU table" \
        --process-obs-metadata taxonomy

biom convert \
        -i $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_abundance_0001_80_manalyst.txt \
        -o $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_abundance_0001_80_manalyst.biom \
        --to-hdf5 \
        --table-type="OTU table" \
        --process-obs-metadata taxonomy

biom convert \
        -i $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_abundance_0001_70_manalyst.txt \
        -o $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_abundance_0001_70_manalyst.biom \
        --to-hdf5 \
        --table-type="OTU table" \
        --process-obs-metadata taxonomy


biom summarize-table \
                -i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_100.biom \
                -o $MYCOUT/otu_picking/rdp/core/biom_summary_core_abundance_0001_100.txt

biom summarize-table \
                -i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_90.biom \
                -o $MYCOUT/otu_picking/rdp/core/biom_summary_core_abundance_0001_90.txt

biom summarize-table \
                -i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_80.biom \
                -o $MYCOUT/otu_picking/rdp/core/biom_summary_core_abundance_0001_80.txt

biom summarize-table \
                -i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_70.biom \
                -o $MYCOUT/otu_picking/rdp/core/biom_summary_core_abundance_0001_70.txt


### 00001 abundance

sed "s/k__Fungi/k__Bacteria/g" \
        $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_100.txt \
        > $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_abundance_00001_100_manalyst.txt

sed "s/k__Fungi/k__Bacteria/g" \
        $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_90.txt \
        > $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_abundance_00001_90_manalyst.txt

sed "s/k__Fungi/k__Bacteria/g" \
        $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_80.txt \
        > $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_abundance_00001_80_manalyst.txt

sed "s/k__Fungi/k__Bacteria/g" \
        $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_70.txt \
        > $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_abundance_00001_70_manalyst.txt

biom convert \
        -i $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_abundance_00001_100_manalyst.txt \
        -o $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_abundance_00001_100_manalyst.biom \
        --to-hdf5 \
        --table-type="OTU table" \
        --process-obs-metadata taxonomy

biom convert \
        -i $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_abundance_00001_90_manalyst.txt \
        -o $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_abundance_00001_90_manalyst.biom \
        --to-hdf5 \
        --table-type="OTU table" \
        --process-obs-metadata taxonomy

biom convert \
        -i $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_abundance_00001_80_manalyst.txt \
        -o $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_abundance_00001_80_manalyst.biom \
        --to-hdf5 \
        --table-type="OTU table" \
        --process-obs-metadata taxonomy

biom convert \
        -i $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_abundance_00001_70_manalyst.txt \
        -o $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_abundance_00001_70_manalyst.biom \
        --to-hdf5 \
        --table-type="OTU table" \
        --process-obs-metadata taxonomy


biom summarize-table \
                -i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_100.biom \
                -o $MYCOUT/otu_picking/rdp/core/biom_summary_core_abundance_00001_100.txt

biom summarize-table \
                -i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_90.biom \
                -o $MYCOUT/otu_picking/rdp/core/biom_summary_core_abundance_00001_90.txt

biom summarize-table \
                -i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_80.biom \
                -o $MYCOUT/otu_picking/rdp/core/biom_summary_core_abundance_00001_80.txt

biom summarize-table \
                -i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_70.biom \
                -o $MYCOUT/otu_picking/rdp/core/biom_summary_core_abundance_00001_70.txt


### Core Singleton

sed "s/k__Fungi/k__Bacteria/g" \
        $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_100.txt \
        > $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_singleton_100_manalyst.txt

sed "s/k__Fungi/k__Bacteria/g" \
        $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_90.txt \
        > $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_singleton_90_manalyst.txt

sed "s/k__Fungi/k__Bacteria/g" \
        $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_80.txt \
        > $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_singleton_80_manalyst.txt

sed "s/k__Fungi/k__Bacteria/g" \
        $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_70.txt \
        > $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_singleton_70_manalyst.txt

biom convert \
        -i $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_singleton_100_manalyst.txt \
        -o $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_singleton_100_manalyst.biom \
        --to-hdf5 \
        --table-type="OTU table" \
        --process-obs-metadata taxonomy

biom convert \
        -i $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_singleton_90_manalyst.txt \
        -o $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_singleton_90_manalyst.biom \
        --to-hdf5 \
        --table-type="OTU table" \
        --process-obs-metadata taxonomy

biom convert \
        -i $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_singleton_80_manalyst.txt \
        -o $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_singleton_80_manalyst.biom \
        --to-hdf5 \
        --table-type="OTU table" \
        --process-obs-metadata taxonomy

biom convert \
        -i $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_singleton_70_manalyst.txt \
        -o $MYCOUT/otu_picking/rdp/core/microbiomeanalyst/merged_otu_table_core_singleton_70_manalyst.biom \
        --to-hdf5 \
        --table-type="OTU table" \
        --process-obs-metadata taxonomy


biom summarize-table \
                -i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_100.biom \
                -o $MYCOUT/otu_picking/rdp/core/biom_summary_core_abundance_00001_100.txt

biom summarize-table \
                -i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_90.biom \
                -o $MYCOUT/otu_picking/rdp/core/biom_summary_core_abundance_00001_90.txt

biom summarize-table \
                -i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_80.biom \
                -o $MYCOUT/otu_picking/rdp/core/biom_summary_core_abundance_00001_80.txt

biom summarize-table \
                -i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_70.biom \
                -o $MYCOUT/otu_picking/rdp/core/biom_summary_core_abundance_00001_70.txt


###############################################################################
# Make OTU network
###############################################################################

echo "###############################"
echo "###############################"
echo "CREATE OTU NETWORK..."
echo "###############################"
echo "###############################"


echo " "
echo " "

mkdir -p $MYCOUT/otu_picking/rdp/networks/core


make_otu_network.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton.biom \
        -m $MAP \
        -o $MYCOUT/otu_picking/rdp/networks/otu_network_singleton

make_otu_network.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered.biom \
        -m $MAP \
        -o $MYCOUT/otu_picking/rdp/networks/otu_network_filtered

make_otu_network.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna.biom \
        -m $MAP \
        -o $MYCOUT/otu_picking/rdp/networks/otu_network

### 0001 abundance

make_otu_network.py \
	-i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001.biom \
	-m $MAP \
	-o $MYCOUT/otu_picking/rdp/networks/otu_network_0001


make_otu_network.py \
        -i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_100.biom \
        -m $MAP \
        -o $MYCOUT/otu_picking/rdp/networks/core/otu_network_core_abundance_0001_100

make_otu_network.py \
        -i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_90.biom \
        -m $MAP \
        -o $MYCOUT/otu_picking/rdp/networks/core/otu_network_core_abundance_0001_90

make_otu_network.py \
        -i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_80.biom \
        -m $MAP \
        -o $MYCOUT/otu_picking/rdp/networks/core/otu_network_core_abundance_0001_80

make_otu_network.py \
        -i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_70.biom \
        -m $MAP \
        -o $MYCOUT/otu_picking/rdp/networks/core/otu_network_core_abundance_0001_70

### 00001 abundance

make_otu_network.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001.biom \
        -m $MAP \
        -o $MYCOUT/otu_picking/rdp/networks/otu_network_00001

make_otu_network.py \
        -i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_100.biom \
        -m $MAP \
        -o $MYCOUT/otu_picking/rdp/networks/core/otu_network_core_abundance_0001_100

make_otu_network.py \
        -i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_90.biom \
        -m $MAP \
        -o $MYCOUT/otu_picking/rdp/networks/core/otu_network_core_abundance_00001_90

make_otu_network.py \
        -i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_80.biom \
        -m $MAP \
        -o $MYCOUT/otu_picking/rdp/networks/core/otu_network_core_abundance_00001_80

make_otu_network.py \
        -i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_70.biom \
        -m $MAP \
        -o $MYCOUT/otu_picking/rdp/networks/core/otu_network_core_abundance_00001_70

### Core Singletons

make_otu_network.py \
        -i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_100.biom \
        -m $MAP \
        -o $MYCOUT/otu_picking/rdp/networks/core/otu_network_core_abundance_0001_100

make_otu_network.py \
        -i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_90.biom \
        -m $MAP \
        -o $MYCOUT/otu_picking/rdp/networks/core/otu_network_core_abundance_00001_90

make_otu_network.py \
        -i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_80.biom \
        -m $MAP \
        -o $MYCOUT/otu_picking/rdp/networks/core/otu_network_core_abundance_00001_80

make_otu_network.py \
        -i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_70.biom \
        -m $MAP \
        -o $MYCOUT/otu_picking/rdp/networks/core/otu_network_core_abundance_00001_70

###############################################################################
# MAKE CORE DIVERSITY
###############################################################################

echo "###############################"
echo "###############################"
echo "STARTING CORE DIVERSITY..."
echo "###############################"
echo "###############################"


echo " "
echo " "


core_diversity_analyses.py \
	-i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton.biom \
	-o $MYCOUT/core_output/singleton/ \
	-m $MAP \
	-e 10000 \
	-p $MYCONFIG \
	--nonphylogenetic_diversity


core_diversity_analyses.py \
	-i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001.biom \
	-o $MYCOUT/core_output/00001/ \
	-m $MAP \
	-e 10000 \
	-p $MYCONFIG \
	--nonphylogenetic_diversity

core_diversity_analyses.py \
        -i $MYCOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001.biom \
        -o $MYCOUT/core_output/0001/ \
        -m $MAP \
        -e 10000 \
        -p $MYCONFIG \
        --nonphylogenetic_diversity


core_diversity_analyses.py \
	-i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_70.biom \
	-o $MYCOUT/core_output/core_0001_70p/ \
	-m $MAP \
	-e 10000 \
	-p $MYCONFIG \
	--nonphylogenetic_diversity

core_diversity_analyses.py \
	-i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_80.biom \
	-o $MYCOUT/core_output/core_0001_80p/ \
	-m $MAP \
	-e 10000 \
	-p $MYCONFIG \
	--nonphylogenetic_diversity

core_diversity_analyses.py \
	-i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_90.biom \
	-o $MYCOUT/core_output/core_0001_90p/ \
	-m $MAP \
	-e 10000 \
	-p $MYCONFIG \
	--nonphylogenetic_diversity

core_diversity_analyses.py \
	-i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_70.biom \
	-o $MYCOUT/core_output/core_00001_70p/ \
	-m $MAP \
	-e 10000 \
	-p $MYCONFIG \
	--nonphylogenetic_diversity

core_diversity_analyses.py \
	-i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_80.biom \
	-o $MYCOUT/core_output/core_00001_80p/ \
	-m $MAP \
	-e 10000 \
	-p $MYCONFIG \
	--nonphylogenetic_diversity

core_diversity_analyses.py \
	-i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_90.biom \
	-o $MYCOUT/core_output/core_00001_90p/ \
	-m $MAP \
	-e 10000 \
	-p $MYCONFIG \
	--nonphylogenetic_diversity

core_diversity_analyses.py \
	-i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_70.biom \
	-o $MYCOUT/core_output/core_singleton_70p/ \
	-m $MAP \
	-e 10000 \
	-p $MYCONFIG \
	--nonphylogenetic_diversity

core_diversity_analyses.py \
	-i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_80.biom \
	-o $MYCOUT/core_output/core_singleton_80p/ \
	-m $MAP \
	-e 10000 \
	-p $MYCONFIG \
	--nonphylogenetic_diversity

core_diversity_analyses.py \
	-i $MYCOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_90.biom \
	-o $MYCOUT/core_output/core_singleton_90p/ \
	-m $MAP \
	-e 10000 \
	-p $MYCONFIG \
	--nonphylogenetic_diversity


###############################################################################
# PROCRUSTES ANALYSIS
###############################################################################

transform_coordinate_matrices.py \
                -i $MYCOUT/core_output/singleton/bdiv_even10000/bray_curtis_pc.txt,$MYCOUT/core_output/singleton/bdiv_even10000/euclidean_pc.txt \
                -r 999 \
                -o $MYCOUT/core_output/singleton/procrustes_results_10000/


transform_coordinate_matrices.py \
                -i $MYCOUT/core_output/00001/bdiv_even10000/bray_curtis_pc.txt,$MYCOUT/core_output/00001/bdiv_even10000/euclidean_pc.txt \
                -r 999 \
                -o $MYCOUT/core_output/00001/procrustes_results_10000/

transform_coordinate_matrices.py \
                -i $MYCOUT/core_output/0001/bdiv_even10000/bray_curtis_pc.txt,$MYCOUT/core_output/0001/bdiv_even10000/euclidean_pc.txt \
                -r 999 \
                -o $MYCOUT/core_output/0001/procrustes_results_10000/

transform_coordinate_matrices.py \
                -i $MYCOUT/core_output/core_0001_90p/bdiv_even10000/bray_curtis_pc.txt,$MYCOUT/core_output/core_0001_90p/bdiv_even10000/euclidean_pc.txt \
                -r 999 \
                -o $MYCOUT/core_output/core_0001_90p/procrustes_results_10000/

transform_coordinate_matrices.py \
                -i $MYCOUT/core_output/core_0001_80p/bdiv_even10000/bray_curtis_pc.txt,$MYCOUT/core_output/core_0001_80p/bdiv_even10000/euclidean_pc.txt \
                -r 999 \
                -o $MYCOUT/core_output/core_0001_80p/procrustes_results_10000/

transform_coordinate_matrices.py \
                -i $MYCOUT/core_output/core_0001_70p/bdiv_even10000/bray_curtis_pc.txt,$MYCOUT/core_output/core_0001_70p/bdiv_even10000/euclidean_pc.txt \
                -r 999 \
                -o $MYCOUT/core_output/core_0001_70p/procrustes_results_10000/

transform_coordinate_matrices.py \
                -i $MYCOUT/core_output/core_00001_90p/bdiv_even10000/bray_curtis_pc.txt,$MYCOUT/core_output/core_00001_90p/bdiv_even10000/euclidean_pc.txt \
                -r 999 \
                -o $MYCOUT/core_output/core_00001_90p/procrustes_results_10000/

transform_coordinate_matrices.py \
                -i $MYCOUT/core_output/core_00001_80p/bdiv_even10000/bray_curtis_pc.txt,$MYCOUT/core_output/core_00001_80p/bdiv_even10000/euclidean_pc.txt \
                -r 999 \
                -o $MYCOUT/core_output/core_00001_80p/procrustes_results_10000/

transform_coordinate_matrices.py \
                -i $MYCOUT/core_output/core_00001_70p/bdiv_even10000/bray_curtis_pc.txt,$MYCOUT/core_output/core_00001_70p/bdiv_even10000/euclidean_pc.txt \
                -r 999 \
                -o $MYCOUT/core_output/core_00001_70p/procrustes_results_10000/

transform_coordinate_matrices.py \
                -i $MYCOUT/core_output/core_singleton_90p/bdiv_even10000/bray_curtis_pc.txt,$MYCOUT/core_output/core_singleton_90p/bdiv_even10000/euclidean_pc.txt \
                -r 999 \
                -o $MYCOUT/core_output/core_singleton_90p/procrustes_results_10000/

transform_coordinate_matrices.py \
                -i $MYCOUT/core_output/core_singleton_80p/bdiv_even10000/bray_curtis_pc.txt,$MYCOUT/core_output/core_singleton_80p/bdiv_even10000/euclidean_pc.txt \
                -r 999 \
                -o $MYCOUT/core_output/core_singleton_80p/procrustes_results_10000/

transform_coordinate_matrices.py \
                -i $MYCOUT/core_output/core_singleton_70p/bdiv_even10000/bray_curtis_pc.txt,$MYCOUT/core_output/core_singleton_70p/bdiv_even10000/euclidean_pc.txt \
                -r 999 \
                -o $MYCOUT/core_output/core_singleton_70p/procrustes_results_10000/

### 00001 Abundance

cd $MYCOUT/core_output/00001
gunzip table_mc10000.biom.gz
gunzip table_even10000.biom.gz


biom summarize-table \
                -i table_mc10000.biom \
                -o biom_summary_table_mc10000.txt

biom summarize-table \
                -i table_even10000.biom \
                -o biom_summary_table_even10000.txt

###############################################################################
# BETA DIVERSITY
###############################################################################

for i in $(cat $WORKDIR/scripts/taxon.txt)

do

beta_diversity.py \
            -i $MYCOUT/core_output/00001/taxa_plots/table_mc10000_sorted.biom \
            -o $MYCOUT/core_output/00001/bdiv_${i}/ \
            --metrics euclidean

mv $MYCOUT/core_output/00001/bdiv_${i}/euclidean_table_mc10000_sorted.txt \
            $MYCOUT/core_output/00001/bdiv_${i}/euclidean_dm.txt

principal_coordinates.py \
            -i $MYCOUT/core_output/00001/bdiv_${i}/euclidean_dm.txt \
            -o $MYCOUT/core_output/00001/bdiv_${i}/euclidean_pc.txt

make_emperor.py \
            -i $MYCOUT/core_output/00001/bdiv_${i}/euclidean_pc.txt \
            -o $MYCOUT/core_output/00001/bdiv_${i}/euclidean_emperor_pcoa_plot/ \
            -m $MAP \
            --ignore_missing_samples

## Bray Curtis

beta_diversity.py \
                -i $MYCOUT/core_output/00001/taxa_plots/table_mc10000_sorted.biom \
                -o $MYCOUT/core_output/00001/bdiv_${i}/ \
                --metrics bray_curtis

mv $MYCOUT/core_output/00001/bdiv_${i}/bray_curtis_table_mc10000_sorted.txt \
                $MYCOUT/core_output/00001/bdiv_${i}/bray_curtis_dm.txt

principal_coordinates.py \
                -i $MYCOUT/core_output/00001/bdiv_${i}/bray_curtis_dm.txt \
                -o $MYCOUT/core_output/00001/bdiv_${i}/bray_curtis_pc.txt

make_emperor.py \
                -i $MYCOUT/core_output/00001/bdiv_${i}/bray_curtis_pc.txt \
                -o $MYCOUT/core_output/00001/bdiv_${i}/bray_curtis_emperor_pcoa_plot/ \
                -m $MAP \
                --ignore_missing_samples


###################################################################
# https://groups.google.com/forum/#!topic/qiime-forum/0TtSBh-0ovg #
###################################################################

###############################################################################
# ALPHA DIVERSITY
###############################################################################

mkdir $MYCOUT/core_output/00001/arare_${i}

cd $MYCOUT/core_output/00001

multiple_rarefactions.py \
            -i $MYCOUT/core_output/00001/table_mc10000.biom \
            -m 10 -x 10000 -s 24 \
            -o $MYCOUT/core_output/00001/arare_${i}/rarefaction/

alpha_diversity.py \
            -i $MYCOUT/core_output/00001/arare_${i}/rarefaction/ \
            -o $MYCOUT/core_output/00001/arare_${i}/alpha_div/ \
            --metrics chao1,observed_otus,goods_coverage,shannon

collate_alpha.py \
            -i $MYCOUT/core_output/00001/arare_${i}/alpha_div/ \
            -o $MYCOUT/core_output/00001/arare_${i}/alpha_div_collated/

rm -r $MYCOUT/core_output/00001/arare_${i}/rarefaction/ \
            $MYCOUT/core_output/00001/arare_${i}/alpha_div/

make_rarefaction_plots.py \
            -i $MYCOUT/core_output/00001/arare_${i}/alpha_div_collated/ \
            -m $MAP \
            -o $MYCOUT/core_output/00001/arare_${i}/alpha_rarefaction_plots/ \
            --generate_average_tables \
            --resolution 1000
            
done

### 0001 Abundance

cd $MYCOUT/core_output/0001
gunzip table_mc10000.biom.gz
gunzip table_even10000.biom.gz

biom summarize-table \
                -i table_mc10000.biom \
                -o biom_summary_table_mc10000.txt

biom summarize-table \
                -i table_even10000.biom \
                -o biom_summary_table_even10000.txt

###############################################################################
# BETA DIVERSITY
###############################################################################

for i in $(cat $WORKDIR/scripts/taxon.txt)

do

beta_diversity.py \
            -i $MYCOUT/core_output/0001/taxa_plots/table_mc10000_sorted.biom \
            -o $MYCOUT/core_output/0001/bdiv_${i}/ \
            --metrics euclidean

mv $MYCOUT/core_output/0001/bdiv_${i}/euclidean_table_mc10000_sorted.txt \
            $MYCOUT/core_output/0001/bdiv_${i}/euclidean_dm.txt

principal_coordinates.py \
            -i $MYCOUT/core_output/0001/bdiv_${i}/euclidean_dm.txt \
            -o $MYCOUT/core_output/0001/bdiv_${i}/euclidean_pc.txt

make_emperor.py \
            -i $MYCOUT/core_output/0001/bdiv_${i}/euclidean_pc.txt \
            -o $MYCOUT/core_output/0001/bdiv_${i}/euclidean_emperor_pcoa_plot/ \
            -m $MAP \
            --ignore_missing_samples

## Bray Curtis

beta_diversity.py \
                -i $MYCOUT/core_output/0001/taxa_plots/table_mc10000_sorted.biom \
                -o $MYCOUT/core_output/0001/bdiv_${i}/ \
                --metrics bray_curtis

mv $MYCOUT/core_output/0001/bdiv_${i}/bray_curtis_table_mc10000_sorted.txt \
                $MYCOUT/core_output/0001/bdiv_${i}/bray_curtis_dm.txt

principal_coordinates.py \
                -i $MYCOUT/core_output/0001/bdiv_${i}/bray_curtis_dm.txt \
                -o $MYCOUT/core_output/0001/bdiv_${i}/bray_curtis_pc.txt

make_emperor.py \
                -i $MYCOUT/core_output/0001/bdiv_${i}/bray_curtis_pc.txt \
                -o $MYCOUT/core_output/0001/bdiv_${i}/bray_curtis_emperor_pcoa_plot/ \
                -m $MAP \
                --ignore_missing_samples

###################################################################
# https://groups.google.com/forum/#!topic/qiime-forum/0TtSBh-0ovg #
###################################################################

###############################################################################
# ALPHA DIVERSITY
###############################################################################


mkdir $MYCOUT/core_output/0001/arare_${i}

cd $MYCOUT/core_output/0001

multiple_rarefactions.py \
            -i $MYCOUT/core_output/0001/table_mc10000.biom \
            -m 10 -x 10000 -s 24 \
            -o $MYCOUT/core_output/0001/arare_${i}/rarefaction/

alpha_diversity.py \
            -i $MYCOUT/core_output/0001/arare_${i}/rarefaction/ \
            -o $MYCOUT/core_output/0001/arare_${i}/alpha_div/ \
            --metrics chao1,observed_otus,goods_coverage,shannon

collate_alpha.py \
            -i $MYCOUT/core_output/0001/arare_${i}/alpha_div/ \
            -o $MYCOUT/core_output/0001/arare_${i}/alpha_div_collated/

rm -r $MYCOUT/core_output/0001/arare_${i}/rarefaction/ \
            $MYCOUT/core_output/0001/arare_${i}/alpha_div/

make_rarefaction_plots.py \
            -i $MYCOUT/core_output/0001/arare_${i}/alpha_div_collated/ \
            -m $MAP \
            -o $MYCOUT/core_output/0001/arare_${i}/alpha_rarefaction_plots/ \
            --generate_average_tables \
            --resolution 1000

done

### SINGLETON

cd $MYCOUT/core_output/singleton
gunzip table_mc10000.biom.gz
gunzip table_even10000.biom.gz

biom summarize-table \
                -i table_mc10000.biom \
                -o biom_summary_table_mc10000.txt

biom summarize-table \
                -i table_even10000.biom \
                -o biom_summary_table_even10000.txt

###############################################################################
# BETA DIVERSITY
###############################################################################

for i in $(cat $WORKDIR/scripts/taxon.txt)

do

beta_diversity.py \
            -i $MYCOUT/core_output/singleton/taxa_plots/table_mc10000_sorted.biom \
            -o $MYCOUT/core_output/singleton/bdiv_${i}/ \
            --metrics euclidean

mv $MYCOUT/core_output/singleton/bdiv_${i}/euclidean_table_mc10000_sorted.txt \
            $MYCOUT/core_output/singleton/bdiv_${i}/euclidean_dm.txt

principal_coordinates.py \
            -i $MYCOUT/core_output/singleton/bdiv_${i}/euclidean_dm.txt \
            -o $MYCOUT/core_output/singleton/bdiv_${i}/euclidean_pc.txt

make_emperor.py \
            -i $MYCOUT/core_output/singleton/bdiv_${i}/euclidean_pc.txt \
            -o $MYCOUT/core_output/singleton/bdiv_${i}/euclidean_emperor_pcoa_plot/ \
            -m $MAP \
            --ignore_missing_samples

## Bray Curtis

beta_diversity.py \
                -i $MYCOUT/core_output/singleton/taxa_plots/table_mc10000_sorted.biom \
                -o $MYCOUT/core_output/singleton/bdiv_${i}/ \
                --metrics bray_curtis

mv $MYCOUT/core_output/singleton/bdiv_${i}/bray_curtis_table_mc10000_sorted.txt \
                $MYCOUT/core_output/singleton/bdiv_${i}/bray_curtis_dm.txt

principal_coordinates.py \
                -i $MYCOUT/core_output/singleton/bdiv_${i}/bray_curtis_dm.txt \
                -o $MYCOUT/core_output/singleton/bdiv_${i}/bray_curtis_pc.txt

make_emperor.py \
                -i $MYCOUT/core_output/singleton/bdiv_${i}/bray_curtis_pc.txt \
                -o $MYCOUT/core_output/singleton/bdiv_${i}/bray_curtis_emperor_pcoa_plot/ \
                -m $MAP \
                --ignore_missing_samples


###################################################################
# https://groups.google.com/forum/#!topic/qiime-forum/0TtSBh-0ovg #
###################################################################

###############################################################################
# ALPHA DIVERSITY
###############################################################################

mkdir $MYCOUT/core_output/singleton/arare_${i}

cd $MYCOUT/core_output/singleton

multiple_rarefactions.py \
            -i $MYCOUT/core_output/singleton/table_mc10000.biom \
            -m 10 -x 10000 -s 24 \
            -o $MYCOUT/core_output/singleton/arare_${i}/rarefaction/

alpha_diversity.py \
            -i $MYCOUT/core_output/singleton/arare_${i}/rarefaction/ \
            -o $MYCOUT/core_output/singleton/arare_${i}/alpha_div/ \
            --metrics chao1,observed_otus,goods_coverage,shannon

collate_alpha.py \
            -i $MYCOUT/core_output/singleton/arare_${i}/alpha_div/ \
            -o $MYCOUT/core_output/singleton/arare_${i}/alpha_div_collated/

rm -r $MYCOUT/core_output/singleton/arare_${i}/rarefaction/ \
            $MYCOUT/core_output/singleton/arare_${i}/alpha_div/

make_rarefaction_plots.py \
            -i $MYCOUT/core_output/singleton/arare_${i}/alpha_div_collated/ \
            -m $MAP \
            -o $MYCOUT/core_output/singleton/arare_${i}/alpha_rarefaction_plots/ \
            --generate_average_tables \
            --resolution 1000
            
done

### CORE 0001 - 70%

cd $MYCOUT/core_output/core_0001_70p
gunzip table_mc10000.biom.gz
gunzip table_even10000.biom.gz

biom summarize-table \
                -i table_mc10000.biom \
                -o biom_summary_table_mc10000.txt

biom summarize-table \
                -i table_even10000.biom \
                -o biom_summary_table_even10000.txt

###############################################################################
# BETA DIVERSITY
###############################################################################

for i in $(cat $WORKDIR/scripts/taxon.txt)

do

beta_diversity.py \
            -i $MYCOUT/core_output/core_0001_70p/taxa_plots/table_mc10000_sorted.biom \
            -o $MYCOUT/core_output/core_0001_70p/bdiv_${i}/ \
            --metrics euclidean

mv $MYCOUT/core_output/core_0001_70p/bdiv_${i}/euclidean_table_mc10000_sorted.txt \
            $MYCOUT/core_output/core_0001_70p/bdiv_${i}/euclidean_dm.txt

principal_coordinates.py \
            -i $MYCOUT/core_output/core_0001_70p/bdiv_${i}/euclidean_dm.txt \
            -o $MYCOUT/core_output/core_0001_70p/bdiv_${i}/euclidean_pc.txt

make_emperor.py \
            -i $MYCOUT/core_output/core_0001_70p/bdiv_${i}/euclidean_pc.txt \
            -o $MYCOUT/core_output/core_0001_70p/bdiv_${i}/euclidean_emperor_pcoa_plot/ \
            -m $MAP \
            --ignore_missing_samples

## Bray Curtis

beta_diversity.py \
                -i $MYCOUT/core_output/core_0001_70p/taxa_plots/table_mc10000_sorted.biom \
                -o $MYCOUT/core_output/core_0001_70p/bdiv_${i}/ \
                --metrics bray_curtis

mv $MYCOUT/core_output/core_0001_70p/bdiv_${i}/bray_curtis_table_mc10000_sorted.txt \
                $MYCOUT/core_output/core_0001_70p/bdiv_${i}/bray_curtis_dm.txt

principal_coordinates.py \
                -i $MYCOUT/core_output/core_0001_70p/bdiv_${i}/bray_curtis_dm.txt \
                -o $MYCOUT/core_output/core_0001_70p/bdiv_${i}/bray_curtis_pc.txt

make_emperor.py \
                -i $MYCOUT/core_output/core_0001_70p/bdiv_${i}/bray_curtis_pc.txt \
                -o $MYCOUT/core_output/core_0001_70p/bdiv_${i}/bray_curtis_emperor_pcoa_plot/ \
                -m $MAP \
                --ignore_missing_samples


###################################################################
# https://groups.google.com/forum/#!topic/qiime-forum/0TtSBh-0ovg #
###################################################################

###############################################################################
# ALPHA DIVERSITY
###############################################################################

mkdir $MYCOUT/core_output/core_0001_70p/arare_${i}

cd $MYCOUT/core_output/core_0001_70p

multiple_rarefactions.py \
            -i $MYCOUT/core_output/core_0001_70p/table_mc10000.biom \
            -m 10 -x 10000 -s 24 \
            -o $MYCOUT/core_output/core_0001_70p/arare_${i}/rarefaction/

alpha_diversity.py \
            -i $MYCOUT/core_output/core_0001_70p/arare_${i}/rarefaction/ \
            -o $MYCOUT/core_output/core_0001_70p/arare_${i}/alpha_div/ \
            --metrics chao1,observed_otus,goods_coverage,shannon

collate_alpha.py \
            -i $MYCOUT/core_output/core_0001_70p/arare_${i}/alpha_div/ \
            -o $MYCOUT/core_output/core_0001_70p/arare_${i}/alpha_div_collated/

rm -r $MYCOUT/core_output/core_0001_70p/arare_${i}/rarefaction/ \
            $MYCOUT/core_output/core_0001_70p/arare_${i}/alpha_div/

make_rarefaction_plots.py \
            -i $MYCOUT/core_output/core_0001_70p/arare_${i}/alpha_div_collated/ \
            -m $MAP \
            -o $MYCOUT/core_output/core_0001_70p/arare_${i}/alpha_rarefaction_plots/ \
            --generate_average_tables \
            --resolution 1000
            
done

### CORE 0001 - 80%

cd $MYCOUT/core_output/core_0001_80p
gunzip table_mc10000.biom.gz
gunzip table_even10000.biom.gz

biom summarize-table \
                -i table_mc10000.biom \
                -o biom_summary_table_mc10000.txt

biom summarize-table \
                -i table_even10000.biom \
                -o biom_summary_table_even10000.txt

###############################################################################
# BETA DIVERSITY
###############################################################################

for i in $(cat $WORKDIR/scripts/taxon.txt)

do

beta_diversity.py \
            -i $MYCOUT/core_output/core_0001_80p/taxa_plots/table_mc10000_sorted.biom \
            -o $MYCOUT/core_output/core_0001_80p/bdiv_${i}/ \
            --metrics euclidean

mv $MYCOUT/core_output/core_0001_80p/bdiv_${i}/euclidean_table_mc10000_sorted.txt \
            $MYCOUT/core_output/core_0001_80p/bdiv_${i}/euclidean_dm.txt

principal_coordinates.py \
            -i $MYCOUT/core_output/core_0001_80p/bdiv_${i}/euclidean_dm.txt \
            -o $MYCOUT/core_output/core_0001_80p/bdiv_${i}/euclidean_pc.txt

make_emperor.py \
            -i $MYCOUT/core_output/core_0001_80p/bdiv_${i}/euclidean_pc.txt \
            -o $MYCOUT/core_output/core_0001_80p/bdiv_${i}/euclidean_emperor_pcoa_plot/ \
            -m $MAP \
            --ignore_missing_samples

## Bray Curtis

beta_diversity.py \
                -i $MYCOUT/core_output/core_0001_80p/taxa_plots/table_mc10000_sorted.biom \
                -o $MYCOUT/core_output/core_0001_80p/bdiv_${i}/ \
                --metrics bray_curtis

mv $MYCOUT/core_output/core_0001_80p/bdiv_${i}/bray_curtis_table_mc10000_sorted.txt \
                $MYCOUT/core_output/core_0001_80p/bdiv_${i}/bray_curtis_dm.txt

principal_coordinates.py \
                -i $MYCOUT/core_output/core_0001_80p/bdiv_${i}/bray_curtis_dm.txt \
                -o $MYCOUT/core_output/core_0001_80p/bdiv_${i}/bray_curtis_pc.txt

make_emperor.py \
                -i $MYCOUT/core_output/core_0001_80p/bdiv_${i}/bray_curtis_pc.txt \
                -o $MYCOUT/core_output/core_0001_80p/bdiv_${i}/bray_curtis_emperor_pcoa_plot/ \
                -m $MAP \
                --ignore_missing_samples


###################################################################
# https://groups.google.com/forum/#!topic/qiime-forum/0TtSBh-0ovg #
###################################################################

###############################################################################
# ALPHA DIVERSITY
###############################################################################

mkdir $MYCOUT/core_output/core_0001_80p/arare_${i}

cd $MYCOUT/core_output/core_0001_80p

multiple_rarefactions.py \
            -i $MYCOUT/core_output/core_0001_80p/table_mc10000.biom \
            -m 10 -x 10000 -s 24 \
            -o $MYCOUT/core_output/core_0001_80p/arare_${i}/rarefaction/

alpha_diversity.py \
            -i $MYCOUT/core_output/core_0001_80p/arare_${i}/rarefaction/ \
            -o $MYCOUT/core_output/core_0001_80p/arare_${i}/alpha_div/ \
            --metrics chao1,observed_otus,goods_coverage,shannon

collate_alpha.py \
            -i $MYCOUT/core_output/core_0001_80p/arare_${i}/alpha_div/ \
            -o $MYCOUT/core_output/core_0001_80p/arare_${i}/alpha_div_collated/

rm -r $MYCOUT/core_output/core_0001_80p/arare_${i}/rarefaction/ \
            $MYCOUT/core_output/core_0001_80p/arare_${i}/alpha_div/

make_rarefaction_plots.py \
            -i $MYCOUT/core_output/core_0001_80p/arare_${i}/alpha_div_collated/ \
            -m $MAP \
            -o $MYCOUT/core_output/core_0001_80p/arare_${i}/alpha_rarefaction_plots/ \
            --generate_average_tables \
            --resolution 1000
            
done

### CORE 0001 - 90%

cd $MYCOUT/core_output/core_0001_90p
gunzip table_mc10000.biom.gz
gunzip table_even10000.biom.gz

biom summarize-table \
                -i table_mc10000.biom \
                -o biom_summary_table_mc10000.txt

biom summarize-table \
                -i table_even10000.biom \
                -o biom_summary_table_even10000.txt

###############################################################################
# BETA DIVERSITY
###############################################################################

for i in $(cat $WORKDIR/scripts/taxon.txt)

do

beta_diversity.py \
            -i $MYCOUT/core_output/core_0001_90p/taxa_plots/table_mc10000_sorted.biom \
            -o $MYCOUT/core_output/core_0001_90p/bdiv_${i}/ \
            --metrics euclidean

mv $MYCOUT/core_output/core_0001_90p/bdiv_${i}/euclidean_table_mc10000_sorted.txt \
            $MYCOUT/core_output/core_0001_90p/bdiv_${i}/euclidean_dm.txt

principal_coordinates.py \
            -i $MYCOUT/core_output/core_0001_90p/bdiv_${i}/euclidean_dm.txt \
            -o $MYCOUT/core_output/core_0001_90p/bdiv_${i}/euclidean_pc.txt

make_emperor.py \
            -i $MYCOUT/core_output/core_0001_90p/bdiv_${i}/euclidean_pc.txt \
            -o $MYCOUT/core_output/core_0001_90p/bdiv_${i}/euclidean_emperor_pcoa_plot/ \
            -m $MAP \
            --ignore_missing_samples

## Bray Curtis

beta_diversity.py \
                -i $MYCOUT/core_output/core_0001_90p/taxa_plots/table_mc10000_sorted.biom \
                -o $MYCOUT/core_output/core_0001_90p/bdiv_${i}/ \
                --metrics bray_curtis

mv $MYCOUT/core_output/core_0001_90p/bdiv_${i}/bray_curtis_table_mc10000_sorted.txt \
                $MYCOUT/core_output/core_0001_90p/bdiv_${i}/bray_curtis_dm.txt

principal_coordinates.py \
                -i $MYCOUT/core_output/core_0001_90p/bdiv_${i}/bray_curtis_dm.txt \
                -o $MYCOUT/core_output/core_0001_90p/bdiv_${i}/bray_curtis_pc.txt

make_emperor.py \
                -i $MYCOUT/core_output/core_0001_90p/bdiv_${i}/bray_curtis_pc.txt \
                -o $MYCOUT/core_output/core_0001_90p/bdiv_${i}/bray_curtis_emperor_pcoa_plot/ \
                -m $MAP \
                --ignore_missing_samples


###################################################################
# https://groups.google.com/forum/#!topic/qiime-forum/0TtSBh-0ovg #
###################################################################

###############################################################################
# ALPHA DIVERSITY
###############################################################################

mkdir $MYCOUT/core_output/core_0001_90p/arare_${i}

cd $MYCOUT/core_output/core_0001_90p

multiple_rarefactions.py \
            -i $MYCOUT/core_output/core_0001_90p/table_mc10000.biom \
            -m 10 -x 10000 -s 24 \
            -o $MYCOUT/core_output/core_0001_90p/arare_${i}/rarefaction/

alpha_diversity.py \
            -i $MYCOUT/core_output/core_0001_90p/arare_${i}/rarefaction/ \
            -o $MYCOUT/core_output/core_0001_90p/arare_${i}/alpha_div/ \
            --metrics chao1,observed_otus,goods_coverage,shannon

collate_alpha.py \
            -i $MYCOUT/core_output/core_0001_90p/arare_${i}/alpha_div/ \
            -o $MYCOUT/core_output/core_0001_90p/arare_${i}/alpha_div_collated/

rm -r $MYCOUT/core_output/core_0001_90p/arare_${i}/rarefaction/ \
            $MYCOUT/core_output/core_0001_90p/arare_${i}/alpha_div/

make_rarefaction_plots.py \
            -i $MYCOUT/core_output/core_0001_90p/arare_${i}/alpha_div_collated/ \
            -m $MAP \
            -o $MYCOUT/core_output/core_0001_90p/arare_${i}/alpha_rarefaction_plots/ \
            --generate_average_tables \
            --resolution 1000
            
done

### CORE 00001 - 70%

cd $MYCOUT/core_output/core_00001_70p
gunzip table_mc10000.biom.gz
gunzip table_even10000.biom.gz

biom summarize-table \
                -i table_mc10000.biom \
                -o biom_summary_table_mc10000.txt

biom summarize-table \
                -i table_even10000.biom \
                -o biom_summary_table_even10000.txt

###############################################################################
# BETA DIVERSITY
###############################################################################

for i in $(cat $WORKDIR/scripts/taxon.txt)

do

beta_diversity.py \
            -i $MYCOUT/core_output/core_00001_70p/taxa_plots/table_mc10000_sorted.biom \
            -o $MYCOUT/core_output/core_00001_70p/bdiv_${i}/ \
            --metrics euclidean

mv $MYCOUT/core_output/core_00001_70p/bdiv_${i}/euclidean_table_mc10000_sorted.txt \
            $MYCOUT/core_output/core_00001_70p/bdiv_${i}/euclidean_dm.txt

principal_coordinates.py \
            -i $MYCOUT/core_output/core_00001_70p/bdiv_${i}/euclidean_dm.txt \
            -o $MYCOUT/core_output/core_00001_70p/bdiv_${i}/euclidean_pc.txt

make_emperor.py \
            -i $MYCOUT/core_output/core_00001_70p/bdiv_${i}/euclidean_pc.txt \
            -o $MYCOUT/core_output/core_00001_70p/bdiv_${i}/euclidean_emperor_pcoa_plot/ \
            -m $MAP \
            --ignore_missing_samples

## Bray Curtis

beta_diversity.py \
                -i $MYCOUT/core_output/core_00001_70p/taxa_plots/table_mc10000_sorted.biom \
                -o $MYCOUT/core_output/core_00001_70p/bdiv_${i}/ \
                --metrics bray_curtis

mv $MYCOUT/core_output/core_00001_70p/bdiv_${i}/bray_curtis_table_mc10000_sorted.txt \
                $MYCOUT/core_output/core_00001_70p/bdiv_${i}/bray_curtis_dm.txt

principal_coordinates.py \
                -i $MYCOUT/core_output/core_00001_70p/bdiv_${i}/bray_curtis_dm.txt \
                -o $MYCOUT/core_output/core_00001_70p/bdiv_${i}/bray_curtis_pc.txt

make_emperor.py \
                -i $MYCOUT/core_output/core_00001_70p/bdiv_${i}/bray_curtis_pc.txt \
                -o $MYCOUT/core_output/core_00001_70p/bdiv_${i}/bray_curtis_emperor_pcoa_plot/ \
                -m $MAP \
                --ignore_missing_samples


###################################################################
# https://groups.google.com/forum/#!topic/qiime-forum/0TtSBh-0ovg #
###################################################################

###############################################################################
# ALPHA DIVERSITY
###############################################################################

mkdir $MYCOUT/core_output/core_00001_70p/arare_${i}

cd $MYCOUT/core_output/core_00001_70p

multiple_rarefactions.py \
            -i $MYCOUT/core_output/core_00001_70p/table_mc10000.biom \
            -m 10 -x 10000 -s 24 \
            -o $MYCOUT/core_output/core_00001_70p/arare_${i}/rarefaction/

alpha_diversity.py \
            -i $MYCOUT/core_output/core_00001_70p/arare_${i}/rarefaction/ \
            -o $MYCOUT/core_output/core_00001_70p/arare_${i}/alpha_div/ \
            --metrics chao1,observed_otus,goods_coverage,shannon

collate_alpha.py \
            -i $MYCOUT/core_output/core_00001_70p/arare_${i}/alpha_div/ \
            -o $MYCOUT/core_output/core_00001_70p/arare_${i}/alpha_div_collated/

rm -r $MYCOUT/core_output/core_00001_70p/arare_${i}/rarefaction/ \
            $MYCOUT/core_output/core_00001_70p/arare_${i}/alpha_div/

make_rarefaction_plots.py \
            -i $MYCOUT/core_output/core_00001_70p/arare_${i}/alpha_div_collated/ \
            -m $MAP \
            -o $MYCOUT/core_output/core_00001_70p/arare_${i}/alpha_rarefaction_plots/ \
            --generate_average_tables \
            --resolution 1000
            
done

### CORE 00001 - 80p

cd $MYCOUT/core_output/core_00001_80p
gunzip table_mc10000.biom.gz
gunzip table_even10000.biom.gz

biom summarize-table \
                -i table_mc10000.biom \
                -o biom_summary_table_mc10000.txt

biom summarize-table \
                -i table_even10000.biom \
                -o biom_summary_table_even10000.txt

###############################################################################
# BETA DIVERSITY
###############################################################################

for i in $(cat $WORKDIR/scripts/taxon.txt)

do

beta_diversity.py \
            -i $MYCOUT/core_output/core_00001_80p/taxa_plots/table_mc10000_sorted.biom \
            -o $MYCOUT/core_output/core_00001_80p/bdiv_${i}/ \
            --metrics euclidean

mv $MYCOUT/core_output/core_00001_80p/bdiv_${i}/euclidean_table_mc10000_sorted.txt \
            $MYCOUT/core_output/core_00001_80p/bdiv_${i}/euclidean_dm.txt

principal_coordinates.py \
            -i $MYCOUT/core_output/core_00001_80p/bdiv_${i}/euclidean_dm.txt \
            -o $MYCOUT/core_output/core_00001_80p/bdiv_${i}/euclidean_pc.txt

make_emperor.py \
            -i $MYCOUT/core_output/core_00001_80p/bdiv_${i}/euclidean_pc.txt \
            -o $MYCOUT/core_output/core_00001_80p/bdiv_${i}/euclidean_emperor_pcoa_plot/ \
            -m $MAP \
            --ignore_missing_samples

## Bray Curtis

beta_diversity.py \
                -i $MYCOUT/core_output/core_00001_80p/taxa_plots/table_mc10000_sorted.biom \
                -o $MYCOUT/core_output/core_00001_80p/bdiv_${i}/ \
                --metrics bray_curtis

mv $MYCOUT/core_output/core_00001_80p/bdiv_${i}/bray_curtis_table_mc10000_sorted.txt \
                $MYCOUT/core_output/core_00001_80p/bdiv_${i}/bray_curtis_dm.txt

principal_coordinates.py \
                -i $MYCOUT/core_output/core_00001_80p/bdiv_${i}/bray_curtis_dm.txt \
                -o $MYCOUT/core_output/core_00001_80p/bdiv_${i}/bray_curtis_pc.txt

make_emperor.py \
                -i $MYCOUT/core_output/core_00001_80p/bdiv_${i}/bray_curtis_pc.txt \
                -o $MYCOUT/core_output/core_00001_80p/bdiv_${i}/bray_curtis_emperor_pcoa_plot/ \
                -m $MAP \
                --ignore_missing_samples


###################################################################
# https://groups.google.com/forum/#!topic/qiime-forum/0TtSBh-0ovg #
###################################################################

###############################################################################
# ALPHA DIVERSITY
###############################################################################

mkdir $MYCOUT/core_output/core_00001_80p/arare_${i}

cd $MYCOUT/core_output/core_00001_80p

multiple_rarefactions.py \
            -i $MYCOUT/core_output/core_00001_80p/table_mc10000.biom \
            -m 10 -x 10000 -s 24 \
            -o $MYCOUT/core_output/core_00001_80p/arare_${i}/rarefaction/

alpha_diversity.py \
            -i $MYCOUT/core_output/core_00001_80p/arare_${i}/rarefaction/ \
            -o $MYCOUT/core_output/core_00001_80p/arare_${i}/alpha_div/ \
            --metrics chao1,observed_otus,goods_coverage,shannon

collate_alpha.py \
            -i $MYCOUT/core_output/core_00001_80p/arare_${i}/alpha_div/ \
            -o $MYCOUT/core_output/core_00001_80p/arare_${i}/alpha_div_collated/

rm -r $MYCOUT/core_output/core_00001_80p/arare_${i}/rarefaction/ \
            $MYCOUT/core_output/core_00001_80p/arare_${i}/alpha_div/

make_rarefaction_plots.py \
            -i $MYCOUT/core_output/core_00001_80p/arare_${i}/alpha_div_collated/ \
            -m $MAP \
            -o $MYCOUT/core_output/core_00001_80p/arare_${i}/alpha_rarefaction_plots/ \
            --generate_average_tables \
            --resolution 1000
            
done

### CORE 00001 - 90p

cd $MYCOUT/core_output/core_00001_90p
gunzip table_mc10000.biom.gz
gunzip table_even10000.biom.gz

biom summarize-table \
                -i table_mc10000.biom \
                -o biom_summary_table_mc10000.txt

biom summarize-table \
                -i table_even10000.biom \
                -o biom_summary_table_even10000.txt

###############################################################################
# BETA DIVERSITY
###############################################################################

for i in $(cat $WORKDIR/scripts/taxon.txt)

do

beta_diversity.py \
            -i $MYCOUT/core_output/core_00001_90p/taxa_plots/table_mc10000_sorted.biom \
            -o $MYCOUT/core_output/core_00001_90p/bdiv_${i}/ \
            --metrics euclidean

mv $MYCOUT/core_output/core_00001_90p/bdiv_${i}/euclidean_table_mc10000_sorted.txt \
            $MYCOUT/core_output/core_00001_90p/bdiv_${i}/euclidean_dm.txt

principal_coordinates.py \
            -i $MYCOUT/core_output/core_00001_90p/bdiv_${i}/euclidean_dm.txt \
            -o $MYCOUT/core_output/core_00001_90p/bdiv_${i}/euclidean_pc.txt

make_emperor.py \
            -i $MYCOUT/core_output/core_00001_90p/bdiv_${i}/euclidean_pc.txt \
            -o $MYCOUT/core_output/core_00001_90p/bdiv_${i}/euclidean_emperor_pcoa_plot/ \
            -m $MAP \
            --ignore_missing_samples

## Bray Curtis

beta_diversity.py \
                -i $MYCOUT/core_output/core_00001_90p/taxa_plots/table_mc10000_sorted.biom \
                -o $MYCOUT/core_output/core_00001_90p/bdiv_${i}/ \
                --metrics bray_curtis

mv $MYCOUT/core_output/core_00001_90p/bdiv_${i}/bray_curtis_table_mc10000_sorted.txt \
                $MYCOUT/core_output/core_00001_90p/bdiv_${i}/bray_curtis_dm.txt

principal_coordinates.py \
                -i $MYCOUT/core_output/core_00001_90p/bdiv_${i}/bray_curtis_dm.txt \
                -o $MYCOUT/core_output/core_00001_90p/bdiv_${i}/bray_curtis_pc.txt

make_emperor.py \
                -i $MYCOUT/core_output/core_00001_90p/bdiv_${i}/bray_curtis_pc.txt \
                -o $MYCOUT/core_output/core_00001_90p/bdiv_${i}/bray_curtis_emperor_pcoa_plot/ \
                -m $MAP \
                --ignore_missing_samples


###################################################################
# https://groups.google.com/forum/#!topic/qiime-forum/0TtSBh-0ovg #
###################################################################

###############################################################################
# ALPHA DIVERSITY
###############################################################################

mkdir $MYCOUT/core_output/core_00001_90p/arare_${i}

cd $MYCOUT/core_output/core_00001_90p

multiple_rarefactions.py \
            -i $MYCOUT/core_output/core_00001_90p/table_mc10000.biom \
            -m 10 -x 10000 -s 24 \
            -o $MYCOUT/core_output/core_00001_90p/arare_${i}/rarefaction/

alpha_diversity.py \
            -i $MYCOUT/core_output/core_00001_90p/arare_${i}/rarefaction/ \
            -o $MYCOUT/core_output/core_00001_90p/arare_${i}/alpha_div/ \
            --metrics chao1,observed_otus,goods_coverage,shannon

collate_alpha.py \
            -i $MYCOUT/core_output/core_00001_90p/arare_${i}/alpha_div/ \
            -o $MYCOUT/core_output/core_00001_90p/arare_${i}/alpha_div_collated/

rm -r $MYCOUT/core_output/core_00001_90p/arare_${i}/rarefaction/ \
            $MYCOUT/core_output/core_00001_90p/arare_${i}/alpha_div/

make_rarefaction_plots.py \
            -i $MYCOUT/core_output/core_00001_90p/arare_${i}/alpha_div_collated/ \
            -m $MAP \
            -o $MYCOUT/core_output/core_00001_90p/arare_${i}/alpha_rarefaction_plots/ \
            --generate_average_tables \
            --resolution 1000
            
done


### Core Singleton

### CORE singleton - 70%

cd $MYCOUT/core_output/core_singleton_70p
gunzip table_mc10000.biom.gz
gunzip table_even10000.biom.gz

biom summarize-table \
                -i table_mc10000.biom \
                -o biom_summary_table_mc10000.txt

biom summarize-table \
                -i table_even10000.biom \
                -o biom_summary_table_even10000.txt

###############################################################################
# BETA DIVERSITY
###############################################################################

for i in $(cat $WORKDIR/scripts/taxon.txt)

do

beta_diversity.py \
            -i $MYCOUT/core_output/core_singleton_70p/taxa_plots/table_mc10000_sorted.biom \
            -o $MYCOUT/core_output/core_singleton_70p/bdiv_${i}/ \
            --metrics euclidean

mv $MYCOUT/core_output/core_singleton_70p/bdiv_${i}/euclidean_table_mc10000_sorted.txt \
            $MYCOUT/core_output/core_singleton_70p/bdiv_${i}/euclidean_dm.txt

principal_coordinates.py \
            -i $MYCOUT/core_output/core_singleton_70p/bdiv_${i}/euclidean_dm.txt \
            -o $MYCOUT/core_output/core_singleton_70p/bdiv_${i}/euclidean_pc.txt

make_emperor.py \
            -i $MYCOUT/core_output/core_singleton_70p/bdiv_${i}/euclidean_pc.txt \
            -o $MYCOUT/core_output/core_singleton_70p/bdiv_${i}/euclidean_emperor_pcoa_plot/ \
            -m $MAP \
            --ignore_missing_samples

## Bray Curtis

beta_diversity.py \
                -i $MYCOUT/core_output/core_singleton_70p/taxa_plots/table_mc10000_sorted.biom \
                -o $MYCOUT/core_output/core_singleton_70p/bdiv_${i}/ \
                --metrics bray_curtis

mv $MYCOUT/core_output/core_singleton_70p/bdiv_${i}/bray_curtis_table_mc10000_sorted.txt \
                $MYCOUT/core_output/core_singleton_70p/bdiv_${i}/bray_curtis_dm.txt

principal_coordinates.py \
                -i $MYCOUT/core_output/core_singleton_70p/bdiv_${i}/bray_curtis_dm.txt \
                -o $MYCOUT/core_output/core_singleton_70p/bdiv_${i}/bray_curtis_pc.txt

make_emperor.py \
                -i $MYCOUT/core_output/core_singleton_70p/bdiv_${i}/bray_curtis_pc.txt \
                -o $MYCOUT/core_output/core_singleton_70p/bdiv_${i}/bray_curtis_emperor_pcoa_plot/ \
                -m $MAP \
                --ignore_missing_samples


###################################################################
# https://groups.google.com/forum/#!topic/qiime-forum/0TtSBh-0ovg #
###################################################################

###############################################################################
# ALPHA DIVERSITY
###############################################################################

mkdir $MYCOUT/core_output/core_singleton_70p/arare_${i}

cd $MYCOUT/core_output/core_singleton_70p

multiple_rarefactions.py \
            -i $MYCOUT/core_output/core_singleton_70p/table_mc10000.biom \
            -m 10 -x 10000 -s 24 \
            -o $MYCOUT/core_output/core_singleton_70p/arare_${i}/rarefaction/

alpha_diversity.py \
            -i $MYCOUT/core_output/core_singleton_70p/arare_${i}/rarefaction/ \
            -o $MYCOUT/core_output/core_singleton_70p/arare_${i}/alpha_div/ \
            --metrics chao1,observed_otus,goods_coverage,shannon

collate_alpha.py \
            -i $MYCOUT/core_output/core_singleton_70p/arare_${i}/alpha_div/ \
            -o $MYCOUT/core_output/core_singleton_70p/arare_${i}/alpha_div_collated/

rm -r $MYCOUT/core_output/core_singleton_70p/arare_${i}/rarefaction/ \
            $MYCOUT/core_output/core_singleton_70p/arare_${i}/alpha_div/

make_rarefaction_plots.py \
            -i $MYCOUT/core_output/core_singleton_70p/arare_${i}/alpha_div_collated/ \
            -m $MAP \
            -o $MYCOUT/core_output/core_singleton_70p/arare_${i}/alpha_rarefaction_plots/ \
            --generate_average_tables \
            --resolution 1000
            
done

### CORE singleton - 80p

cd $MYCOUT/core_output/core_singleton_80p
gunzip table_mc10000.biom.gz
gunzip table_even10000.biom.gz

biom summarize-table \
                -i table_mc10000.biom \
                -o biom_summary_table_mc10000.txt

biom summarize-table \
                -i table_even10000.biom \
                -o biom_summary_table_even10000.txt

###############################################################################
# BETA DIVERSITY
###############################################################################

for i in $(cat $WORKDIR/scripts/taxon.txt)

do

beta_diversity.py \
            -i $MYCOUT/core_output/core_singleton_80p/taxa_plots/table_mc10000_sorted.biom \
            -o $MYCOUT/core_output/core_singleton_80p/bdiv_${i}/ \
            --metrics euclidean

mv $MYCOUT/core_output/core_singleton_80p/bdiv_${i}/euclidean_table_mc10000_sorted.txt \
            $MYCOUT/core_output/core_singleton_80p/bdiv_${i}/euclidean_dm.txt

principal_coordinates.py \
            -i $MYCOUT/core_output/core_singleton_80p/bdiv_${i}/euclidean_dm.txt \
            -o $MYCOUT/core_output/core_singleton_80p/bdiv_${i}/euclidean_pc.txt

make_emperor.py \
            -i $MYCOUT/core_output/core_singleton_80p/bdiv_${i}/euclidean_pc.txt \
            -o $MYCOUT/core_output/core_singleton_80p/bdiv_${i}/euclidean_emperor_pcoa_plot/ \
            -m $MAP \
            --ignore_missing_samples

## Bray Curtis

beta_diversity.py \
                -i $MYCOUT/core_output/core_singleton_80p/taxa_plots/table_mc10000_sorted.biom \
                -o $MYCOUT/core_output/core_singleton_80p/bdiv_${i}/ \
                --metrics bray_curtis

mv $MYCOUT/core_output/core_singleton_80p/bdiv_${i}/bray_curtis_table_mc10000_sorted.txt \
                $MYCOUT/core_output/core_singleton_80p/bdiv_${i}/bray_curtis_dm.txt

principal_coordinates.py \
                -i $MYCOUT/core_output/core_singleton_80p/bdiv_${i}/bray_curtis_dm.txt \
                -o $MYCOUT/core_output/core_singleton_80p/bdiv_${i}/bray_curtis_pc.txt

make_emperor.py \
                -i $MYCOUT/core_output/core_singleton_80p/bdiv_${i}/bray_curtis_pc.txt \
                -o $MYCOUT/core_output/core_singleton_80p/bdiv_${i}/bray_curtis_emperor_pcoa_plot/ \
                -m $MAP \
                --ignore_missing_samples


###################################################################
# https://groups.google.com/forum/#!topic/qiime-forum/0TtSBh-0ovg #
###################################################################

###############################################################################
# ALPHA DIVERSITY
###############################################################################

mkdir $MYCOUT/core_output/core_singleton_80p/arare_${i}

cd $MYCOUT/core_output/core_singleton_80p

multiple_rarefactions.py \
            -i $MYCOUT/core_output/core_singleton_80p/table_mc10000.biom \
            -m 10 -x 10000 -s 24 \
            -o $MYCOUT/core_output/core_singleton_80p/arare_${i}/rarefaction/

alpha_diversity.py \
            -i $MYCOUT/core_output/core_singleton_80p/arare_${i}/rarefaction/ \
            -o $MYCOUT/core_output/core_singleton_80p/arare_${i}/alpha_div/ \
            --metrics chao1,observed_otus,goods_coverage,shannon

collate_alpha.py \
            -i $MYCOUT/core_output/core_singleton_80p/arare_${i}/alpha_div/ \
            -o $MYCOUT/core_output/core_singleton_80p/arare_${i}/alpha_div_collated/

rm -r $MYCOUT/core_output/core_singleton_80p/arare_${i}/rarefaction/ \
            $MYCOUT/core_output/core_singleton_80p/arare_${i}/alpha_div/

make_rarefaction_plots.py \
            -i $MYCOUT/core_output/core_singleton_80p/arare_${i}/alpha_div_collated/ \
            -m $MAP \
            -o $MYCOUT/core_output/core_singleton_80p/arare_${i}/alpha_rarefaction_plots/ \
            --generate_average_tables \
            --resolution 1000
            
done

### CORE singleton - 90p

cd $MYCOUT/core_output/core_singleton_90p
gunzip table_mc10000.biom.gz
gunzip table_even10000.biom.gz

biom summarize-table \
                -i table_mc10000.biom \
                -o biom_summary_table_mc10000.txt

biom summarize-table \
                -i table_even10000.biom \
                -o biom_summary_table_even10000.txt

###############################################################################
# BETA DIVERSITY
###############################################################################

for i in $(cat $WORKDIR/scripts/taxon.txt)

do

beta_diversity.py \
            -i $MYCOUT/core_output/core_singleton_90p/taxa_plots/table_mc10000_sorted.biom \
            -o $MYCOUT/core_output/core_singleton_90p/bdiv_${i}/ \
            --metrics euclidean

mv $MYCOUT/core_output/core_singleton_90p/bdiv_${i}/euclidean_table_mc10000_sorted.txt \
            $MYCOUT/core_output/core_singleton_90p/bdiv_${i}/euclidean_dm.txt

principal_coordinates.py \
            -i $MYCOUT/core_output/core_singleton_90p/bdiv_${i}/euclidean_dm.txt \
            -o $MYCOUT/core_output/core_singleton_90p/bdiv_${i}/euclidean_pc.txt

make_emperor.py \
            -i $MYCOUT/core_output/core_singleton_90p/bdiv_${i}/euclidean_pc.txt \
            -o $MYCOUT/core_output/core_singleton_90p/bdiv_${i}/euclidean_emperor_pcoa_plot/ \
            -m $MAP \
            --ignore_missing_samples

## Bray Curtis

beta_diversity.py \
                -i $MYCOUT/core_output/core_singleton_90p/taxa_plots/table_mc10000_sorted.biom \
                -o $MYCOUT/core_output/core_singleton_90p/bdiv_${i}/ \
                --metrics bray_curtis

mv $MYCOUT/core_output/core_singleton_90p/bdiv_${i}/bray_curtis_table_mc10000_sorted.txt \
                $MYCOUT/core_output/core_singleton_90p/bdiv_${i}/bray_curtis_dm.txt

principal_coordinates.py \
                -i $MYCOUT/core_output/core_singleton_90p/bdiv_${i}/bray_curtis_dm.txt \
                -o $MYCOUT/core_output/core_singleton_90p/bdiv_${i}/bray_curtis_pc.txt

make_emperor.py \
                -i $MYCOUT/core_output/core_singleton_90p/bdiv_${i}/bray_curtis_pc.txt \
                -o $MYCOUT/core_output/core_singleton_90p/bdiv_${i}/bray_curtis_emperor_pcoa_plot/ \
                -m $MAP \
                --ignore_missing_samples


###################################################################
# https://groups.google.com/forum/#!topic/qiime-forum/0TtSBh-0ovg #
###################################################################

###############################################################################
# ALPHA DIVERSITY
###############################################################################

mkdir $MYCOUT/core_output/core_singleton_90p/arare_${i}

cd $MYCOUT/core_output/core_singleton_90p

multiple_rarefactions.py \
            -i $MYCOUT/core_output/core_singleton_90p/table_mc10000.biom \
            -m 10 -x 10000 -s 24 \
            -o $MYCOUT/core_output/core_singleton_90p/arare_${i}/rarefaction/

alpha_diversity.py \
            -i $MYCOUT/core_output/core_singleton_90p/arare_${i}/rarefaction/ \
            -o $MYCOUT/core_output/core_singleton_90p/arare_${i}/alpha_div/ \
            --metrics chao1,observed_otus,goods_coverage,shannon

collate_alpha.py \
            -i $MYCOUT/core_output/core_singleton_90p/arare_${i}/alpha_div/ \
            -o $MYCOUT/core_output/core_singleton_90p/arare_${i}/alpha_div_collated/

rm -r $MYCOUT/core_output/core_singleton_90p/arare_${i}/rarefaction/ \
            $MYCOUT/core_output/core_singleton_90p/arare_${i}/alpha_div/

make_rarefaction_plots.py \
            -i $MYCOUT/core_output/core_singleton_90p/arare_${i}/alpha_div_collated/ \
            -m $MAP \
            -o $MYCOUT/core_output/core_singleton_90p/arare_${i}/alpha_rarefaction_plots/ \
            --generate_average_tables \
            --resolution 1000
            
done


echo "###############################"
echo "###############################"
echo "COMPARE CATEGORIES..."
echo "###############################"
echo "###############################"


echo " "
echo " "

mkdir $MYCOUT/core_output/compare_categories

for i in $(cat $WORKDIR/scripts/taxon.txt)

do


### 00001 Abundance

compare_categories.py \
            --method adonis \
            -i $MYCOUT/core_output/00001/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/adonis_out/00001/euclidean/ \
            -n 999


compare_categories.py \
            --method adonis \
            -i $MYCOUT/core_output/00001/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/adonis_out/00001/bray_curtis/ \
            -n 999


compare_categories.py \
            --method anosim \
            -i $MYCOUT/core_output/00001/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/anosim_out/00001/euclidean/ \
            -n 999


compare_categories.py \
            --method anosim \
            -i $MYCOUT/core_output/00001/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/anosim_out/00001/bray_curtis/ \
            -n 999

compare_categories.py \
            --method permanova \
            -i $MYCOUT/core_output/00001/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permanova_out/00001/euclidean/ \
            -n 999


compare_categories.py \
            --method permanova \
            -i $MYCOUT/core_output/00001/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permanova_out/00001/bray_curtis/ \
            -n 999


compare_categories.py \
            --method permdisp \
            -i $MYCOUT/core_output/00001/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permdisp_out/00001/euclidean/ \
            -n 999


compare_categories.py \
            --method permdisp \
            -i $MYCOUT/core_output/00001/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permdisp_out/00001/bray_curtis/ \
            -n 999


compare_categories.py \
            --method dbrda \
            -i $MYCOUT/core_output/00001/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/dbrda_out/00001/euclidean/ \
            -n 999


compare_categories.py \
            --method dbrda \
            -i $MYCOUT/core_output/00001/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/dbrda_out/00001/bray_curtis/ \
            -n 999


compare_categories.py \
            --method mrpp \
            -i $MYCOUT/core_output/00001/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/mrpp_out/00001/euclidean/ \
            -n 999


compare_categories.py \
            --method mrpp \
            -i $MYCOUT/core_output/00001/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/mrpp_out/00001/bray_curtis/ \
            -n 999

#### 0001 Abundance

compare_categories.py \
            --method adonis \
            -i $MYCOUT/core_output/0001/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/adonis_out/0001/euclidean/ \
            -n 999


compare_categories.py \
            --method adonis \
            -i $MYCOUT/core_output/0001/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/adonis_out/0001/bray_curtis/ \
            -n 999


compare_categories.py \
            --method anosim \
            -i $MYCOUT/core_output/0001/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/anosim_out/0001/euclidean/ \
            -n 999


compare_categories.py \
            --method anosim \
            -i $MYCOUT/core_output/0001/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/anosim_out/0001/bray_curtis/ \
            -n 999

compare_categories.py \
            --method permanova \
            -i $MYCOUT/core_output/0001/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permanova_out/0001/euclidean/ \
            -n 999


compare_categories.py \
            --method permanova \
            -i $MYCOUT/core_output/0001/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permanova_out/0001/bray_curtis/ \
            -n 999


compare_categories.py \
            --method permdisp \
            -i $MYCOUT/core_output/0001/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permdisp_out/0001/euclidean/ \
            -n 999


compare_categories.py \
            --method permdisp \
            -i $MYCOUT/core_output/0001/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permdisp_out/0001/bray_curtis/ \
            -n 999


compare_categories.py \
            --method dbrda \
            -i $MYCOUT/core_output/0001/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/dbrda_out/0001/euclidean/ \
            -n 999


compare_categories.py \
            --method dbrda \
            -i $MYCOUT/core_output/0001/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/dbrda_out/0001/bray_curtis/ \
            -n 999


compare_categories.py \
            --method mrpp \
            -i $MYCOUT/core_output/0001/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/mrpp_out/0001/euclidean/ \
            -n 999


compare_categories.py \
            --method mrpp \
            -i $MYCOUT/core_output/0001/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/mrpp_out/0001/bray_curtis/ \
            -n 999

### Singleton

compare_categories.py \
            --method adonis \
            -i $MYCOUT/core_output/singleton/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/adonis_out/singleton/euclidean/ \
            -n 999


compare_categories.py \
            --method adonis \
            -i $MYCOUT/core_output/singleton/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/adonis_out/singleton/bray_curtis/ \
            -n 999


compare_categories.py \
            --method anosim \
            -i $MYCOUT/core_output/singleton/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/anosim_out/singleton/euclidean/ \
            -n 999


compare_categories.py \
            --method anosim \
            -i $MYCOUT/core_output/singleton/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/anosim_out/singleton/bray_curtis/ \
            -n 999

compare_categories.py \
            --method permanova \
            -i $MYCOUT/core_output/singleton/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permanova_out/singleton/euclidean/ \
            -n 999


compare_categories.py \
            --method permanova \
            -i $MYCOUT/core_output/singleton/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permanova_out/singleton/bray_curtis/ \
            -n 999


compare_categories.py \
            --method permdisp \
            -i $MYCOUT/core_output/singleton/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permdisp_out/singleton/euclidean/ \
            -n 999


compare_categories.py \
            --method permdisp \
            -i $MYCOUT/core_output/singleton/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permdisp_out/singleton/bray_curtis/ \
            -n 999


compare_categories.py \
            --method dbrda \
            -i $MYCOUT/core_output/singleton/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/dbrda_out/singleton/euclidean/ \
            -n 999


compare_categories.py \
            --method dbrda \
            -i $MYCOUT/core_output/singleton/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/dbrda_out/singleton/bray_curtis/ \
            -n 999


compare_categories.py \
            --method mrpp \
            -i $MYCOUT/core_output/singleton/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/mrpp_out/singleton/euclidean/ \
            -n 999


compare_categories.py \
            --method mrpp \
            -i $MYCOUT/core_output/singleton/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/mrpp_out/singleton/bray_curtis/ \
            -n 999

### CORE 00001 - 70%

compare_categories.py \
            --method adonis \
            -i $MYCOUT/core_output/core_00001_70p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/adonis_out/core_00001_70p/euclidean/ \
            -n 999


compare_categories.py \
            --method adonis \
            -i $MYCOUT/core_output/core_00001_70p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/adonis_out/core_00001_70p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method anosim \
            -i $MYCOUT/core_output/core_00001_70p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/anosim_out/core_00001_70p/euclidean/ \
            -n 999


compare_categories.py \
            --method anosim \
            -i $MYCOUT/core_output/core_00001_70p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/anosim_out/core_00001_70p/bray_curtis/ \
            -n 999

compare_categories.py \
            --method permanova \
            -i $MYCOUT/core_output/core_00001_70p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permanova_out/core_00001_70p/euclidean/ \
            -n 999


compare_categories.py \
            --method permanova \
            -i $MYCOUT/core_output/core_00001_70p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permanova_out/core_00001_70p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method permdisp \
            -i $MYCOUT/core_output/core_00001_70p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permdisp_out/core_00001_70p/euclidean/ \
            -n 999


compare_categories.py \
            --method permdisp \
            -i $MYCOUT/core_output/core_00001_70p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permdisp_out/core_00001_70p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method dbrda \
            -i $MYCOUT/core_output/core_00001_70p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/dbrda_out/core_00001_70p/euclidean/ \
            -n 999


compare_categories.py \
            --method dbrda \
            -i $MYCOUT/core_output/core_00001_70p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/dbrda_out/core_00001_70p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method mrpp \
            -i $MYCOUT/core_output/core_00001_70p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/mrpp_out/core_00001_70p/euclidean/ \
            -n 999


compare_categories.py \
            --method mrpp \
            -i $MYCOUT/core_output/core_00001_70p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/mrpp_out/core_00001_70p/bray_curtis/ \
            -n 999

### CORE 00001 - 80%

compare_categories.py \
            --method adonis \
            -i $MYCOUT/core_output/core_00001_80p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/adonis_out/core_00001_80p/euclidean/ \
            -n 999


compare_categories.py \
            --method adonis \
            -i $MYCOUT/core_output/core_00001_80p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/adonis_out/core_00001_80p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method anosim \
            -i $MYCOUT/core_output/core_00001_80p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/anosim_out/core_00001_80p/euclidean/ \
            -n 999


compare_categories.py \
            --method anosim \
            -i $MYCOUT/core_output/core_00001_80p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/anosim_out/core_00001_80p/bray_curtis/ \
            -n 999

compare_categories.py \
            --method permanova \
            -i $MYCOUT/core_output/core_00001_80p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permanova_out/core_00001_80p/euclidean/ \
            -n 999


compare_categories.py \
            --method permanova \
            -i $MYCOUT/core_output/core_00001_80p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permanova_out/core_00001_80p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method permdisp \
            -i $MYCOUT/core_output/core_00001_80p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permdisp_out/core_00001_80p/euclidean/ \
            -n 999


compare_categories.py \
            --method permdisp \
            -i $MYCOUT/core_output/core_00001_80p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permdisp_out/core_00001_80p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method dbrda \
            -i $MYCOUT/core_output/core_00001_80p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/dbrda_out/core_00001_80p/euclidean/ \
            -n 999


compare_categories.py \
            --method dbrda \
            -i $MYCOUT/core_output/core_00001_80p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/dbrda_out/core_00001_80p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method mrpp \
            -i $MYCOUT/core_output/core_00001_80p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/mrpp_out/core_00001_80p/euclidean/ \
            -n 999


compare_categories.py \
            --method mrpp \
            -i $MYCOUT/core_output/core_00001_80p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/mrpp_out/core_00001_80p/bray_curtis/ \
            -n 999

### CORE 00001 - 90%

compare_categories.py \
            --method adonis \
            -i $MYCOUT/core_output/core_00001_90p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/adonis_out/core_00001_90p/euclidean/ \
            -n 999


compare_categories.py \
            --method adonis \
            -i $MYCOUT/core_output/core_00001_90p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/adonis_out/core_00001_90p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method anosim \
            -i $MYCOUT/core_output/core_00001_90p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/anosim_out/core_00001_90p/euclidean/ \
            -n 999


compare_categories.py \
            --method anosim \
            -i $MYCOUT/core_output/core_00001_90p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/anosim_out/core_00001_90p/bray_curtis/ \
            -n 999

compare_categories.py \
            --method permanova \
            -i $MYCOUT/core_output/core_00001_90p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permanova_out/core_00001_90p/euclidean/ \
            -n 999


compare_categories.py \
            --method permanova \
            -i $MYCOUT/core_output/core_00001_90p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permanova_out/core_00001_90p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method permdisp \
            -i $MYCOUT/core_output/core_00001_90p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permdisp_out/core_00001_90p/euclidean/ \
            -n 999


compare_categories.py \
            --method permdisp \
            -i $MYCOUT/core_output/core_00001_90p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permdisp_out/core_00001_90p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method dbrda \
            -i $MYCOUT/core_output/core_00001_90p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/dbrda_out/core_00001_90p/euclidean/ \
            -n 999


compare_categories.py \
            --method dbrda \
            -i $MYCOUT/core_output/core_00001_90p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/dbrda_out/core_00001_90p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method mrpp \
            -i $MYCOUT/core_output/core_00001_90p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/mrpp_out/core_00001_90p/euclidean/ \
            -n 999


compare_categories.py \
            --method mrpp \
            -i $MYCOUT/core_output/core_00001_90p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/mrpp_out/core_00001_90p/bray_curtis/ \
            -n 999


### CORE 0001 - 70%

compare_categories.py \
            --method adonis \
            -i $MYCOUT/core_output/core_0001_70p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/adonis_out/core_0001_70p/euclidean/ \
            -n 999


compare_categories.py \
            --method adonis \
            -i $MYCOUT/core_output/core_0001_70p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/adonis_out/core_0001_70p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method anosim \
            -i $MYCOUT/core_output/core_0001_70p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/anosim_out/core_0001_70p/euclidean/ \
            -n 999


compare_categories.py \
            --method anosim \
            -i $MYCOUT/core_output/core_0001_70p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/anosim_out/core_0001_70p/bray_curtis/ \
            -n 999

compare_categories.py \
            --method permanova \
            -i $MYCOUT/core_output/core_0001_70p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permanova_out/core_0001_70p/euclidean/ \
            -n 999


compare_categories.py \
            --method permanova \
            -i $MYCOUT/core_output/core_0001_70p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permanova_out/core_0001_70p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method permdisp \
            -i $MYCOUT/core_output/core_0001_70p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permdisp_out/core_0001_70p/euclidean/ \
            -n 999


compare_categories.py \
            --method permdisp \
            -i $MYCOUT/core_output/core_0001_70p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permdisp_out/core_0001_70p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method dbrda \
            -i $MYCOUT/core_output/core_0001_70p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/dbrda_out/core_0001_70p/euclidean/ \
            -n 999


compare_categories.py \
            --method dbrda \
            -i $MYCOUT/core_output/core_0001_70p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/dbrda_out/core_0001_70p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method mrpp \
            -i $MYCOUT/core_output/core_0001_70p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/mrpp_out/core_0001_70p/euclidean/ \
            -n 999


compare_categories.py \
            --method mrpp \
            -i $MYCOUT/core_output/core_0001_70p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/mrpp_out/core_0001_70p/bray_curtis/ \
            -n 999

### CORE 0001 - 80%

compare_categories.py \
            --method adonis \
            -i $MYCOUT/core_output/core_0001_80p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/adonis_out/core_0001_80p/euclidean/ \
            -n 999


compare_categories.py \
            --method adonis \
            -i $MYCOUT/core_output/core_0001_80p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/adonis_out/core_0001_80p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method anosim \
            -i $MYCOUT/core_output/core_0001_80p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/anosim_out/core_0001_80p/euclidean/ \
            -n 999


compare_categories.py \
            --method anosim \
            -i $MYCOUT/core_output/core_0001_80p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/anosim_out/core_0001_80p/bray_curtis/ \
            -n 999

compare_categories.py \
            --method permanova \
            -i $MYCOUT/core_output/core_0001_80p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permanova_out/core_0001_80p/euclidean/ \
            -n 999


compare_categories.py \
            --method permanova \
            -i $MYCOUT/core_output/core_0001_80p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permanova_out/core_0001_80p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method permdisp \
            -i $MYCOUT/core_output/core_0001_80p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permdisp_out/core_0001_80p/euclidean/ \
            -n 999


compare_categories.py \
            --method permdisp \
            -i $MYCOUT/core_output/core_0001_80p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permdisp_out/core_0001_80p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method dbrda \
            -i $MYCOUT/core_output/core_0001_80p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/dbrda_out/core_0001_80p/euclidean/ \
            -n 999


compare_categories.py \
            --method dbrda \
            -i $MYCOUT/core_output/core_0001_80p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/dbrda_out/core_0001_80p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method mrpp \
            -i $MYCOUT/core_output/core_0001_80p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/mrpp_out/core_0001_80p/euclidean/ \
            -n 999


compare_categories.py \
            --method mrpp \
            -i $MYCOUT/core_output/core_0001_80p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/mrpp_out/core_0001_80p/bray_curtis/ \
            -n 999

### CORE 0001 - 90%

compare_categories.py \
            --method adonis \
            -i $MYCOUT/core_output/core_0001_90p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/adonis_out/core_0001_90p/euclidean/ \
            -n 999


compare_categories.py \
            --method adonis \
            -i $MYCOUT/core_output/core_0001_90p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/adonis_out/core_0001_90p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method anosim \
            -i $MYCOUT/core_output/core_0001_90p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/anosim_out/core_0001_90p/euclidean/ \
            -n 999


compare_categories.py \
            --method anosim \
            -i $MYCOUT/core_output/core_0001_90p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/anosim_out/core_0001_90p/bray_curtis/ \
            -n 999

compare_categories.py \
            --method permanova \
            -i $MYCOUT/core_output/core_0001_90p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permanova_out/core_0001_90p/euclidean/ \
            -n 999


compare_categories.py \
            --method permanova \
            -i $MYCOUT/core_output/core_0001_90p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permanova_out/core_0001_90p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method permdisp \
            -i $MYCOUT/core_output/core_0001_90p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permdisp_out/core_0001_90p/euclidean/ \
            -n 999


compare_categories.py \
            --method permdisp \
            -i $MYCOUT/core_output/core_0001_90p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permdisp_out/core_0001_90p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method dbrda \
            -i $MYCOUT/core_output/core_0001_90p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/dbrda_out/core_0001_90p/euclidean/ \
            -n 999


compare_categories.py \
            --method dbrda \
            -i $MYCOUT/core_output/core_0001_90p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/dbrda_out/core_0001_90p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method mrpp \
            -i $MYCOUT/core_output/core_0001_90p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/mrpp_out/core_0001_90p/euclidean/ \
            -n 999


compare_categories.py \
            --method mrpp \
            -i $MYCOUT/core_output/core_0001_90p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/mrpp_out/core_0001_90p/bray_curtis/ \
            -n 999

### Core Singleton

### CORE singleton - 70%

compare_categories.py \
            --method adonis \
            -i $MYCOUT/core_output/core_singleton_70p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/adonis_out/core_singleton_70p/euclidean/ \
            -n 999


compare_categories.py \
            --method adonis \
            -i $MYCOUT/core_output/core_singleton_70p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/adonis_out/core_singleton_70p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method anosim \
            -i $MYCOUT/core_output/core_singleton_70p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/anosim_out/core_singleton_70p/euclidean/ \
            -n 999


compare_categories.py \
            --method anosim \
            -i $MYCOUT/core_output/core_singleton_70p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/anosim_out/core_singleton_70p/bray_curtis/ \
            -n 999

compare_categories.py \
            --method permanova \
            -i $MYCOUT/core_output/core_singleton_70p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permanova_out/core_singleton_70p/euclidean/ \
            -n 999


compare_categories.py \
            --method permanova \
            -i $MYCOUT/core_output/core_singleton_70p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permanova_out/core_singleton_70p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method permdisp \
            -i $MYCOUT/core_output/core_singleton_70p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permdisp_out/core_singleton_70p/euclidean/ \
            -n 999


compare_categories.py \
            --method permdisp \
            -i $MYCOUT/core_output/core_singleton_70p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permdisp_out/core_singleton_70p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method dbrda \
            -i $MYCOUT/core_output/core_singleton_70p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/dbrda_out/core_singleton_70p/euclidean/ \
            -n 999


compare_categories.py \
            --method dbrda \
            -i $MYCOUT/core_output/core_singleton_70p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/dbrda_out/core_singleton_70p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method mrpp \
            -i $MYCOUT/core_output/core_singleton_70p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/mrpp_out/core_singleton_70p/euclidean/ \
            -n 999


compare_categories.py \
            --method mrpp \
            -i $MYCOUT/core_output/core_singleton_70p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/mrpp_out/core_singleton_70p/bray_curtis/ \
            -n 999

### CORE singleton - 80%

compare_categories.py \
            --method adonis \
            -i $MYCOUT/core_output/core_singleton_80p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/adonis_out/core_singleton_80p/euclidean/ \
            -n 999


compare_categories.py \
            --method adonis \
            -i $MYCOUT/core_output/core_singleton_80p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/adonis_out/core_singleton_80p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method anosim \
            -i $MYCOUT/core_output/core_singleton_80p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/anosim_out/core_singleton_80p/euclidean/ \
            -n 999


compare_categories.py \
            --method anosim \
            -i $MYCOUT/core_output/core_singleton_80p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/anosim_out/core_singleton_80p/bray_curtis/ \
            -n 999

compare_categories.py \
            --method permanova \
            -i $MYCOUT/core_output/core_singleton_80p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permanova_out/core_singleton_80p/euclidean/ \
            -n 999


compare_categories.py \
            --method permanova \
            -i $MYCOUT/core_output/core_singleton_80p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permanova_out/core_singleton_80p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method permdisp \
            -i $MYCOUT/core_output/core_singleton_80p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permdisp_out/core_singleton_80p/euclidean/ \
            -n 999


compare_categories.py \
            --method permdisp \
            -i $MYCOUT/core_output/core_singleton_80p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permdisp_out/core_singleton_80p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method dbrda \
            -i $MYCOUT/core_output/core_singleton_80p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/dbrda_out/core_singleton_80p/euclidean/ \
            -n 999


compare_categories.py \
            --method dbrda \
            -i $MYCOUT/core_output/core_singleton_80p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/dbrda_out/core_singleton_80p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method mrpp \
            -i $MYCOUT/core_output/core_singleton_80p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/mrpp_out/core_singleton_80p/euclidean/ \
            -n 999


compare_categories.py \
            --method mrpp \
            -i $MYCOUT/core_output/core_singleton_80p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/mrpp_out/core_singleton_80p/bray_curtis/ \
            -n 999

### CORE singleton - 90%

compare_categories.py \
            --method adonis \
            -i $MYCOUT/core_output/core_singleton_90p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/adonis_out/core_singleton_90p/euclidean/ \
            -n 999


compare_categories.py \
            --method adonis \
            -i $MYCOUT/core_output/core_singleton_90p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/adonis_out/core_singleton_90p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method anosim \
            -i $MYCOUT/core_output/core_singleton_90p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/anosim_out/core_singleton_90p/euclidean/ \
            -n 999


compare_categories.py \
            --method anosim \
            -i $MYCOUT/core_output/core_singleton_90p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/anosim_out/core_singleton_90p/bray_curtis/ \
            -n 999

compare_categories.py \
            --method permanova \
            -i $MYCOUT/core_output/core_singleton_90p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permanova_out/core_singleton_90p/euclidean/ \
            -n 999


compare_categories.py \
            --method permanova \
            -i $MYCOUT/core_output/core_singleton_90p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permanova_out/core_singleton_90p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method permdisp \
            -i $MYCOUT/core_output/core_singleton_90p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permdisp_out/core_singleton_90p/euclidean/ \
            -n 999


compare_categories.py \
            --method permdisp \
            -i $MYCOUT/core_output/core_singleton_90p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/permdisp_out/core_singleton_90p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method dbrda \
            -i $MYCOUT/core_output/core_singleton_90p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/dbrda_out/core_singleton_90p/euclidean/ \
            -n 999


compare_categories.py \
            --method dbrda \
            -i $MYCOUT/core_output/core_singleton_90p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/dbrda_out/core_singleton_90p/bray_curtis/ \
            -n 999


compare_categories.py \
            --method mrpp \
            -i $MYCOUT/core_output/core_singleton_90p/bdiv_${i}/euclidean_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/mrpp_out/core_singleton_90p/euclidean/ \
            -n 999


compare_categories.py \
            --method mrpp \
            -i $MYCOUT/core_output/core_singleton_90p/bdiv_${i}/bray_curtis_dm.txt \
            -m $MAP \
            -c ALL \
            -o $MYCOUT/core_output/compare_categories/mrpp_out/core_singleton_90p/bray_curtis/ \
            -n 999

done

echo "###############################"
echo "###############################"
echo "FINISH HIM!"
echo "###############################"
echo "###############################"


echo " "
echo " "
