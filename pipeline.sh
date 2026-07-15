#!/bin/bash

### Metagenomic Hybrid Assembling, Correction and Binning pipeline (MetaAHCB)
### version 1.1
### 2024/11/6 edition
#####################################################################################

scripts_DIR="/dellfsqd1/ST_OCEAN/ST_OCEAN/USRS/c-microbio/xiajun/TGS_assembly_polish_binning/20240627_test/pipeline_test/pipeline/scripts/"

### Modify the directory and parameters

### Input file directory
### Please reformat and rename TGS reads file as {sample_name}.fq.gz
### Please reformat and rename paired-end NGS reads file as {sample_name}_1.fq.gz and {sample_name}_2.fq.gz
### "sample_name" of TGS and NGS must be corresponded for each sample. 
### For example, name of a sample is "NPS_001", then in the TGS_DIR folder, the TGS reads file should be named as "NPS_001.fq.gz"; in the NGS_DIR folder, the paired end NGS reads file should be named as "NPS_001_1.fq.gz" and "NPS_001_2.fq.gz". 

TGS_DIR="/dellfsqd1/ST_OCEAN/ST_OCEAN/USRS/c-microbio/xiajun/TGS_assembly_polish_binning/20240627_test/pipeline_test/TGS_seq/"
NGS_DIR="/dellfsqd1/ST_OCEAN/ST_OCEAN/USRS/c-microbio/xiajun/TGS_assembly_polish_binning/20240627_test/pipeline_test/NGS_seq/"

### Output file directory
### If it does not exist, this pipeline will build automatically
OUTPUT_DIR="/dellfsqd1/ST_OCEAN/ST_OCEAN/USRS/c-microbio/xiajun/TGS_assembly_polish_binning/20240627_test/pipeline_test/pipeline_output_20240819/"

### Set parameters

### Number of threads (default=16)
### If the core dump (for example, core.12345) appear in the folder of this pipeline, please run "Step2: Machine learning" individually with the unfinished samples with lower threads (for example, num_proc=1 or 2) in the qsub command. 
THREADS="16"
### TGS reads minimum length after filtering (default=1000)
TGS_reads_min_len="1000"
### k-mer size for Bayesian Gaussian clustering (default=3)
KMER="3"
### PCA dimension reduction to this value (default=16)
pca_dimension_n="16"
### Bayesian Gaussian will assign TGS reads to this number of clusters (default=20)
bgm_cluster_n="20"
### the high confident overlap (HCO) depend on this value (1.0 - 5.0, default=3.0)
quickmerge_HCO="3.0"
### the minimum seed contig length when quickmerge do contig merging based on overlaps (default=5000)
quickmerge_SEED_LEN="5000"
### the minimum overlap length by quickmerge (default=1000)
quickmerge_MERGE_LEN="1000"
### the nucleotide identity when mmseqs2 deduplicate contigs (0-1, default=0.99)
dedup_ident="0.99"
### the coverage on shorter sequence when mmseqs2 deduplicate contigs (0-1, default=0.9)
dedup_cov="0.9"
### the minimum length of short reads assembled contig megahit will output (default=1000)
sr_contig_min_len="1000"
### the minimum length of contig input into metabat2 binning process (default=1500)
bin_min_contig="1500"
### the minimum bin length output from metabat2 binning process (default=50000)
bin_min_len="50000"

### Software directory and activate
seqkit_DIR="/dellfsqd3/ST_OCEAN/USER/chenjianwei/Software/"
gc_nmer_DIR="/dellfsqd1/ST_OCEAN/ST_OCEAN/USRS/xumengyang/software/Symbiont-Screener/sources/"
flye_DIR="/dellfsqd1/ST_OCEAN/ST_OCEAN/USRS/xumengyang/software/anaconda3/bin/"
flye_env="/dellfsqd1/ST_OCEAN/ST_OCEAN/USRS/xumengyang/software/anaconda2/bin/"
mmseqs_DIR="/dellfsqd2/ST_OCEAN/USER/zhouchanghao/software/miniconda3/envs/EukMS_run/bin/"
minimap2_DIR="/dellfsqd3/ST_OCEAN/USER/chenjianwei/Software/minimap2-2.26/"
HERO_DIR="/dellfsqd3/ST_OCEAN/USER/chenjianwei/Software/HERO/bin/"
#####################################################################################


#####################################################################################
### Pipeline start

### Build output directory if not exist. 
if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p "$OUTPUT_DIR"
fi

### Build log file, print pipeline info.
timestamp=$(date +%Y_%m_%d_%H_%M)
log_file="${OUTPUT_DIR}log_${timestamp}.txt"
touch $log_file

echo " " >> $log_file
echo "#####################################################" >> $log_file
echo "Metagenomic Hybrid Assembling, Correction and Binning pipeline (MetaAHCB)" >> $log_file
echo "version 1.1" >> $log_file
echo "2024/11/6 edition" >> $log_file
echo "Main program by J.X., theoretical and technical support by J.C., M.X., Y.S., Y.Q." >> $log_file
echo "#####################################################" >> $log_file
echo " " >> $log_file
echo "pipeline start..." >> $log_file
echo "#####################################################" >> $log_file
echo " " >> $log_file
echo "Start time: $(date +%Y/%m/%d"   "%H:%M)" >> $log_file
echo "#####################################################" >> $log_file
echo "      " >> $log_file

### Count number, load names of samples.
echo "loading samples..." >> $log_file
sample_num=0
SAMPLE_LIST=()
for READS_FILE in ${TGS_DIR}*.fq.gz;do
  IND_SAMPLE=$(basename $READS_FILE .fq.gz)
  SAMPLE_LIST+=($IND_SAMPLE)
  sample_num=$((sample_num+1))
done
echo "${sample_num} samples loaded in total" >> $log_file
echo "      " >> $log_file

### !!! DO NOT modify the part above. 
#####################################################################################
### Add "#" in front of each line of the part you want to skip. 


### NGS paired end reads quality control by fastp.
echo "#####################################################" >> $log_file
echo "Step00:short reads QC..." >> $log_file
echo "Current time: $(date +%Y/%m/%d"   "%H:%M)" >> $log_file
mkdir ${OUTPUT_DIR}00_passed_NGS_reads
source /dellfsqd2/ST_OCEAN/USER/sunying6/program/anaconda3.2023.03/bin/activate && conda activate fastp.0.23.4

for SAMPLE in "${SAMPLE_LIST[@]}"; do
  fastp \
    -i ${NGS_DIR}${SAMPLE}_1.fq.gz \
    -I ${NGS_DIR}${SAMPLE}_2.fq.gz \
    -o ${OUTPUT_DIR}00_passed_NGS_reads/${SAMPLE}_1_fastp.fq.gz \
    -O ${OUTPUT_DIR}00_passed_NGS_reads/${SAMPLE}_2_fastp.fq.gz \
    --thread $THREADS
done
echo "fastp performed" >> $log_file
echo "Current time: $(date +%Y/%m/%d"   "%H:%M)" >> $log_file
echo "      " >> $log_file


### TGS long reads length filter, and get k-mer frequency. 
echo "#####################################################" >> $log_file
echo "Step01: TGS reads gc_nmer processing..." >> $log_file
echo "Current time: $(date +%Y/%m/%d"   "%H:%M)" >> $log_file
echo "filtering < ${TGS_reads_min_len} TGS reads..." >> $log_file
mkdir ${OUTPUT_DIR}00_passed_TGS_reads
for SAMPLE in "${SAMPLE_LIST[@]}"; do
  ${seqkit_DIR}seqkit seq \
    -m ${TGS_reads_min_len} \
    -o ${OUTPUT_DIR}00_passed_TGS_reads/${SAMPLE}_m${TGS_reads_min_len}.fq.gz \
    ${TGS_DIR}${SAMPLE}.fq.gz
  gunzip -c ${OUTPUT_DIR}00_passed_TGS_reads/${SAMPLE}_m${TGS_reads_min_len}.fq.gz > ${OUTPUT_DIR}00_passed_TGS_reads/${SAMPLE}_m${TGS_reads_min_len}.fastq
done

mkdir ${OUTPUT_DIR}01_gc_nmer
i=1
for SAMPLE in "${SAMPLE_LIST[@]}"; do
  echo "  processing ${SAMPLE} (${i}/${sample_num})..." >> $log_file
 (${gc_nmer_DIR}gc_nmer \
    --thread $THREADS \
    --format fastq \
    --kmer ${KMER} \
    --read ${OUTPUT_DIR}00_passed_TGS_reads/${SAMPLE}_m${TGS_reads_min_len}.fq.gz ) > ${OUTPUT_DIR}01_gc_nmer/${SAMPLE}_gc_nmer_k${KMER}_out.txt
  i=$((i+1))
done
echo "gc_nmer finished." >> $log_file
echo "Current time: $(date +%Y/%m/%d"   "%H:%M)" >> $log_file
echo "      " >> $log_file


### Machine learning: TGS long reads cluster by Bayesian Gaussian mixture model, and dimension reduction by PCA.
### If the core dump (for example, core.12345) appear in the folder of this pipeline, please run this step individually with the unfinished samples with lower threads (for example, num_proc=1 or 2) in the qsub command. 
echo "#####################################################" >> $log_file
echo "Step02: TGS reads Bayesian Gaussian mixture model clustering..." >> $log_file
echo "Current time: $(date +%Y/%m/%d"   "%H:%M)" >> $log_file
mkdir ${OUTPUT_DIR}02_bgm_cluster
i=1
source /dellfsqd2/ST_OCEAN/USER/zhouchanghao/software/miniconda3/bin/activate && conda activate python-3.8

for SAMPLE in "${SAMPLE_LIST[@]}"; do
  echo "  processing ${SAMPLE} (${i}/${sample_num})..." >> $log_file
  python3.8 ${scripts_DIR}bgm_cluster.py \
    -i ${OUTPUT_DIR}01_gc_nmer/${SAMPLE}_gc_nmer_k${KMER}_out.txt \
    -o ${OUTPUT_DIR}02_bgm_cluster/${SAMPLE}_cluster.txt \
    -x ${pca_dimension_n} \
    -n ${bgm_cluster_n}
  i=$((i+1))
done
echo "bgm cluster finished." >> $log_file
echo "Current time: $(date +%Y/%m/%d"   "%H:%M)" >> $log_file
echo "      " >> $log_file


### Assign the TGS long reads into clusters based on the result of machine learning.
echo "#####################################################" >> $log_file
echo "Step03: assigning reads into clusters..." >> $log_file
echo "Current time: $(date +%Y/%m/%d"   "%H:%M)" >> $log_file
mkdir ${OUTPUT_DIR}03_clustered_reads
source /dellfsqd2/ST_OCEAN/USER/zhouchanghao/software/miniconda3/bin/activate && conda activate python-3.8

for SAMPLE in "${SAMPLE_LIST[@]}"; do
  mkdir ${OUTPUT_DIR}03_clustered_reads/${SAMPLE}
  python3.8 ${scripts_DIR}assign_reads_clusters.py \
    -i ${OUTPUT_DIR}02_bgm_cluster/${SAMPLE}_cluster.txt \
    -s ${OUTPUT_DIR}00_passed_TGS_reads/${SAMPLE}_m${TGS_reads_min_len}.fq.gz \
    -c ${bgm_cluster_n} \
    -t $THREADS \
    -o ${OUTPUT_DIR}03_clustered_reads/${SAMPLE}/${SAMPLE}
done
echo "TGS reads assigned by clusters." >> $log_file
echo "Current time: $(date +%Y/%m/%d"   "%H:%M)" >> $log_file
echo "      " >> $log_file


### Assemble clustered TGS long reads into contigs by metaflye.
echo "#####################################################" >> $log_file
echo "Step04: assemble TGS reads using metaflye..." >> $log_file
echo "Current time: $(date +%Y/%m/%d"   "%H:%M)" >> $log_file
mkdir ${OUTPUT_DIR}04_metaflye_assembly
i=1
export PATH=$flye_env:$PATH
for SAMPLE in "${SAMPLE_LIST[@]}"; do
  echo "  assembling ${SAMPLE} (${i}/${sample_num}) clustered TGS reads..." >> $log_file
  mkdir ${OUTPUT_DIR}04_metaflye_assembly/${SAMPLE}
  cd ${OUTPUT_DIR}03_clustered_reads/${SAMPLE}
  for CLUSTERED_READS in ${SAMPLE}*.fq.gz;do
    mkdir ${OUTPUT_DIR}04_metaflye_assembly/${SAMPLE}/${CLUSTERED_READS%.fq.gz}
    ${flye_DIR}flye \
      --nano-raw ${CLUSTERED_READS} \
      --meta \
      --out-dir ${OUTPUT_DIR}04_metaflye_assembly/${SAMPLE}/${CLUSTERED_READS%.fq.gz}/ \
      --threads $THREADS
  done
  i=$((i+1))
done
echo "TGS reads assembled." >> $log_file
echo "Current time: $(date +%Y/%m/%d"   "%H:%M)" >> $log_file
echo "      " >> $log_file


### Merge the clustered contigs for each sample by quickmerge, and deduplicate by mmseqs2. 
echo "#####################################################" >> $log_file
echo "Step05: merge clustered assembly" >> $log_file
echo "Current time: $(date +%Y/%m/%d"   "%H:%M)" >> $log_file
mkdir ${OUTPUT_DIR}05_quickmerge_dedup
i=1
source /dellfsqd2/ST_OCEAN/USER/zhouchanghao/software/miniconda3/bin/activate && conda activate quickmerge

for SAMPLE in "${SAMPLE_LIST[@]}"; do
  echo "  merging ${SAMPLE} (${i}/${sample_num}) contigs..." >> $log_file
  mkdir ${OUTPUT_DIR}05_quickmerge_dedup/${SAMPLE}
  cd ${OUTPUT_DIR}05_quickmerge_dedup/${SAMPLE}
  touch ${OUTPUT_DIR}05_quickmerge_dedup/${SAMPLE}_contigs_quickmerge.fasta
  NUM_LIST=()
  TOTAL=0
  for folder in ${OUTPUT_DIR}04_metaflye_assembly/${SAMPLE}/${SAMPLE}_cluster*/;do
    if [ -f ${folder}assembly.fasta ]; then
      num=$(basename "$folder" | sed -e "s/${SAMPLE}_cluster//")
      NUM_LIST+=($num)
      TOTAL=$((TOTAL+1))
      cat ${OUTPUT_DIR}04_metaflye_assembly/${SAMPLE}/${SAMPLE}_cluster${num}/assembly.fasta > ${OUTPUT_DIR}05_quickmerge_dedup/${SAMPLE}/merged_${SAMPLE}_cluster${num}_quickmerge000.fasta
    fi
  done

  j=0

  while [ $j -lt $TOTAL ]; do
    x=0
    while [ $x -lt $TOTAL ];do
      if [ $x -lt 10 ]; then
        code_x="00$x"
      else
        code_x="0$x"
      fi
      y=$((x+1))
      if [ $y -lt 10 ]; then
        code_y="00$y"
      else
        code_y="0$y"
      fi
      echo "merged_${SAMPLE}_cluster${NUM_LIST[$j]}_quickmerge${code_x}.fa"
      echo "${SAMPLE}_cluster${NUM_LIST[$x]}"
      echo "${SAMPLE}_cluster${NUM_LIST[$j]}_quickmerge${code_y}"

      merge_wrapper.py \
        ${OUTPUT_DIR}05_quickmerge_dedup/${SAMPLE}/merged_${SAMPLE}_cluster${NUM_LIST[$j]}_quickmerge${code_x}.fasta \
        ${OUTPUT_DIR}04_metaflye_assembly/${SAMPLE}/${SAMPLE}_cluster${NUM_LIST[$x]}/assembly.fasta \
        -hco ${quickmerge_HCO} \
        -c 1.5 \
        -l ${quickmerge_SEED_LEN} \
        -ml ${quickmerge_MERGE_LEN} \
        -pre ${SAMPLE}_cluster${NUM_LIST[$j]}_quickmerge${code_y}
      x=$((x+1))
    done
    j=$((j+1))
  done
  cat merged_${SAMPLE}_cluster*_quickmerge${code_y}.fasta>>${OUTPUT_DIR}05_quickmerge_dedup/${SAMPLE}_contigs_quickmerge.fasta
  awk -v sample="${SAMPLE}" '/^>/ {print ">" sample "_quickmerge_contig_" sprintf("%05d", ++i); next} {print}' ${OUTPUT_DIR}05_quickmerge_dedup/${SAMPLE}_contigs_quickmerge.fasta > ${OUTPUT_DIR}05_quickmerge_dedup/${SAMPLE}_contigs_quickmerge_renamed.fasta
  ${mmseqs_DIR}mmseqs easy-cluster \
    ${OUTPUT_DIR}05_quickmerge_dedup/${SAMPLE}_contigs_quickmerge_renamed.fasta \
    ${OUTPUT_DIR}05_quickmerge_dedup/${SAMPLE}_contigs_quickmerge_mmseq_res \
    ${OUTPUT_DIR}05_quickmerge_dedup/${SAMPLE}_contigs_quickmerge_mmseq_tmp \
    --threads $THREADS \
    --min-seq-id ${dedup_ident} \
    --cov-mode 1 \
    -c ${dedup_cov}
  ${seqkit_DIR}seqkit seq \
    -m 1000 \
    ${OUTPUT_DIR}05_quickmerge_dedup/${SAMPLE}_contigs_quickmerge_mmseq_res_rep_seq.fasta \
    -o ${OUTPUT_DIR}05_quickmerge_dedup/${SAMPLE}_contigs_quickmerge_mmseq_m1000.fasta
  i=$((i+1))
done
echo "clustered contig merged and deduplicated." >> $log_file
echo "Current time: $(date +%Y/%m/%d"   "%H:%M)" >> $log_file
echo "      " >> $log_file


### Hybrid correction for merged contigs by NGS short reads and TGS long reads, using fmlrc2 and HERO. 
echo "#####################################################" >> $log_file
echo "Step06: long contigs correcting by short reads and long reads..." >> $log_file
echo "Current time: $(date +%Y/%m/%d"   "%H:%M)" >> $log_file
mkdir ${OUTPUT_DIR}06_correction
i=1
for SAMPLE in "${SAMPLE_LIST[@]}"; do
  echo "  correcting ${SAMPLE} (${i}/${sample_num})" >> $log_file
  mkdir ${OUTPUT_DIR}06_correction/${SAMPLE}
  ${minimap2_DIR}minimap2 \
    -ax sr \
    -t $THREADS \
    ${OUTPUT_DIR}05_quickmerge_dedup/${SAMPLE}_contigs_quickmerge_mmseq_m1000.fasta \
    ${OUTPUT_DIR}00_passed_NGS_reads/${SAMPLE}_1_fastp.fq.gz \
    ${OUTPUT_DIR}00_passed_NGS_reads/${SAMPLE}_2_fastp.fq.gz \
    -o ${OUTPUT_DIR}06_correction/${SAMPLE}/${SAMPLE}_sr_lrcontigs.sam
  
  source /dellfsqd2/ST_OCEAN/USER/sunying6/program/anaconda3/bin/activate && conda activate samtools

  samtools view -Sb ${OUTPUT_DIR}06_correction/${SAMPLE}/${SAMPLE}_sr_lrcontigs.sam > ${OUTPUT_DIR}06_correction/${SAMPLE}/${SAMPLE}_sr_lrcontigs.bam
  samtools sort --threads $THREADS ${OUTPUT_DIR}06_correction/${SAMPLE}/${SAMPLE}_sr_lrcontigs.bam > ${OUTPUT_DIR}06_correction/${SAMPLE}/${SAMPLE}_sr_lrcontigs_sorted.bam
  samtools index ${OUTPUT_DIR}06_correction/${SAMPLE}/${SAMPLE}_sr_lrcontigs_sorted.bam
  samtools fastq ${OUTPUT_DIR}06_correction/${SAMPLE}/${SAMPLE}_sr_lrcontigs_sorted.bam -o ${OUTPUT_DIR}06_correction/${SAMPLE}/${SAMPLE}_sr_lrcontigs_sorted.fastq
  
  source /dellfsqd2/ST_OCEAN/USER/zhouchanghao/software/miniconda3/bin/activate && conda activate hero

  cat ${OUTPUT_DIR}06_correction/${SAMPLE}/${SAMPLE}_sr_lrcontigs_sorted.fastq | awk 'NR % 4 == 2' | tr NT TN | \
    ropebwt2 -LR | \
    tr NT TN | \
    fmlrc2-convert ${OUTPUT_DIR}06_correction/${SAMPLE}/${SAMPLE}_sr_lrcontigs_comp_msbwt.npy

  fmlrc2 \
    -t $THREADS \
    ${OUTPUT_DIR}06_correction/${SAMPLE}/${SAMPLE}_sr_lrcontigs_comp_msbwt.npy \
    ${OUTPUT_DIR}05_quickmerge_dedup/${SAMPLE}_contigs_quickmerge_mmseq_res_rep_seq.fasta \
    ${OUTPUT_DIR}06_correction/${SAMPLE}/${SAMPLE}_contigs_fmlrc1.fasta

  fmlrc2 \
    -t $THREADS \
    ${OUTPUT_DIR}06_correction/${SAMPLE}/${SAMPLE}_sr_lrcontigs_comp_msbwt.npy \
    ${OUTPUT_DIR}06_correction/${SAMPLE}/${SAMPLE}_contigs_fmlrc1.fasta \
    ${OUTPUT_DIR}06_correction/${SAMPLE}/${SAMPLE}_contigs_fmlrc2.fasta

  fmlrc2 \
    -t $THREADS \
    ${OUTPUT_DIR}06_correction/${SAMPLE}/${SAMPLE}_sr_lrcontigs_comp_msbwt.npy \
    ${OUTPUT_DIR}06_correction/${SAMPLE}/${SAMPLE}_contigs_fmlrc2.fasta \
    ${OUTPUT_DIR}06_correction/${SAMPLE}/${SAMPLE}_contigs_fmlrc3.fasta
  
  mkdir ${OUTPUT_DIR}06_correction/${SAMPLE}/HERO1/
  cd ${OUTPUT_DIR}06_correction/${SAMPLE}/HERO1/
  
  python3 ${HERO_DIR}HERO.py \
    -r ${OUTPUT_DIR}00_passed_TGS_reads/${SAMPLE}_m${TGS_reads_min_len}.fastq \
    -lc ${OUTPUT_DIR}06_correction/${SAMPLE}/${SAMPLE}_contigs_fmlrc3.fasta \
    -pl ont \
    -s 500 \
    -i 1 \
    -t ${THREADS}
  
  cat ${OUTPUT_DIR}06_correction/${SAMPLE}/HERO1/sub*/*polished.fa  ${OUTPUT_DIR}06_correction/${SAMPLE}/HERO1/polished.fa > ${OUTPUT_DIR}06_correction/${SAMPLE}/${SAMPLE}_contigs_HERO_corrected_1.fasta
  
  #This part will consume large memory resource (e.g. >2TB). It is recommended to be turned off. 

  #mkdir ${OUTPUT_DIR}06_correction/${SAMPLE}/HERO2/
  #cd ${OUTPUT_DIR}06_correction/${SAMPLE}/HERO2/
 
  #python ${HERO_DIR}HERO.py \
  #  -r ${OUTPUT_DIR}06_correction/${SAMPLE}/${SAMPLE}_sr_lrcontigs_sorted.fastq \
  #  -lc ${OUTPUT_DIR}06_correction/${SAMPLE}/${SAMPLE}_contigs_HERO_corrected_1.fasta \
  #  -pl ont \
  #  -s 500 \
  #  -i 1 \
  #  -t ${THREADS}

  #cat ${OUTPUT_DIR}06_correction/${SAMPLE}/HERO2/sub*/*polished.fa  ${OUTPUT_DIR}06_correction/${SAMPLE}/HERO2/polished.fa > ${OUTPUT_DIR}06_correction/${SAMPLE}/${SAMPLE}_contigs_HERO_corrected_2.fasta
  
  if [-f ${OUTPUT_DIR}06_correction/${SAMPLE}/${SAMPLE}_contigs_HERO_corrected_2.fasta];then
    ${seqkit_DIR}seqkit seq \
      -m 1000 \
      ${OUTPUT_DIR}06_correction/${SAMPLE}/${SAMPLE}_contigs_HERO_corrected_2.fasta \
      -o ${OUTPUT_DIR}06_correction/${SAMPLE}/${SAMPLE}_contigs_HERO_corrected_2_m1000.fasta
  else
    ${seqkit_DIR}seqkit seq \
      -m 1000 \
      ${OUTPUT_DIR}06_correction/${SAMPLE}/${SAMPLE}_contigs_HERO_corrected_1.fasta \
      -o ${OUTPUT_DIR}06_correction/${SAMPLE}/${SAMPLE}_contigs_HERO_corrected_2_m1000.fasta
  fi
  
  i=$((i+1))
done
echo "long contigs corrected by short reads and long reads." >> $log_file
echo "Current time: $(date +%Y/%m/%d"   "%H:%M)" >> $log_file
echo "      " >> $log_file


### Assemble NGS short reads by megahit. 
echo "#####################################################" >> $log_file
echo "Step07: assemble NGS reads using megahit..." >> $log_file
echo "Current time: $(date +%Y/%m/%d"   "%H:%M)" >> $log_file
mkdir ${OUTPUT_DIR}07_megahit_assembly
i=1
source /dellfsqd2/ST_OCEAN/USER/sunying6/program/anaconda3/bin/activate && conda activate megahit

for SAMPLE in "${SAMPLE_LIST[@]}"; do
  echo "  assembling ${SAMPLE} (${i}/${sample_num}) NGS reads..." >> $log_file
  megahit \
    -1 ${OUTPUT_DIR}00_passed_NGS_reads/${SAMPLE}_1_fastp.fq.gz \
    -2 ${OUTPUT_DIR}00_passed_NGS_reads/${SAMPLE}_2_fastp.fq.gz \
    -o ${OUTPUT_DIR}07_megahit_assembly/${SAMPLE}/ \
    --presets meta-large \
    --min-contig-len ${sr_contig_min_len} \
    --continue \
    -m 0.9 \
    -t $THREADS
  i=$((i+1))
done
echo "NGS reads assembled."  >> $log_file
echo "Current time: $(date +%Y/%m/%d"   "%H:%M)"  >> $log_file
echo "      "  >> $log_file


### Merge hybrid-corrected contigs and short reads assembled contigs by quickmerge, deduplicate by mmseqs2. 
echo "#####################################################" >> $log_file
echo "Step08: merge contigs assembled by long reads and short reads..." >> $log_file
echo "Current time: $(date +%Y/%m/%d"   "%H:%M)" >> $log_file
mkdir ${OUTPUT_DIR}08_quickmerge_dedup_lc_sc
cd ${OUTPUT_DIR}08_quickmerge_dedup_lc_sc
i=1
source /dellfsqd2/ST_OCEAN/USER/zhouchanghao/software/miniconda3/bin/activate && conda activate quickmerge

for SAMPLE in "${SAMPLE_LIST[@]}"; do
  echo "  merging ${SAMPLE} (${i}/${sample_num}) lr/sr contigs..."  >> $log_file
  merge_wrapper.py \
    ${OUTPUT_DIR}06_correction/${SAMPLE}/${SAMPLE}_contigs_HERO_corrected_2_m1000.fasta \
    ${OUTPUT_DIR}07_megahit_assembly/${SAMPLE}/final.contigs.fa \
    -hco ${quickmerge_HCO} \
    -c 1.5 \
    -l ${quickmerge_SEED_LEN} \
    -ml ${quickmerge_MERGE_LEN} \
    -pre ${SAMPLE}_lc_sc_quickmerge
  awk -v sample="${SAMPLE}" '/^>/ {print ">" sample "_lc_sc_quickmerge_contig_" sprintf("%05d", ++i); next} {print}' \
    ${OUTPUT_DIR}08_quickmerge_dedup_lc_sc/merged_${SAMPLE}_lc_sc_quickmerge.fasta \
    > ${OUTPUT_DIR}08_quickmerge_dedup_lc_sc/${SAMPLE}_lc_sc_quickmerge_renamed.fasta
  ${mmseqs_DIR}mmseqs easy-cluster \
    ${OUTPUT_DIR}08_quickmerge_dedup_lc_sc/${SAMPLE}_lc_sc_quickmerge_renamed.fasta \
    ${OUTPUT_DIR}08_quickmerge_dedup_lc_sc/${SAMPLE}_lc_sc_quickmerge_mmseqs_res \
    ${OUTPUT_DIR}08_quickmerge_dedup_lc_sc/${SAMPLE}_lc_sc_quickmerge_mmseqs_tmp \
    --threads $THREADS \
    --min-seq-id ${dedup_ident} \
    --cov-mode 1 \
    -c ${dedup_cov}
  i=$((i+1))
done
echo "lr/sr contigs merged and deduplicated." >> $log_file
echo "         " >> $log_file
echo "High quality contigs acquired!!!" >> $log_file
echo "Current time: $(date +%Y/%m/%d"   "%H:%M)" >> $log_file
echo "      " >> $log_file


### Binning the merged contigs by MetaBat2, based on the depth of short reads. 
echo "#####################################################" >> $log_file
echo "Step09: binning..." >> $log_file
echo "Current time: $(date +%Y/%m/%d"   "%H:%M)" >> $log_file
mkdir ${OUTPUT_DIR}09_metabat2_bin
mkdir ${OUTPUT_DIR}09_metabat2_bin/mapping/
mkdir ${OUTPUT_DIR}09_metabat2_bin/depth/
mkdir ${OUTPUT_DIR}09_metabat2_bin/bins/
mkdir ${OUTPUT_DIR}09_metabat2_bin/dRep/
i=1

for SAMPLE_i in "${SAMPLE_LIST[@]}"; do
  echo "  binning ${SAMPLE_i} (${i}/${sample_num}) contigs..." >> $log_file
  for SAMPLE_j in "${SAMPLE_LIST[@]}"; do
    ${minimap2_DIR}minimap2 \
      -ax sr \
      -t $THREADS \
      ${OUTPUT_DIR}08_quickmerge_dedup_lc_sc/${SAMPLE_i}_lc_sc_quickmerge_mmseqs_res_rep_seq.fasta \
      ${OUTPUT_DIR}00_passed_NGS_reads/${SAMPLE_j}_1_fastp.fq.gz \
      ${OUTPUT_DIR}00_passed_NGS_reads/${SAMPLE_j}_1_fastp.fq.gz \
      -o ${OUTPUT_DIR}09_metabat2_bin/mapping/${SAMPLE_j}_sr_${SAMPLE_i}_contigs.sam

    source /dellfsqd2/ST_OCEAN/USER/sunying6/program/anaconda3/bin/activate && conda activate samtools
    samtools view -Sb ${OUTPUT_DIR}09_metabat2_bin/mapping/${SAMPLE_j}_sr_${SAMPLE_i}_contigs.sam > ${OUTPUT_DIR}09_metabat2_bin/mapping/${SAMPLE_j}_sr_${SAMPLE_i}_contigs.bam
    samtools sort --threads $THREADS ${OUTPUT_DIR}09_metabat2_bin/mapping/${SAMPLE_j}_sr_${SAMPLE_i}_contigs.bam > ${OUTPUT_DIR}09_metabat2_bin/mapping/${SAMPLE_j}_sr_${SAMPLE_i}_contigs_sorted.bam
    samtools index ${OUTPUT_DIR}09_metabat2_bin/mapping/${SAMPLE_j}_sr_${SAMPLE_i}_contigs_sorted.bam

    source /dellfsqd2/ST_OCEAN/USER/sunying6/program/anaconda3.2023.03/bin/activate && conda activate metabat2.2.15
    jgi_summarize_bam_contig_depths \
      --outputDepth ${OUTPUT_DIR}09_metabat2_bin/depth/${SAMPLE_j}_sr_${SAMPLE_i}_contigs_depth.txt \
      ${OUTPUT_DIR}09_metabat2_bin/mapping/${SAMPLE_j}_sr_${SAMPLE_i}_contigs_sorted.bam
  done
  
  source /dellfsqd2/ST_OCEAN/USER/sunying6/program/anaconda3.2023.03/bin/activate && conda activate metabat2.2.15
  cd ${OUTPUT_DIR}09_metabat2_bin/depth/

  depth_files=(*_sr_${SAMPLE_i}_contigs_depth.txt)
  file_count=${#depth_files[@]}

  if [ "$file_count" -eq 1 ]; then
    cp "${depth_files[0]}" "${OUTPUT_DIR}09_metabat2_bin/depth/${SAMPLE_i}_depth.txt"
  elif [ "$file_count" -gt 1 ]; then
    merge_depths.pl *_sr_${SAMPLE_i}_contigs_depth.txt > ${OUTPUT_DIR}09_metabat2_bin/depth/${SAMPLE_i}_depth.txt
  else
    echo "No matching depth files found for ${SAMPLE_i}"
  fi

  mkdir ${OUTPUT_DIR}09_metabat2_bin/bins/${SAMPLE_i}_metabat2bin/

  metabat2 \
    -i ${OUTPUT_DIR}08_quickmerge_dedup_lc_sc/${SAMPLE_i}_lc_sc_quickmerge_mmseqs_res_rep_seq.fasta \
    -a ${OUTPUT_DIR}09_metabat2_bin/depth/${SAMPLE_i}_depth.txt \
    -o ${OUTPUT_DIR}09_metabat2_bin/bins/${SAMPLE_i}_metabat2bin/${SAMPLE_i}_lc_metabat2bin \
    -m ${bin_min_contig} \
    -s ${bin_min_len} \
    --minS 60 \
    -t $THREADS

  file_counter=1
  cd ${OUTPUT_DIR}09_metabat2_bin/bins/${SAMPLE_i}_metabat2bin/

  for file in *.fa; do
    new_filename=$(printf "metabat2bin_lc_%05d.fasta" "$file_counter")
    mv "$file" ${new_filename}
    file_counter=$((file_counter + 1))
  done

  i=$((i+1))
done


###Binning with short reads assembled contigs which are not merged on long reads assembled contigs. 

for SAMPLE in "${SAMPLE_LIST[@]}"; do
  cat ${OUTPUT_DIR}08_quickmerge_dedup_lc_sc/aln_summary_${SAMPLE}_lc_sc_quickmerge.tsv | grep -v "^REF" | cut -f 1 | sed 's/_flag.*$//g' > ${OUTPUT_DIR}08_quickmerge_dedup_lc_sc/${SAMPLE}_pattern.txt
  ${seqkit_DIR}seqkit grep \
    -f ${OUTPUT_DIR}08_quickmerge_dedup_lc_sc/${SAMPLE}_pattern.txt \
    -v ${OUTPUT_DIR}07_megahit_assembly/${SAMPLE}/final.contigs.fa \
    -o ${OUTPUT_DIR}08_quickmerge_dedup_lc_sc/${SAMPLE}_sc_no_merged.fa
done

mkdir ${OUTPUT_DIR}08_sc_no_merged
cd ${OUTPUT_DIR}08_merged_dedup
for SAMPLE in "${SAMPLE_LIST[@]}"; do
  awk -v sample="${SAMPLE}" '/^>/ {print ">" sample "_no_merged_contig_" sprintf("%05d", ++i); next} {print}' \
    ${OUTPUT_DIR}08_quickmerge_dedup_lc_sc/${SAMPLE}_sc_no_merged.fa \
    > ${OUTPUT_DIR}08_sc_no_merged/${SAMPLE}_sc_no_merged_renamed.fasta
done

mkdir ${OUTPUT_DIR}09_metabat2_bin_sc_no_merged
mkdir ${OUTPUT_DIR}09_metabat2_bin_sc_no_merged/mapping/
mkdir ${OUTPUT_DIR}09_metabat2_bin_sc_no_merged/depth/
mkdir ${OUTPUT_DIR}09_metabat2_bin_sc_no_merged/bins/
mkdir ${OUTPUT_DIR}09_metabat2_bin_sc_no_merged/dRep/

for SAMPLE_i in "${SAMPLE_LIST[@]}"; do
  for SAMPLE_j in "${SAMPLE_LIST[@]}"; do
    ${minimap2_DIR}minimap2 \
      -ax sr \
      -t $THREADS \
      ${OUTPUT_DIR}08_sc_no_merged/${SAMPLE_i}_sc_no_merged_renamed.fasta \
      ${OUTPUT_DIR}00_passed_NGS_reads/${SAMPLE_j}_1_fastp.fq.gz \
      ${OUTPUT_DIR}00_passed_NGS_reads/${SAMPLE_j}_1_fastp.fq.gz \
      -o ${OUTPUT_DIR}09_metabat2_bin_sc_no_merged/mapping/${SAMPLE_j}_sr_${SAMPLE_i}_contigs.sam
    source /dellfsqd2/ST_OCEAN/USER/sunying6/program/anaconda3/bin/activate && conda activate samtools
    samtools view -Sb ${OUTPUT_DIR}09_metabat2_bin_sc_no_merged/mapping/${SAMPLE_j}_sr_${SAMPLE_i}_contigs.sam > ${OUTPUT_DIR}09_metabat2_bin_sc_no_merged/mapping/${SAMPLE_j}_sr_${SAMPLE_i}_contigs.bam
    samtools sort --threads $THREADS ${OUTPUT_DIR}09_metabat2_bin_sc_no_merged/mapping/${SAMPLE_j}_sr_${SAMPLE_i}_contigs.bam > ${OUTPUT_DIR}09_metabat2_bin_sc_no_merged/mapping/${SAMPLE_j}_sr_${SAMPLE_i}_contigs_sorted.bam
    samtools index ${OUTPUT_DIR}09_metabat2_bin_sc_no_merged/mapping/${SAMPLE_j}_sr_${SAMPLE_i}_contigs_sorted.bam

    source /dellfsqd2/ST_OCEAN/USER/sunying6/program/anaconda3.2023.03/bin/activate && conda activate metabat2.2.15
    jgi_summarize_bam_contig_depths \
      --outputDepth ${OUTPUT_DIR}09_metabat2_bin_sc_no_merged/depth/${SAMPLE_j}_sr_${SAMPLE_i}_contigs_depth.txt \
      ${OUTPUT_DIR}09_metabat2_bin_sc_no_merged/mapping/${SAMPLE_j}_sr_${SAMPLE_i}_contigs_sorted.bam
  done

  source /dellfsqd2/ST_OCEAN/USER/sunying6/program/anaconda3.2023.03/bin/activate && conda activate metabat2.2.15
  cd ${OUTPUT_DIR}09_metabat2_bin_sc_no_merged/depth/

  depth_files=(*_sr_${SAMPLE_i}_contigs_depth.txt)
  file_count=${#depth_files[@]}

  if [ "$file_count" -eq 1 ]; then
    cp "${depth_files[0]}" "${OUTPUT_DIR}09_metabat2_bin_sc_no_merged/depth/${SAMPLE_i}_depth.txt"
  elif [ "$file_count" -gt 1 ]; then
    merge_depths.pl *_sr_${SAMPLE_i}_contigs_depth.txt > ${OUTPUT_DIR}09_metabat2_bin_sc_no_merged/depth/${SAMPLE_i}_depth.txt
  else
    echo "No matching depth files found for ${SAMPLE_i}"
  fi

  mkdir ${OUTPUT_DIR}09_metabat2_bin_sc_no_merged/bins/${SAMPLE_i}_metabat2bin/
  metabat2 \
    -i ${OUTPUT_DIR}08_sc_no_merged/${SAMPLE_i}_sc_no_merged_renamed.fasta \
    -a ${OUTPUT_DIR}09_metabat2_bin_sc_no_merged/depth/${SAMPLE_i}_depth.txt \
    -o ${OUTPUT_DIR}09_metabat2_bin_sc_no_merged/bins/${SAMPLE_i}_metabat2bin/${SAMPLE_i}_sc_metabat2bin \
    -m ${bin_min_contig} \
    -s ${bin_min_len} \
    --minS 60 \
    -t $THREADS

  file_counter=1
  cd ${OUTPUT_DIR}09_metabat2_bin_sc_no_merged/bins/${SAMPLE_i}_metabat2bin/

  for file in *.fa; do
    new_filename=$(printf "metabat2bin_sc_no_merged_%05d.fasta" "$file_counter")
    mv "$file" ${new_filename}
    file_counter=$((file_counter + 1))
  done

  i=$((i+1))
done

echo "binning finished." >> $log_file
echo "Current time: $(date +%Y/%m/%d"   "%H:%M)" >> $log_file
echo "      " >> $log_file
echo "#####################################################" >> $log_file


### Dereplicate bins, checkM2, GTDB
echo "#####################################################" >> $log_file
echo "Step10: dereplicate bins..." >> $log_file
echo "Current time: $(date +%Y/%m/%d"   "%H:%M)" >> $log_file
mkdir ${OUTPUT_DIR}10_dRep_bins
mkdir ${OUTPUT_DIR}10_dRep_bins/raw_bins

for SAMPLE in "${SAMPLE_LIST[@]}"; do
  cd ${OUTPUT_DIR}09_metabat2_bin/bins/${SAMPLE}_sr_metabat2bin/
  for FILE in *.fasta;do
    cp $FILE ${OUTPUT_DIR}10_dRep_bins/raw_bins/${SAMPLE}_${FILE}
  done
  cd ${OUTPUT_DIR}09_metabat2_bin_sc_no_merged/bins/${SAMPLE}_sr_metabat2bin/
  for FILE in *.fasta;do
    cp $FILE ${OUTPUT_DIR}10_dRep_bins/raw_bins/${SAMPLE}_${FILE}
  done
done

source /dellfsqd2/ST_OCEAN/USER/sunying6/program/anaconda3.2023.03/bin/activate && conda activate drep.3.4.5

dRep dereplicate \
  -g ${OUTPUT_DIR}10_dRep_bins/raw_bins/*.fasta \
  -sa 0.99 \
  --ignoreGenomeQuality \
  --skip_plots \
  ${OUTPUT_DIR}10_dRep_bins/dRep_99/

source /dellfsqd2/ST_OCEAN/USER/sunying6/program/anaconda3.2023.03/bin/activate
conda activate gtdbtk.2.3.2

time gtdbtk classify_wf \
  --genome_dir ${OUTPUT_DIR}10_dRep_bins/dRep_99/dereplicated_genomes \
  --out_dir ${OUTPUT_DIR}10_dRep_bins/bins_dRep99_GTDB \
  -x fasta \
  --cpus $THREADS \
  --scratch_dir ${OUTPUT_DIR}10_dRep_bins_all_sep/bins_dRep99_gtdbtk \
  --mash_db gtdbtk.mash

cd ${OUTPUT_DIR}10_dRep_bins/

source /dellfsqd2/ST_OCEAN/USER/sunying6/program/anaconda3.2023.03/bin/activate
conda activate checkm2

time checkm2 predict \
  --threads $THREADS \
  -x fasta \
  --input dRep_99/dereplicated_genomes \
  --output-directory bins_dRep99_checkm2 \
  --database_path /dellfsqd3/ST_OCEAN/USER/chenjianwei/Database/checkm2/uniref100.KO.1.dmnd

echo "bins dereplicated. GTDB, checkM2 done." >> $log_file
echo "Current time: $(date +%Y/%m/%d"   "%H:%M)" >> $log_file
echo "      " >> $log_file
echo "#####################################################" >> $log_file


## Pipeline end
echo "done." >> $log_file
echo "Current time: $(date +%Y/%m/%d"   "%H:%M)" >> $log_file
echo "      " >> $log_file

