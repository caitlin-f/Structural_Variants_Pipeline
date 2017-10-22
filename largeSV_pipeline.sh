#!/bin/bash

# Caitlin Falconer - 22/10/2017
# https://github.com/caitlin-f/Structural_Variants_Pipeline

### Implements various tools including BreakDancer, Crest, Delly and Pindel to find large structural variations
# Dependencies:
# bcftools, BreakDancer, BWA, Crest, Delly, Picard, Pindel, Samtools, BLAT (gfServer)
# Scripts collate_output.py, dellytree.py and largeSVtree.R required in current working directory or in path

# Requires:
# Indexed reference: samtools faidx reference.fa and bwa index reference.fa
# Reference in 2bit format: faToTwoBit reference.fa reference.2bit
# If script manually stopped during processing, gfServer will need stop message to quit

# Usage: final_pipeline.sh -f Path/To/REFERENCE.fa -b Path/To/REFERENCE.2bit -o Path/To/OutDir -s SAMPLEA_1 SAMBLEA_2 SAMPLEB_1 SAMPLEB_2 ....
# Path to reference files and output directory must given as absolute
# -B : Perform BreakDancer
# -C : Perform Crest
# -D : Perform Delly
# -P int : Perform Pindel, specify library insert size

while getopts ":f:b:o:s:BCDP:" flag; do
  case "${flag}" in
    f ) REF="${OPTARG}" ;;
    b ) BIT="${OPTARG}" ;;
    o ) OUTDIR="${OPTARG}" ;;
    s ) SAMPLES=("${OPTARG}") # Save each arg after -s until end or next - into SAMPLES
      until [[ $(eval "echo \${$OPTIND}") =~ ^-.* ]] || [ -z $(eval "echo \${$OPTIND}") ]; do
              SAMPLES+=($(eval "echo \${$OPTIND}"))
              OPTIND=$((OPTIND + 1))
            done 
            ;;
    B ) BREAK=True ;;
    C ) CREST=True ;;
    D ) DELLY=True ;;
    P ) PINDEL=True ; INSERT="${OPTARG}" ;;
  esac
done

# Error checking to ensure all sample files are paired
if [ $((${#SAMPLES[@]} %2)) -eq 1 ] ; then
  echo "Error: Not all files are paired"
  exit
fi

# Remove last '/' if has been put entered in OUTDIR by user
OUTDIR=${OUTDIR%\/}

# Split input sample files into array for each pair
for i in "${!SAMPLES[@]}" # for each index in SAMPLE array
do
  ((i%2==0)) && SAMPLE_1+=(${SAMPLES[i]})
  ((i%2==1)) && SAMPLE_2+=(${SAMPLES[i]})
done

# Make output directories
mkdir -p ${OUTDIR}/1_Mapping
if [ ${BREAK} ] ; then mkdir -p ${OUTDIR}/2_LargeSVs/BreakDancer ; fi
if [ ${CREST} ] ; then mkdir -p ${OUTDIR}/2_LargeSVs/Crest ; fi
if [ ${DELLY} ] ; then mkdir -p ${OUTDIR}/2_LargeSVs/Delly ; fi
if [ ${PINDEL} ] ; then mkdir -p ${OUTDIR}/2_LargeSVs/Pindel/vcf ; fi

# File for full collated results
touch ${OUTDIR}/2_LargeSVs/all_data.txt

# CREST Requires absolute file paths as CREST.pl changes to a temporary directory
DIR=$(pwd)

# Due to system integrity protection (SIP) in mac when invoke .sh, all DYLD_* environment variables are purged
# Specifying paths required for Delly
DELLY=$(which delly)
DLPATH=${DELLY%\/*}
export DYLD_FALLBACK_LIBRARY_PATH=$DYLD_FALLBACK_LIBRARY_PATH:${DLPATH}/modular-boost/stage/lib/:${DLPATH}/htslib/


echo -e "\nReference file: ${REF}"
echo "2bit file: ${BIT}"
echo -e "Output directory: ${OUTDIR}\n"


# Run bwa in its own loop, possible improvement in speed by keeping FMD-index of reference in memory
for i in "${!SAMPLE_1[@]}"
do
  # Extract just filename portion
  FILE=${SAMPLE_1[i]}
  FILENAME=${FILE##*/}
  FILENAME=${FILENAME%_*}
  echo -e "\nRunning BWA on ${FILENAME}...\n"

  ### bwa mem - map reads to reference, run on 4 threads (-t 4)
  # bam files sorted, necessary for breakdancer config files (bam2cfg.pl)
  bwa mem -t 4 ${REF} ${SAMPLE_1[i]} ${SAMPLE_2[i]} | samtools sort | samtools view -b \
  > ${OUTDIR}/1_Mapping/${FILENAME}.noRG.bam
done

# Add RG header to bam files (necessary for BreakDancer) and index
for i in "${!SAMPLE_1[@]}"
do
  FILE=${SAMPLE_1[i]}
  FILENAME=${FILE##*/}
  FILENAME=${FILENAME%_*}
  ### picard - add RG header to bam files
  java -jar /Users/Caitlin/Applications/picard.jar AddOrReplaceReadGroups \
  I=${OUTDIR}/1_Mapping/${FILENAME}.noRG.bam \
  O=${OUTDIR}/1_Mapping/${FILENAME}.bam \
  RGID=${FILENAME} RGLB=lib${FILENAME} RGPL=illumina RGPU=unit${FILENAME} RGSM=SM${FILENAME}

  ### index bam file
  samtools index ${OUTDIR}/1_Mapping/${FILENAME}.bam
done

echo -e "\nMapping finished...\n"

# CREST requires BLAT server (gfServer) -canStop arg so that quit message will stop server
gfServer start -canStop localhost 6666 ${BIT} &

for i in "${!SAMPLE_1[@]}"
do
  FILE=${SAMPLE_1[i]}
  FILENAME=${FILE##*/}
  FILENAME=${FILENAME%_*}

  ### BreakDancer
  # Requires config files, use bam2cfg.pl
  # -d : output fastq files of reads covering SV locations
  if [ ${BREAK} ] ; then
    echo -e "\nRunning BreakDancer on ${FILENAME}...\n"
    # make config files
    bam2cfg.pl ${OUTDIR}/1_Mapping/${FILENAME}.bam \
    > ${OUTDIR}/2_LargeSVs/BreakDancer/${FILENAME}.bd.conf
    
    breakdancer-max -d ${OUTDIR}/2_LargeSVs/BreakDancer/${FILENAME}.bd \
    ${OUTDIR}/2_LargeSVs/BreakDancer/${FILENAME}.bd.conf \
    > ${OUTDIR}/2_LargeSVs/BreakDancer/${FILENAME}.bd.out
  fi

  ### Pindel
  # Requires tab separated config file: Path/To/sample.bam  insertsize  samplename
  # -l : report long insertions
  # -M 10 : change minimum support for event to 10 (default 1 returning false positives)
  # -x 3 : change maximum size of structural variant to detect to 2048bp (default 2 = 512bp)
  # -q : detect distant duplications (DD) - requires a separate run
  if [ ${PINDEL} ] ; then
    echo -e "\nRunning Pindel on ${FILENAME}...\n"
    mkdir ${OUTDIR}/2_LargeSVs/Pindel/${FILENAME}_Results/
    
    #make config file
    printf "${OUTDIR}/1_Mapping/${FILENAME}.bam\t${INSERT}\t${FILENAME}" > ${OUTDIR}/2_LargeSVs/Pindel/${FILENAME}_Results/${FILENAME}.pd.conf

    pindel -l -M 10 -x 3 -f ${REF} -i ${OUTDIR}/2_LargeSVs/Pindel/${FILENAME}_Results/${FILENAME}.pd.conf \
    -o ${OUTDIR}/2_LargeSVs/Pindel/${FILENAME}_Results/${FILENAME}.pd.out
    pindel -q -M 10 -x 3 -f ${REF} -i ${OUTDIR}/2_LargeSVs/Pindel/${FILENAME}_Results/${FILENAME}.pd.conf \
    -o ${OUTDIR}/2_LargeSVs/Pindel/${FILENAME}_Results/${FILENAME}.pd.q.out

    # Convert files to vcf format using pindel2vcf, delete empty files (where none of SV type found)
    for SV in BP CloseEndMapped D INT_final INV LI RP SI TD
    do
      if [[ $(wc -c "${OUTDIR}/2_LargeSVs/Pindel/${FILENAME}_Results/${FILENAME}.pd.out_${SV}" | awk '{print $1}') != 0 ]]; # If file is not empty
      then
        pindel2vcf -p ${OUTDIR}/2_LargeSVs/Pindel/${FILENAME}_Results/${FILENAME}.pd.out_${SV} -r ${REF} \
        -R ${REF} -d 20170000 -v ${OUTDIR}/2_LargeSVs/Pindel/vcf/${FILENAME}_${SV}.vcf
      else
        rm ${OUTDIR}/2_LargeSVs/Pindel/${FILENAME}_Results/${FILENAME}.pd.out_${SV} # delete empty files
      fi
    done

    for SV in BP CloseEndMapped D DD INV LI SI TD
    do
      if [[ $(wc -c "${OUTDIR}/2_LargeSVs/Pindel/${FILENAME}_Results/${FILENAME}.pd.q.out_${SV}" | awk '{print $1}') != 0 ]];
      then
        pindel2vcf -p ${OUTDIR}/2_LargeSVs/Pindel/${FILENAME}_Results/${FILENAME}.pd.q.out_${SV} -r ${REF} \
        -R ${REF} -d 20170000 -v ${OUTDIR}/2_LargeSVs/Pindel/vcf/${FILENAME}_${SV}.q.vcf
      else 
        rm ${OUTDIR}/2_LargeSVs/Pindel/${FILENAME}_Results/${FILENAME}.pd.q.out_${SV}
      fi    
    done
  fi

  ### Crest
  # Requires sorted and indexed bam files
  # Get soft-clipping positions, use extractSClip.pl, outputs filename.cover
  # Required parameters: 
  # -d (input bam) -f (soft clip cover file) --ref_genome (ref in fa format)
  # -t (2bit ref used for blat server) --blatserver (server name)
  # --blatport (blat server port)
  if [ ${CREST} ] ; then
    echo -e "\nRunning Crest on ${FILENAME}...\n"

    cd ${OUTDIR}/2_LargeSVs/Crest
    extractSClip.pl -i ${OUTDIR}/1_Mapping/${FILENAME}.bam --ref_genome ${REF}
    cd ${DIR}

    CREST.pl -f ${OUTDIR}/2_LargeSVs/Crest/${FILENAME}.bam.cover \
    -d ${OUTDIR}/1_Mapping/${FILENAME}.bam \
    -o ${OUTDIR}/2_LargeSVs/Crest/ \
    --ref_genome ${REF} -t ${BIT} \
    --blatserver localhost --blatport 6666
  fi

  ### Delly
  # Requires sorted indexed bam files
  # Calls variants separately -t {arg} (DEL, DUP, INV, BND, INS), default DEL
  # -n : no small InDel calling
  # -g : genome reference files
  # -o : output file
  if [ ${DELLY} ] ; then
    echo -e "\nRunning Delly on ${FILENAME}...\n"
    for SV in DEL INS DUP BND
    do
      delly call -n -t ${SV} -g ${REF} -o ${OUTDIR}/2_LargeSVs/Delly/${FILENAME}.${SV}.bcf \
      ${OUTDIR}/1_Mapping/${FILENAME}.bam
    done

    # merge output bcf files
    bcftools merge --force-samples \
    ${OUTDIR}/2_LargeSVs/Delly/${FILENAME}.DEL.bcf ${OUTDIR}/2_LargeSVs/Delly/${FILENAME}.INS.bcf \
    ${OUTDIR}/2_LargeSVs/Delly/${FILENAME}.DUP.bcf ${OUTDIR}/2_LargeSVs/Delly/${FILENAME}.BND.bcf \
    -O v > ${OUTDIR}/2_LargeSVs/Delly/${FILENAME}.vcf

    # Collate results
    echo -e "\nCollating results...\n"
    python3 collate_output.py ${OUTDIR}/2_LargeSVs ${FILENAME} ${REF##*/}

    sort -k2 -n ${OUTDIR}/2_LargeSVs/${FILENAME}.txt >> ${OUTDIR}/2_LargeSVs/all_data.txt
    echo "" >> ${OUTDIR}/2_LargeSVs/all_data.txt
  fi

done

# Send stop message to quit gfServer
gfServer stop localhost 6666

echo -e "\nAnalysis complete, building tree...\n"
# Build tree
python3 dellytree.py ${OUTDIR}/2_LargeSVs
Rscript largeSVtree.R ${OUTDIR}/2_LargeSVs

echo -e "\nFinished"