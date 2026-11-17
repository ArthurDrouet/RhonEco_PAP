#!/bin/bash
#SBATCH --job-name=Preprocessing_cutadapt_PAP
#SBATCH --mem=8G
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --partition=normal
#SBATCH --output=/beegfs/home/adrouet/RHONECO/MB_PAP_LTS_MultiTaxa/Script/2.Traitement_de_données/1_Preprocessing/Log_Preprocessing.out
#SBATCH --error=/beegfs/home/adrouet/RHONECO/MB_PAP_LTS_MultiTaxa/Script/2.Traitement_de_données/1_Preprocessing/Log_Preprocessing.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arthur.drouet@inrae.fr

set -euo pipefail
shopt -s nullglob

#### PRETRAITEMENT PIPELINE ####
# Date :   11/2026
# Lieu :   Lyon, France
# Projet : Piège à sédiments métabaroding
# Caractèrisation de diverisité multitaxonomique (Diatomées, Insectes, Poissons) dans les sédiments capturé par des pièges à sédiments dans le Rhône.

# Auteur : Arthur DROUET
# Equipe : Johan PANSU (LEHNA, CNRS), Arthur DROUET (RIVERLY, INRAE)
# Titre :  Preprocess Script

# Objectif :
# Utiliser Cutadapt pour retirer les amorces des reads et effectuer un contrôle qualité avec FastQC et MultiQC.

# Steps of Preprocessing :
# 1. Paths and working directory
# 2. Select Primer & Tag
# 3. Demultiplexing
# 4. RemovePrimer
# 5. Number of reads
# 6. FastQC
# 7. MultiQC

#####################################
###### CutAdapt version : 3.5  ######
###### FastQC version : 0.11.9 ######
###### MultiQC version : 1.33  ######
#####################################

### 1. Paths and working directory ####
echo "- [STEP 1] : Paths and working directory - "

WORKDIR="/beegfs/home/adrouet/RHONECO/MB_PAP_LTS_MultiTaxa"
## - Choix du jeu de donneé - ##
DATA="Raw_miniData" # Options : Raw_data, Raw_miniData
# Raw_data = donnée
# Raw_miniData = petit jeu de donnée pour les tests

INDIR="$WORKDIR/Data/$DATA"
echo "Data selctionné : $INDIR "

OUTDIR="$WORKDIR/Results/2.Traitement_de_données/$DATA/1_Preprocessing"
echo "Chemin d'écriture : $OUTDIR "

## - Quality checks (FastQC and MultiQC): choisir "YES" or "NO" - ##
RUN_QUALITY_CHECK="YES"
echo " Quality Check activé ? : $RUN_QUALITY_CHECK "

# Log directory
LOGDIR="$OUTDIR/logs"

# Demultiplexing directory :
DEMULTDIR="$OUTDIR/Demultiplexing"

# Remouve primer directory :
RMPRIMERDIR="$OUTDIR/RemovePrimer"
COUNTDIR="$OUTDIR/NumberReads"

# Build directory :
mkdir -p "$OUTDIR" "$LOGDIR" "$RMPRIMERDIR" "$DEMULTDIR" "$COUNTDIR"

# Quality check directory :
case "$RUN_QUALITY_CHECK" in
  YES)
    FASTQCDIR="$OUTDIR/FastQC"
    MULTIQCDIR="$OUTDIR/MultiQC"
    mkdir -p "$FASTQCDIR" "$MULTIQCDIR"
    ;;
  NO)
    echo "Quality check désactivé (/FastQC et /MultiQC seront ignorés)."
    ;;
  *)
    echo "Erreur : RUN_QUALITY_CHECK doit être défini à YES ou NO." >&2
    exit 1
    ;;
esac

# Environnement directory
ENVDIR="/beegfs/home/adrouet/Env"

# Environnement logiciel
CUTADAPT_SIF="$ENVDIR/cutadapt_3.5.sif"
FASTQC_SIF="$ENVDIR/fastqc_0.11.9.sif"
MULTIQC_SIF="$ENVDIR/multiqc_v1.33.sif"

# Fonction counts :
count_pairs_gz () {
  if [ -s "$1" ]; then
    zcat "$1" | awk 'END{print NR/4}'
  else
    echo 0
  fi
}


### Function : sample_name_from_r1 ####
sample_name_from_r1() {
    local R1="$1"
    local BN
    BN="$(basename "$R1")"
    if [[ "$BN" =~ ^(.+)_S[0-9]+_L[0-9]+_R1_001\.fastq\.gz$ ]]; then
        echo "${BASH_REMATCH[1]}"
    else
        echo "ERROR : format de fichier R1 non reconnu : $BN" >&2
        return 1
    fi
}


### 2. Primer & Tags####
echo "- [STEP 2] : Select Primer & Tag-"

##############################################
echo " For Diatomées"
## Primer selection for Diatomées
echo " -- Primer selection for Diatomées -- "

FORWARD_PRIMER_DIA="ATTCCAGCTCCAATAGCGTA"
REVERSE_PRIMER_DIA="CTTTRAACRCKCTRATTTYTTCAC"

echo "FORWARD PRIMER : $FORWARD_PRIMER_DIA"
echo "REVERSE PRIMER : $REVERSE_PRIMER_DIA"

## Tag selection for Diatomées
echo " -- Tag selection for Diatomées -- "

TAG_DIA_1_R1="GTCTCGTC"
TAG_DIA_1_R2="GTCTCGTC" #-> "CTGCTCTG"
TAG_DIA_2_R1="GAACGACC"
TAG_DIA_2_R2="GAACGACC" #-> "CCAGCAAG"
TAG_DIA_3_R1="CTACTTAC"
TAG_DIA_3_R2="CTACTTAC" #-> "CATTCATC"

echo "TAG 1 : $TAG_DIA_1_R1"
echo "TAG 2 : $TAG_DIA_2_R1"
echo "TAG 3 : $TAG_DIA_3_R1"

##############################################
echo " For Insectes"

## Primer selection for Insectes
echo " -- Primer selection for Insectes -- "

FORWARD_PRIMER_INSE="RGACGAGAAGACCCTATARA"
REVERSE_PRIMER_INSE="ACGCTGTTATCCCTAARGTA"

echo "FORWARD PRIMER : $FORWARD_PRIMER_INSE"
echo "REVERSE PRIMER : $REVERSE_PRIMER_INSE"

## Tag selection for Insectes
echo " -- Tag selection for Insectes -- "

TAG_INSE_1_R1="AAGGACAC"
TAG_INSE_1_R2="AAGGACAC" #-> "GTGTCCTT"
TAG_INSE_2_R1="GATCGCAT"
TAG_INSE_2_R2="GATCGCAT" #-> "TACGCTAG"
TAG_INSE_3_R1="CCAGTAGC"
TAG_INSE_3_R2="CCAGTAGC" #-> "CGATGACC"

echo "TAG 1 : $TAG_INSE_1_R1 "
echo "TAG 2 : $TAG_INSE_2_R1"
echo "TAG 3 : $TAG_INSE_3_R1"

##############################################
echo " For Fish"

## Primer selection for Fish
echo " -- Primer selection for Fish -- "

FORWARD_PRIMER_FISH="ACACCGCCCGTCACTCT"
REVERSE_PRIMER_FISH="CTTCCGGTACACTTACCATG"

echo "FORWARD PRIMER : $FORWARD_PRIMER_FISH"
echo "REVERSE PRIMER : $REVERSE_PRIMER_FISH"

## Tag selection for Fish
echo " -- Tag selection for Fish -- "

# pas de tag 1 pour les poissons
TAG_FISH_1_R1=TCGTTCAC
TAG_FISH_1_R2=TTCGTGCT #-> "TCGTGCTT"
TAG_FISH_2_R1="TCGTTCAC"
TAG_FISH_2_R2="TCGTTCAC" #-> "CACTTGCT"
TAG_FISH_3_R1="TTCGTGCT"
TAG_FISH_3_R2="TTCGTGCT" #-> "TCGTGCTT"

echo "TAG 1 : $TAG_FISH_1_R1 - $TAG_FISH_1_R2"
echo "TAG 2 : $TAG_FISH_2_R1"
echo "TAG 3 : $TAG_FISH_3_R1"


### 3. Demultiplexing  ####
echo "- [STEP 3] : Demultiplexing   -"

# Pour chaque échantillon, on démultiplexe les réplicats (tags) par taxon.
# Pour les échantillons RH : DIA, INSE et FISH.
# Pour les autres échantillons : DIA uniquement.
#
# Le démultiplexage est séquentiel :
# A : lectures correspondant au tag A / lectures restantes
# B : recherche du tag B parmi les lectures restantes
# C : recherche du tag C parmi les lectures restantes
# Les lectures restantes après C correspondent aux NoTag.

# 1. Récupérer les noms des échantillons :
R1_FILES=("$INDIR"/*_R1_001.fastq.gz)

for R1 in "${R1_FILES[@]}"; do

  SAMPLE="$(sample_name_from_r1 "$R1")"
  R2="${R1%_R1_001.fastq.gz}_R2_001.fastq.gz"

  echo "Demultiplexing for sample : $SAMPLE"

  # Vérifier la présence du fichier R2
  if [[ ! -f "$R2" ]]; then
    echo "WARNING : R2 missing for $SAMPLE" >&2
    continue
  fi

  if [[ "$SAMPLE" == RH* ]]; then

    # Demultiplexing for RH samples
    ListeTaxon=("DIA" "INSE" "FISH")

    for TAX in "${ListeTaxon[@]}"; do

      echo "Demultiplexing for taxon : $TAX"

      # Récupérer les séquences des tags correspondant au taxon
      TAG_R1_A_NAME="TAG_${TAX}_1_R1"
      TAG_R2_A_NAME="TAG_${TAX}_1_R2"
      TAG_R1_B_NAME="TAG_${TAX}_2_R1"
      TAG_R2_B_NAME="TAG_${TAX}_2_R2"
      TAG_R1_C_NAME="TAG_${TAX}_3_R1"
      TAG_R2_C_NAME="TAG_${TAX}_3_R2"

      TAG_R1_A="${!TAG_R1_A_NAME}"
      TAG_R2_A="${!TAG_R2_A_NAME}"
      TAG_R1_B="${!TAG_R1_B_NAME}"
      TAG_R2_B="${!TAG_R2_B_NAME}"
      TAG_R1_C="${!TAG_R1_C_NAME}"
      TAG_R2_C="${!TAG_R2_C_NAME}"

      # Réplicat A
      echo "Demultiplexing replicate A for taxon : $TAX for sample : $SAMPLE"

      apptainer exec --bind /beegfs:/beegfs "$CUTADAPT_SIF" \
        cutadapt \
          -j 8 \
          -e 0 \
          -g "$TAG_R1_A" \
          -G "$TAG_R2_A" \
          -o "$DEMULTDIR/${SAMPLE}_${TAX}_A_R1_Demultiplexed.fastq.gz" \
          -p "$DEMULTDIR/${SAMPLE}_${TAX}_A_R2_Demultiplexed.fastq.gz" \
          --untrimmed-output "$DEMULTDIR/${SAMPLE}_${TAX}_NoTag_A_R1.fastq.gz" \
          --untrimmed-paired-output "$DEMULTDIR/${SAMPLE}_${TAX}_NoTag_A_R2.fastq.gz" \
          "$R1" "$R2" \
          >> "$LOGDIR/${SAMPLE}_${TAX}_A_Demultiplexing.log" 2>&1

      # Réplicat B
      echo "Demultiplexing replicate B for taxon : $TAX for sample : $SAMPLE"

      apptainer exec --bind /beegfs:/beegfs "$CUTADAPT_SIF" \
        cutadapt \
          -j 8 \
          -e 0 \
          -g "$TAG_R1_B" \
          -G "$TAG_R2_B" \
          -o "$DEMULTDIR/${SAMPLE}_${TAX}_B_R1_Demultiplexed.fastq.gz" \
          -p "$DEMULTDIR/${SAMPLE}_${TAX}_B_R2_Demultiplexed.fastq.gz" \
          --untrimmed-output "$DEMULTDIR/${SAMPLE}_${TAX}_NoTag_B_R1.fastq.gz" \
          --untrimmed-paired-output "$DEMULTDIR/${SAMPLE}_${TAX}_NoTag_B_R2.fastq.gz" \
          "$DEMULTDIR/${SAMPLE}_${TAX}_NoTag_A_R1.fastq.gz" "$DEMULTDIR/${SAMPLE}_${TAX}_NoTag_A_R2.fastq.gz" \
          >> "$LOGDIR/${SAMPLE}_${TAX}_B_Demultiplexing.log" 2>&1

      # Réplicat C
      echo "Demultiplexing replicate C for taxon : $TAX for sample : $SAMPLE"

      apptainer exec --bind /beegfs:/beegfs "$CUTADAPT_SIF" \
        cutadapt \
          -j 8 \
          -e 0 \
          -g "$TAG_R1_C" \
          -G "$TAG_R2_C" \
          -o "$DEMULTDIR/${SAMPLE}_${TAX}_C_R1_Demultiplexed.fastq.gz" \
          -p "$DEMULTDIR/${SAMPLE}_${TAX}_C_R2_Demultiplexed.fastq.gz" \
          --untrimmed-output "$DEMULTDIR/${SAMPLE}_${TAX}_R1_NoTag.fastq.gz" \
          --untrimmed-paired-output "$DEMULTDIR/${SAMPLE}_${TAX}_R2_NoTag.fastq.gz" \
          "$DEMULTDIR/${SAMPLE}_${TAX}_NoTag_B_R1.fastq.gz" "$DEMULTDIR/${SAMPLE}_${TAX}_NoTag_B_R2.fastq.gz" \
          >> "$LOGDIR/${SAMPLE}_${TAX}_C_Demultiplexing.log" 2>&1

    done

  else

    # Demultiplexing for non-RH samples : DIA only

    # Réplicat A
    echo "Demultiplexing replicate A of the sample : $SAMPLE"

    apptainer exec --bind /beegfs:/beegfs "$CUTADAPT_SIF" \
      cutadapt \
        -j 8 \
        -e 0 \
        -g "$TAG_DIA_1_R1" \
        -G "$TAG_DIA_1_R2" \
        -o "$DEMULTDIR/${SAMPLE}_A_R1_Demultiplexed.fastq.gz" \
        -p "$DEMULTDIR/${SAMPLE}_A_R2_Demultiplexed.fastq.gz" \
        --untrimmed-output "$DEMULTDIR/${SAMPLE}_NoTag_A_R1.fastq.gz" \
        --untrimmed-paired-output "$DEMULTDIR/${SAMPLE}_NoTag_A_R2.fastq.gz" \
        "$R1" "$R2" \
        >> "$LOGDIR/${SAMPLE}_A_Demultiplexing.log" 2>&1

    # Réplicat B
    echo "Demultiplexing replicate B of the sample : $SAMPLE"

    apptainer exec --bind /beegfs:/beegfs "$CUTADAPT_SIF" \
      cutadapt \
        -j 8 \
        -e 0 \
        -g "$TAG_DIA_2_R1" \
        -G "$TAG_DIA_2_R2" \
        -o "$DEMULTDIR/${SAMPLE}_B_R1_Demultiplexed.fastq.gz" \
        -p "$DEMULTDIR/${SAMPLE}_B_R2_Demultiplexed.fastq.gz" \
        --untrimmed-output "$DEMULTDIR/${SAMPLE}_NoTag_B_R1.fastq.gz" \
        --untrimmed-paired-output "$DEMULTDIR/${SAMPLE}_NoTag_B_R2.fastq.gz" \
        "$DEMULTDIR/${SAMPLE}_NoTag_A_R1.fastq.gz" "$DEMULTDIR/${SAMPLE}_NoTag_A_R2.fastq.gz" \
        >> "$LOGDIR/${SAMPLE}_B_Demultiplexing.log" 2>&1

    # Réplicat C
    echo "Demultiplexing replicate C of the sample : $SAMPLE"

    apptainer exec --bind /beegfs:/beegfs "$CUTADAPT_SIF" \
      cutadapt \
        -j 8 \
        -e 0 \
        -g "$TAG_DIA_3_R1" \
        -G "$TAG_DIA_3_R2" \
        -o "$DEMULTDIR/${SAMPLE}_C_R1_Demultiplexed.fastq.gz" \
        -p "$DEMULTDIR/${SAMPLE}_C_R2_Demultiplexed.fastq.gz" \
        --untrimmed-output "$DEMULTDIR/${SAMPLE}_NoTag_R1.fastq.gz" \
        --untrimmed-paired-output "$DEMULTDIR/${SAMPLE}_NoTag_R2.fastq.gz" \
        "$DEMULTDIR/${SAMPLE}_NoTag_B_R1.fastq.gz" "$DEMULTDIR/${SAMPLE}_NoTag_B_R2.fastq.gz" \
        >> "$LOGDIR/${SAMPLE}_C_Demultiplexing.log" 2>&1

  fi

done

### 4. RemovePrimer ####
echo " - [STEP 4] : Remove Primer - "

# Récupérer les fichiers R1 démultiplexés
R1_FILES=("$DEMULTDIR"/*_R1_Demultiplexed.fastq.gz)

for R1 in "${R1_FILES[@]}"; do
    # Nom du fichier
    BN="$(basename "$R1")"
    # Récupérer le fichier R2 correspondant
    R2="${R1%_R1_Demultiplexed.fastq.gz}_R2_Demultiplexed.fastq.gz"

    # Préfixe du fichier (conserve taxon et réplicat)
    PREFIX="${BN%_R1_Demultiplexed.fastq.gz}"

    echo "Remove Primer for : $PREFIX"
    # Vérifier la présence du R2
    if [[ ! -f "$R2" ]]; then
        echo "WARNING : R2 missing for $PREFIX" >&2
        continue
    fi

    if [[ "$BN" == RH* ]]; then
        # Identifier le taxon dans le nom du fichier
        TAX="$(echo "$BN" | sed -E 's/^RH.*_(DIA|INSE|FISH)_[ABC]_R1_Demultiplexed.fastq.gz$/\1/')"
        # Vérifier que le taxon a été identifié
        if [[ "$TAX" != "DIA" && "$TAX" != "INSE" && "$TAX" != "FISH" ]]; then
            echo "ERROR : Taxon non reconnu pour $BN" >&2
            exit 1
        fi

        # Sélectionner les amorces adaptées
        forward_name="FORWARD_PRIMER_${TAX}"
        reverse_name="REVERSE_PRIMER_${TAX}"
        FORWARD_PRIMER="${!forward_name}"
        REVERSE_PRIMER="${!reverse_name}"

    else
        # Pour les autres échantillons : Diatomées
        FORWARD_PRIMER="$FORWARD_PRIMER_DIA"
        REVERSE_PRIMER="$REVERSE_PRIMER_DIA"

    fi

    # Remove Primer
    apptainer exec --bind /beegfs:/beegfs "$CUTADAPT_SIF" \
      cutadapt \
        -j 8 \
        -e 2 \
        -g "$FORWARD_PRIMER" \
        -G "$REVERSE_PRIMER" \
        --discard-untrimmed \
        -o "$RMPRIMERDIR/${PREFIX}_R1_RemovePrimer.fastq.gz" \
        -p "$RMPRIMERDIR/${PREFIX}_R2_RemovePrimer.fastq.gz" \
        "$R1" "$R2" \
        >> "$LOGDIR/${PREFIX}_RemovePrimer.log" 2>&1

done

### 5. Number of reads ####
echo "- [STEP 5] : Number of reads -"

COUNTFILE="$COUNTDIR/NumberReads.tsv"
SUMMARYFILE="$COUNTDIR/Reads_Summary.tsv"

# Autoriser les tableaux de fichiers vides si aucun fichier ne correspond au motif
shopt -s nullglob

# Récupérer les fichiers bruts
RAW_R1_FILES=("$INDIR"/*_R1_001.fastq.gz)

if (( ${#RAW_R1_FILES[@]} == 0 )); then
  echo "Aucun fichier R1 brut trouvé dans $INDIR" >&2
  exit 1
fi

echo -e "Stage\tSample\tTaxon\tReplicate\tR1_reads\tR2_reads" > "$COUNTFILE"

# Pour chaque échantillon
for RAW_R1 in "${RAW_R1_FILES[@]}"; do

  SAMPLE="$(sample_name_from_r1 "$RAW_R1")"
  RAW_R2="${RAW_R1%_R1_001.fastq.gz}_R2_001.fastq.gz"

  if [[ ! -f "$RAW_R2" ]]; then
    echo "Fichier R2 manquant : $RAW_R2" >&2
    exit 1
  fi

  echo "Counting reads for sample : $SAMPLE"

    # 1. Raw
    echo -e "Raw\t${SAMPLE}\tNA\tNA\t$(count_pairs_gz "$RAW_R1")\t$(count_pairs_gz "$RAW_R2")" >> "$COUNTFILE"

    # 2. Demultiplexed : A, B et C
    for REP in A B C; do

        if [[ "$SAMPLE" == RH* ]]; then
            PREFIX="${SAMPLE}_${TAX}_${REP}"
        else
            PREFIX="${SAMPLE}_${REP}"
        fi

        DEMULT_R1="$DEMULTDIR/${PREFIX}_R1_Demultiplexed.fastq.gz"
        DEMULT_R2="$DEMULTDIR/${PREFIX}_R2_Demultiplexed.fastq.gz"

        echo -e "Demultiplexed\t${SAMPLE}\t${TAX}\t${REP}\t$(count_pairs_gz "$DEMULT_R1")\t$(count_pairs_gz "$DEMULT_R2")" >> "$COUNTFILE"

        # RemovePrimer
        REMOVE_R1="$RMPRIMERDIR/${PREFIX}_R1_RemovePrimer.fastq.gz"
        REMOVE_R2="$RMPRIMERDIR/${PREFIX}_R2_RemovePrimer.fastq.gz"

        echo -e "RemovePrimer\t${SAMPLE}\t${TAX}\t${REP}\t$(count_pairs_gz "$REMOVE_R1")\t$(count_pairs_gz "$REMOVE_R2")" >> "$COUNTFILE"

    done

    # 3. NoTag final uniquement
    if [[ "$SAMPLE" == RH* ]]; then
        NOTAG_R1="$DEMULTDIR/${SAMPLE}_${TAX}_R1_NoTag.fastq.gz"
        NOTAG_R2="$DEMULTDIR/${SAMPLE}_${TAX}_R2_NoTag.fastq.gz"
    else
        NOTAG_R1="$DEMULTDIR/${SAMPLE}_NoTag_R1.fastq.gz"
        NOTAG_R2="$DEMULTDIR/${SAMPLE}_NoTag_R2.fastq.gz"
    fi

    echo -e "NoTag\t${SAMPLE}\t${TAX}\tNA\t$(count_pairs_gz "$NOTAG_R1")\t$(count_pairs_gz "$NOTAG_R2")" >> "$COUNTFILE"

    done

    # 4. NoTag final (lectures non attribuées après A, B et C)
    if [[ "$SAMPLE" == RH* ]]; then
      NOTAG_R1="$DEMULTDIR/${SAMPLE}_${TAX}_R1_NoTag.fastq.gz"
      NOTAG_R2="$DEMULTDIR/${SAMPLE}_${TAX}_R2_NoTag.fastq.gz"
    else
      NOTAG_R1="$DEMULTDIR/${SAMPLE}_NoTag_R1.fastq.gz"
      NOTAG_R2="$DEMULTDIR/${SAMPLE}_NoTag_R2.fastq.gz"
    fi

    if [[ ! -f "$NOTAG_R1" || ! -f "$NOTAG_R2" ]]; then
      echo "Fichier NoTag final manquant pour ${SAMPLE} - ${TAX}." >&2
      exit 1
    fi

    echo -e "NoTag\t${SAMPLE}\t${TAX}\tNA\t$(count_pairs_gz "$NOTAG_R1")\t$(count_pairs_gz "$NOTAG_R2")" >> "$COUNTFILE"


# Totaux globaux par étape
echo -e "Stage\tTotal_R1\tTotal_R2" > "$SUMMARYFILE"

for STAGE in Raw Demultiplexed RemovePrimer NoTag; do

  SUM_R1=$(awk -F'\t' -v s="$STAGE" '$1==s {sum+=$5} END{print sum+0}' "$COUNTFILE")
  SUM_R2=$(awk -F'\t' -v s="$STAGE" '$1==s {sum+=$6} END{print sum+0}' "$COUNTFILE")

  echo -e "${STAGE}\t${SUM_R1}\t${SUM_R2}" >> "$SUMMARYFILE"

done

echo "Nombre de reads par échantillon (fichier : $COUNTFILE) :"
cat "$COUNTFILE"

echo "Totaux par étape (fichier : $SUMMARYFILE) :"
cat "$SUMMARYFILE"


### 6. FastQC ####

if [[ "$RUN_QUALITY_CHECK" == "YES" ]]; then

  echo "- [STEP 6] : FastQC -"

  # Autoriser les tableaux de fichiers vides
  shopt -s nullglob

  # Créer les dossiers de sortie
  mkdir -p "$OUTDIR/FastQC" "$OUTDIR/MultiQC"

  # Récupérer les fichiers R1 après retrait des amorces
  REMOVE_R1_FILES=("$RMPRIMERDIR"/*_R1_RemovePrimer.fastq.gz)

  if (( ${#REMOVE_R1_FILES[@]} == 0 )); then

    echo "Aucun fichier RemovePrimer à analyser : FastQC et MultiQC ignorés."

  else

    FASTQC_HAS_REPORTS=0
    : > "$LOGDIR/FastQC.log"

    for REMOVE_R1 in "${REMOVE_R1_FILES[@]}"; do

      REMOVE_R2="${REMOVE_R1%_R1_RemovePrimer.fastq.gz}_R2_RemovePrimer.fastq.gz"

      if [[ ! -f "$REMOVE_R2" ]]; then
        echo "Fichier R2 RemovePrimer manquant : $REMOVE_R2" >&2
        exit 1
      fi

      # Analyser R1 et R2 séparément
      for FASTQ_FILE in "$REMOVE_R1" "$REMOVE_R2"; do

        READ_COUNT="$(count_pairs_gz "$FASTQ_FILE")"

        if [[ ! "$READ_COUNT" =~ ^[0-9]+$ ]]; then
          echo "Erreur : nombre de reads invalide pour $FASTQ_FILE : $READ_COUNT" >&2
          exit 1
        fi

        if (( READ_COUNT == 0 )); then
          echo "FastQC ignoré pour $FASTQ_FILE (0 reads après retrait des amorces)."
          continue
        fi

        FASTQ_BASENAME="${FASTQ_FILE##*/}"

        echo "FastQC : $FASTQ_BASENAME"

        # /work évite le problème d'encodage du chemin contenant un accent
        if ! apptainer exec --bind "$OUTDIR:/work" "$FASTQC_SIF" \
          fastqc -t 8 -o /work/FastQC "/work/RemovePrimer/$FASTQ_BASENAME" \
          >> "$LOGDIR/FastQC.log" 2>&1; then

          echo "ERREUR : FastQC a échoué pour $FASTQ_BASENAME" >&2
          exit 1
        fi

        FASTQC_HAS_REPORTS=1

      done
    done

    ### 7. MultiQC ####

    echo "- [STEP 7] : MultiQC -"

    if (( FASTQC_HAS_REPORTS )); then

      if ! apptainer exec --bind "$OUTDIR:/work" "$MULTIQC_SIF" \
        multiqc /work/FastQC -o /work/MultiQC \
        > "$LOGDIR/MultiQC.log" 2>&1; then

        echo "ERREUR : MultiQC a échoué." >&2
        exit 1
      fi

    else

      echo "Aucun rapport FastQC à regrouper avec MultiQC."

    fi

  fi

else

  echo "- [STEP 6 & 7] : Quality Control -> Ignored"

fi
