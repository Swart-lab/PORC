#!/bin/sh

#Run as: ./PORC.sh $PATH_TO_PORC_SCRIPTS HMMER_PROCESSES GENOME_FILE_PREFIX GENOME_FILE_SUFFIX PFAM_PATH
#with genome file, e.g. mygenome.fasta, specified as two separate arguments,
#i.e. mygenome fasta
#PFAM_PATH is the directory that contains all the Pfam-A.hmm derived database
#files produced by HMMER's hmmpress command

#Example command:
#./PORC.sh /path_to_my_porc 32 a_genome fa ../db

PORC_PATH=$1
HMMER_CPU=$2
GENOME=$3.$4
CDS=$3.cds
PEP=$3.pep
WEBLOGO_MATRIX=$3.mat
HMMER_OUT=$PEP.hmmer_out
PORC_OUT=$3.porc
PFAM_PATH=$5
PFAM_DB=$PFAM_PATH/Pfam-A.hmm

SIX_FRAME_TRANSLATE='$PORC_PATH/six_frame_pep_and_cds.py $GENOME $CDS $PEP'
echo $SIX_FRAME_TRANSLATE
eval $SIX_FRAME_TRANSLATE

HMMER_COMMAND='hmmsearch --cpu $HMMER_CPU -o $HMMER_OUT $PFAM_DB $PEP'
echo $HMMER_COMMAND
eval $HMMER_COMMAND

MAIN_COMMAND='$PORC_PATH/porc_cod_usage.py $CDS $HMMER_OUT $WEBLOGO_MATRIX > $PORC_OUT'
echo $MAIN_COMMAND
eval $MAIN_COMMAND

WEBLOGO_COMMAND='weblogo -n 64 --scale-width no -c chemistry -U probability -A protein -F pdf < $WEBLOGO_MATRIX > $GENOME.weblogo.pdf'
echo $WEBLOGO_COMMAND
eval $WEBLOGO_COMMAND
