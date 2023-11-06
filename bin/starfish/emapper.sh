#!/bin/bash
#SBATCH -A naiss2023-22-643 
#SBATCH -p core
#SBATCH -n 6
#SBATCH -t 0-00:20:00
#SBATCH --job-name=emapper
#SBATCH --output=/home/adrianf/systbio-project-folder/Pvar_starfish/emapper/emapper_%A.out

threads=8
run_name=$1
PROTEOME=$2
STRAIN=$3

if [[ $USER == "adrianf" ]]; then
  module load bioinfo-tools eggNOG-mapper
  ROOT_DIR=/home/adrianf/systbio-project-folder
else
  ROOT_DIR=/mnt/sda/johannesson_lab/adrian
fi

cd $ROOT_DIR || exit

WORKDIR=$ROOT_DIR/$run_name
EGGNOG_DIR=$WORKDIR/emapper
OUTDIR=$EGGNOG_DIR/$STRAIN
INDIR=$EGGNOG_DIR/input

mkdir -p "$OUTDIR"

# run eggnog mapper
# -i FASTA_FILE         Input FASTA file containing query sequences (proteins
                      # by default; see --itype and --translate). Required
                      # unless -m no_search. (default: None)

# -m MODE
                      # how input queries will be searched against eggNOG sequences. Default is -m diamond.

# --dmnd_algo [auto, 0, 1, ctg]
                      # Diamond's --algo option. Use --dmnd_algo ctg to perform a faster search when the input set of sequences is small. Default is Diamond's default. (since version 2.1.6).

# --data_dir DIR
#                       Specify a path to the eggNOG-mapper databases. By default, data/ folder or the one specified by the EGGNOG_DATA_DIR environment variable. The annotation DBs must be always present in this directory, whereas HMMER, Diamond and MMseqs2 DBs could be in a different directory when using -d/--database, --dmnd_db and --mmseqs_db, respectively.

# -d HMMER_DB_PREFIX, --database HMMER_DB_PREFIX
#                       specify the target database for sequence searches.
#                       Choose among: euk,bact,arch, or a database loaded in a
#                       server, db.hmm:host:port (see hmm_server.py) (default:
#                       None)

# --usemem              Use this option to allocate the whole database (-d) in
#                       memory using hmmpgmd. If --dbtype hmm, the database
#                       must be a hmmpress-ed database. If --dbtype seqdb, the
#                       database must be a HMMER-format database created with
#                       esl-reformat. Database will be unloaded after
#                       execution. Note that this only works for HMMER based
#                       searches. To load the eggnog-mapper annotation DB into
#                       memory use --dbmem. (default: False)

# emapper.py --override -m mmseqs --output_dir "$OUTDIR" -i "$PROTEOME" -d fuNOG --scratch_dir "$SNIC_TMP" --usemem --cpu "$threads"

if [[ "$PROTEOME" == *".cds.fasta" ]]; then
  emapper.py --override -o "$STRAIN" --output_dir "$OUTDIR" --itype CDS --translate -i "$PROTEOME" --scratch_dir "$SNIC_TMP" --cpu "$threads"
elif [[ "$PROTEOME" == *".faa" ]]; then
  emapper.py --override -o "$STRAIN" --output_dir "$OUTDIR" -i "$PROTEOME" --scratch_dir "$SNIC_TMP" --cpu "$threads"
else
  emapper.py --override -o "$STRAIN" --output_dir "$OUTDIR" --itype genome -i "$PROTEOME" --scratch_dir "$SNIC_TMP" --cpu "$threads"
fi