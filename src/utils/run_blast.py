from Bio.Blast.Applications import (
    NcbiblastnCommandline,
    NcbiblastpCommandline,
    NcbitblastnCommandline,
)
import tempfile
import pandas as pd


def run_blast(
    seq_type=None,
    blast_type=None,
    tmp_fasta=None,
    input_eval=None,
    threads=None,
    stitch=False,
):
    # Define paths to databases
    db_list = {
        "ship": {"nucl": "project-vol/Starships/ships/fna/blastdb/concatenated.fa"},
        "gene": {
            "tyr": {
                "prot": "project-vol/Starships/captain/tyr/faa/blastdb/concatenated.faa"
            },
            "fre": {
                "prot": "project-vol/Starships/cargo/fre/faa/blastdb/fre.mycoDB.faa",
                "nucl": "project-vol/Starships/cargo/fre/fna/blastdb/fre.fa",
            },
            "nlr": {
                "prot": "project-vol/Starships/cargo/nlr/faa/blastdb/nlr.mycoDB.faa",
                "nucl": "project-vol/Starships/cargo/nlr/fna/blastdb/nlr.fa",
            },
            "DUF3723": {
                "prot": "project-vol/Starships/cargo/duf3723/faa/blastdb/duf3723.mycoDB.faa",
                "nucl": "project-vol/Starships/cargo/duf3723/fna/blastdb/duf3723.fa",
            },
            "plp": {
                "prot": "project-vol/Starships/cargo/plp/faa/blastdb/plp.mycoDB.faa",
                "nucl": "project-vol/Starships/cargo/plp/fna/blastdb/plp.fa",
            },
        },
    }

    # Determine blast program and database based on sequence type and blast type
    if seq_type == "nucl":
        if blast_type == "ship":
            blast_program = NcbitblastnCommandline
            blastdb = db_list["ship"]["nucl"]
        else:
            blast_program = NcbiblastnCommandline
            blastdb = db_list["gene"][blast_type]["nucl"]
    else:
        if blast_type == "ship":
            blast_program = NcbitblastnCommandline
            blastdb = db_list["ship"]["nucl"]
        else:
            blast_program = NcbiblastpCommandline
            blastdb = db_list["gene"][blast_type]["prot"]

    if len(blastdb) != 1:
        raise ValueError("Issue accessing BLAST database")

    # Perform the BLAST search
    blast_tmp = tempfile.NamedTemporaryFile(suffix=".blast").name
    blast_cline = blast_program(
        query=tmp_fasta,
        db=blastdb,
        evalue=input_eval,
        out=blast_tmp,
        outfmt=6,
        num_threads=threads,
    )
    stdout, stderr = blast_cline()

    # Read the BLAST results into a DataFrame
    blast_results = pd.read_csv(
        blast_tmp,
        sep="\t",
        names=[
            "query_id",
            "hit_IDs",
            "pident",
            "aln_length",
            "mismatches",
            "gaps",
            "query_start",
            "query_end",
            "subject_start",
            "subject_end",
            "evalue",
            "bitscore",
            "query_seq",
            "subject_seq",
        ],
    )

    # Optionally stitch BLAST results
    if stitch:
        stitched_blast_tmp = tempfile.NamedTemporaryFile(suffix=".stitch").name
        stitch_blast_cmd = (
            f"python bin/BLASTstitcher.py -i {blast_tmp} -o {stitched_blast_tmp}"
        )
        subprocess.run(stitch_blast_cmd, shell=True)
        ship_blast_out = stitched_blast_tmp

    return blast_results
