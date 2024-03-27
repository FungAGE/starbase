import subprocess
import tempfile
import pandas as pd


def run_hmmer(
    seq_type=None, tmp_hmmer=None, input_eval=None, tmp_fasta=None, threads=4
):
    hmmer_program = "hmmsearch"
    hmmer_db = "/home/project-vol/Starships/captain/tyr/hmm/YRsuperfams.p1-512.hmm"

    # Run HMMER search
    hmmer_cmd = f"{hmmer_program} -o {tmp_hmmer} --cpu {threads} --domE {input_eval} {hmmer_db} {tmp_fasta}"
    subprocess.run(hmmer_cmd, shell=True)

    # Parse HMMER output
    tmp_hmmer_parsed = tempfile.NamedTemporaryFile(suffix=".hmmer.parsed.txt").name
    parse_cmd = (
        f"python bin/hmm.py --hmmer_output_file {tmp_hmmer} --parsed {tmp_hmmer_parsed}"
    )
    subprocess.run(parse_cmd, shell=True)

    # Read parsed output into DataFrame
    hmmer_results = pd.read_csv(tmp_hmmer_parsed, sep="\t")
    return hmmer_results
