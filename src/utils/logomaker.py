# standard imports
import io
import base64
import tempfile
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import logomaker as lm

from Bio import AlignIO, SeqIO
from Bio.Align.Applications import ClustalwCommandline
from io import StringIO


def make_logo(seqs):  # Create a temporary file to hold the input sequences
    temp_file = tempfile.NamedTemporaryFile(suffix=".fa", delete=False)
    with open(temp_file.name, "w") as file:
        for idx, seq in enumerate(seqs):
            header = f">seq{idx + 1}"
            file.write(f"{header}\n{seq}\n")

    clustalw_cline = ClustalwCommandline("clustalw", infile=temp_file.name)
    clustalw_cline()
    print(temp_file.name)

    aligner = AlignIO.read(temp_file.name, "fasta")

    aligned_sequences = []
    for record in aligner:
        aligned_sequences.append(record.seq)

    fig, ax = plt.subplots()
    counts_mat = lm.alignment_to_matrix(
        sequences=aligned_sequences, to_type="counts", characters_to_ignore=".-X"
    )

    lm.Logo(counts_mat, ax=ax)

    # Save the figure to a BytesIO object
    buf = io.BytesIO()
    fig.savefig(buf, format="png")
    buf.seek(0)

    # Encode the image to base64
    img_base64 = base64.b64encode(buf.getvalue()).decode("utf-8")

    plt.close(fig)  # Close the figure to free up resources
    return img_base64
