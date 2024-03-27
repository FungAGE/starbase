from Bio import SeqIO
import re

def clean_lines(query_list):
    cleaned_seq = ""
    nucl_char = set("ATGC")
    prot_char = set("ARNDBCEQZGHILKMFPSTWYV")

    with open(query_list["tmp_file"], "r") as f:
        # Read the sequence records from the FASTA file
        records = SeqIO.parse(f, "fasta")

        # Concatenate sequence lines, removing non-letter characters
        for record in records:
            for line in record.seq:
                cleaned_seq += re.sub("[^A-Za-z]", "", str(line))

    # Count characters for deciding type later
    nucl_count = sum(cleaned_seq.count(base) for base in nucl_char)
    prot_count = sum(cleaned_seq.count(aa) for aa in prot_char)

    # Guess if sequence is nucleotide or protein
    if prot_count >= (0.1 * len(cleaned_seq)):
        query_type = "prot"
        print("Query is protein sequence")
    else:
        query_type = "nucl"
        print("Query is nucleotide sequence")

    return {
        **query_list,
        "query_type": query_type,
        "cleaned_query": cleaned_seq
    }