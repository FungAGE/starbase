from Bio import SeqIO

# Path to the HMMER output file
hmmer_output_file = "/home/adrian/Systematics/output.txt"

# Path to the query sequence file in FASTA format
query_sequence_file = "/home/adrian/Systematics/Starship_Database/Starships/genes/cap_tyr/YRsuperfamRefs.faa.split/YRsuperfamRefs.part_aaoarx1_319010.faa"

# Create a dictionary to store sequences
sequences = {}

# Parse the HMMER results
with open(hmmer_output_file, "r") as hmmer_output:
    for line in hmmer_output:
        if line.startswith("#"):
            continue
        fields = line.strip().split()
        hit_name = fields[3]
        start = int(fields[19])
        end = int(fields[20])
        # Load the query sequence
        query_sequences = SeqIO.to_dict(SeqIO.parse(query_sequence_file, "fasta"))
        if hit_name in query_sequences:
            hit_sequence = query_sequences[hit_name]
            gene_sequence = hit_sequence[start - 1 : end]  # Adjust for 0-based indexing
            sequences[hit_name] = gene_sequence
            print(sequences[hit_name])

# Now, 'sequences' contains the extracted gene sequences
