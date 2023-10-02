from Bio import AlignIO
import csv

# Open the Stockholm file for reading
input_file = "tmp/alignment.sto"

# Create a CSV file for writing
output_file = "tmp/alignment.csv"
csv_columns = ["ID", "Sequence"]

# Create a list to store data
data_list = []

# Parse the Stockholm file and extract sequences
with open(input_file, "r") as stockholm_file:
    for alignment in AlignIO.parse(stockholm_file, "stockholm"):
        for record in alignment:
            data_list.append({"ID": record.id, "Sequence": str(record.seq)})

# Write the extracted data to a CSV file
with open(output_file, "w", newline="") as csv_file:
    writer = csv.DictWriter(csv_file, fieldnames=csv_columns)
    writer.writeheader()
    for data in data_list:
        writer.writerow(data)
