import argparse
import os
import pandas as pd

# Create a command-line argument parser
parser = argparse.ArgumentParser(description="Process TSV and BED files")

# Add an argument for the bed_file
parser.add_argument("bed_file", help="Path to the BED file")

# Parse the command-line arguments
args = parser.parse_args()

# make coordinates in gff relative to the element, rather than the chr. Use columns 4 and 5 in `mycodb.final.starships.feat` to correct. and remember, gff is 1-indexed
# Load the TSV file into a DataFrame
file_path = "Starships/ships/mycodb/output/mycodb.final.starships.feat"
df = pd.read_csv(file_path, sep="\t")

# Initialize a dictionary to store data on element lengths
begin_dict = {}
end_dict = {}
for index, row in df.iterrows():
    id_value = row["starshipID"]
    # Calculate the length of elements
    start = row["elementBegin"]
    end = row["elementEnd"]
    begin_dict[id_value] = start
    end_dict[id_value] = end

# Initialize a dictionary to store data for each ID
data_by_id = {}

with open(args.bed_file, "r") as bed_file:
    for line in bed_file:
        fields = line.strip().split("\t")
        contig = fields[0]
        name = fields[3]
        feature = fields[4]
        strand = fields[5]
        element = fields[6]  # 7th column with comma-separated items

        if feature == "flank":
            type = "region"
        else:
            type = "gene"
        score = "."

        # Split the 7th column by commas
        elements = element.split(",")

        for item in elements:
            start = abs(int(fields[1]) - int(begin_dict[item])) + 1
            end = abs(int(fields[2]) - int(end_dict[item]))

            output = f"{contig}\tbed\t{type}\t{start}\t{end}\t{score}\t{strand}\t{feature}\tID={name};Description={item};parent={contig}"

            name_key = item

            # Check if the ID already exists in the dictionary
            if name_key in data_by_id:
                data_by_id[name_key].append(output)
            else:
                data_by_id[name_key] = [output]

# Create a directory to store the output files
# output_directory = args.bed_file.strip(".bed") + "_split"
output_directory = "SQL/data/gff/mycodb"
os.makedirs(output_directory, exist_ok=True)

# Create separate GFF files for each ID
for name_key, data in data_by_id.items():
    output_file_name = f"{output_directory}/{name_key}.gff"
    with open(output_file_name, "w") as output_file:
        for line in data:
            output_file.write(line + "\n")
