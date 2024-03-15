import os
import sqlite3

# Connect to SQLite database (or create a new one if it doesn't exist)
conn = sqlite3.connect("SQL/starbase.sqlite")
cursor = conn.cursor()

# Create table if it doesn't exist
cursor.execute(
    """
    CREATE TABLE IF NOT EXISTS checksums (
        id INTEGER PRIMARY KEY,
        checksum TEXT,
        checksum_file_path TEXT,
        data_type TEXT,
        type TEXT,
        gene_type TEXT
    )
    """
)

# list of directories to search through
directories = (
    "SQL/data/tree/cargo",
    "SQL/data/tree/captain",
    "SQL/data/fna/cargo/nlr",
    "SQL/data/fna/cargo/fre",
    "SQL/data/fna/cargo/plp",
    "SQL/data/fna/cargo/duf3723",
    "SQL/data/fna/captain/tyr/mycodb",
    "SQL/data/fna/ships/manual-annotations",
    "SQL/data/fna/ships/mycodb",
    "metadata/ships/starfish/gff/mycodb",
    "SQL/data/hmm/cargo",
    "SQL/data/hmm/captain",
    "SQL/data/faa/cargo/nlr/mycodb",
    "SQL/data/faa/cargo/fre/mycodb",
    "SQL/data/faa/cargo/plp/mycodb",
    "SQL/data/faa/cargo/duf3723/mycodb",
    "SQL/data/faa/captain/tyr/mycodb",
    "SQL/data/faa/captain/tyr/manual",
)

# for directory in directories:
#     data_type = directory.split("/")[2]

#     if "manual" in directory:
#         db_source = "manual-annotations"
#     elif "mycodb" or "mtdb" in directory:
#         db_source = "mycodb"
#     else:
#         db_source = ""

#     if "cargo" in directory:
#         type = "cargo"
#     elif "captain" in directory:
#         type = "captain"
#     elif "ship" in directory:
#         type = "ship"
#     else:
#         type = ""

#     if (type == "captain" or "cargo") and len(directory.split("/")) > 4:
#         gene_type = directory.split("/")[4]

#     # Create table if it doesn't exist
#     cursor.execute(
#         f"""
#         CREATE TABLE IF NOT EXISTS {data_type} (
#             id INTEGER PRIMARY KEY,
#             name TEXT,
#             data_type TEXT,
#             db_source TEXT,
#             type TEXT,
#             gene_type TEXT,
#             file_path TEXT
#         )
#     """
#     )

#     # Iterate through files in the directory
#     for file_name in os.listdir(directory):
#         name = os.path.splitext(file_name)[0]
#         file_path = os.path.join(directory, file_name)
#         # Insert record into the fasta table
#         cursor.execute(
#             f"""
#             INSERT INTO {data_type} (name, data_type, db_source, type, gene_type, file_path)
#             VALUES (?, ?, ?, ?, ?, ?)
#             """,
#             (name, data_type, db_source, type, gene_type, file_path),
#         )

checksum_files = (
    "SQL/data/fna/cargo/fna.checksums.txt",
    "SQL/data/fna/ships/fna.checksums.txt",
    "SQL/data/faa/cargo/nlr/faa.checksums.txt",
    "SQL/data/faa/cargo/fre/faa.checksums.txt",
    "SQL/data/faa/cargo/plp/faa.checksums.txt",
    "SQL/data/faa/cargo/duf3723/faa.checksums.txt",
    "SQL/data/faa/captain/tyr/faa.checksums.txt",
)
for checksum_file in checksum_files:
    if checksum_file.endswith(".txt"):
        print(checksum_file)
        data_type = checksum_file.split("/")[2]
        type = checksum_file.split("/")[3]
        if type == "ships":
            gene_type = ""
        else:
            gene_type = checksum_file.split("/")[4]
        with open(checksum_file, "r") as file:
            lines = file.readlines()
        # Insert data into the table
        for line in lines:
            parts = line.strip().split("\t")
            checksum, checksum_file_path = parts[0], parts[1]

            insert_query = """
                INSERT INTO checksums (checksum, checksum_file_path, data_type, type, gene_type)
                VALUES (?, ?, ?, ?, ?);
            """
            cursor.execute(
                insert_query,
                (
                    checksum,
                    checksum_file_path,
                    data_type,
                    type,
                    gene_type,
                ),
            )


# Commit changes and close the connection
conn.commit()
conn.close()
