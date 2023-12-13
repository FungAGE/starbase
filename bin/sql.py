import os
import sqlite3

# Connect to SQLite database (or create a new one if it doesn't exist)
conn = sqlite3.connect("SQL/starbase.sqlite")
cursor = conn.cursor()

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
    "SQL/data/gff/mycodb",
    "SQL/data/hmm/cargo",
    "SQL/data/hmm/captain",
    "SQL/data/faa/cargo/nlr/mycodb",
    "SQL/data/faa/cargo/fre/mycodb",
    "SQL/data/faa/cargo/plp/mycodb",
    "SQL/data/faa/cargo/duf3723/mycodb",
    "SQL/data/faa/captain/tyr/mycodb",
    "SQL/data/faa/captain/tyr/manual",
)

for directory in directories:
    data_type = directory.split("/")[2]

    if "manual" in directory:
        db_source = "manual-annotations"
    elif "mycodb" or "mtdb" in directory:
        db_source = "mycodb"
    else:
        db_source = ""

    if "cargo" in directory:
        type = "cargo"
    elif "captain" in directory:
        type = "captain"
    elif "ship" in directory:
        type = "ship"
    else:
        type = ""

    if (type == "captain" or "cargo") and len(directory.split("/")) > 4:
        gene_type = directory.split("/")[4]

    # Create table if it doesn't exist
    cursor.execute(
        f"""
        CREATE TABLE IF NOT EXISTS {data_type} (
            id INTEGER PRIMARY KEY,
            name TEXT,
            data_type TEXT,
            db_source TEXT,
            type TEXT,
            gene_type TEXT,
            file_path TEXT
        )
    """
    )

    # Iterate through files in the directory
    for file_name in os.listdir(directory):
        file_path = os.path.join(directory, file_name)
        name = os.path.splitext(file_name)[0]

        # Insert record into the fasta table
        cursor.execute(
            f"""
            INSERT INTO {data_type} (name, data_type, db_source, type, gene_type, file_path)
            VALUES (?, ?, ?, ?, ?, ?)
            """,
            (name, data_type, db_source, type, gene_type, file_path),
        )

# Commit changes and close the connection
conn.commit()
conn.close()
