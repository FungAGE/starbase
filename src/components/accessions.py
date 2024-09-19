import sqlite3

# Connect to the SQLite database
conn = sqlite3.connect("database_folder/starbase.sqlite")
cursor = conn.cursor()

# Fetch all rows from the table (assuming 'accession' is the primary key or unique identifier)
cursor.execute("SELECT accession FROM ships")
accessions = cursor.fetchall()

# Generate and update accession tags
for acc in accessions:
    # Extract the single value from the tuple
    accession = acc[0]
    
    # Generate the accession tag with zero-padded numbers
    accession_tag = f"SBS{str(accession).zfill(6)}"
    
    # Print for debugging purposes
    print(f"Updating accession {accession} with tag {accession_tag}")
    
    # Update the table with the generated accession tag
    cursor.execute("UPDATE ships SET accession_tag = ? WHERE accession = ?", (accession_tag, accession))

# Commit the transaction to save changes
conn.commit()

# Close the connection
conn.close()