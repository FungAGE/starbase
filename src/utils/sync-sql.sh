#!/bin/bash

# MariaDB credentials
DB_USER="username"
DB_PASS="password"
DB_NAME="your_database_name"
TABLES="table1 table2"

# Export the schema and data
mysqldump -u $DB_USER -p$DB_PASS --no-data --databases $DB_NAME --tables $TABLES > schema.sql
mysqldump -u $DB_USER -p$DB_PASS --databases $DB_NAME --tables $TABLES > data.sql

# Convert schema.sql to SQLite compatible schema
# (you could add sed or awk commands here to automate conversions)

# Create the SQLite database and import the modified schema and data
sqlite3 your_database.sqlite < modified_schema.sql
sqlite3 your_database.sqlite < modified_data.sql

# Run the script to sync your MariaDB database with the SQLite subset.

# 6. Docker Integration

# If you want to automate this within a Docker environment, you can create a Dockerfile and entrypoint script that runs the above commands when the container is built or started.
# Here’s a simple example of a Dockerfile:

# FROM python:3.9-slim

# # Install MariaDB client
# RUN apt-get update && apt-get install -y mariadb-client sqlite3

# # Copy the script into the container
# COPY sync_databases.sh /sync_databases.sh

# # Make the script executable
# RUN chmod +x /sync_databases.sh

# # Set the entrypoint
# ENTRYPOINT ["/sync_databases.sh"]

# # Make sure to pass the necessary environment variables (e.g., credentials) when you run the container:

# docker run --rm -e DB_USER=user -e DB_PASS=pass my_database_sync_image