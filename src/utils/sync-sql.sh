#!/bin/bash

# Enhanced MySQL/MariaDB to SQLite Sync Script
# This script provides utilities for bidirectional database synchronization

set -e  # Exit on any error

# Default configuration - can be overridden by environment variables
DB_HOST="${DB_HOST:-localhost}"
DB_PORT="${DB_PORT:-3306}"
DB_USER="${DB_USER:-username}"
DB_PASS="${DB_PASS:-password}"
DB_NAME="${DB_NAME:-your_database_name}"
SQLITE_DB="${SQLITE_DB:-converted_database.sqlite}"
TABLES="${TABLES:-}"  # Empty means all tables

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Logging functions
log_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Check if required tools are installed
check_dependencies() {
    local missing_tools=()

    if ! command -v mysqldump &> /dev/null; then
        missing_tools+=("mysqldump")
    fi

    if ! command -v mysql &> /dev/null; then
        missing_tools+=("mysql")
    fi

    if ! command -v sqlite3 &> /dev/null; then
        missing_tools+=("sqlite3")
    fi

    if [ ${#missing_tools[@]} -ne 0 ]; then
        log_error "Missing required tools: ${missing_tools[*]}"
        log_error "Please install MySQL/MariaDB client and SQLite3"
        exit 1
    fi
}

# Test MySQL connection
test_mysql_connection() {
    log_info "Testing MySQL connection..."
    if ! mysql -h"$DB_HOST" -P"$DB_PORT" -u"$DB_USER" -p"$DB_PASS" -e "SELECT 1;" &> /dev/null; then
        log_error "Cannot connect to MySQL database"
        exit 1
    fi
    log_info "MySQL connection successful"
}

# Export MySQL schema
export_mysql_schema() {
    local output_file="$1"
    local tables_option=""

    if [ -n "$TABLES" ]; then
        tables_option="--tables $TABLES"
    fi

    log_info "Exporting MySQL schema..."
    mysqldump -h"$DB_HOST" -P"$DB_PORT" -u"$DB_USER" -p"$DB_PASS" \
        --no-data --databases "$DB_NAME" $tables_option > "$output_file"

    if [ $? -eq 0 ]; then
        log_info "Schema exported to $output_file"
    else
        log_error "Failed to export schema"
        exit 1
    fi
}

# Export MySQL data
export_mysql_data() {
    local output_file="$1"
    local tables_option=""

    if [ -n "$TABLES" ]; then
        tables_option="--tables $TABLES"
    fi

    log_info "Exporting MySQL data..."
    mysqldump -h"$DB_HOST" -P"$DB_PORT" -u"$DB_USER" -p"$DB_PASS" \
        --databases "$DB_NAME" $tables_option > "$output_file"

    if [ $? -eq 0 ]; then
        log_info "Data exported to $output_file"
    else
        log_error "Failed to export data"
        exit 1
    fi
}

# Convert MySQL schema to SQLite compatible format
convert_schema_for_sqlite() {
    local input_file="$1"
    local output_file="$2"

    log_info "Converting schema for SQLite compatibility..."

    # Create a temporary file for the converted schema
    cat "$input_file" | \
    # Remove MySQL-specific statements
    grep -v '^--' | \
    grep -v '^/\*' | \
    grep -v '^\*' | \
    grep -v '^LOCK TABLES' | \
    grep -v '^UNLOCK TABLES' | \
    grep -v '^DROP TABLE' | \
    grep -v '^SET ' | \
    grep -v '^USE ' | \
    # Convert AUTO_INCREMENT to AUTOINCREMENT
    sed 's/AUTO_INCREMENT/AUTOINCREMENT/g' | \
    # Convert ENGINE=InnoDB and similar
    sed 's/ENGINE=[^ ]* //g' | \
    sed 's/DEFAULT CHARSET=[^ ]* //g' | \
    sed 's/COLLATE=[^ ]* //g' | \
    # Convert MySQL types to SQLite types
    sed 's/TINYINT(1)/INTEGER/g' | \
    sed 's/INT([^)]*)/INTEGER/g' | \
    sed 's/BIGINT([^)]*)/INTEGER/g' | \
    sed 's/VARCHAR([^)]*)/TEXT/g' | \
    sed 's/CHAR([^)]*)/TEXT/g' | \
    sed 's/DATETIME/TEXT/g' | \
    sed 's/TIMESTAMP/TEXT/g' | \
    sed 's/DECIMAL([^)]*)/REAL/g' | \
    sed 's/DOUBLE/REAL/g' | \
    sed 's/FLOAT/REAL/g' | \
    # Remove backticks around table names in CREATE statements
    sed 's/CREATE TABLE `\([^`]*\)`/CREATE TABLE \1/g' | \
    # Clean up empty lines and trailing commas
    sed '/^$/d' | \
    sed 's/,$//' > "$output_file"

    log_info "Schema conversion completed: $output_file"
}

# Convert MySQL data dump to SQLite compatible format
convert_data_for_sqlite() {
    local input_file="$1"
    local output_file="$2"

    log_info "Converting data for SQLite compatibility..."

    # Extract INSERT statements and clean them up
    grep '^INSERT' "$input_file" | \
    # Remove backticks around table names
    sed 's/INSERT INTO `\([^`]*\)`/INSERT INTO \1/g' | \
    # Remove backticks around column names
    sed 's/`\([^`]*\)`/\1/g' > "$output_file"

    log_info "Data conversion completed: $output_file"
}

# Create SQLite database from converted schema
create_sqlite_db() {
    local schema_file="$1"
    local data_file="$2"
    local sqlite_db="$3"

    log_info "Creating SQLite database..."

    # Remove existing database if it exists
    if [ -f "$sqlite_db" ]; then
        log_warn "Removing existing SQLite database: $sqlite_db"
        rm "$sqlite_db"
    fi

    # Create database and import schema
    sqlite3 "$sqlite_db" < "$schema_file"

    if [ -f "$data_file" ] && [ -s "$data_file" ]; then
        log_info "Importing data into SQLite database..."
        sqlite3 "$sqlite_db" < "$data_file"
    else
        log_info "No data to import"
    fi

    log_info "SQLite database created: $sqlite_db"
}

# Main sync function
sync_mysql_to_sqlite() {
    log_info "Starting MySQL to SQLite synchronization..."

    check_dependencies
    test_mysql_connection

    local timestamp=$(date +"%Y%m%d_%H%M%S")
    local temp_dir="/tmp/mysql_sqlite_sync_$timestamp"

    mkdir -p "$temp_dir"

    local mysql_schema="$temp_dir/mysql_schema.sql"
    local mysql_data="$temp_dir/mysql_data.sql"
    local sqlite_schema="$temp_dir/sqlite_schema.sql"
    local sqlite_data="$temp_dir/sqlite_data.sql"

    # Export from MySQL
    export_mysql_schema "$mysql_schema"
    export_mysql_data "$mysql_data"

    # Convert for SQLite
    convert_schema_for_sqlite "$mysql_schema" "$sqlite_schema"
    convert_data_for_sqlite "$mysql_data" "$sqlite_data"

    # Create SQLite database
    create_sqlite_db "$sqlite_schema" "$sqlite_data" "$SQLITE_DB"

    # Cleanup
    rm -rf "$temp_dir"

    log_info "Synchronization completed successfully!"
    log_info "SQLite database created: $SQLITE_DB"
}

# Show usage information
show_usage() {
    cat << EOF
Enhanced MySQL/MariaDB to SQLite Sync Script

USAGE:
    $0 [OPTIONS]

ENVIRONMENT VARIABLES:
    DB_HOST        MySQL host (default: localhost)
    DB_PORT        MySQL port (default: 3306)
    DB_USER        MySQL username (default: username)
    DB_PASS        MySQL password (default: password)
    DB_NAME        MySQL database name (default: your_database_name)
    SQLITE_DB      Output SQLite database file (default: converted_database.sqlite)
    TABLES         Specific tables to sync (default: all tables)

EXAMPLES:
    # Using environment variables
    export DB_USER=myuser
    export DB_PASS=mypass
    export DB_NAME=mydb
    $0

    # Sync specific tables only
    TABLES="users products" DB_USER=myuser DB_PASS=mypass DB_NAME=mydb $0

    # Custom output file
    SQLITE_DB=my_database.sqlite DB_USER=myuser DB_PASS=mypass DB_NAME=mydb $0

EOF
}

# Main script logic
case "${1:-sync}" in
    "sync"|"-s"|"--sync")
        sync_mysql_to_sqlite
        ;;
    "help"|"-h"|"--help")
        show_usage
        ;;
    *)
        log_error "Unknown option: $1"
        show_usage
        exit 1
        ;;
esac