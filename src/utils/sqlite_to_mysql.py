#!/usr/bin/env python3
"""
SQLite to MySQL Converter

This script converts a SQLite database to MySQL format for multi-user access.
It handles schema conversion, data migration, and relationship preservation.

Usage:
    python sqlite_to_mysql.py --sqlite-db path/to/sqlite.db --mysql-config config.json

Requirements:
    - pymysql
    - sqlalchemy
    - mysql-connector-python or similar MySQL driver
"""

import argparse
import json
import logging
import os
import sys
from pathlib import Path
from typing import Dict, Any, List, Optional

import sqlalchemy as sa
from sqlalchemy import create_engine, MetaData, Table, Column, ForeignKey
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, relationship
from sqlalchemy.sql import text

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class DatabaseConverter:
    """Handles conversion between SQLite and MySQL databases."""

    def __init__(self, sqlite_path: str, mysql_config: Dict[str, Any]):
        self.sqlite_path = sqlite_path
        self.mysql_config = mysql_config

        # SQLite connection
        self.sqlite_engine = create_engine(f'sqlite:///{sqlite_path}')

        # MySQL connection
        mysql_url = self._build_mysql_url()
        self.mysql_engine = create_engine(mysql_url, echo=False)

        # Metadata for schema operations
        self.sqlite_metadata = MetaData()
        self.mysql_metadata = MetaData()

    def _build_mysql_url(self) -> str:
        """Build MySQL connection URL from config."""
        config = self.mysql_config
        return (f"mysql+pymysql://{config['user']}:{config['password']}"
                f"@{config['host']}:{config.get('port', 3306)}/{config['database']}")

    def _get_sqlite_tables(self) -> List[str]:
        """Get all table names from SQLite database."""
        with self.sqlite_engine.connect() as conn:
            result = conn.execute(text("SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'"))
            return [row[0] for row in result.fetchall()]

    def _convert_column_type(self, sqlite_type: str, column_name: str, is_primary_key: bool = False) -> str:
        """Convert SQLite column types to MySQL equivalents."""
        sqlite_type = sqlite_type.upper()
        column_name_lower = column_name.lower()

        # Handle common SQLite to MySQL type conversions
        type_mapping = {
            'INTEGER': 'INT',
            'BIGINT': 'BIGINT',
            'SMALLINT': 'SMALLINT',
            'TINYINT': 'TINYINT',
            'VARCHAR': 'VARCHAR(255)',
            'TEXT': 'TEXT',
            'BLOB': 'BLOB',
            'REAL': 'DOUBLE',
            'FLOAT': 'DOUBLE',
            'DOUBLE': 'DOUBLE',
            'DECIMAL': 'DECIMAL',
            'BOOLEAN': 'TINYINT(1)',
            'DATE': 'DATE',
            'DATETIME': 'DATETIME',
            'TIMESTAMP': 'TIMESTAMP',
            'TIME': 'TIME',
        }

        # Handle VARCHAR with specific lengths
        if sqlite_type.startswith('VARCHAR('):
            return sqlite_type  # Keep as is

        # Handle auto-increment: only for primary keys that are integers
        if 'AUTOINCREMENT' in sqlite_type and is_primary_key:
            return 'INT AUTO_INCREMENT'

        # Special handling for sequence-related columns that need LONGTEXT
        sequence_columns = ['sequence', 'seq', 'dna', 'rna', 'nucleotide', 'protein']
        if sqlite_type == 'TEXT' and any(seq_term in column_name_lower for seq_term in sequence_columns):
            return 'LONGTEXT'

        # Special handling for abstract and other long text fields
        long_text_columns = ['abstract', 'description', 'content', 'body', 'comment', 'notes',
                           'details', 'source', 'evidence', 'curated_status', 'starshipID']
        if sqlite_type == 'TEXT' and any(long_term in column_name_lower for long_term in long_text_columns):
            return 'LONGTEXT'

        # Special handling for TIMESTAMP vs DATETIME
        if sqlite_type == 'TIMESTAMP':
            return 'DATETIME'

        return type_mapping.get(sqlite_type, 'TEXT')

    def _get_table_schema(self, table_name: str) -> Dict[str, Any]:
        """Get table schema information from SQLite."""
        with self.sqlite_engine.connect() as conn:
            # Get column information
            columns_result = conn.execute(text(f"PRAGMA table_info({table_name})"))
            columns = columns_result.fetchall()

            # Get foreign key information
            fk_result = conn.execute(text(f"PRAGMA foreign_key_list({table_name})"))
            foreign_keys = fk_result.fetchall()

            # Get indexes
            index_result = conn.execute(text(f"PRAGMA index_list({table_name})"))
            indexes = index_result.fetchall()

            return {
                'columns': columns,
                'foreign_keys': foreign_keys,
                'indexes': indexes
            }

    def _create_mysql_table(self, table_name: str, schema_info: Dict[str, Any]) -> None:
        """Create MySQL table with proper schema."""
        columns = schema_info['columns']
        foreign_keys = schema_info['foreign_keys']

        # Build column definitions
        column_defs = []

        for col in columns:
            col_name = col[1]
            col_type = col[2]
            not_null = col[3]
            default_value = col[4]
            is_pk = col[5]

            # Convert type
            mysql_type = self._convert_column_type(col_type, col_name, is_pk)

            # Build column definition
            col_def = f"`{col_name}` {mysql_type}"

            if not_null:
                col_def += " NOT NULL"

            if default_value is not None:
                if isinstance(default_value, str):
                    col_def += f" DEFAULT '{default_value}'"
                else:
                    col_def += f" DEFAULT {default_value}"

            if is_pk:
                col_def += " PRIMARY KEY"

            column_defs.append(col_def)

        # Add foreign keys
        for fk in foreign_keys:
            fk_col = fk[3]
            ref_table = fk[2]
            ref_col = fk[4]
            column_defs.append(f"FOREIGN KEY (`{fk_col}`) REFERENCES `{ref_table}`(`{ref_col}`)")

        # Special handling for cleanup_issues tables with problematic defaults
        if table_name in ['cleanup_issues_old', 'cleanup_issues']:
            # Fix the default value syntax for status column
            column_defs = [col.replace("'OPEN''", "'OPEN'") for col in column_defs]
            column_defs = [col.replace("''OPEN''", "'OPEN'") for col in column_defs]

        # Create table
        create_sql = f"CREATE TABLE IF NOT EXISTS `{table_name}` ({', '.join(column_defs)})"

        with self.mysql_engine.connect() as conn:
            try:
                conn.execute(text(create_sql))
                # Note: commit is automatic in SQLAlchemy 2.0 context manager
                logger.info(f"Created table: {table_name}")
            except Exception as e:
                if "foreign key constraint" in str(e).lower():
                    # Try creating table without foreign keys first
                    logger.warning(f"Foreign key constraint failed for {table_name}, creating without constraints: {e}")
                    # Remove foreign key definitions and try again
                    simple_defs = [col for col in column_defs if "FOREIGN KEY" not in col]
                    simple_sql = f"CREATE TABLE IF NOT EXISTS `{table_name}` ({', '.join(simple_defs)})"
                    conn.execute(text(simple_sql))
                    logger.info(f"Created table {table_name} without foreign key constraints")
                else:
                    raise

    def _migrate_table_data(self, table_name: str) -> None:
        """Migrate data from SQLite table to MySQL table."""
        # Get row count
        with self.sqlite_engine.connect() as conn:
            count_result = conn.execute(text(f"SELECT COUNT(*) FROM `{table_name}`"))
            row_count = count_result.fetchone()[0]

        if row_count == 0:
            logger.info(f"No data to migrate for table: {table_name}")
            return

        logger.info(f"Migrating {row_count} rows from table: {table_name}")

        # Get all data from SQLite
        with self.sqlite_engine.connect() as conn:
            result = conn.execute(text(f"SELECT * FROM `{table_name}`"))
            rows = result.fetchall()
            column_names = result.keys()

        if not rows:
            return

        # Insert data into MySQL
        with self.mysql_engine.connect() as conn:
            # Prepare insert statement with named parameters
            param_names = [f'p{i}' for i in range(len(column_names))]
            placeholders = ', '.join([f':{p}' for p in param_names])
            columns_str = ', '.join([f'`{col}`' for col in column_names])
            insert_sql = f"INSERT INTO `{table_name}` ({columns_str}) VALUES ({placeholders})"

            # Insert in batches to handle large tables
            batch_size = 1000
            for i in range(0, len(rows), batch_size):
                batch = rows[i:i + batch_size]

                # Insert each row individually to handle data type issues
                for row in batch:
                    try:
                        # Create parameter dictionary
                        params = {f'p{i}': val for i, val in enumerate(row)}

                        # Execute the insert with named parameters
                        conn.execute(text(insert_sql), params)

                    except Exception as row_error:
                        logger.warning(f"Error inserting row: {row_error}, row data: {row}")
                        logger.warning(f"Row types: {[type(val).__name__ for val in row]}")
                        # Try alternative approach with explicit type conversion
                        try:
                            # Convert all values to strings first, let MySQL handle conversion
                            string_params = {f'p{i}': str(val) if val is not None else None for i, val in enumerate(row)}
                            conn.execute(text(insert_sql), string_params)
                        except Exception as fallback_error:
                            logger.error(f"Fallback insert also failed: {fallback_error}")
                            # Continue with other rows
                            continue

                # Note: commit is automatic in SQLAlchemy 2.0 context manager
                logger.info(f"Migrated batch {i//batch_size + 1} for table: {table_name}")

    def convert(self) -> None:
        """Main conversion method."""
        logger.info("Starting SQLite to MySQL conversion...")

        # Create MySQL database if it doesn't exist
        with self.mysql_engine.connect() as conn:
            conn.execute(text(f"CREATE DATABASE IF NOT EXISTS `{self.mysql_config['database']}`"))
            conn.execute(text(f"USE `{self.mysql_config['database']}`"))
            # Note: commit is automatic in SQLAlchemy 2.0 context manager

        # Get all tables
        tables = self._get_sqlite_tables()
        logger.info(f"Found {len(tables)} tables to convert: {tables}")

        # Define table creation order to handle foreign key dependencies
        # Tables without foreign keys first, then tables with foreign keys
        priority_tables = [
            'taxonomy', 'papers', 'family_names', 'genome', 'accessions', 'ships',
            'captains', 'starship_features', 'gff', 'joined_ships', 'paper_family_association',
            'haplotype_names', 'navis_names', 'sequence_reference', 'NAR_S13',
            'cleanup_issues_old', 'cleanup_issues', 'alembic_version'
        ]

        # Separate tables into priority and remaining
        ordered_tables = []
        remaining_tables = []

        for table in tables:
            if table in priority_tables:
                ordered_tables.append(table)
            else:
                remaining_tables.append(table)

        # Sort priority tables by dependency order
        ordered_tables.sort(key=lambda x: priority_tables.index(x))
        final_table_order = ordered_tables + remaining_tables

        logger.info(f"Processing tables in order: {final_table_order}")

        # Convert each table in dependency order
        for table_name in final_table_order:
            try:
                logger.info(f"Converting table: {table_name}")

                # Get schema info
                schema_info = self._get_table_schema(table_name)

                # Create MySQL table
                self._create_mysql_table(table_name, schema_info)

                # Migrate data
                self._migrate_table_data(table_name)

                logger.info(f"Successfully converted table: {table_name}")

            except Exception as e:
                logger.error(f"Error converting table {table_name}: {str(e)}")
                # Continue with other tables instead of failing completely
                continue

        logger.info("SQLite to MySQL conversion completed successfully!")

    def verify_conversion(self) -> None:
        """Verify that the conversion was successful by comparing row counts."""
        logger.info("Verifying conversion...")

        tables = self._get_sqlite_tables()

        for table_name in tables:
            # Get SQLite row count
            with self.sqlite_engine.connect() as conn:
                sqlite_count = conn.execute(text(f"SELECT COUNT(*) FROM `{table_name}`")).fetchone()[0]

            # Get MySQL row count
            with self.mysql_engine.connect() as conn:
                mysql_count = conn.execute(text(f"SELECT COUNT(*) FROM `{table_name}`")).fetchone()[0]

            if sqlite_count != mysql_count:
                logger.warning(f"Row count mismatch for {table_name}: SQLite={sqlite_count}, MySQL={mysql_count}")
            else:
                logger.info(f"Row count verified for {table_name}: {sqlite_count} rows")


def load_config(config_path: str) -> Dict[str, Any]:
    """Load MySQL configuration from JSON file."""
    with open(config_path, 'r') as f:
        return json.load(f)


def main():
    parser = argparse.ArgumentParser(description='Convert SQLite database to MySQL')
    parser.add_argument('--sqlite-db', required=True, help='Path to SQLite database file')
    parser.add_argument('--mysql-config', required=True, help='Path to MySQL configuration JSON file')
    parser.add_argument('--verify', action='store_true', help='Verify conversion after completion')

    args = parser.parse_args()

    # Check if SQLite file exists
    if not os.path.exists(args.sqlite_db):
        logger.error(f"SQLite database file not found: {args.sqlite_db}")
        sys.exit(1)

    # Check if config file exists
    if not os.path.exists(args.mysql_config):
        logger.error(f"MySQL config file not found: {args.mysql_config}")
        sys.exit(1)

    # Load configuration
    try:
        mysql_config = load_config(args.mysql_config)
    except Exception as e:
        logger.error(f"Error loading MySQL config: {str(e)}")
        sys.exit(1)

    # Perform conversion
    try:
        converter = DatabaseConverter(args.sqlite_db, mysql_config)
        converter.convert()

        if args.verify:
            converter.verify_conversion()

        logger.info("Conversion completed successfully!")

    except Exception as e:
        logger.error(f"Conversion failed: {str(e)}")
        sys.exit(1)


if __name__ == '__main__':
    main()
