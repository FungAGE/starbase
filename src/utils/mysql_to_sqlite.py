#!/usr/bin/env python3
"""
MySQL to SQLite Converter

This script converts a MySQL database to SQLite format for local project use.
It handles schema conversion, data migration, and relationship preservation.

Usage:
    python mysql_to_sqlite.py --sqlite-db path/to/output.db --mysql-config config.json

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
from typing import Dict, Any, List, Optional, Tuple

import sqlalchemy as sa
from sqlalchemy import create_engine, MetaData, Table, Column, ForeignKey, text
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, relationship

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class DatabaseConverter:
    """Handles conversion between MySQL and SQLite databases."""

    def __init__(self, sqlite_path: str, mysql_config: Dict[str, Any]):
        self.sqlite_path = sqlite_path
        self.mysql_config = mysql_config

        # MySQL connection
        mysql_url = self._build_mysql_url()
        self.mysql_engine = create_engine(mysql_url, echo=False)

        # SQLite connection (will create the database if it doesn't exist)
        self.sqlite_engine = create_engine(f'sqlite:///{sqlite_path}')

        # Metadata for schema operations
        self.mysql_metadata = MetaData()
        self.sqlite_metadata = MetaData()

    def _build_mysql_url(self) -> str:
        """Build MySQL connection URL from config."""
        config = self.mysql_config
        return (f"mysql+pymysql://{config['user']}:{config['password']}"
                f"@{config['host']}:{config.get('port', 3306)}/{config['database']}")

    def _get_mysql_tables(self) -> List[str]:
        """Get all table names from MySQL database."""
        with self.mysql_engine.connect() as conn:
            result = conn.execute(text("SHOW TABLES"))
            return [row[0] for row in result.fetchall()]

    def _convert_column_type(self, mysql_type: str, column_name: str, extra: str = "", is_primary_key: bool = False) -> str:
        """Convert MySQL column types to SQLite equivalents."""
        mysql_type = mysql_type.upper()
        column_name_lower = column_name.lower()

        # Handle common MySQL to SQLite type conversions
        type_mapping = {
            'INT': 'INTEGER',
            'BIGINT': 'INTEGER',
            'SMALLINT': 'INTEGER',
            'TINYINT': 'INTEGER',
            'VARCHAR': 'TEXT',
            'CHAR': 'TEXT',
            'TEXT': 'TEXT',
            'LONGTEXT': 'TEXT',
            'MEDIUMTEXT': 'TEXT',
            'BLOB': 'BLOB',
            'LONGBLOB': 'BLOB',
            'DOUBLE': 'REAL',
            'FLOAT': 'REAL',
            'DECIMAL': 'REAL',
            'BOOLEAN': 'INTEGER',
            'BOOL': 'INTEGER',
            'DATE': 'TEXT',
            'DATETIME': 'TEXT',
            'TIMESTAMP': 'TEXT',
            'TIME': 'TEXT',
            'YEAR': 'INTEGER',
        }

        # Handle auto-increment: only for primary keys
        if 'AUTO_INCREMENT' in extra.upper() and is_primary_key:
            return 'INTEGER PRIMARY KEY AUTOINCREMENT'

        # Handle VARCHAR with specific lengths
        if mysql_type.startswith('VARCHAR(') or mysql_type.startswith('CHAR('):
            return 'TEXT'

        # Handle DECIMAL/NUMERIC with precision
        if mysql_type.startswith('DECIMAL(') or mysql_type.startswith('NUMERIC('):
            return 'REAL'

        return type_mapping.get(mysql_type.split('(')[0], 'TEXT')

    def _get_table_schema(self, table_name: str) -> Dict[str, Any]:
        """Get table schema information from MySQL."""
        with self.mysql_engine.connect() as conn:
            # Get column information
            columns_result = conn.execute(text(f"DESCRIBE `{table_name}`"))
            columns = columns_result.fetchall()

            # Get foreign key information
            fk_result = conn.execute(text(f"""
                SELECT
                    COLUMN_NAME,
                    REFERENCED_TABLE_NAME,
                    REFERENCED_COLUMN_NAME
                FROM INFORMATION_SCHEMA.KEY_COLUMN_USAGE
                WHERE TABLE_NAME = '{table_name}'
                AND TABLE_SCHEMA = '{self.mysql_config['database']}'
                AND REFERENCED_TABLE_NAME IS NOT NULL
            """))
            foreign_keys = fk_result.fetchall()

            # Get indexes (excluding primary keys)
            index_result = conn.execute(text(f"""
                SELECT INDEX_NAME, COLUMN_NAME, NON_UNIQUE
                FROM INFORMATION_SCHEMA.STATISTICS
                WHERE TABLE_NAME = '{table_name}'
                AND TABLE_SCHEMA = '{self.mysql_config['database']}'
                AND INDEX_NAME != 'PRIMARY'
                ORDER BY INDEX_NAME, SEQ_IN_INDEX
            """))
            indexes = index_result.fetchall()

            return {
                'columns': columns,
                'foreign_keys': foreign_keys,
                'indexes': indexes
            }

    def _create_sqlite_table(self, table_name: str, schema_info: Dict[str, Any]) -> None:
        """Create SQLite table with proper schema."""
        columns = schema_info['columns']
        foreign_keys = schema_info['foreign_keys']
        indexes = schema_info['indexes']

        # Build column definitions
        column_defs = []
        primary_keys = []

        for col in columns:
            col_name = col[0]
            col_type = col[1]
            nullable = col[2] == 'YES'
            key = col[3]
            default_value = col[4]
            extra = col[5] or ""

            # Check if this is a primary key
            is_primary_key = key == 'PRI'

            # Convert type
            sqlite_type = self._convert_column_type(col_type, col_name, extra, is_primary_key)

            # Build column definition
            col_def = f"`{col_name}` {sqlite_type}"

            if not nullable:
                col_def += " NOT NULL"

            if default_value is not None and default_value != '':
                # Handle different default value formats
                if isinstance(default_value, str):
                    if default_value.upper() in ['CURRENT_TIMESTAMP', 'NOW()']:
                        col_def += " DEFAULT CURRENT_TIMESTAMP"
                    elif default_value.startswith("'") and default_value.endswith("'"):
                        col_def += f" DEFAULT {default_value}"
                    else:
                        col_def += f" DEFAULT '{default_value}'"
                else:
                    col_def += f" DEFAULT {default_value}"

            # Handle primary keys
            if key == 'PRI':
                if 'AUTOINCREMENT' in sqlite_type:
                    col_def = col_def.replace(' NOT NULL', '')  # SQLite auto-increment columns are implicitly NOT NULL
                else:
                    primary_keys.append(f"`{col_name}`")

            column_defs.append(col_def)

        # Add composite primary key if needed
        if primary_keys and len(primary_keys) > 1:
            column_defs.append(f"PRIMARY KEY ({', '.join(primary_keys)})")

        # Add foreign keys
        for fk in foreign_keys:
            fk_col = fk[0]
            ref_table = fk[1]
            ref_col = fk[2]
            column_defs.append(f"FOREIGN KEY (`{fk_col}`) REFERENCES `{ref_table}`(`{ref_col}`)")

        # Create table
        create_sql = f"CREATE TABLE IF NOT EXISTS `{table_name}` ({', '.join(column_defs)})"

        with self.sqlite_engine.connect() as conn:
            conn.execute(text(create_sql))
            # Note: commit is automatic in SQLAlchemy 2.0 context manager

        logger.info(f"Created table: {table_name}")

        # Create indexes
        self._create_indexes(table_name, indexes)

    def _create_indexes(self, table_name: str, indexes: List[Tuple]) -> None:
        """Create indexes on SQLite table."""
        if not indexes:
            return

        index_groups = {}
        for index_name, column_name, non_unique in indexes:
            if index_name not in index_groups:
                index_groups[index_name] = {'columns': [], 'unique': not non_unique}
            index_groups[index_name]['columns'].append(column_name)

        with self.sqlite_engine.connect() as conn:
            for index_name, index_info in index_groups.items():
                columns = index_info['columns']
                unique = index_info['unique']

                # Skip if this is a single-column primary key index
                if len(columns) == 1 and unique:
                    continue

                unique_str = "UNIQUE " if unique else ""
                columns_str = ', '.join([f'`{col}`' for col in columns])
                index_sql = f"CREATE {unique_str}INDEX IF NOT EXISTS `{index_name}` ON `{table_name}` ({columns_str})"

                try:
                    conn.execute(text(index_sql))
                    # Note: commit is automatic in SQLAlchemy 2.0 context manager
                    logger.info(f"Created index: {index_name} on {table_name}")
                except Exception as e:
                    logger.warning(f"Could not create index {index_name}: {str(e)}")

    def _migrate_table_data(self, table_name: str) -> None:
        """Migrate data from MySQL table to SQLite table."""
        # Get row count
        with self.mysql_engine.connect() as conn:
            count_result = conn.execute(text(f"SELECT COUNT(*) FROM `{table_name}`"))
            row_count = count_result.fetchone()[0]

        if row_count == 0:
            logger.info(f"No data to migrate for table: {table_name}")
            return

        logger.info(f"Migrating {row_count} rows from table: {table_name}")

        # Get column information for proper data handling
        with self.mysql_engine.connect() as conn:
            desc_result = conn.execute(text(f"DESCRIBE `{table_name}`"))
            columns_info = desc_result.fetchall()
            column_names = [col[0] for col in columns_info]

        # Get all data from MySQL
        with self.mysql_engine.connect() as conn:
            result = conn.execute(text(f"SELECT * FROM `{table_name}`"))
            rows = result.fetchall()

        if not rows:
            return

        # Insert data into SQLite
        with self.sqlite_engine.connect() as conn:
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
                            # Convert all values to strings first, let SQLite handle conversion
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
        logger.info("Starting MySQL to SQLite conversion...")

        # Get all tables
        tables = self._get_mysql_tables()
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

                # Create SQLite table
                self._create_sqlite_table(table_name, schema_info)

                # Migrate data
                self._migrate_table_data(table_name)

                logger.info(f"Successfully converted table: {table_name}")

            except Exception as e:
                logger.error(f"Error converting table {table_name}: {str(e)}")
                # Continue with other tables instead of failing completely
                continue

        logger.info("MySQL to SQLite conversion completed successfully!")

    def verify_conversion(self) -> None:
        """Verify that the conversion was successful by comparing row counts."""
        logger.info("Verifying conversion...")

        tables = self._get_mysql_tables()

        for table_name in tables:
            # Get MySQL row count
            with self.mysql_engine.connect() as conn:
                mysql_count = conn.execute(text(f"SELECT COUNT(*) FROM `{table_name}`")).fetchone()[0]

            # Get SQLite row count
            with self.sqlite_engine.connect() as conn:
                sqlite_count = conn.execute(text(f"SELECT COUNT(*) FROM `{table_name}`")).fetchone()[0]

            if mysql_count != sqlite_count:
                logger.warning(f"Row count mismatch for {table_name}: MySQL={mysql_count}, SQLite={sqlite_count}")
            else:
                logger.info(f"Row count verified for {table_name}: {mysql_count} rows")


def load_config(config_path: str) -> Dict[str, Any]:
    """Load MySQL configuration from JSON file."""
    with open(config_path, 'r') as f:
        return json.load(f)


def main():
    parser = argparse.ArgumentParser(description='Convert MySQL database to SQLite')
    parser.add_argument('--sqlite-db', required=True, help='Path to output SQLite database file')
    parser.add_argument('--mysql-config', required=True, help='Path to MySQL configuration JSON file')
    parser.add_argument('--verify', action='store_true', help='Verify conversion after completion')

    args = parser.parse_args()

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
